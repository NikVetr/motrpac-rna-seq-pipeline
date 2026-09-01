#!/usr/bin/env python3
"""Propagate UMI-tools representatives to every matching STAR alignment."""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator


ACGT_UMI = re.compile(r"^[ACGT]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representatives-bam", required=True)
    parser.add_argument("--genomic-bam", type=Path, required=True)
    parser.add_argument("--transcriptome-bam", type=Path, required=True)
    parser.add_argument("--genomic-output", type=Path, required=True)
    parser.add_argument("--transcriptome-output", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--umi-length", type=int, required=True)
    parser.add_argument("--representation", choices=("rx_v1",), required=True)
    parser.add_argument("--container", required=True)
    parser.add_argument("--batch-records", type=int, default=20_000)
    return parser.parse_args()


def batched(values: list[str], size: int = 500) -> Iterator[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def umi_from_record(record, representation: str, umi_length: int) -> str:
    if representation != "rx_v1":
        raise ValueError(f"unsupported UMI representation: {representation}")
    if not record.has_tag("RX"):
        raise ValueError(f"representative lacks RX tag: {record.query_name}")
    umi = str(record.get_tag("RX")).upper()
    if len(umi) != umi_length or not ACGT_UMI.fullmatch(umi):
        raise ValueError(
            f"representative has invalid {umi_length}-base UMI {umi!r}: "
            f"{record.query_name}"
        )
    return umi


def valid_qname(query_name: str | None) -> str:
    if not query_name or any(character.isspace() for character in query_name):
        raise ValueError(f"invalid representative QNAME: {query_name!r}")
    return query_name


def configure_database(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=MEMORY;
        PRAGMA locking_mode=EXCLUSIVE;
        CREATE TABLE qnames (
            qname TEXT PRIMARY KEY,
            selected INTEGER NOT NULL DEFAULT 0 CHECK (selected IN (0, 1)),
            representative_records INTEGER NOT NULL DEFAULT 0,
            representative_r1 INTEGER NOT NULL DEFAULT 0,
            representative_r2 INTEGER NOT NULL DEFAULT 0,
            representative_umi TEXT,
            representative_umi_mismatch INTEGER NOT NULL DEFAULT 0,
            genomic_records INTEGER NOT NULL DEFAULT 0,
            transcriptome_records INTEGER NOT NULL DEFAULT 0
        ) WITHOUT ROWID;
        """
    )


def collect_representatives(
    connection: sqlite3.Connection,
    representatives_bam: str,
    representation: str,
    umi_length: int,
) -> dict[str, int]:
    import pysam

    counters: Counter[str] = Counter()
    with pysam.AlignmentFile(representatives_bam, "rb") as source:
        with connection:
            for record in source.fetch(until_eof=True):
                counters["representative_alignment_records"] += 1
                qname = valid_qname(record.query_name)
                if record.is_secondary or record.is_supplementary:
                    raise ValueError(
                        f"representative is secondary/supplementary: {qname}"
                    )
                if (
                    not record.is_paired
                    or record.is_unmapped
                    or record.mate_is_unmapped
                    or not record.is_proper_pair
                ):
                    raise ValueError(
                        f"representative is not a mapped proper pair: {qname}"
                    )
                if record.is_read1 == record.is_read2:
                    raise ValueError(
                        f"representative lacks exactly one mate flag: {qname}"
                    )
                umi = umi_from_record(record, representation, umi_length)
                connection.execute(
                    """
                    INSERT INTO qnames (
                        qname, selected, representative_records,
                        representative_r1, representative_r2, representative_umi
                    ) VALUES (?, 1, 1, ?, ?, ?)
                    ON CONFLICT(qname) DO UPDATE SET
                        representative_records = representative_records + 1,
                        representative_r1 = representative_r1 + excluded.representative_r1,
                        representative_r2 = representative_r2 + excluded.representative_r2,
                        representative_umi_mismatch =
                            representative_umi_mismatch OR
                            representative_umi != excluded.representative_umi
                    """,
                    (qname, int(record.is_read1), int(record.is_read2), umi),
                )

    selected = connection.execute(
        "SELECT COUNT(*) FROM qnames WHERE selected = 1"
    ).fetchone()[0]
    if selected == 0:
        raise ValueError("UMI-tools emitted no representative templates")
    invalid = connection.execute(
        """
        SELECT qname, representative_records, representative_r1,
               representative_r2, representative_umi_mismatch
        FROM qnames
        WHERE selected = 1 AND (
            representative_records != 2 OR representative_r1 != 1 OR
            representative_r2 != 1 OR representative_umi_mismatch != 0
        )
        ORDER BY qname LIMIT 5
        """
    ).fetchall()
    if invalid:
        raise ValueError(
            "duplicate or invalid representative template records: " + repr(invalid)
        )
    if counters["representative_alignment_records"] != 2 * selected:
        raise ValueError("representative record/template denominator mismatch")
    counters["selected_representative_qnames"] = selected
    return dict(counters)


def selected_qnames(
    connection: sqlite3.Connection, query_names: Iterable[str]
) -> set[str]:
    names = sorted(set(query_names))
    selected: set[str] = set()
    for batch in batched(names):
        placeholders = ",".join("?" for _ in batch)
        selected.update(
            row[0]
            for row in connection.execute(
                f"SELECT qname FROM qnames WHERE selected = 1 "
                f"AND qname IN ({placeholders})",
                batch,
            )
        )
    return selected


def record_semantics(record) -> tuple[str, ...]:
    labels: list[str] = []
    if record.is_secondary:
        labels.append("secondary")
    if record.is_supplementary:
        labels.append("supplementary")
    if record.is_reverse:
        labels.append("reverse")
    if record.is_read1:
        labels.append("r1")
    if record.is_read2:
        labels.append("r2")
    if record.cigartuples and any(
        operation == 3 for operation, _ in record.cigartuples
    ):
        labels.append("spliced")
    return tuple(labels)


def filter_stream(
    connection: sqlite3.Connection,
    source_path: Path,
    output_path: Path,
    stream: str,
    batch_records: int,
) -> dict[str, int]:
    import pysam

    if stream not in {"genomic", "transcriptome"}:
        raise ValueError(f"invalid stream: {stream}")
    column = f"{stream}_records"
    counters: Counter[str] = Counter()

    def process(records: list) -> None:
        if not records:
            return
        names = Counter(record.query_name for record in records)
        if None in names or "" in names:
            raise ValueError(f"{stream} BAM contains an empty QNAME")
        with connection:
            connection.executemany(
                f"""
                INSERT INTO qnames(qname, {column}) VALUES (?, ?)
                ON CONFLICT(qname) DO UPDATE SET
                    {column} = {column} + excluded.{column}
                """,
                names.items(),
            )
        selected = selected_qnames(connection, names)
        for record in records:
            counters["source_alignment_records"] += 1
            if record.query_name in selected:
                output.write(record)
                counters["kept_alignment_records"] += 1
                for label in record_semantics(record):
                    counters[f"kept_{label}_records"] += 1
            else:
                counters["dropped_alignment_records"] += 1

    with pysam.AlignmentFile(str(source_path), "rb") as source:
        with pysam.AlignmentFile(str(output_path), "wb", template=source) as output:
            records: list = []
            for record in source.fetch(until_eof=True):
                records.append(record)
                if len(records) == batch_records:
                    process(records)
                    records = []
            process(records)
    if counters["source_alignment_records"] == 0:
        raise ValueError(f"{stream} BAM contains no alignment records")
    if (
        counters["kept_alignment_records"] + counters["dropped_alignment_records"]
        != counters["source_alignment_records"]
    ):
        raise ValueError(f"{stream} kept/dropped alignment denominators disagree")
    quickcheck = pysam.quickcheck("-v", str(output_path))
    if quickcheck:
        raise ValueError(f"invalid filtered {stream} BAM: {quickcheck}")
    return dict(counters)


def scalar(connection: sqlite3.Connection, query: str) -> int:
    value = connection.execute(query).fetchone()[0]
    if not isinstance(value, int):
        raise ValueError(f"noninteger SQLite metric from: {query}")
    return value


def qname_metrics(connection: sqlite3.Connection) -> dict[str, int]:
    metrics = {
        "source_union_unique_qnames": scalar(connection, "SELECT COUNT(*) FROM qnames"),
        "genomic_source_unique_qnames": scalar(
            connection, "SELECT COUNT(*) FROM qnames WHERE genomic_records > 0"
        ),
        "transcriptome_source_unique_qnames": scalar(
            connection, "SELECT COUNT(*) FROM qnames WHERE transcriptome_records > 0"
        ),
        "selected_representative_qnames": scalar(
            connection, "SELECT COUNT(*) FROM qnames WHERE selected = 1"
        ),
        "selected_genomic_present_qnames": scalar(
            connection,
            "SELECT COUNT(*) FROM qnames WHERE selected = 1 AND genomic_records > 0",
        ),
        "selected_transcriptome_present_qnames": scalar(
            connection,
            "SELECT COUNT(*) FROM qnames WHERE selected = 1 AND transcriptome_records > 0",
        ),
        "transcriptome_qnames_absent_from_genomic_source": scalar(
            connection,
            "SELECT COUNT(*) FROM qnames WHERE transcriptome_records > 0 "
            "AND genomic_records = 0",
        ),
        "selected_genomic_alignment_records": scalar(
            connection,
            "SELECT COALESCE(SUM(genomic_records), 0) FROM qnames WHERE selected = 1",
        ),
        "selected_transcriptome_alignment_records": scalar(
            connection,
            "SELECT COALESCE(SUM(transcriptome_records), 0) FROM qnames "
            "WHERE selected = 1",
        ),
        "selected_qname_logical_bytes": scalar(
            connection,
            "SELECT COALESCE(SUM(LENGTH(CAST(qname AS BLOB))), 0) "
            "FROM qnames WHERE selected = 1",
        ),
        "source_union_qname_logical_bytes": scalar(
            connection,
            "SELECT COALESCE(SUM(LENGTH(CAST(qname AS BLOB))), 0) FROM qnames",
        ),
    }
    metrics["selected_genomic_absent_qnames"] = (
        metrics["selected_representative_qnames"]
        - metrics["selected_genomic_present_qnames"]
    )
    metrics["selected_transcriptome_absent_qnames"] = (
        metrics["selected_representative_qnames"]
        - metrics["selected_transcriptome_present_qnames"]
    )
    return metrics


def validate_denominators(
    qnames: dict[str, int], genomic: dict[str, int], transcriptome: dict[str, int]
) -> None:
    selected = qnames["selected_representative_qnames"]
    if qnames["selected_genomic_present_qnames"] != selected:
        raise ValueError(
            f"{qnames['selected_genomic_absent_qnames']} representative QNAMEs are "
            "absent from the original genomic BAM"
        )
    if (
        qnames["selected_transcriptome_present_qnames"]
        + qnames["selected_transcriptome_absent_qnames"]
        != selected
    ):
        raise ValueError("transcriptome representative present/absent counts disagree")
    if qnames["transcriptome_qnames_absent_from_genomic_source"]:
        raise ValueError(
            "STAR transcriptome BAM contains QNAMEs absent from its genomic BAM: "
            f"{qnames['transcriptome_qnames_absent_from_genomic_source']}"
        )
    if genomic["kept_alignment_records"] != qnames[
        "selected_genomic_alignment_records"
    ]:
        raise ValueError("not every selected genomic alignment was retained")
    if transcriptome["kept_alignment_records"] != qnames[
        "selected_transcriptome_alignment_records"
    ]:
        raise ValueError("not every selected transcriptome alignment was retained")


def main() -> int:
    import pysam

    args = parse_args()
    if args.umi_length < 1 or args.batch_records < 1:
        raise ValueError("UMI length and batch size must be positive")
    for path in (args.genomic_bam, args.transcriptome_bam):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"input BAM does not exist or is empty: {path}")
    for path in (
        args.genomic_output,
        args.transcriptome_output,
        args.metrics,
        args.database,
    ):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite output: {path}")

    connection = sqlite3.connect(args.database)
    try:
        configure_database(connection)
        representatives = collect_representatives(
            connection,
            args.representatives_bam,
            args.representation,
            args.umi_length,
        )
        genomic = filter_stream(
            connection,
            args.genomic_bam,
            args.genomic_output,
            "genomic",
            args.batch_records,
        )
        transcriptome = filter_stream(
            connection,
            args.transcriptome_bam,
            args.transcriptome_output,
            "transcriptome",
            args.batch_records,
        )
        qnames = qname_metrics(connection)
        validate_denominators(qnames, genomic, transcriptome)
        page_size = scalar(connection, "PRAGMA page_size")
        page_count = scalar(connection, "PRAGMA page_count")
        connection.commit()
    finally:
        connection.close()

    metrics = {
        "schema_version": 1,
        "acceptance_class": "scientific-truth",
        "status": "production-shadow",
        "algorithm": "propagate_genomic_umi_representative_qnames_v1",
        "container": args.container,
        "runtime_versions": {
            "pysam": pysam.__version__,
            "sqlite": sqlite3.sqlite_version,
        },
        "contract": {
            "family_source": (
                "UMI-tools primary mapped proper-pair genomic alignments"
            ),
            "genomic": (
                "retain every original alignment for each representative QNAME"
            ),
            "transcriptome": (
                "retain every STAR transcript alternative for each representative "
                "QNAME present; measure transcript-incompatible representatives as absent"
            ),
            "ordering": "preserve original BAM header and alignment-record order",
            "qname_index": (
                "exact task-local SQLite TEXT primary key; never uploaded"
            ),
        },
        "representatives": representatives,
        "qnames": qnames,
        "genomic": genomic,
        "transcriptome": transcriptome,
        "task_local_storage": {
            "sqlite_database_bytes": args.database.stat().st_size,
            "sqlite_allocated_bytes": page_size * page_count,
            "filtered_genomic_bam_bytes": args.genomic_output.stat().st_size,
            "filtered_transcriptome_bam_bytes": args.transcriptome_output.stat().st_size,
        },
    }
    args.metrics.write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.database.unlink()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FileExistsError,
        FileNotFoundError,
        OSError,
        sqlite3.Error,
        ValueError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
