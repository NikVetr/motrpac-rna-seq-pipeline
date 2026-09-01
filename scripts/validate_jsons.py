"""
Usage:
python3 scripts/validate_jsons.py first.json second.json
"""
import argparse
import json
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description="script to compare and validate json files")
    parser.add_argument("infile1", type=str, help="Input json filename")
    parser.add_argument("infile2", type=str, help="Name of the second input json")
    args = parser.parse_args(argv)

    try:
        with open(args.infile1, encoding="utf-8") as json_file1:
            file_a = json.load(json_file1)

        with open(args.infile2, encoding="utf-8") as json_file2:
            file_b = json.load(json_file2)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        print("JSON comparison failed: {}".format(exc), file=sys.stderr)
        return 2

    a, b = json.dumps(file_a, sort_keys=True), json.dumps(file_b, sort_keys=True)

    print("Validation results")
    if a == b:
        print("Two jsons are identical")
        return 0
    print("Two jsons don't match")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
