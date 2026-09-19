#ifndef RSEM_CONVERGENCE_H
#define RSEM_CONVERGENCE_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "GroupInfo.h"
#include "Transcripts.h"

// Counts are evaluated at each theta, including RSEM's existing final E-step.
inline void writeConvergence(const char* prefix, const char* reference,
        Transcripts& transcripts, int round, double tolerance,
        const std::vector<double>& previous, const std::vector<double>& current,
        const std::vector<double>& previousCounts, const double* counts) {
    std::ofstream transcriptFile, geneFile;
    transcriptFile.exceptions(std::ios::failbit | std::ios::badbit);
    geneFile.exceptions(std::ios::failbit | std::ios::badbit);
    transcriptFile.open(std::string(prefix) + ".rsem_convergence.tsv");
    geneFile.open(std::string(prefix) + ".rsem_gene_convergence.tsv");
    transcriptFile << std::setprecision(17);
    geneFile << std::setprecision(17);
    const char* values = "previous_theta\tfinal_theta\trelative_change\t"
        "expected_count_previous_theta\texpected_count_final_theta\texpected_count_change";
    transcriptFile << "component\ttranscript_id\tgene_id\titeration\t" << values << '\n';
    geneFile << "gene_id\titeration\tflagged_transcripts\t" << values
             << "\tsum_absolute_isoform_count_change\n";

    const auto flagged = [&](int i) {
        return previous[i] >= 1e-7 &&
            std::fabs(current[i] - previous[i]) / previous[i] >= tolerance;
    };
    const auto writeValues = [](std::ostream& out, double before, double after,
                               double beforeCount, double afterCount) {
        out << before << '\t' << after << '\t' << std::fabs(after - before) / before
            << '\t' << beforeCount << '\t' << afterCount << '\t' << afterCount - beforeCount;
    };
    for (int i = 0; i <= transcripts.getM(); ++i) {
        if (!flagged(i)) continue;
        transcriptFile << (i == 0 ? "background" : "transcript") << '\t'
            << (i == 0 ? "." : transcripts.getTranscriptAt(i).getTranscriptID()) << '\t'
            << (i == 0 ? "." : transcripts.getTranscriptAt(i).getGeneID()) << '\t' << round << '\t';
        writeValues(transcriptFile, previous[i], current[i], previousCounts[i], counts[i]);
        transcriptFile << '\n';
    }

    GroupInfo groups;
    groups.load((std::string(reference) + ".grp").c_str());
    for (int g = 0; g < groups.getm(); ++g) {
        int affected = 0;
        double before = 0, after = 0, beforeCount = 0, afterCount = 0, absoluteChange = 0;
        for (int i = groups.spAt(g); i < groups.spAt(g + 1); ++i) {
            affected += flagged(i);
            before += previous[i];
            after += current[i];
            beforeCount += previousCounts[i];
            afterCount += counts[i];
            absoluteChange += std::fabs(counts[i] - previousCounts[i]);
        }
        if (!affected) continue;
        geneFile << transcripts.getTranscriptAt(groups.spAt(g)).getGeneID() << '\t'
                 << round << '\t' << affected << '\t';
        writeValues(geneFile, before, after, beforeCount, afterCount);
        geneFile << '\t' << absoluteChange << '\n';
    }
    transcriptFile.close();
    geneFile.close();
}

#endif
