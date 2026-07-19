// Lightweight replacement for Metabuli's LocalParameters / mmseqs Parameters.
// Holds only the options the three benchmark tools actually read, plus a small
// hand-rolled argument parser. No mmseqs Command/Parameters framework involved.
#ifndef BENCHTOOLS_PARAMETERS_H
#define BENCHTOOLS_PARAMETERS_H

#include <string>
#include <vector>

struct Parameters {
    // Positional arguments (input file paths), in order.
    std::vector<std::string> filenames;

    // --- grade ---
    int readIdCol = 1;          // --read-id-col : 0-based column with the read id
    int taxidCol = 2;           // --tax-id-col  : 0-based column with the taxon id
    int scoreCol = 0;           // --score-col   : 0-based column with the score
    int verbosity = 2;          // --verbosity   : 3 prints per-read decisions
    bool skipSecondary = false; // --skip-secondary
    bool topHitOnly = false;    // --top-hit-only
    std::string testRank;       // --rank        : comma-separated ranks (default set below)
    std::string printColumns;   // --print-cols  : comma-separated column indices to dump for TP/FP/FN

    // --- shared ---
    int threads = 1;            // --threads
    std::string testType = "gtdb"; // --test-type

    // --- split ---
    // Positional args: <assemblyList> <taxonomyDir> <acc2taxid> <outputPrefix>
    unsigned int randomSeed = 0;   // --seed
    bool skipValidation = false;   // --skip-validation : skip exclusion/inclusion checks

    // --- sample-queries ---
    int sampleNumber = 0;          // --number : total sample size across categories
    std::string sampleRatio;       // --ratio  : 6 comma-separated category weights
};

// Parse argv (starting at `start`) into `par`. Recognized flags are consumed;
// everything else is collected into par.filenames in order. Returns false and
// fills `err` on an unknown flag or a flag missing its value.
bool parseArguments(Parameters &par, int argc, const char **argv, int start, std::string &err);

#endif // BENCHTOOLS_PARAMETERS_H
