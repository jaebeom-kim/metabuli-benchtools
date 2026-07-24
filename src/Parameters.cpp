#include "Parameters.h"

#include <cstdlib>
#include <string>

namespace {

// Split "--flag=value" into name and inline value. Returns true if an inline
// value was present.
bool splitInlineValue(const std::string &arg, std::string &name, std::string &value) {
    size_t eq = arg.find('=');
    if (eq == std::string::npos) {
        name = arg;
        return false;
    }
    name = arg.substr(0, eq);
    value = arg.substr(eq + 1);
    return true;
}

} // namespace

bool parseArguments(Parameters &par, int argc, const char **argv, int start, std::string &err) {
    for (int i = start; i < argc; ++i) {
        std::string arg = argv[i];

        // Positional argument.
        if (arg.size() < 2 || arg[0] != '-') {
            par.filenames.push_back(arg);
            continue;
        }

        std::string name, inlineValue;
        bool hasInline = splitInlineValue(arg, name, inlineValue);

        // Boolean flags (no value).
        if (name == "--skip-secondary") { par.skipSecondary = true; continue; }
        if (name == "--top-hit-only")   { par.topHitOnly = true; continue; }
        if (name == "--skip-validation") { par.skipValidation = true; continue; }
        if (name == "--score-summary")  { par.scoreSummary = true; continue; }
        if (name == "--score-hist")     { par.scoreHist = true; continue; }

        // Value flags: fetch the value from inline (--flag=value) or the next token.
        auto nextValue = [&](std::string &out) -> bool {
            if (hasInline) { out = inlineValue; return true; }
            if (i + 1 >= argc) { err = "Missing value for " + name; return false; }
            out = argv[++i];
            return true;
        };

        std::string value;
        if (name == "--read-id-col") { if (!nextValue(value)) return false; par.readIdCol = std::atoi(value.c_str()); }
        else if (name == "--tax-id-col") { if (!nextValue(value)) return false; par.taxidCol = std::atoi(value.c_str()); }
        else if (name == "--score-col") { if (!nextValue(value)) return false; par.scoreCol = std::atoi(value.c_str()); }
        else if (name == "--score-bins") { if (!nextValue(value)) return false; par.scoreBins = std::atoi(value.c_str()); }
        else if (name == "--verbosity") { if (!nextValue(value)) return false; par.verbosity = std::atoi(value.c_str()); }
        else if (name == "--rank") { if (!nextValue(value)) return false; par.testRank = value; }
        else if (name == "--print-cols") { if (!nextValue(value)) return false; par.printColumns = value; }
        else if (name == "--threads") { if (!nextValue(value)) return false; par.threads = std::atoi(value.c_str()); }
        else if (name == "--test-type") { if (!nextValue(value)) return false; par.testType = value; }
        else if (name == "--seed") { if (!nextValue(value)) return false; par.randomSeed = (unsigned int) std::strtoul(value.c_str(), nullptr, 10); }
        else if (name == "--number") { if (!nextValue(value)) return false; par.sampleNumber = std::atoi(value.c_str()); }
        else if (name == "--ratio") { if (!nextValue(value)) return false; par.sampleRatio = value; }
        else if (name == "--min-abundance") { if (!nextValue(value)) return false; par.minAbundance = std::atof(value.c_str()); }
        else if (name == "--filter") { if (!nextValue(value)) return false; par.filterPercent = std::atof(value.c_str()); }
        else { err = "Unknown option: " + name; return false; }
    }
    return true;
}
