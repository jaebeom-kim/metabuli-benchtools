#include "Parameters.h"

#include <iostream>
#include <string>
#include <vector>

// Tool entry points (defined in src/tools/*.cpp).
int grade(const Parameters &par);
int split(const Parameters &par);
int sampleQueries(const Parameters &par);

namespace {

struct Tool {
    std::string name;
    int minPositionals;
    std::string usage;
    int (*run)(const Parameters &);
};

const std::vector<Tool> TOOLS = {
    {"grade", 3,
     "grade <classificationList> <mappingList> <taxonomyDir> [options]\n"
     "    Score classification results against answer sheets.\n"
     "    <classificationList>  file listing one classification-result path per line\n"
     "    <mappingList>         file listing one accession->taxid answer-sheet path per line (aligned to the above)\n"
     "    <taxonomyDir>         directory with names.dmp, nodes.dmp, merged.dmp\n"
     "  Options:\n"
     "    --test-type STR       gtdb | gtdb-amgsim | cami | cami-long | cami-euk | hiv | hiv-ex | over | kapk (default: gtdb)\n"
     "    --rank STR            comma-separated ranks to score (default: class,order,family,genus,species)\n"
     "    --read-id-col INT     0-based column holding the read id (default: 1)\n"
     "    --tax-id-col INT      0-based column holding the taxon id (default: 2)\n"
     "    --score-col INT       0-based column holding the score (default: 0)\n"
     "    --print-cols STR      comma-separated column indices to dump into .tp/.fp/.fn files\n"
     "    --skip-secondary      count only the first classification per read (GTDB only)\n"
     "    --verbosity INT       3 prints per-read decisions (default: 2)\n"
     "    --threads INT         number of OpenMP threads (default: 1)",
     grade},

    {"split", 4,
     "split <assemblyList> <taxonomyDir> <acc2taxid> <outputPrefix> [options]\n"
     "    Build exclusion + inclusion benchmark sets from an assembly list.\n"
     "    Emits the family/genus/species/subspecies exclusion sets, the database\n"
     "    assembly list, and the species/subspecies inclusion query pairs.\n"
     "    <assemblyList>        file listing one assembly accession per line\n"
     "    <taxonomyDir>         directory with names.dmp, nodes.dmp, merged.dmp\n"
     "    <acc2taxid>           assembly accession -> taxid mapping file\n"
     "    <outputPrefix>        output file prefix\n"
     "  Options:\n"
     "    --seed INT            random seed (default: 0)\n"
     "    --skip-validation     skip the exclusion/inclusion validation checks\n"
     "  Outputs: <outputPrefix>.database, <outputPrefix>.query.tsv, <outputPrefix>.summary",
     split},

    {"sample-queries", 2,
     "sample-queries <queryTsv> <outPrefix> [options]\n"
     "    Sample a diversity-maximizing subset of query genomes from a\n"
     "    split .query.tsv manifest (round-robin over SubjectTaxID).\n"
     "    <queryTsv>            the .query.tsv manifest from split\n"
     "    <outPrefix>           output prefix (required)\n"
     "  Options:\n"
     "    --number INT          total sample size across categories (required)\n"
     "    --ratio a,b,c,d,e,f   6 category weights: familyExcl,genusExcl,speciesExcl,\n"
     "                          subspeciesExcl,speciesIncl,subspeciesIncl (default: all 1);\n"
     "                          the inclusion weights count pairs (both members kept)\n"
     "    --seed INT            random seed (default: 0)\n"
     "  Outputs: <outPrefix>.query.tsv, <outPrefix>.summary",
     sampleQueries},
};

void printUsage(const std::string &program) {
    std::cerr << "metabuli-benchtools -- standalone benchmark & grading tools\n\n";
    std::cerr << "Usage: " << program << " <command> [args...]\n\n";
    std::cerr << "Commands:\n";
    for (const auto &tool : TOOLS) {
        std::cerr << "  " << tool.name << "\n";
    }
    std::cerr << "\nRun '" << program << " <command> --help' for command-specific usage.\n";
}

} // namespace

int main(int argc, const char **argv) {
    const std::string program = "benchtools";

    if (argc < 2) {
        printUsage(program);
        return 1;
    }

    std::string command = argv[1];
    if (command == "-h" || command == "--help" || command == "help") {
        printUsage(program);
        return 0;
    }

    for (const auto &tool : TOOLS) {
        if (tool.name != command) {
            continue;
        }

        // Command-specific help.
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "-h" || a == "--help") {
                std::cerr << "Usage: " << program << " " << tool.usage << "\n";
                return 0;
            }
        }

        Parameters par;
        std::string err;
        if (!parseArguments(par, argc, argv, 2, err)) {
            std::cerr << "Error: " << err << "\n\n";
            std::cerr << "Usage: " << program << " " << tool.usage << "\n";
            return 1;
        }

        if ((int) par.filenames.size() < tool.minPositionals) {
            std::cerr << "Error: " << tool.name << " expects at least " << tool.minPositionals
                      << " positional argument(s), got " << par.filenames.size() << ".\n\n";
            std::cerr << "Usage: " << program << " " << tool.usage << "\n";
            return 1;
        }

        return tool.run(par);
    }

    std::cerr << "Error: unknown command '" << command << "'\n\n";
    printUsage(program);
    return 1;
}
