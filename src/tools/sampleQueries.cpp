// sample-queries: draw a diversity-maximizing subset of query genomes from the
// <assemblyList>.query.tsv manifest produced by makeBenchmarkSet.
//
// --number sets the total sample size and --ratio (6 comma-separated weights, in
// the order familyExcl,genusExcl,speciesExcl,subspeciesExcl,speciesIncl,
// subspeciesIncl) splits it across categories. Because makeBenchmarkSet excludes
// exactly one taxon per parent, every SubjectTaxID group in a category sits under
// a distinct parent, so "maximize diversity across taxa" reduces to round-robin
// over SubjectTaxID groups: we take one genome from each group before taking a
// second from any, which spreads the picks over as many parents as possible.
//
// The two inclusion categories are pairs: their weight/allocation counts *pairs*,
// and both members of a chosen pair are always emitted together (a pair is one
// unit). If a category has fewer units than its target, we take all of them and
// warn; the total may then fall short of --number.
//
// Given a required output prefix, it writes <prefix>.query.tsv (the sampled rows,
// same columns as the input) and <prefix>.summary (the input path, sampling
// parameters, the per-category sampling table, and any capping warnings).
#include "Parameters.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace {

const int NCAT = 6;
const char *CAT_NAMES[NCAT] = {
    "familyExclusion", "genusExclusion", "speciesExclusion",
    "subspeciesExclusion", "speciesInclusionPair", "subspeciesInclusionPair"};
const bool CAT_IS_PAIR[NCAT] = {false, false, false, false, true, true};

int catIndex(const string &c) {
    for (int i = 0; i < NCAT; i++) {
        if (c == CAT_NAMES[i]) return i;
    }
    return -1;
}

vector<string> splitTab(const string &s) {
    vector<string> out;
    size_t start = 0;
    while (true) {
        size_t pos = s.find('\t', start);
        if (pos == string::npos) { out.push_back(s.substr(start)); break; }
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

} // namespace

int sampleQueries(const Parameters &par) {
    const string inputPath = par.filenames[0];
    const string prefix = par.filenames[1];
    const string outputPath = prefix + ".query.tsv";
    const string summaryPath = prefix + ".summary";

    if (outputPath == inputPath) {
        cerr << "Error: output " << outputPath << " would overwrite the input; "
             << "choose a different prefix." << endl;
        return 1;
    }
    if (par.sampleNumber <= 0) {
        cerr << "Error: --number must be given and greater than 0." << endl;
        return 1;
    }

    // Parse --ratio (default: equal weights).
    vector<long> ratio(NCAT, 1);
    if (!par.sampleRatio.empty()) {
        vector<string> toks;
        stringstream ss(par.sampleRatio);
        string t;
        while (getline(ss, t, ',')) toks.push_back(t);
        if ((int) toks.size() != NCAT) {
            cerr << "Error: --ratio needs " << NCAT << " comma-separated values "
                 << "(familyExcl,genusExcl,speciesExcl,subspeciesExcl,"
                 << "speciesIncl,subspeciesIncl)." << endl;
            return 1;
        }
        for (int i = 0; i < NCAT; i++) {
            ratio[i] = atol(toks[i].c_str());
            if (ratio[i] < 0) { cerr << "Error: --ratio values must be >= 0." << endl; return 1; }
        }
    }
    long ratioSum = accumulate(ratio.begin(), ratio.end(), 0L);
    if (ratioSum <= 0) { cerr << "Error: --ratio sums to 0." << endl; return 1; }

    // Read the manifest, grouping each category's rows by SubjectTaxID (col 4).
    ifstream in(inputPath);
    if (!in.is_open()) { cerr << "Error: cannot open " << inputPath << endl; return 1; }
    string header;
    if (!getline(in, header)) { cerr << "Error: " << inputPath << " is empty." << endl; return 1; }

    vector<map<long, vector<string>>> byCat(NCAT); // [cat][subject] -> member lines
    string line;
    long skipped = 0;
    while (getline(in, line)) {
        if (line.empty()) continue;
        vector<string> f = splitTab(line);
        if (f.size() < 5) { skipped++; continue; }
        int ci = catIndex(f[1]);
        if (ci < 0) { skipped++; continue; }
        long subject = atol(f[4].c_str());
        byCat[ci][subject].push_back(line);
    }
    in.close();

    // Build per-category groups of "units". A unit is the list of output lines
    // taken together: one genome for exclusion categories, a whole pair for
    // inclusion categories. Each group corresponds to one SubjectTaxID (parent).
    vector<vector<vector<vector<string>>>> groups(NCAT); // [cat][group][unit][line]
    vector<long> available(NCAT, 0);
    for (int ci = 0; ci < NCAT; ci++) {
        for (auto &kv : byCat[ci]) {
            vector<vector<string>> g;
            if (CAT_IS_PAIR[ci]) {
                g.push_back(kv.second); // the pair is a single unit
            } else {
                for (auto &ln : kv.second) g.push_back(vector<string>{ln});
            }
            available[ci] += (long) g.size();
            groups[ci].push_back(std::move(g));
        }
    }

    // Split --number across categories by --ratio (largest-remainder so the
    // targets sum exactly to --number before capping).
    const long N = par.sampleNumber;
    vector<long> target(NCAT, 0);
    {
        vector<double> rem(NCAT);
        long allocated = 0;
        for (int i = 0; i < NCAT; i++) {
            double raw = (double) N * ratio[i] / ratioSum;
            long fl = (long) floor(raw);
            target[i] = fl;
            rem[i] = raw - fl;
            allocated += fl;
        }
        long leftover = N - allocated;
        vector<int> idx(NCAT);
        iota(idx.begin(), idx.end(), 0);
        sort(idx.begin(), idx.end(), [&](int a, int b) {
            if (rem[a] != rem[b]) return rem[a] > rem[b];
            return a < b;
        });
        for (long k = 0; k < leftover; k++) target[idx[k]]++;
    }

    // Sample each category: shuffle group order and members, then round-robin
    // across groups so distinct parents are covered before any parent repeats.
    mt19937 rng(par.randomSeed);
    vector<vector<vector<string>>> chosen(NCAT); // [cat] -> chosen units
    vector<long> sampled(NCAT, 0);
    vector<string> warnings;
    for (int ci = 0; ci < NCAT; ci++) {
        auto &G = groups[ci];
        shuffle(G.begin(), G.end(), rng);
        for (auto &g : G) shuffle(g.begin(), g.end(), rng);

        vector<vector<string> *> ordered;
        size_t maxUnits = 0;
        for (auto &g : G) maxUnits = max(maxUnits, g.size());
        for (size_t k = 0; k < maxUnits; k++) {
            for (auto &g : G) {
                if (k < g.size()) ordered.push_back(&g[k]);
            }
        }

        long want = target[ci];
        long take = min<long>(want, (long) ordered.size());
        for (long i = 0; i < take; i++) chosen[ci].push_back(*ordered[i]);
        sampled[ci] = take;

        if (take < want) {
            const char *unit = CAT_IS_PAIR[ci] ? " pair(s)" : " genome(s)";
            warnings.push_back(string(CAT_NAMES[ci]) + ": requested " + to_string(want) +
                               unit + " but only " + to_string(ordered.size()) +
                               " available; taking all.");
        }
    }

    // Write the sampled rows (same columns as the manifest).
    ofstream out(outputPath);
    if (!out.is_open()) { cerr << "Error: cannot write " << outputPath << endl; return 1; }
    out << header << "\n";
    long rowCount = 0;
    for (int ci = 0; ci < NCAT; ci++) {
        for (auto &unit : chosen[ci]) {
            for (auto &ln : unit) { out << ln << "\n"; rowCount++; }
        }
    }
    out.close();

    // Write the summary: provenance, per-category sampling table, and any
    // capping warnings.
    ofstream sum(summaryPath);
    if (!sum.is_open()) { cerr << "Error: cannot write " << summaryPath << endl; return 1; }
    sum << "# sample-queries summary\n";
    sum << "input\t" << inputPath << "\n";
    sum << "number\t" << par.sampleNumber << "\n";
    sum << "ratio\t" << (par.sampleRatio.empty() ? "1,1,1,1,1,1" : par.sampleRatio) << "\n";
    sum << "seed\t" << par.randomSeed << "\n";
    sum << "sampled_rows\t" << rowCount << "\n";
    if (skipped > 0) sum << "ignored_rows\t" << skipped << "\n";
    sum << "\n# per-category sampling\n";
    sum << "Category\tUnit\tAvailable\tRequested\tSampled\n";
    for (int ci = 0; ci < NCAT; ci++) {
        sum << CAT_NAMES[ci] << "\t" << (CAT_IS_PAIR[ci] ? "pair" : "genome")
            << "\t" << available[ci] << "\t" << target[ci] << "\t" << sampled[ci] << "\n";
    }
    sum << "\n# capping warnings\n";
    if (warnings.empty()) {
        sum << "(none)\n";
    } else {
        for (auto &w : warnings) sum << w << "\n";
    }
    sum.close();

    // Echo warnings to stderr and a short confirmation to stdout.
    for (auto &w : warnings) cerr << "Warning: " << w << endl;
    cout << "Sampled " << rowCount << " query row(s) -> " << outputPath << endl;
    cout << "Summary -> " << summaryPath << endl;
    return 0;
}
