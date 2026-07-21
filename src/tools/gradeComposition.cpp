// grade-composition: evaluate abundance / profiling estimates against MGSIM
// ground-truth abundances, for multiple tools over multiple communities.
//
// <profileList> is a 5-column TSV: group (tool name), community, profile path,
// taxid column, read-count column. The last two are 1-based column indices into
// that report (so different tools can use different layouts). The `community`
// keys into the MGSIM truth's `Community` column. Each profile is graded against
// its community, per rank; results are then grouped by `group` and aggregated
// across that group's communities into a mean and standard deviation.
//
// MGSIM `communities` writes the truth (columns: Community, Taxon, Perc_rel_abund,
// Rank); Taxon is the assembly accession (benchtools_to_mgsim.py labels genomes by
// accession). Two flavours exist: *_abund.txt (cell / genome-copy fraction) and
// *_wAbund.txt (DNA-pool fraction) -- pass whichever matches the profiler.
//
// The estimate is a report with a taxID column and a *per-taxon* read-count column
// (e.g. Metabuli's taxon_count, not the cumulative clade_count). Each taxon's
// directly-assigned reads are rolled up to the target rank. The truth is compared
// against the
// *ideal (best-possible) profile*: <queryTsv> (split/sample-queries .query.tsv)
// gives each genome its taxid (QueryTaxID) and ExpectedRank -- the deepest rank
// reportable given the database -- so a held-out (exclusion) genome only counts at
// ranks at or coarser than its ExpectedRank, never as a phantom species.
//
// Per (group, rank), the summary reports the mean and sample SD across communities
// of: L1 distance, Bray-Curtis dissimilarity (= L1/2), Purity (precision), and
// Completeness (recall). A taxon counts as predicted-present when its estimated
// fraction exceeds --min-abundance.
#include "TaxonomyWrapper.h"
#include "Parameters.h"
#include "Util.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

namespace {

// Split on tab without collapsing empty fields (unlike Util::split).
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

// Per query genome: its own taxid and the deepest rank reportable given the
// database (from split/sample-queries .query.tsv).
struct QueryMeta {
    TaxID taxid;
    string expectedRank;
};

// Load accession -> {QueryTaxID, ExpectedRank} from a .query.tsv manifest
// (columns: Accession, Category, ExpectedRank, QueryTaxID, ...), with a
// version-stripped fallback key.
unordered_map<string, QueryMeta> loadQueryMeta(const string &path, bool &ok) {
    unordered_map<string, QueryMeta> m;
    ifstream in(path);
    if (!in.is_open()) { ok = false; return m; }
    string line;
    bool first = true;
    while (getline(in, line)) {
        if (line.empty()) continue;
        vector<string> f = splitTab(line);
        if (f.size() < 4) continue;
        if (first) { first = false; if (f[0] == "Accession") continue; } // skip header
        char *end = nullptr;
        long taxid = strtol(f[3].c_str(), &end, 10);
        if (end == f[3].c_str()) continue;
        m[f[0]] = QueryMeta{(TaxID) taxid, f[2]};
    }
    unordered_map<string, QueryMeta> noVersion;
    for (auto &it : m) {
        size_t dot = it.first.find('.');
        if (dot != string::npos) noVersion[it.first.substr(0, dot)] = it.second;
    }
    for (auto &it : noVersion) m.emplace(it.first, it.second);
    ok = true;
    return m;
}

// The taxon's ancestor exactly at `rank`, or 0 if its lineage has no node at that
// rank (e.g. the taxon is coarser than `rank`). getTaxIdAtRank alone returns the
// taxon itself when it sits above `rank`, so we confirm the result's rank matches.
TaxID ancestorAtRank(TaxonomyWrapper &tax, TaxID t, const string &rank) {
    TaxID r = tax.getTaxIdAtRank(t, rank);
    if (r == 0) return 0;
    const TaxonNode *node = tax.taxonNode(r);
    if (node == nullptr || rank != tax.getString(node->rankIdx)) return 0;
    return r;
}

struct CompMetrics {
    double l1 = 0, bray = 0, purity = 0, completeness = 0;
};

// Grade one profile against one community's truth at one rank.
CompMetrics gradeOne(TaxonomyWrapper &tax,
                     unordered_map<string, QueryMeta> &queryMeta,
                     const vector<pair<string, double>> &truthRows,
                     const unordered_map<TaxID, double> &taxonReads,
                     const string &rank, double minAbund,
                     set<string> &missingAcc) {
    const int rankIdx = tax.findRankIndex(rank);

    // Ideal-profile truth: a genome counts only if reportable at this rank (its
    // ExpectedRank is at least as deep), rolled up to its ancestor at the rank.
    map<TaxID, double> truthAtRank;
    double truthTotal = 0;
    for (auto &tp : truthRows) {
        auto it = queryMeta.find(tp.first);
        if (it == queryMeta.end()) {
            size_t dot = tp.first.find('.');
            if (dot != string::npos) it = queryMeta.find(tp.first.substr(0, dot));
        }
        if (it == queryMeta.end()) { missingAcc.insert(tp.first); continue; }
        int expIdx = tax.findRankIndex(it->second.expectedRank);
        if (expIdx < 0 || rankIdx < expIdx) continue;
        TaxID r = ancestorAtRank(tax, it->second.taxid, rank);
        if (r == 0) continue;
        truthAtRank[r] += tp.second;
        truthTotal += tp.second;
    }

    // Estimate: roll each taxon's directly-assigned reads up to the rank.
    map<TaxID, double> estAtRank;
    double estTotal = 0;
    for (auto &kv : taxonReads) {
        TaxID r = ancestorAtRank(tax, kv.first, rank);
        if (r == 0) continue;
        estAtRank[r] += kv.second;
        estTotal += kv.second;
    }

    set<TaxID> uni;
    for (auto &kv : truthAtRank) uni.insert(kv.first);
    for (auto &kv : estAtRank) uni.insert(kv.first);

    double l1 = 0;
    int tp = 0, fp = 0, fn = 0;
    for (TaxID t : uni) {
        double tv = truthTotal > 0 && truthAtRank.count(t) ? truthAtRank[t] / truthTotal : 0.0;
        double ev = estTotal   > 0 && estAtRank.count(t)   ? estAtRank[t]   / estTotal   : 0.0;
        l1 += fabs(ev - tv);
        bool inTruth = tv > 0;
        bool inEst = ev > minAbund;
        if (inTruth && inEst) tp++;
        else if (inEst) fp++;
        else if (inTruth) fn++;
    }

    CompMetrics m;
    m.l1 = l1;
    m.bray = l1 / 2.0;
    m.purity = (tp + fp) > 0 ? (double) tp / (tp + fp) : 0.0;
    m.completeness = (tp + fn) > 0 ? (double) tp / (tp + fn) : 0.0;
    return m;
}

// Mean and sample standard deviation of a series.
pair<double, double> meanSD(const vector<double> &xs) {
    if (xs.empty()) return {0.0, 0.0};
    double sum = 0;
    for (double x : xs) sum += x;
    double mean = sum / xs.size();
    if (xs.size() < 2) return {mean, 0.0};
    double ss = 0;
    for (double x : xs) ss += (x - mean) * (x - mean);
    return {mean, sqrt(ss / (xs.size() - 1))};
}

} // namespace

int gradeComposition(const Parameters &par) {
    const string profileList = par.filenames[0];
    const string truthFile   = par.filenames[1];
    const string queryTsv    = par.filenames[2];
    const string taxDir      = par.filenames[3];

    vector<string> ranks = par.testRank.empty()
        ? vector<string>{"species", "genus"}
        : Util::split(par.testRank, ",");
    const double minAbund = par.minAbundance;

    TaxonomyWrapper taxonomy(taxDir + "/names.dmp", taxDir + "/nodes.dmp",
                             taxDir + "/merged.dmp", false);

    bool ok;
    unordered_map<string, QueryMeta> queryMeta = loadQueryMeta(queryTsv, ok);
    if (!ok) { cerr << "Error: cannot open query.tsv file " << queryTsv << endl; return 1; }

    // MGSIM truth: Community -> list of (accession, Perc_rel_abund).
    unordered_map<string, vector<pair<string, double>>> truthByComm;
    {
        ifstream in(truthFile);
        if (!in.is_open()) { cerr << "Error: cannot open truth file " << truthFile << endl; return 1; }
        string line;
        bool first = true;
        while (getline(in, line)) {
            if (line.empty()) continue;
            vector<string> f = splitTab(line);
            if (f.size() < 3) continue;
            if (first) { first = false; if (f[0] == "Community") continue; } // skip header
            truthByComm[f[0]].push_back({f[1], atof(f[2].c_str())});
        }
    }
    if (truthByComm.empty()) { cerr << "Error: no truth rows parsed from " << truthFile << endl; return 1; }

    // Parse the 5-column profile list: group, community, profile path,
    // taxid column, read-count column (the last two are 1-based indices into the
    // report, so different tools can use different report layouts).
    struct Job { string group, community, path; int taxidCol, countCol; };
    vector<Job> jobs;
    vector<string> groupOrder;
    set<string> groupSeen;
    {
        ifstream in(profileList);
        if (!in.is_open()) { cerr << "Error: cannot open profile list " << profileList << endl; return 1; }
        string line;
        bool first = true;
        while (getline(in, line)) {
            if (line.empty()) continue;
            vector<string> f = splitTab(line);
            if (f.size() < 5) continue;
            if (first) { first = false; if (f[0] == "group" || f[0] == "Group") continue; } // optional header
            int taxidCol = atoi(f[3].c_str()) - 1; // 1-based -> 0-based
            int countCol = atoi(f[4].c_str()) - 1;
            if (taxidCol < 0 || countCol < 0) {
                cerr << "Error: taxid/count columns must be 1-based indices (>= 1): " << line << endl;
                return 1;
            }
            jobs.push_back(Job{f[0], f[1], f[2], taxidCol, countCol});
            if (!groupSeen.count(f[0])) { groupSeen.insert(f[0]); groupOrder.push_back(f[0]); }
        }
    }
    if (jobs.empty()) { cerr << "Error: no rows in profile list " << profileList << endl; return 1; }

    // group -> rank -> per-community metrics
    unordered_map<string, unordered_map<string, vector<CompMetrics>>> data;
    set<string> missingAcc;
    set<string> missingComm;

    for (auto &job : jobs) {
        auto tc = truthByComm.find(job.community);
        if (tc == truthByComm.end()) { missingComm.insert(job.community); continue; }

        // Parse the report using this row's taxid/count columns; roll per-taxon
        // (directly-assigned) reads up to the target rank later.
        unordered_map<TaxID, double> taxonReads;
        {
            ifstream in(job.path);
            if (!in.is_open()) { cerr << "Error: cannot open profile " << job.path << endl; return 1; }
            const int taxidCol = job.taxidCol, countCol = job.countCol;
            string line;
            while (getline(in, line)) {
                if (line.empty() || line[0] == '#') continue; // skip blank / header lines
                vector<string> f = splitTab(line);
                if ((int) f.size() <= max(taxidCol, countCol)) continue;
                char *end = nullptr;
                long taxid = strtol(f[taxidCol].c_str(), &end, 10);
                if (end == f[taxidCol].c_str()) continue; // non-numeric taxID (header)
                double reads = atof(f[countCol].c_str());
                if (taxid <= 0 || reads <= 0) continue;
                taxonReads[(TaxID) taxid] += reads;
            }
        }

        for (const string &rank : ranks) {
            data[job.group][rank].push_back(
                gradeOne(taxonomy, queryMeta, tc->second, taxonReads, rank, minAbund, missingAcc));
        }
    }

    // Per-group summary: mean + sample SD across communities, per rank.
    printf("Group\tRank\tN\tL1_mean\tL1_sd\tBrayCurtis_mean\tBrayCurtis_sd"
           "\tPurity_mean\tPurity_sd\tCompleteness_mean\tCompleteness_sd\n");
    for (const string &group : groupOrder) {
        for (const string &rank : ranks) {
            const vector<CompMetrics> &v = data[group][rank];
            if (v.empty()) continue;
            vector<double> l1, bray, pur, comp;
            for (const CompMetrics &m : v) {
                l1.push_back(m.l1); bray.push_back(m.bray);
                pur.push_back(m.purity); comp.push_back(m.completeness);
            }
            auto ml1 = meanSD(l1), mbray = meanSD(bray), mpur = meanSD(pur), mcomp = meanSD(comp);
            printf("%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                   group.c_str(), rank.c_str(), (int) v.size(),
                   ml1.first, ml1.second, mbray.first, mbray.second,
                   mpur.first, mpur.second, mcomp.first, mcomp.second);
        }
    }

    if (!missingComm.empty()) {
        cerr << "Warning: " << missingComm.size() << " community id(s) in the profile list "
             << "not found in the truth; those rows were skipped." << endl;
    }
    if (!missingAcc.empty()) {
        cerr << "Warning: " << missingAcc.size() << " truth accession(s) not found in "
             << queryTsv << "; their abundance was ignored." << endl;
    }
    return 0;
}
