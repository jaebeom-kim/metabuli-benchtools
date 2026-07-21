// grade-classification: the grouped, multi-tool version of `grade`.
//
// <classificationList> is a 3-column TSV: group (tool), community, classification
// result path. Each result is graded per read against a single shared answer
// sheet (accession<TAB>taxid) using the same per-read logic as `grade` (parse the
// read id -> source accession via --test-type, look up the true taxid, compare
// predicted vs. true at each rank). Results are grouped by `group` and aggregated
// across communities into a mean and sample standard deviation of precision,
// sensitivity, and F1, per rank -- mirroring grade-composition.
#include "TaxonomyWrapper.h"
#include "Parameters.h"
#include "Util.h"
#include "gradeCommon.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

namespace {

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

// Parse a classification read id into the answer-sheet key (source accession),
// matching grade's per-test-type handling. Returns false if the id cannot be
// parsed (the read is then skipped).
bool parseReadKey(const string &readIdField, const string &testType, string &keyOut) {
    static const regex accWithVersion("(GC[AF]_[0-9]+\\.[0-9]+)");
    static const regex accNoVersion("(GC[AF]_[0-9]+)");
    smatch m;
    string id = readIdField;
    if (testType == "gtdb") {
        id = regex_search(id, m, accWithVersion) ? string(m[0]) : string();
        size_t pos = id.find('.');
        if (pos != string::npos) id = id.substr(0, pos);
    } else if (testType == "gtdb-amgsim") {
        if (!regex_search(id, m, accNoVersion)) {
            cerr << "Cannot parse GTDB-aMGSIM read ID: " << id << endl;
            return false;
        }
        id = m[0];
    } else if (testType == "hiv" || testType == "hiv-ex") {
        id = id.substr(0, id.find('_'));
    } else if (testType == "cami" || testType == "cami-long" || testType == "cami-euk") {
        id = id.substr(0, id.find('/'));
    } else if (testType == "over") {
        id = regex_search(id, m, accWithVersion) ? string(m[0]) : string();
    } else if (testType == "kapk") {
        size_t start = id.find("----");
        if (start == string::npos) { cerr << "Cannot parse KapK read ID: " << id << endl; return false; }
        size_t accStart = start + 4;
        size_t end = id.find("__", accStart);
        if (end == string::npos) { cerr << "Cannot parse KapK read ID: " << id << endl; return false; }
        id = id.substr(accStart, end - accStart);
    }
    keyOut = id;
    return true;
}

// Grade one classification file against the answer sheet; returns rank -> counts.
unordered_map<string, CountAtRank> gradeOneFile(const string &path,
                                                const unordered_map<string, int> &assacc2taxid,
                                                TaxonomyWrapper &tax,
                                                const vector<string> &ranks,
                                                const Parameters &par) {
    unordered_map<string, CountAtRank> counts;
    for (const string &r : ranks) counts[r] = CountAtRank{}; // zero-init

    ifstream in(path);
    if (!in.is_open()) { cerr << "Error: cannot open classification " << path << endl; return counts; }

    string line;
    unordered_map<string, int> observed; // for --skip-secondary
    const int idCol = par.readIdCol, taxCol = par.taxidCol;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        vector<string> fields = splitTab(line);
        if ((int) fields.size() <= max(idCol, taxCol)) continue;
        if (fields[taxCol].empty() || !isdigit((unsigned char) fields[taxCol][0])) continue;

        string fullId = fields[idCol];
        string key;
        if (!parseReadKey(fullId, par.testType, key)) continue;
        int classInt = stoi(fields[taxCol]);

        if (par.skipSecondary) {
            fullId = fullId.substr(0, fullId.find('/'));
            if (observed.find(fullId) == observed.end()) {
                if (classInt != 0) observed[fullId] = 1;
                else { observed[fullId] = 0; continue; }
            } else if (observed[fullId] == 1) {
                continue;
            }
        }

        int trueTaxid = 0;
        auto it = assacc2taxid.find(key);
        if (it != assacc2taxid.end()) trueTaxid = it->second;

        for (const string &rank : ranks) {
            if (par.testType == "over") {
                compareTaxon_overclassification(classInt, trueTaxid, tax, counts[rank], rank);
            } else if (par.testType == "hiv-ex") {
                compareTaxon_hivExclusion(classInt, 11676, counts[rank]);
            } else if (par.testType == "cami-euk") {
                compareTaxonAtRank_CAMI_euk(classInt, trueTaxid, tax, counts[rank], rank);
            } else {
                compareTaxonAtRank_CAMI(classInt, trueTaxid, tax, counts[rank], rank);
            }
        }
    }
    for (const string &r : ranks) counts[r].calculate();
    return counts;
}

double clean(float x) { return std::isfinite(x) ? (double) x : 0.0; }

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

struct PRF { double precision, sensitivity, f1; };

} // namespace

int gradeClassification(const Parameters &par) {
    const string classificationList = par.filenames[0];
    const string mappingFile        = par.filenames[1];
    const string taxDir             = par.filenames[2];

    vector<string> ranks = par.testRank.empty()
        ? vector<string>{"class", "order", "family", "genus", "species"}
        : Util::split(par.testRank, ",");

    TaxonomyWrapper taxonomy(taxDir + "/names.dmp", taxDir + "/nodes.dmp",
                             taxDir + "/merged.dmp", false);

    // Single shared answer sheet: accession -> taxid (version-stripped, as grade).
    unordered_map<string, int> assacc2taxid;
    {
        ifstream m(mappingFile);
        if (!m.is_open()) { cerr << "Error: cannot open mapping file " << mappingFile << endl; return 1; }
        string key, value;
        while (getline(m, key, '\t')) {
            if (!getline(m, value, '\n')) break;
            if (par.testType != "kapk") {
                size_t pos = key.find('.');
                if (pos != string::npos) key = key.substr(0, pos);
            }
            try { assacc2taxid[key] = stoi(value); } catch (...) { /* skip header/garbage */ }
        }
    }
    if (assacc2taxid.empty()) { cerr << "Error: no entries parsed from " << mappingFile << endl; return 1; }

    // Parse the 3-column list: group, community, classification path.
    struct Job { string group, community, path; };
    vector<Job> jobs;
    vector<string> groupOrder;
    set<string> groupSeen;
    {
        ifstream in(classificationList);
        if (!in.is_open()) { cerr << "Error: cannot open " << classificationList << endl; return 1; }
        string line;
        bool first = true;
        while (getline(in, line)) {
            if (line.empty()) continue;
            vector<string> f = splitTab(line);
            if (f.size() < 3) continue;
            if (first) { first = false; if (f[0] == "group" || f[0] == "Group") continue; } // optional header
            jobs.push_back(Job{f[0], f[1], f[2]});
            if (!groupSeen.count(f[0])) { groupSeen.insert(f[0]); groupOrder.push_back(f[0]); }
        }
    }
    if (jobs.empty()) { cerr << "Error: no rows in " << classificationList << endl; return 1; }

    // group -> rank -> per-community PRF
    unordered_map<string, unordered_map<string, vector<PRF>>> data;
    for (const Job &job : jobs) {
        unordered_map<string, CountAtRank> counts = gradeOneFile(job.path, assacc2taxid, taxonomy, ranks, par);
        for (const string &rank : ranks) {
            const CountAtRank &c = counts[rank];
            data[job.group][rank].push_back({clean(c.precision), clean(c.sensitivity), clean(c.f1)});
        }
    }

    printf("Group\tRank\tN\tPrecision_mean\tPrecision_sd\tSensitivity_mean\tSensitivity_sd\tF1_mean\tF1_sd\n");
    for (const string &group : groupOrder) {
        for (const string &rank : ranks) {
            const vector<PRF> &v = data[group][rank];
            if (v.empty()) continue;
            vector<double> pr, se, f1;
            for (const PRF &m : v) { pr.push_back(m.precision); se.push_back(m.sensitivity); f1.push_back(m.f1); }
            auto mp = meanSD(pr), ms = meanSD(se), mf = meanSD(f1);
            printf("%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                   group.c_str(), rank.c_str(), (int) v.size(),
                   mp.first, mp.second, ms.first, ms.second, mf.first, mf.second);
        }
    }
    return 0;
}
