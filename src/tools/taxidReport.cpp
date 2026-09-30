// taxid-report: turn a plain list of NCBI taxonomy IDs into a Kraken/Metabuli-
// style hierarchical report. Each input line contributes one "read" at its
// taxon; counts are rolled up the taxonomy and emitted as an indented tree:
//
//   #clade_proportion  clade_count  taxon_count  rank  taxID  name
//
// Input: one taxid per line (the first whitespace-delimited token of the line is
// used; blank lines and lines starting with '#' are skipped; a taxid of 0 counts
// as "unclassified"). Positional args:
//   <taxidList> <taxonomyDir> [outReport]   (writes to stdout if outReport omitted)

#include "Parameters.h"
#include "TaxonomyWrapper.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

int taxidReport(const Parameters &par) {
    const string taxidListFile = par.filenames[0];
    const string taxonomyDir   = par.filenames[1];
    const string outFile       = (par.filenames.size() >= 3) ? par.filenames[2] : "";

    // Load taxonomy. Its progress messages go to stdout (Debug::INFO); redirect
    // stdout to stderr during construction so a piped report stays clean.
    const string names  = taxonomyDir + "/names.dmp";
    const string nodes  = taxonomyDir + "/nodes.dmp";
    const string merged = taxonomyDir + "/merged.dmp";
    std::streambuf * coutBuf = std::cout.rdbuf(std::cerr.rdbuf());
    TaxonomyWrapper taxonomy(names, nodes, merged, false);
    std::cout.rdbuf(coutBuf);

    // Count taxids. taxid 0 -> unclassified; taxids not in the taxonomy are
    // skipped (reported to stderr) and excluded from the total.
    ifstream in(taxidListFile);
    if (!in.is_open()) {
        cerr << "Cannot open taxid list: " << taxidListFile << endl;
        return 1;
    }
    unordered_map<TaxID, long> taxonCount; // key 0 = unclassified
    long total = 0;
    long unresolved = 0;
    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t s = line.find_first_not_of(" \t\r\n");
        if (s == string::npos) continue;
        size_t e = line.find_first_of(" \t\r\n", s);
        const string tok = line.substr(s, e == string::npos ? string::npos : e - s);
        TaxID taxid;
        try { taxid = (TaxID) stol(tok); } catch (...) { continue; }
        if (taxid < 0) continue;
        if (taxid == 0) { taxonCount[0]++; total++; continue; }
        if (!taxonomy.nodeExists(taxid)) { unresolved++; continue; }
        // Resolve to the canonical taxon (maps merged/old ids to their node).
        const TaxonNode * node = taxonomy.taxonNode(taxid, false);
        if (node == nullptr) { unresolved++; continue; }
        taxonCount[node->taxId]++;
        total++;
    }
    in.close();

    if (unresolved > 0) {
        cerr << "Warning: " << unresolved << " taxid(s) not found in the taxonomy were skipped" << endl;
    }

    // Roll classified taxon counts up the lineage into clade counts.
    unordered_map<TaxID, long> cladeCount;
    for (const auto & kv : taxonCount) {
        if (kv.first == 0) continue;
        const long c = kv.second;
        cladeCount[kv.first] += c;
        const TaxonNode * node = taxonomy.taxonNode(kv.first);
        while (node->parentTaxId != node->taxId && taxonomy.nodeExists(node->parentTaxId)) {
            node = taxonomy.taxonNode(node->parentTaxId);
            cladeCount[node->taxId] += c;
        }
    }

    // Output
    ofstream ofs;
    if (!outFile.empty()) {
        ofs.open(outFile);
        if (!ofs.is_open()) {
            cerr << "Cannot open output file: " << outFile << endl;
            return 1;
        }
    }
    ostream & out = outFile.empty() ? cout : ofs;

    out << "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname\n";
    out << std::fixed << std::setprecision(4);
    if (total == 0) {
        cerr << "No taxids counted." << endl;
        return 0;
    }

    // Unclassified line first (Kraken convention).
    auto uit = taxonCount.find(0);
    if (uit != taxonCount.end() && uit->second > 0) {
        out << 100.0 * (double) uit->second / (double) total << "\t" << uit->second << "\t"
            << uit->second << "\tno rank\t0\tunclassified\n";
    }

    // Depth-first from the root (taxID 1). Push children in ascending clade count
    // so the largest is popped (and printed) first.
    const unordered_map<TaxID, vector<TaxID>> parentToChildren = taxonomy.getParentToChildren();
    vector<pair<TaxID, int>> stack;
    stack.emplace_back(1, 0);
    while (!stack.empty()) {
        const TaxID tid = stack.back().first;
        const int depth = stack.back().second;
        stack.pop_back();
        auto ccIt = cladeCount.find(tid);
        if (ccIt == cladeCount.end() || ccIt->second <= 0) continue;
        const long cc = ccIt->second;
        long tc = 0;
        auto tcIt = taxonCount.find(tid);
        if (tcIt != taxonCount.end()) tc = tcIt->second;
        const TaxonNode * node = taxonomy.taxonNode(tid, false);
        const char * rank = node ? taxonomy.getString(node->rankIdx) : "no rank";
        const char * name = node ? taxonomy.getString(node->nameIdx) : "";
        out << 100.0 * (double) cc / (double) total << "\t" << cc << "\t" << tc << "\t"
            << rank << "\t" << tid << "\t" << string(2 * depth, ' ') << name << "\n";
        auto chIt = parentToChildren.find(tid);
        if (chIt != parentToChildren.end()) {
            vector<TaxID> kids;
            for (const TaxID k : chIt->second) {
                auto it = cladeCount.find(k);
                if (it != cladeCount.end() && it->second > 0) kids.push_back(k);
            }
            sort(kids.begin(), kids.end(), [&](TaxID a, TaxID b) {
                return cladeCount[a] < cladeCount[b];
            });
            for (const TaxID k : kids) stack.emplace_back(k, depth + 1);
        }
    }

    if (!outFile.empty()) ofs.close();
    return 0;
}
