#include "TaxonomyWrapper.h"
#include "Parameters.h"
#include "Assembly.h"
#include "Debug.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <cstring>
#include <cstdio>
#include <cstdint>

using namespace std;

// One query genome in the benchmark, with the answer-key metadata needed to
// grade it. All query sets (exclusion + inclusion) are emitted as rows of this
// record into a single <assemblyList>.query.tsv manifest, so downstream tools
// filter by `category` instead of reading a dozen separate files.
struct QueryRecord {
    std::string accession;    // the query assembly accession
    std::string category;     // familyExclusion | genusExclusion | speciesExclusion |
                              // subspeciesExclusion | speciesInclusionPair |
                              // subspeciesInclusionPair
    std::string expectedRank; // deepest rank a correct call can reach given the DB
    TaxID       queryTaxid;   // taxid of the query assembly itself
    TaxID       subjectTaxid; // the taxon that defines the case (excluded/shared);
                              // paired rows share this value
    std::string subjectRank;  // rank of subjectTaxid
    int         subjectSize;  // number of direct members of subjectTaxid
};

// Append the subspecies/species inclusion query pairs to `queries`. Uses its own
// RNG streams (srand(4), mt19937(0)) that run independently of the exclusion
// stage. The caller passes a database-only assembly graph (species2assembly /
// genus2species built from the .database genomes), so every pair member is
// guaranteed to be in the database.
static int appendInclusionQueries(
    std::vector<QueryRecord> & queries,
    std::unordered_map<std::string, TaxID> & observedAcc2taxid,
    std::unordered_map<TaxID, std::vector<Assembly>> & species2assembly,
    std::unordered_map<TaxID, std::vector<TaxID>> & genus2species);

// Input-tree counts for the summary report (filled by each benchmark builder).
struct BenchmarkSummary {
    long totalGenomes = 0; // total input assemblies
    long nOrders = 0;
    long nFamilies = 0;
    long nGenera = 0;
    long nSpecies = 0;
};

// Write the three benchmark outputs: <prefix>.database (one accession per line,
// the genomes to build the classifier DB from), <prefix>.query.tsv (the unified
// query manifest), and <prefix>.summary (totals + per-rank exclusion/inclusion
// counts).
static void writeBenchmarkOutputs(
    const std::string & prefix,
    const BenchmarkSummary & summary,
    const std::vector<std::string> & databaseAssemblies,
    const std::vector<QueryRecord> & queries);

int split(const Parameters &par) {
    std::srand(par.randomSeed);
    string assemblyList  = par.filenames[0];
    string taxonomyPath  = par.filenames[1];
    string acc2taxidFile = par.filenames[2];
    string outputPrefix  = par.filenames[3];

    vector<QueryRecord> queries;

    // Load taxonomy
    TaxonomyWrapper taxonomy(taxonomyPath + "/names.dmp",
                             taxonomyPath + "/nodes.dmp",
                             taxonomyPath + "/merged.dmp",
                             false);

    // Load the assembly accession -> taxid mapping. Also index by version-stripped
    // accession so lookups work whether or not the assembly list carries a version
    // suffix.
    cout << "Loading accession to taxid mapping...";
    std::unordered_map<std::string, TaxID> acc2taxid;
    {
        FILE *handle = fopen(acc2taxidFile.c_str(), "r");
        if (handle == NULL) {
            cerr << "Error: could not open acc2taxid file " << acc2taxidFile << endl;
            return 1;
        }
        char buffer[4096];
        int fileTaxid;
        while (fscanf(handle, "%s\t%d", buffer, &fileTaxid) == 2) {
            acc2taxid[buffer] = (TaxID) fileTaxid;
        }
        fclose(handle);
        std::unordered_map<std::string, TaxID> noVersion;
        for (auto &it : acc2taxid) {
            size_t dot = it.first.find(".");
            if (dot != std::string::npos) {
                noVersion[it.first.substr(0, dot)] = it.second;
            }
        }
        for (auto &it : noVersion) {
            acc2taxid.emplace(it.first, it.second);
        }
    }
    cout << "done." << endl;


    cout << "Loading assembly list...";
    vector<string> totalAssemblyAccessions;
    vector<Assembly> assemblies;
    ifstream assemblyListFile(assemblyList);
    if (!assemblyListFile.is_open()) {
        cerr << "Error: could not open assembly list file " << assemblyList << endl;
        return 1;
    }

    string assemblyAccession;
    TaxID taxid;
    cout << "Making observedAcc2taxid map...";
    std::unordered_map<std::string, TaxID> observedAcc2taxid;
    while (getline(assemblyListFile, assemblyAccession)) {
        string assAccNoVersion = assemblyAccession.substr(0, assemblyAccession.find("."));
        
        // Check if a different version of the same assembly has already been observed
        if (observedAcc2taxid.find(assemblyAccession) != observedAcc2taxid.end()) {
            cout << "Warning: assembly " << assemblyAccession << " has already been observed" << endl;
        }
        
        // Get the taxonomy ID of the current assembly
        if (acc2taxid.find(assemblyAccession) != acc2taxid.end()) {
            taxid = acc2taxid[assemblyAccession];
        } else if (acc2taxid.find(assAccNoVersion) != acc2taxid.end()) {
            taxid = acc2taxid[assAccNoVersion];
        } else {
            cerr << "Error: accession " << assemblyAccession << " " << assAccNoVersion << " not found in --acc2taxid" << endl;
            return 1;
        }
        observedAcc2taxid[assemblyAccession] = taxid;

        // Record the assembly
        totalAssemblyAccessions.push_back(assemblyAccession);
        assemblies.emplace_back(assemblyAccession);
        assemblies.back().taxid = taxid;
        assemblies.back().speciesId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "species");
        assemblies.back().genusId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "genus");
        assemblies.back().familyId = taxonomy.getTaxIdAtRank(assemblies.back().taxid, "family");
    }
    cout << "done." << endl;


    unordered_map<TaxID, vector<Assembly>> species2assembly;
    for (auto &assembly : assemblies) {
        species2assembly[assembly.speciesId].push_back(assembly);
    }

    unordered_map<TaxID, vector<TaxID>> genus2species;
    for (auto &species : species2assembly) {
        TaxID genusId = taxonomy.getTaxIdAtRank(species.first, "genus");
        genus2species[genusId].push_back(species.first);
    }

    unordered_map<TaxID, vector<TaxID>> family2genus;
    for (auto &genus : genus2species) {
        TaxID familyId = taxonomy.getTaxIdAtRank(genus.first, "family");
        family2genus[familyId].push_back(genus.first);
    }

    cout << "Number of total families: " << family2genus.size() << endl;


    unordered_map<TaxID, vector<TaxID>> order2family;
    for (auto &family : family2genus) {
        TaxID orderId = taxonomy.getTaxIdAtRank(family.first, "order");
        order2family[orderId].push_back(family.first);
    }


    vector<string> totalExcludedAssemblies;
    vector<string> currentExcludedAssemblies;

    vector<TaxID> excludedFamilies;
    vector<TaxID> excludedGenera;
    vector<TaxID> excludedSpecies;

    // *** Family Exclusion ***
    // 1. Find orders with multiple families
    vector<TaxID> orderWithMultipleFamilies;
    for (auto &order : order2family) {
        if (order.second.size() > 1) {
            orderWithMultipleFamilies.push_back(order.first);
        }
    }
    // 2. Randomly choose n orders with multiple families without replacement
    vector<TaxID> selectedOrders;
    int n = orderWithMultipleFamilies.size() / 3;
    cout << "Excluding " << n << " families" << endl;
    for (int i = 0; i < n; i++) {
        int idx = rand() % orderWithMultipleFamilies.size();
        selectedOrders.push_back(orderWithMultipleFamilies[idx]);
        orderWithMultipleFamilies.erase(orderWithMultipleFamilies.begin() + idx);
    }
    // 3. Randomly choose one family from each selected order; its whole subtree
    //    is dropped from the DB and one member becomes the family-exclusion query.
    vector<string> excludedFamilyAssemblies;
    for (auto &order : selectedOrders) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % order2family[order].size();
        TaxID excludingFamily = order2family[order][idx];
        excludedFamilies.push_back(excludingFamily);
        // Exclude the selected family
        for (size_t i = 0; i < family2genus[excludingFamily].size(); i++) {
            TaxID genus = family2genus[excludingFamily][i];
            excludedGenera.push_back(genus);
            for (size_t j = 0; j < genus2species[genus].size(); j++) {
                TaxID species = genus2species[genus][j];
                excludedSpecies.push_back(species);
                for (size_t k = 0; k < species2assembly[species].size(); k++) {
                    const string & assemblyName = species2assembly[species][k].name;
                    totalExcludedAssemblies.push_back(assemblyName);
                    excludedFamilyAssemblies.push_back(assemblyName);
                    currentExcludedAssemblies.push_back(assemblyName);
                 }
            }
        }
        // Report every excluded member genome as a family-exclusion query candidate.
        for (const string & member : currentExcludedAssemblies) {
            queries.push_back({member, "familyExclusion", "order", observedAcc2taxid[member],
                               excludingFamily, "family", (int) family2genus[excludingFamily].size()});
        }
    }


    // *** Genus Exclusion ***
    // 1. Find families with multiple genera
    vector<TaxID> familyWithMultipleGenera;

    size_t genusNumBeforeExclusion = 0;
    for (auto &family : family2genus) {
        if (find(excludedFamilies.begin(), excludedFamilies.end(), family.first) != excludedFamilies.end()) {
            continue;
        }
        genusNumBeforeExclusion += family.second.size();
        if (family.second.size() > 1) {
            familyWithMultipleGenera.push_back(family.first);
        }
    }
    cout << "Number of genera before exclusion: " << genusNumBeforeExclusion << endl;
    // 2. Randomly choose n families with multiple genera without replacement
    n = familyWithMultipleGenera.size() / 3;
    cout << "Excluding " << n << " genera" << endl;
    vector<TaxID> selectedFamilies;
    for (int i = 0; i < n; i++) {
        int idx = rand() % familyWithMultipleGenera.size();
        selectedFamilies.push_back(familyWithMultipleGenera[idx]);
        familyWithMultipleGenera.erase(familyWithMultipleGenera.begin() + idx);
    }
    // 3. Randomly choose one genus from each selected family; drop it from the DB
    //    and take one member as the genus-exclusion query.
    vector<string> excludedGenusAssemblies;
    for (auto &family : selectedFamilies) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % family2genus[family].size();
        TaxID excludingGenus = family2genus[family][idx];
        excludedGenera.push_back(excludingGenus);
        for (size_t i = 0; i < genus2species[excludingGenus].size(); i++) {
            excludedSpecies.push_back(genus2species[excludingGenus][i]);
            for (size_t j = 0; j < species2assembly[genus2species[excludingGenus][i]].size(); j++) {
                totalExcludedAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j].name);
                excludedGenusAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j].name);
                currentExcludedAssemblies.push_back(species2assembly[genus2species[excludingGenus][i]][j].name);
            }
        }
        // Report every excluded member genome as a genus-exclusion query candidate.
        for (const string & member : currentExcludedAssemblies) {
            queries.push_back({member, "genusExclusion", "family", observedAcc2taxid[member],
                               excludingGenus, "genus", (int) genus2species[excludingGenus].size()});
        }
    }


    // *** Species Exclusion ***
    // 1. Find genera with multiple species
    vector<TaxID> genusWithMultipleSpecies;
    for (auto &genus : genus2species) {
        if (genus.second.size() > 1) {
            if (find(excludedGenera.begin(), excludedGenera.end(), genus.first) != excludedGenera.end()) {
                continue;
            }
            genusWithMultipleSpecies.push_back(genus.first);
        }
    }

    cout << "Number of genera with multiple species: " << genusWithMultipleSpecies.size() << endl;
    // Select n genera with multiple species
    n = int(genusWithMultipleSpecies.size() / 3);
    cout << "Excluding " << n << " species" << endl;
    vector<TaxID> selectedGenera;
    for (int i = 0; i < n; i++) {
        int idx = rand() % genusWithMultipleSpecies.size();
        selectedGenera.push_back(genusWithMultipleSpecies[idx]);
        genusWithMultipleSpecies.erase(genusWithMultipleSpecies.begin() + idx);
    }

    // For each selected genus, randomly choose one species to exclude; drop it
    // from the DB and take one member as the species-exclusion query.
    vector<string> excludedSpeciesAssemblies;
    for (auto &genus : selectedGenera) {
        currentExcludedAssemblies.clear();
        int random = rand();
        int idx = random % genus2species[genus].size();
        TaxID excludingSpecies = genus2species[genus][idx];
        excludedSpecies.push_back(excludingSpecies);
        for (size_t i = 0; i < species2assembly[excludingSpecies].size(); i++) {
            currentExcludedAssemblies.push_back(species2assembly[excludingSpecies][i].name);
            totalExcludedAssemblies.push_back(species2assembly[excludingSpecies][i].name);
            excludedSpeciesAssemblies.push_back(species2assembly[excludingSpecies][i].name);
        }
        // Report every excluded member genome as a species-exclusion query candidate.
        for (const string & member : currentExcludedAssemblies) {
            queries.push_back({member, "speciesExclusion", "genus", observedAcc2taxid[member],
                               excludingSpecies, "species", (int) species2assembly[excludingSpecies].size()});
        }
    }

    for (auto &excludedGenus : excludedGenera) {
        for (auto & species : genus2species[excludedGenus]) {
            excludedSpecies.push_back(species);
        }
    }

    // *** Subspecies Exclusion and Inclusion ***

    // 1. Find species with multiple assemblies
    vector<TaxID> speciesWithMultipleAssemblies;
    for (auto &species : species2assembly) {
        if (species.second.size() > 1) {
            if (find(excludedSpecies.begin(), excludedSpecies.end(), species.first) != excludedSpecies.end()) {
                continue;
            }
            speciesWithMultipleAssemblies.push_back(species.first);
        }
    }
  
    // 2. Select n species with multiple assemblies
    n = int(speciesWithMultipleAssemblies.size()/3);
    cout << "Excluding " << n << " assemblies" << endl;
    vector<TaxID> selectedSpecies;
    for (int i = 0; i < n; i++) {
        int idx = rand() % speciesWithMultipleAssemblies.size();
        selectedSpecies.push_back(speciesWithMultipleAssemblies[idx]);
        speciesWithMultipleAssemblies.erase(speciesWithMultipleAssemblies.begin() + idx);
    }

    // 3. For each selected species, randomly choose one assembly to exclude; its
    //    species stays in the DB via siblings, so it is the subspecies-exclusion query.
    vector<string> excludedSubspeciesAssemblies;
    for (auto &species : selectedSpecies) {
        int idx = rand() % species2assembly[species].size();
        const string & query = species2assembly[species][idx].name;
        totalExcludedAssemblies.push_back(query);
        excludedSubspeciesAssemblies.push_back(query);
        queries.push_back({query, "subspeciesExclusion", "species", observedAcc2taxid[query],
                           species, "species", (int) species2assembly[species].size()});
    }

    cout << "Total excluded assemblies: " << totalExcludedAssemblies.size() << endl;

    // Compute the database assembly list (everything not excluded).
    vector<string> databaseAssemblies;
    for (auto &assembly : totalAssemblyAccessions) {
        if (find(totalExcludedAssemblies.begin(), totalExcludedAssemblies.end(), assembly) == totalExcludedAssemblies.end()) {
            databaseAssemblies.push_back(assembly);
        }
    }



    // Validate the database assembly list
    if (par.skipValidation) {
        cout << "Skipping validation (--skip-validation)." << endl;
    } else {

    // Validate Family Exclusion
    cout << "Validating excluded family..." << endl;
    for (size_t i = 0; i < excludedFamilyAssemblies.size(); i++) {
        TaxID excludedTaxid = observedAcc2taxid[excludedFamilyAssemblies[i]];
        int orderCount = 0;
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(excludedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);
            if (rank == "order") {
                orderCount++;
                continue;
            } 
            
            if (strcmp(taxonomy.getString(lcaTaxon->nameIdx), "root") == 0) {
                continue;
            }

            if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("order")) {    
                cout << "Error: " << excludedFamilyAssemblies[i] << " is not a valid family exclusion. LCA is below order rank." << endl;
                cout << "Name: " << taxonomy.getString(lcaTaxon->nameIdx) << " TaxID: " << lcaTaxon->taxId << endl;
                cout << "Rank: " << rank << " RankIdx: " << lcaTaxon->rankIdx << endl;
                cout << excludedFamilyAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
                return 1;
            }
        }
        if (orderCount == 0) {
            cout << "Error: " << excludedFamilyAssemblies[i] << " is not a valid exclusion. No order rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of excluded families complete." << endl;

    // Validate genus exclusions
    cout << "Validating excluded genera..." << endl;
    for (size_t i = 0; i < excludedGenusAssemblies.size(); i++) {
        TaxID excludedTaxid = observedAcc2taxid[excludedGenusAssemblies[i]];
        // There must be at least one Family rank LCA
        // There must not be any LCA below Family rank
        int familyCount = 0;
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(excludedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);
            if (rank == "family") {
                familyCount++;
                continue;
            } 

            if (strcmp(taxonomy.getString(lcaTaxon->nameIdx), "root") == 0) {
                continue;
            }

            if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("family")) {    
                cout << "Error: " << excludedFamilyAssemblies[i] << " is not a valid genus exclusion. LCA is below family rank." << endl;
                cout << "Name: " << taxonomy.getString(lcaTaxon->nameIdx) << " TaxID: " << lcaTaxon->taxId << endl;
                cout << "Rank: " << rank << " RankIdx: " << lcaTaxon->rankIdx << endl;
                cout << excludedFamilyAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
                return 1;
            }            
            // else if (taxonomy.findRankIndex(rank) == -1) {
            //     cout << "Error: LCA rank index is -1." << endl;
            //     return 1;
            // } else if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("family")) {
            //     cout << "Error: " << excludedGenusAssemblies[i] << " is not a valid exclusion. LCA is below Family rank." << endl;
            //     cout << excludedGenusAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
            //     return 1;
            // }
        }
        if (familyCount == 0) {
            cout << "Error: " << excludedGenusAssemblies[i] << " is not a valid exclusion. No Family rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of excluded genera complete." << endl;

    // Validate species exclusions
    cout << "Validating excluded species..." << endl;
    for (size_t i = 0; i < excludedSpeciesAssemblies.size(); i++) {
        TaxID excludedTaxid = observedAcc2taxid[excludedSpeciesAssemblies[i]];
        // There must be at least one Genus rank LCA
        // There must not be any LCA below Genus rank
        int genusCount = 0;
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(excludedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);

            if (rank == "genus") {
                genusCount++;
                continue;
            } 

            if (strcmp(taxonomy.getString(lcaTaxon->nameIdx), "root") == 0) {
                continue;
            }

            if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("genus")) {    
                cout << "Error: " << excludedFamilyAssemblies[i] << " is not a valid species exclusion. LCA is below genus rank." << endl;
                cout << "Name: " << taxonomy.getString(lcaTaxon->nameIdx) << " TaxID: " << lcaTaxon->taxId << endl;
                cout << "Rank: " << rank << " RankIdx: " << lcaTaxon->rankIdx << endl;
                cout << excludedFamilyAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
                return 1;
            }            

            // if (rank == "genus") {
            //     genusCount++;
            // } else if (taxonomy.findRankIndex(rank) == -1) {
            //     cout << "Error: LCA rank index is -1." << endl;
            //     return 1;
            // } else if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("genus")) {
            //     cout << "Error: " << excludedSpeciesAssemblies[i] << " is not a valid exclusion. LCA is below Genus rank." << endl;
            //     cout << excludedSpeciesAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
            //     return 1;
            // }
        }
        if (genusCount == 0) {
            cout << "Error: " << excludedSpeciesAssemblies[i] << " is not a valid exclusion. No Genus rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of excluded species complete." << endl;

    // Validate subspecies exclusions
    cout << "Validating excluded subspecies..." << endl;
    for (size_t i = 0; i < excludedSubspeciesAssemblies.size(); i++) {
        TaxID excludedTaxid = observedAcc2taxid[excludedSubspeciesAssemblies[i]];
        // There must be at least one Species rank LCA
        // There must not be any LCA below Species rank
        int speciesCount = 0;
        for (size_t j = 0; j < databaseAssemblies.size(); j++) {
            TaxID databaseTaxid = observedAcc2taxid[databaseAssemblies[j]];
            const TaxonNode * lcaTaxon = taxonomy.taxonNode(taxonomy.LCA(excludedTaxid, databaseTaxid));
            const string & rank = taxonomy.getString(lcaTaxon->rankIdx);
            if (rank == "species") {
                speciesCount++;
                continue;
            } 

            if (strcmp(taxonomy.getString(lcaTaxon->nameIdx), "root") == 0) {
                continue;
            }

            if (taxonomy.findRankIndex(rank) < taxonomy.findRankIndex("species")) {    
                cout << "Error: " << excludedFamilyAssemblies[i] << " is not a valid exclusion. LCA is below species rank." << endl;
                cout << "Name: " << taxonomy.getString(lcaTaxon->nameIdx) << " TaxID: " << lcaTaxon->taxId << endl;
                cout << "Rank: " << rank << " RankIdx: " << lcaTaxon->rankIdx << endl;
                cout << excludedFamilyAssemblies[i] << " " << databaseAssemblies[j] << " " << rank << endl;
                return 1;
            } 
        }
        if (speciesCount == 0) {
            cout << "Error: " << excludedSubspeciesAssemblies[i] << " is not a valid exclusion. No Species rank LCA." << endl;
            return 1;
        }
    }
    cout << "Validation of excluded subspecies complete." << endl;
    } // end validation (--skip-validation)

    // Inclusion query pairs (subspecies + species). These must be genomes that
    // are actually in the database, so build a database-only assembly graph
    // (excluding every held-out genome) and draw the pairs from it.
    std::unordered_set<std::string> databaseSet(databaseAssemblies.begin(), databaseAssemblies.end());
    unordered_map<TaxID, vector<Assembly>> dbSpecies2assembly;
    for (auto &assembly : assemblies) {
        if (databaseSet.count(assembly.name)) {
            dbSpecies2assembly[assembly.speciesId].push_back(assembly);
        }
    }
    unordered_map<TaxID, vector<TaxID>> dbGenus2species;
    for (auto &sp : dbSpecies2assembly) {
        TaxID genusId = taxonomy.getTaxIdAtRank(sp.first, "genus");
        dbGenus2species[genusId].push_back(sp.first);
    }
    if (int rc = appendInclusionQueries(queries, observedAcc2taxid, dbSpecies2assembly, dbGenus2species)) {
        return rc;
    }

    BenchmarkSummary summary;
    summary.totalGenomes = (long) totalAssemblyAccessions.size();
    summary.nOrders = (long) order2family.size();
    summary.nFamilies = (long) family2genus.size();
    summary.nGenera = (long) genus2species.size();
    summary.nSpecies = (long) species2assembly.size();

    writeBenchmarkOutputs(outputPrefix, summary, databaseAssemblies, queries);
    return 0;
}

static int appendInclusionQueries(
    std::vector<QueryRecord> & queries,
    std::unordered_map<std::string, TaxID> & observedAcc2taxid,
    std::unordered_map<TaxID, std::vector<Assembly>> & species2assembly,
    std::unordered_map<TaxID, std::vector<TaxID>> & genus2species)
{
    std::srand(4);

    /* Subspecies inclusion test: two subspecies (assemblies) of the same species. */
    cout << "Making query sets for subspecies inclusion test..." << endl;

    // Find species with multiple assemblies(subspecies)
    vector<TaxID> speciesWithMultipleAssemblies;
    for (auto &species : species2assembly) {
        if (species.second.size() > 1) {
            speciesWithMultipleAssemblies.push_back(species.first);
        }
    }
    cout << "Found " << speciesWithMultipleAssemblies.size() << " species with multiple assemblies. A random eigth will be used." << endl;

    std::mt19937 rng1(0);
    std::shuffle(speciesWithMultipleAssemblies.begin(), speciesWithMultipleAssemblies.end(), rng1);
    std::size_t eigth = speciesWithMultipleAssemblies.size() / 8; // floor if odd
    std::vector<TaxID> selectedSpecies(speciesWithMultipleAssemblies.begin(),
                                  speciesWithMultipleAssemblies.begin() + eigth);

    // Select two assemblies per species with multiple assemblies
    for (auto &species : selectedSpecies) {
        if (species2assembly[species].size() < 2) {
            cerr << "Error: species " << species << " has less than 2 assemblies." << endl;
            return 1;
        }
        int idx1 = rand() % species2assembly[species].size();
        int idx2 = rand() % species2assembly[species].size();
        while (idx2 == idx1) {
            idx2 = rand() % species2assembly[species].size();
        }
        int size = (int) species2assembly[species].size();
        for (int idx : {idx1, idx2}) {
            const string & query = species2assembly[species][idx].name;
            queries.push_back({query, "subspeciesInclusionPair", "species",
                               observedAcc2taxid[query], species, "species", size});
        }
    }

    /* Species inclusion test: two species of the same genus. */
    cout << "Making query sets for species inclusion test..." << endl;

    // Find genera with multiple species
    vector<TaxID> genusWithMultipleSpecies;
    for (auto &genus : genus2species) {
        if (genus.second.size() > 1) {
            genusWithMultipleSpecies.push_back(genus.first);
        }
    }
    cout << "Found " << genusWithMultipleSpecies.size() << " genera with multiple species. A random quater will be used." << endl;

    std::mt19937 rng(0);
    std::shuffle(genusWithMultipleSpecies.begin(), genusWithMultipleSpecies.end(), rng);
    std::size_t quater2 = genusWithMultipleSpecies.size() / 4; // floor if odd
    std::vector<TaxID> selectedGenera(genusWithMultipleSpecies.begin(),
                                  genusWithMultipleSpecies.begin() + quater2);

    // Select two species per genus with multiple species
    for (auto &genus : selectedGenera) {
        if (genus2species[genus].size() < 2) {
            cerr << "Error: genus " << genus << " has less than 2 species." << endl;
            return 1;
        }
        int idx1 = rand() % genus2species[genus].size();
        int idx2 = rand() % genus2species[genus].size();
        while (idx2 == idx1) {
            idx2 = rand() % genus2species[genus].size();
        }
        int species1 = genus2species[genus][idx1];
        int species2 = genus2species[genus][idx2];
        if (species1 == species2) {
            cerr << "Error: selected species are the same for genus " << genus << endl;
            return 1;
        }
        // Choose one assembly from each species
        idx1 = rand() % species2assembly[species1].size();
        idx2 = rand() % species2assembly[species2].size();
        int size = (int) genus2species[genus].size();
        const string & q1 = species2assembly[species1][idx1].name;
        const string & q2 = species2assembly[species2][idx2].name;
        // ExpectedRank is species: the query genome's own species is in the
        // database. The shared genus (SubjectTaxID/SubjectRank) is the same-genus
        // distractor the test adds, not the expected classification rank.
        queries.push_back({q1, "speciesInclusionPair", "species", observedAcc2taxid[q1], genus, "genus", size});
        queries.push_back({q2, "speciesInclusionPair", "species", observedAcc2taxid[q2], genus, "genus", size});
    }

    cout << "Query sets for inclusion tests created." << endl;
    return 0;
}

static void writeBenchmarkOutputs(
    const std::string & prefix,
    const BenchmarkSummary & summary,
    const std::vector<std::string> & databaseAssemblies,
    const std::vector<QueryRecord> & queries)
{
    string databaseFile = prefix + ".database";
    string queryFile = prefix + ".query.tsv";
    string summaryFile = prefix + ".summary";

    ofstream db(databaseFile);
    for (auto &assembly : databaseAssemblies) {
        db << assembly << "\n";
    }
    db.close();

    ofstream q(queryFile);
    q << "Accession\tCategory\tExpectedRank\tQueryTaxID\tSubjectTaxID\tSubjectRank\tSubjectSize\n";
    for (auto &r : queries) {
        q << r.accession << "\t" << r.category << "\t" << r.expectedRank << "\t"
          << r.queryTaxid << "\t" << r.subjectTaxid << "\t" << r.subjectRank << "\t"
          << r.subjectSize << "\n";
    }
    q.close();

    // Aggregate the manifest per category: query-genome rows and distinct
    // subjects (excluded taxa, or pairs for the inclusion categories).
    unordered_map<string, long> rows;
    unordered_map<string, unordered_set<TaxID>> subjects;
    for (auto &r : queries) {
        rows[r.category]++;
        subjects[r.category].insert(r.subjectTaxid);
    }
    auto nRows = [&](const string &c) { return rows.count(c) ? rows[c] : 0L; };
    auto nSubj = [&](const string &c) { return subjects.count(c) ? (long) subjects[c].size() : 0L; };

    long databaseGenomes = (long) databaseAssemblies.size();
    long excludedGenomes = summary.totalGenomes - databaseGenomes;

    ofstream s(summaryFile);
    s << "# totals\n";
    s << "total_genomes\t" << summary.totalGenomes << "\n";
    s << "database_genomes\t" << databaseGenomes << "\n";
    s << "excluded_genomes\t" << excludedGenomes << "\n";
    s << "n_orders\t" << summary.nOrders << "\n";
    s << "n_families\t" << summary.nFamilies << "\n";
    s << "n_genera\t" << summary.nGenera << "\n";
    s << "n_species\t" << summary.nSpecies << "\n";
    s << "\n# per-rank exclusion / inclusion\n";
    s << "Rank\tExcludedTaxa\tExclQueryGenomes\tInclusionPairs\tInclQueryGenomes\n";
    s << "family\t"     << nSubj("familyExclusion")     << "\t" << nRows("familyExclusion")     << "\t-\t-\n";
    s << "genus\t"      << nSubj("genusExclusion")      << "\t" << nRows("genusExclusion")      << "\t-\t-\n";
    s << "species\t"    << nSubj("speciesExclusion")    << "\t" << nRows("speciesExclusion")    << "\t"
      << nSubj("speciesInclusionPair")    << "\t" << nRows("speciesInclusionPair")    << "\n";
    s << "subspecies\t" << nSubj("subspeciesExclusion") << "\t" << nRows("subspeciesExclusion") << "\t"
      << nSubj("subspeciesInclusionPair") << "\t" << nRows("subspeciesInclusionPair") << "\n";
    s.close();

    cout << "Wrote " << databaseGenomes << " database genome(s) to " << databaseFile << endl;
    cout << "Wrote " << queries.size() << " query genome(s) to " << queryFile << endl;
    cout << "Wrote summary to " << summaryFile << endl;
}
