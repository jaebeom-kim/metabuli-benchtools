// Extracted from Metabuli's src/commons/common.h (only the Assembly struct is
// used by the benchmark tools).
#ifndef BENCHTOOLS_ASSEMBLY_H
#define BENCHTOOLS_ASSEMBLY_H

#include <string>
#include <iostream>
#include "NcbiTaxonomy.h" // for TaxID

struct Assembly {
    std::string name;
    TaxID taxid;
    TaxID speciesId;
    TaxID genusId;
    TaxID familyId;
    TaxID orderId;

    Assembly(std::string name) : name(name), taxid(0), speciesId(0), genusId(0), familyId(0), orderId(0) {}
    Assembly() : name(""), taxid(0), speciesId(0), genusId(0), familyId(0), orderId(0) {}

    void print() const {
        std::cout << "Assembly: " << name << ", TaxID: " << taxid
                  << ", SpeciesID: " << speciesId
                  << ", GenusID: " << genusId
                  << ", FamilyID: " << familyId
                  << ", OrderID: " << orderId << std::endl;
    }
};

#endif // BENCHTOOLS_ASSEMBLY_H
