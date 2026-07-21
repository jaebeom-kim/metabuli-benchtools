// Shared per-read grading pieces used by both `grade` and `grade-classification`.
// The comparison functions are defined in grade.cpp; both tools link them.
#ifndef BENCHTOOLS_GRADE_COMMON_H
#define BENCHTOOLS_GRADE_COMMON_H

#include "TaxonomyWrapper.h"
#include <string>

// Tally of read outcomes at one rank, plus the derived scores.
struct CountAtRank {
    int total;
    int FP;
    int TP;
    int FN;
    float precision;
    float sensitivity;
    float f1;
    void calculate() {
        precision = (float)TP / (float)(TP + FP);
        sensitivity = (float)TP / (float)(total);
        f1 = 2 * precision * sensitivity / (precision + sensitivity);
    }
};

// Score one read (shot = predicted taxid, target = true taxid) at `rank`,
// updating `count`. Returns 'O' (TP), 'X' (FP), or 'N' (FN). Defined in grade.cpp.
char compareTaxonAtRank_CAMI(TaxID shot, TaxID target, const TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                             const std::string & rank);
char compareTaxonAtRank_CAMI_euk(TaxID shot, TaxID target, TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                                 const std::string & rank);
char compareTaxon_overclassification(TaxID shot, TaxID target, TaxonomyWrapper & ncbiTaxonomy, CountAtRank & count,
                                     const std::string & rank);
char compareTaxon_hivExclusion(TaxID shot, TaxID target, CountAtRank & count);

#endif // BENCHTOOLS_GRADE_COMMON_H
