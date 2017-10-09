/*
 *  DWARF - Genomic Analysis Software
 *
 *  (c) Dr Mark Daniel Preston, LSHTM, 2011
 *  www.markdpreston.com
 *
 */
#include "analysis.h"
#include "subject.h"
#include "maths.h"
#include "statistics.h"
#include "random.h"

void CAnalysisCAlpha::statistics(TRVector& rmStatistics, const bool pbShuffle) {
    unsigned int    i, j, liCount, liTotal;
    TReal           lrP, lrNumerator, lrDenominator;
    TIVector        lmDosageI, lmPhenotypes;
    TRVector        lmDosageR, lmCases, lmTotal, lmCounts;

    lrP = 1.0 * miCases / miSubjects;

    lmCases.resize(miSNPs+1);
    lmTotal.resize(miSNPs+1);
    lmCases = 0;
    lmTotal = 0;

    //  Perform a random shuffle of the phenotypes.
    if (pbShuffle) { CRandom::shuffle(mmPhenotypes); }

    //  Count the occurences in the cases and in total of each variant.
    liCount = miSubjects;
    for (i = 0; i < miSNPs; i++) {
        lmDosageR.resize(miSubjects);
        lmPhenotypes.resize(miSubjects);
        lmPhenotypes = mmPhenotypes;
        switch (moPopulation->getType()) {
            case cePopulationDosages:    lmDosageR = CDosage::convertToDosage(moPopulation->getDosages()->getSNP(i)); break;
            case cePopulationGenotypes:  lmDosageR = CDosage::convertToDosage(moPopulation->getGenotypes()->getSNP(i)); break;
            case cePopulationHaplotypes: lmDosageR = CDosage::convertToDosage(CHaplotype::transform(moPopulation->getHaplotypes()->getSNP(i))); break;
            default: msError << "Population type not implemented"; error();
        }
        if (mbCluster) {
            liCount = cluster(lmDosageR, lmPhenotypes);
        }
        lmDosageI.resize(liCount);
        lmDosageI = CDosage::convertToIntegerDosage(lmDosageR);
        lmCases(i) = sum((lmPhenotypes == 2) * lmDosageI);
        lmTotal(i) = sum(lmDosageI);
    }

    //  Handle singletons by pooling them, as they have no variance of their own.
    lmCases(miSNPs) = sum((lmTotal == 1) * lmCases);
    lmTotal(miSNPs) = sum(lmTotal == 1);

    //  Calculate the statistic, via the numerator and denominator separately.
    lrNumerator   = 0.0;
    lrDenominator = 0.0;
    liTotal = max(lmTotal);

    //  Zeros ignored and singletons pooled above.
    for (i = 2; i <= liTotal; i++) {
        liCount = sum(lmTotal == i);
        if (0 < liCount) {
            for (j = 0; j <= i; j++) {
                lrNumerator   += 1.0 * (pow(j - i * lrP,  (TReal) 2.0) - i * lrP * (1 - lrP)) * sum(lmCases == j && lmTotal == i);
                lrDenominator += 1.0 * liCount * pow(pow(j - i * lrP, (TReal) 2.0) - i * lrP * (1 - lrP), (TReal) 2.0) * CMaths::binomial(i,j) * pow(lrP,(TReal) j) * pow(1 - lrP,(TReal) i - j);
            }
        }
    }
    lrDenominator = sqrt(lrDenominator);

    if (1e-12 < lrDenominator) {
        rmStatistics(0) = -lrNumerator / lrDenominator;
    } else {
        rmStatistics(0) = 1;
    }
}

void CAnalysisCAlpha::asymptotics (TRVector& rmAsymptotics) {
    statistics(mmStatistics, false);
    rmAsymptotics(0) = CStatistics::pvalueNormal(mmStatistics(0));
}

void CAnalysisCAlpha::test() {
//    cout << setw(24) << msCommand << mmDataLine(0) << " " << (lrNumerator / lrDenominator) << endl;
}
