/*
 *  DWARF - genomic analysis software
 *
 *  software.markdpreston.com/dwarf
 *
 *  (c) Mark Daniel Preston 2011-
 *
 *  When using this software, commercially or academically, please
 *  get in contact and reference:
 *
 *  XXXX
 *
 *
 *  This software was written, in large part, at the London School of
 *  Hygiene and Tropical Medicine, UK for the XXX project funded by
 *  YYYY.
 *
 *
 *  This file is part of DWARF.
 *
 *  DWARF is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DWARF is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DWARF.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "analysis.h"
#include "process.h"
#include "random.h"
#include "maths.h"
#include "statistics.h"
#include "utility.h"
#include <exception>
#include "contrib_kbac.h"

void CAnalysisKBAC::run (const vector<string>& paParameters) {
    unsigned int  i, j, liCount = 0;
    int           lbQuiet, liPermutations, liSNPs, liSubjects, liSided;
    double        lrAlpha, lrMAFUpper, lrP, *lmPhenotypes, *lmMAF, *lmData;
    KbacTest     *loKBAC;

    if (! preRun(paParameters)) { return; }

//cout << miPermutations << endl;

    //  Get the phenotype/SNP data.
    liCount = 0;
    lmData       = new double[miSNPs * miSubjects];
    lmMAF        = new double[miSNPs];
    for (i = 0; i < miSubjects; i++) {
        for (j = 0; j < miSNPs; j++) {
            switch (moPopulation->getType()) {
                case cePopulationGenotypes:  lmData[liCount] = (double) moPopulation->getGenotypes()->getDosage(i,j,false); break;
                case cePopulationHaplotypes: lmData[liCount] = (double) moPopulation->getHaplotypes()->getDosage(i,j); break;
                default: msError << "Population type not implemented"; error();
            }
            lmMAF[j] += lmData[liCount] / miSubjects / 2;
            liCount++;
        }
    }
    lmPhenotypes = new double[miSubjects];
    for (i = 0; i < miSubjects; i++) {
        lmPhenotypes[i] = moPopulation->getPhenotype(i) - 1.0;
    }

//for (i = 0; i < 10; i++) { cout << lmPhenotypes[i] << " "; } cout << endl;
//for (i = 0; i < 10; i++) { for (j = 0; j < 10; j++) { cout << lmData[i*miSNPs+j] << " "; } cout << endl; } cout << endl;

    //  Run the KBAC code.
    liPermutations  = miPermutations;
    lbQuiet         = true;
    lrAlpha         = 2.0;
    lrMAFUpper      = 1.0;
    liSNPs          = miSNPs;
    liSubjects      = miSubjects;
    lrP             = 9.0;  //  Must be > 1.0 otherwise contrib_kbac thinks that the first iteration has a p-value!
    liSided         = 1;
    loKBAC = new KbacTest(&liPermutations, &lbQuiet, &lrAlpha, &lrMAFUpper, lmData, lmPhenotypes, lmMAF, &liSNPs, &liSubjects);
    loKBAC->calcKbacP(&lrP, &liSided);

//cout << lrP << endl;

    //  Expose the result.
    mmDataLine(0) = lrP;

    postRun();
}

void CAnalysisKBAC::test () {
#ifdef INCLUDE_R
    TReal       lrP;

    lrP = getKBAC();

    cout << setw(24) << "Test KBAC";
    if (abs(mmDataLine(0) - lrP) < 0.005) {
        cout << "Pass" << endl;
//        cerr << left;
//        cerr << setw(24) << "Calculated" << mmDataLine(0) << endl;
//        cerr << setw(24) << "Reference" << lrP << endl;
    } else {
        cout << "Fail" << endl;
        cerr << left;
        cerr << setw(24) << "     Calculated" << mmDataLine(0) << endl;
        cerr << setw(24) << "     Reference" << lrP << endl;
    }

#endif
}

TReal CAnalysisKBAC::getKBAC() {
    if (0 == soProcess->miR.rinside || 0 == soProcess->miR.rcpp || 0 == soProcess->miR.kbac) {
        msWarning << "R libraries required for testing: RInside, Rcpp, KBAC";
        warning();
        return -1;
    }

    TRMatrix    lmKBACData;
//    string      lsPermutations;
    TReal       rrP;

    //  Get the phenotype/SNP data.
    lmKBACData.resize(miSubjects,miSNPs+1);
    lmKBACData(CUtility::getAll(),0) = moPopulation->getPhenotypes() - 1;
    switch (moPopulation->getType()) {
        case cePopulationGenotypes:  lmKBACData(CUtility::getAll(),blitz::Range(1,miSNPs)) = moPopulation->getGenotypes()->getDosage(false); break;
        case cePopulationHaplotypes: lmKBACData(CUtility::getAll(),blitz::Range(1,miSNPs)) = moPopulation->getHaplotypes()->getDosage(); break;
        default: msError << "Population type not implemented"; error();
    }

    //  Call out to R to analyse.
    CInterfaceR::setMatrix(string("x"), lmKBACData);
//    lsPermutations = boost::lexical_cast<std::string>(miPermutations);
    rrP = CInterfaceR::runCommandR("KbacTest(x,alpha=2,num.permutation=10000);");
    CInterfaceR::runCommand("remove(x);");
    return rrP;
}

/*
TReal CAnalysisKBAC::kbacMonteCarlo(const vector<TCount>& mmCounts, const TReal prT) {
    unsigned int    i, j, liCounts, liCases;
    blitz::Array<TBinomialD*,1> lmDistributions;
    blitz::Array<TBinomialG*,1> lmGenerators;
    TReal           lrProbability;
    vector<TCount>  lmCounts;
    TReal           rrP;

    liCounts = mmCounts.size();
    lmCounts = mmCounts;

    lrProbability = 1.0 * miCases / miSubjects;
    lmDistributions.resize(liCounts);
    lmGenerators.resize(liCounts);
    for (i = 0; i < liCounts; i++) {
        lmDistributions(i) = new TBinomialD(mmCounts[i].miTotal, lrProbability);
        lmGenerators(i)    = new TBinomialG(CRandom::soBoostRandom, *lmDistributions(i));
    }

    rrP = 0.0;
    for (i = 0; i < miPermutations; i++) {
        for (j = 1; j < liCounts; j++) {
            liCases = (*lmGenerators(j))();
            lmCounts[j].miControls  = lmCounts[j].miTotal - liCases;
            lmCounts[0].miControls += lmCounts[j].miTotal - liCases;
            lmCounts[j].miCases     = liCases;
            lmCounts[0].miCases    += liCases;
        }
        rrP += (kbacStatistic(lmCounts) >= prT ? 1.0 : 0.0);
    }
    rrP = (rrP + 1.0) / (miPermutations + 1);

    for (i = 0; i < liCounts; i++) {
        delete lmDistributions(i);
        delete lmGenerators(i);
    }

    return rrP;
}
*/
