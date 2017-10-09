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
#include "statistics.h"

CAnalysis::CAnalysis () {
    clean();
}

CAnalysis::~CAnalysis () {
    clean();
}

void CAnalysis::clean () {
    miPermutations = 0;
    miDataEntry    = 0;
    miDataSize     = 0;
    miDataFields   = -1;
    mrAlpha        = 0.05;
    mbCluster      = false;
    mmClusters.resize(0);
    mmPermutationLine.resize(0);
    mmDataLine.resize(0);
    mmPermutations.resize(0,0);
    mmData         = NULL;
}

bool CAnalysis::preRun(const vector<string>& paParameters) {
    string  lsPermutations;

    //  Common pre-run actions.
    CCommand::preRun(paParameters);

    //  Check, extract and use parameters
    if (2 > paParameters.size()) {
        help();
        return false;
    }

    //  Get command line parameters.
    miDataSize     = CUtility::stringToInt(paParameters[0], 1);
    miDataEntry    = CUtility::stringToInt(paParameters[1], 1);
    miPermutations = CUtility::getParameterInteger(paParameters,"permutations",0);
    msData         = CUtility::getParameterString(paParameters,"data","");
    mbCluster      = CUtility::getParameterString(paParameters, "cluster", "false") != "false";
    mmClusters.resize(miSubjects);

    //  Initialise the data.
    lsPermutations = miPermutations > 0 ? "P" : "";
    mmData         = soProcess->getData(msCommand + msData + lsPermutations);
    mrAlpha        = CUtility::getParameterReal(paParameters,"alpha",0.05);
    miDataFields   = miStatistics * (mbSNPs ? miSNPs : 1);
    if (1 == miDataEntry) {
        mmStatistics.resize(miDataFields);
        mmPermutationLine.resize(miDataFields);
        mmDataLine.resize(miDataFields);
        mmPermutations.resize(miPermutations,miDataFields);
        mmData->resize(miDataSize,miDataFields);
        (*mmData) = 0;
    }

    mmDataLine = 1;
    mmStatistics = 1;
    mmPermutations = 0;

    return true;
}

void CAnalysis::run(const vector<string>& paParameters) {
    unsigned int    i;

    if (! preRun(paParameters)) { return; }

    if (miPermutations > 0) {
        statistics(mmStatistics, false);
        for (i = 0; i < miPermutations; i++) {
            statistics(mmPermutationLine, true);
            mmPermutations((int) i,CUtility::getAll()) = mmPermutationLine;
        }
        for (i = 0; i < miStatistics; i++) {
            mmDataLine(i) = (TReal) blitz::sum(mmPermutations(CUtility::getAll(),i) <= mmStatistics(i)) / miPermutations;
        }
    } else {
        asymptotics(mmDataLine);
    }

    postRun();
}

//  Store and output final (mean) values
bool CAnalysis::postRun (const bool pbSuppress) {
    bool        rbLast = false;
    int         i, j;
    TRVector    lmMean, lmVariance;

    CCommand::postRun();

    if (-1 != miDataFields) {
        (*mmData)(miDataEntry-1,CUtility::getAll()) = mmDataLine;
        rbLast = (miDataEntry == miDataSize);
        if (rbLast && ! pbSuppress) {
            (*out) << msCommand << " " << msData << endl << left;
            lmMean.resize(mmData->cols());
            lmVariance.resize(mmData->cols());
            lmMean = CStatistics::mean(*mmData);
            lmVariance = CStatistics::variance(*mmData);
            if (mbSNPs) {
                (*out) << setw(10) << "SNP" << setw(10) << "Mean P" << setw(10) << "SD" << "P<" << mrAlpha << endl;
                for (i = 0; i < (int)miSNPs; i++) {
                    for (j = 0; j < (int)miDataFields / (int)miSNPs; j++) {
                        (*out) << setw(10) << moPopulation->getSNPInfo(i).getName();
                        (*out) << setw(10) << lmMean(i*miDataFields/miSNPs+j);
                        (*out) << setw(10) << sqrt(lmVariance(i*miDataFields/miSNPs+j));
                        (*out) << 1.0 * sum(abs((*mmData)(CUtility::getAll(),i)) < mrAlpha) / miDataSize << endl;
                    }
                }
            } else {
                (*out) << setw(10) << "Mean "; for (i = 0; i < lmMean.rows(); i++) { (*out) << setw(10) << lmMean(i); } (*out) << endl;
                (*out) << setw(10) << "SD   "; for (i = 0; i < lmVariance.rows(); i++) { (*out) << setw(10) << lmVariance(i); } (*out) << endl;
                (*out) << setw(10) << "0.05 "; for (i = 0; i < mmData->cols(); i++) { (*out) << setw(10) << 1.0 * sum(abs((*mmData)(CUtility::getAll(),i)) < mrAlpha) / miDataSize; } (*out) << endl;
            }
        }
    }

    return rbLast;
}

unsigned int CAnalysis::cluster(TRVector& pmDosage, TIVector& pmPhenotypes) {
    unsigned int    riCount;
    TRMatrix        lmDosage;
    lmDosage.resize(pmDosage.rows(),1);
    lmDosage(CUtility::getAll(),0) = pmDosage;
    riCount = cluster(lmDosage, pmPhenotypes);
    pmDosage.resize(riCount);
    pmDosage = lmDosage(0,CUtility::getAll());
    return riCount;
}

unsigned int CAnalysis::cluster(TRMatrix& pmDosage, TIVector& pmPhenotypes) {
    unsigned int        i, liCluster, riCount;
    blitz::firstIndex   loIndex1;

    mmClusters.resize(miSubjects);
    mmClusters = loIndex1;
    for (i = 0; i < miSubjects; i++) {
        liCluster = first(mmPedigrees == mmPedigrees(i) && mmPhenotypes == mmPhenotypes(i));
        if (liCluster != i) {
            mmClusters(i) = -1;
            pmDosage(liCluster,CUtility::getAll()) += pmDosage(i,CUtility::getAll());
            pmDosage(i,CUtility::getAll()) = 0;
        }
    }

    riCount = 0;
    for (i = 0; i < miSubjects; i++) {
        if (mmClusters(i) != -1) {
            pmDosage(riCount,CUtility::getAll()) = pmDosage(i,CUtility::getAll());
            pmPhenotypes(riCount) = mmPhenotypes(i);
            riCount++;
        }
    }
    pmDosage.resizeAndPreserve(riCount, miSNPs);
    pmPhenotypes.resizeAndPreserve(riCount);

    return riCount;
}

unsigned int CAnalysis::cluster(ETRMatrix& pmDosage, TIVector& pmPhenotypes) {
    unsigned int        i, liCluster, riCount;
    blitz::firstIndex   loIndex1;

    mmClusters = loIndex1;
    for (i = 0; i < miSubjects; i++) {
        liCluster = first(mmPedigrees == mmPedigrees(i) && mmPhenotypes == mmPhenotypes(i));
        if (liCluster != i) {
            mmClusters(i) = -1;
            pmDosage.row(liCluster) += pmDosage.row(i);
            pmDosage.row(i).setZero();
        }
    }
    riCount = 0;
    for (i = 0; i < miSubjects; i++) {
        if (mmClusters(i) != -1) {
            pmDosage.row(riCount) = pmDosage.row(i);
            pmPhenotypes(riCount) = mmPhenotypes(i);
            riCount++;
        }
    }
    pmDosage.conservativeResize(riCount, miSNPs);
    pmPhenotypes.resizeAndPreserve(riCount);

    return riCount;
}
