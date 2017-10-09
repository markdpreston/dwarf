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
#include "analysis_score.h"
#include "process.h"
#include "statistics.h"
#include "random.h"
#include "command_data.h"

CAnalysisScore::CAnalysisScore () {
    moSolver1 = NULL;
    moSolver2 = NULL;
}

CAnalysisScore::~CAnalysisScore () {
   clean();
}

void CAnalysisScore::clean() {
    CAnalysis::clean();
    mmP.resize(0);
    mmCorrelation.resize(0,0);
    mmVBlitz.resize(0,0);
    mmU.resize(0);
    mmV.resize(0,0);
    mmEvals.resize(0);
    mmTemp.resize(0,0);
    mmX.resize(0,0);
    mmY.resize(0,1);
    mmMeanX.resize(0,0);
    mmAffected.resize(0);
    mmBminusC.resize(0,0);
    mmMean.resize(0,0);
    if (NULL != moSolver1) { delete moSolver1; }
    if (NULL != moSolver2) { delete moSolver2; }
    moSolver1 = NULL;
    moSolver2 = NULL;
#ifdef INCLUDE_R
//    CInterfaceR::soR.assign(1.0,"mmCorrelation");
//    CInterfaceR::removeData("mmCorrelation");
#endif
}

bool CAnalysisScore::preRun(const vector<string>& paParameters) {
    bool    rbSuccess;

    rbSuccess = CAnalysis::preRun(paParameters);
    msMethod = CUtility::getParameterString(paParameters, "method", "unrelated");
    miAffected = moPopulation->countAffected();

    if (rbSuccess && (int) miSNPs != mmP.rows()) {
        //  General/calculation.
        mmZeros.resize(miSNPs);
        mmP.resize(miSNPs);
        mmCorrelation.resize(miSNPs,miSNPs);
        mmVBlitz.resize(miSNPs,miSNPs);
        mmU.resize(miSNPs);
        mmV.resize(miSNPs,miSNPs);
        mmEvals.resize(miSNPs);
        mmTemp.resize(miSNPs,miSNPs);
        if (NULL != moSolver1) { delete moSolver1; }
        if (NULL != moSolver2) { delete moSolver2; }
        moSolver1 = new Eigen::SelfAdjointEigenSolver<ETRMatrix>(miSNPs);
        moSolver2 = new Eigen::ComplexEigenSolver<ETRMatrix>(miSNPs);
    }
    if (rbSuccess && msMethod == "unrelated" && ((int) miSNPs != mmMeanX.cols() || (int) miSubjects != mmMeanX.rows())) {
        //  Unrelateds.
        mmClusters.resize(miSubjects);
        mmX.resize(miSubjects,miSNPs);
        mmY.resize(miSubjects,1);
        mmMeanX.resize(miSubjects,miSNPs);
    }
    if (rbSuccess && msMethod == "family" && ((int) miSNPs != mmBminusC.cols() || (int) miAffected != mmBminusC.rows())) {
        //  Family.
        mmAffected.resize(miAffected);
        mmBminusC.resize(miAffected,miSNPs);
        mmMean.resize(miAffected,miSNPs);
    }

    if (0 == soProcess->miR.rinside || 0 == soProcess->miR.rcpp || 0 == soProcess->miR.mvtnorm) {
        msWarning << "R libraries required: RInside, Rcpp, mvtnorm";
        warning();
    }

    return rbSuccess;
}

void CAnalysisScore::statistics(TRVector& rmStatistics, const bool pbShuffle) {
    if (pbShuffle) { CRandom::shuffle(mmPhenotypes); }
    calculation(rmStatistics, mmDataLine);
}

void CAnalysisScore::asymptotics(TRVector& rmAsymptotics) {
    calculation(rmAsymptotics, mmDataLine);
}

void CAnalysisScore::calculation(TRVector& rmStatistics, TRVector& rmAsymptotics) {
    unsigned int        i, j, liZeros, liSNPs, liCount1, liCount2;
    TReal               lrMeanY;
    blitz::firstIndex   loIndex1;
    ESNPCode            leSNP, leFather, leMother;
    CSubject           *loSubject;
    THetero             laBC;

    //  Need to know the size to declare these Eigen variables.
    mmU.setZero(miSNPs);
    mmV.setZero(miSNPs,miSNPs);
    mmEvals.setZero(miSNPs);
    mmTemp.setZero(miSNPs,miSNPs);

    if ("unrelated" == msMethod) {
        //  Resize, if not needed then Eigen will not bother changing anything!
        mmMeanX.resize(miSubjects, miSNPs);
        mmX.resize(miSubjects, miSNPs);
        mmY.resize(miSubjects);
        //  Calculate the U and V matrices in Eigen format.
        switch (mePopulation) {
            case cePopulationGenotypes:  mmX = CUtility::makeEigen(moPopulation->getGenotypes()->getDosage(false)); break;
            case cePopulationHaplotypes: mmX = CUtility::makeEigen(moPopulation->getHaplotypes()->getDosage());     break;
            case cePopulationDosages:    mmX = CUtility::makeEigen(moPopulation->getDosages()->getDosage());        break;
            default: msError << "Not supported population type"; error(); break;
        }
        if (mbCluster) {
            liCount1 = cluster(mmX, mmPhenotypes);
            mmMeanX.resize(liCount1,miSNPs);
            mmY.resize(liCount1);
            mmY = CUtility::makeEigen(mmPhenotypes);
        } else {
            liCount1 = miSubjects;
            mmY = CUtility::makeEigen(mmPhenotypes);
        }
        //  Pre-solving calculations.
        mmMeanX.row(0) = mmX.colwise().sum() / liCount1;
        mmMeanX = mmMeanX.row(0).colwise().replicate(liCount1);
        mmY.array() -= 1;
        //  Change X to (X - X-bar) and Y to (Y Y-bar) for score and covariance calculation.
        lrMeanY = mmY.sum() / liCount1;
        mmY.array() -= lrMeanY;
        mmU = mmX.transpose() * mmY;
        mmX = mmX - mmMeanX;
        mmV = lrMeanY * (1 - lrMeanY) * mmX.transpose() * mmX;
    } else if ("family" == msMethod) {
        //  Calculate the U and V matrices in Eigen format.
        mmAffected = moPopulation->getAffected();
        mmBminusC.resize(miAffected,miSNPs);
        mmBminusC.setZero(miAffected,miSNPs);
        liCount1 = 0;
        for (i = 0; i < miAffected; i++) {
            loSubject = mmAffected(i);
            if (2 != loSubject->countParents()) {
                continue;
            }
            for (j = 0; j < miSNPs; j++) {
                switch (moPopulation->getType()) {
                    case cePopulationGenotypes: {
                        leSNP    = moPopulation->getGenotypes()->get(loSubject->getIdNumber(),j);
                        leFather = moPopulation->getGenotypes()->get(loSubject->getFather()->getIdNumber(),j);
                        leMother = moPopulation->getGenotypes()->get(loSubject->getMother()->getIdNumber(),j);
                        break;
                    }
                    case cePopulationHaplotypes: {
                        leSNP    = CHaplotype::transform(moPopulation->getHaplotypes()->get(loSubject->getIdNumber(),j));
                        leFather = CHaplotype::transform(moPopulation->getHaplotypes()->get(loSubject->getFather()->getIdNumber(),j));
                        leMother = CHaplotype::transform(moPopulation->getHaplotypes()->get(loSubject->getMother()->getIdNumber(),j));
                        break;
                    }
                    default: msError << "Population type not implemented"; error();
                }
                laBC = score(leSNP, leFather, leMother);
                mmBminusC(liCount1,j) = 0.5 * laBC.miB - 0.5 * laBC.miC;
            }
            liCount1++;
        }
        mmBminusC.conservativeResize(liCount1,miSNPs);
        mmMean.resize(liCount1,miSNPs);
        //  Calculate the score.
        mmU = mmBminusC.colwise().sum();
        //  Change B-C to (B-C - \bar{B-C}) for covariance calculation.
        mmMean.row(0) = mmU / liCount1;
        mmMean = mmMean.row(0).colwise().replicate(liCount1);
        mmBminusC = mmBminusC - mmMean;
        //  Calculated the score covariance matrix.
        mmV = mmBminusC.transpose() * mmBminusC;
    } else {
        msError << "Unrecognised score test: " << msMethod;
        error();
    }

    mmZeros = 0;
    for (i = 0; i < miSNPs; i++) {
        if (abs(mmV.diagonal()(i)) < 1e-8) {
            mmZeros(i) = 1;
        }
    }
    liZeros = (unsigned int) sum(mmZeros);

    //  Use the existing data.
    if (liZeros == miSNPs) {
        mmDataLine(0) = 1.0;
        mmDataLine(1) = 1.0;
        mmDataLine(2) = 1.0;
        mmDataLine(3) = 1.0;
        mmDataLine(4) = 1.0;
    } else if (0 == liZeros) {
        pvalues(rmStatistics, rmAsymptotics, miSNPs, mmU, mmV, mmP, mmCorrelation, mmVBlitz, mmEvals, mmTemp);
    } else {
        TRMatrix    lmCorrelation, lmVBlitz;
        ETRVector   lmU, lmEvals, lmP;
        ETRMatrix   lmV, lmTemp;
        liSNPs = miSNPs - liZeros;
        lmU.resize(liSNPs);
        lmV.resize(liSNPs,liSNPs);
        lmP.resize(liSNPs);
        lmCorrelation.resize(liSNPs,liSNPs);
        lmVBlitz.resize(liSNPs,liSNPs);
        lmEvals.resize(liSNPs);
        lmTemp.resize(liSNPs,liSNPs);
        liCount1 = 0;
        for (i = 0; i < miSNPs; i++) {
            if (0 == mmZeros(i)) {
                lmU(liCount1) = mmU(i);
                liCount2 = 0;
                for (j = 0; j < miSNPs; j++) {
                    if (0 == mmZeros(j)) {
                        lmV(liCount1,liCount2) = mmV(i,j);
                        liCount2++;
                    }
                }
                liCount1++;
            }
        }
        pvalues(rmStatistics, rmAsymptotics, liSNPs, lmU, lmV, lmP, lmCorrelation, lmVBlitz, lmEvals, lmTemp);
    }
}

void CAnalysisScore::pvalues(TRVector& rmStatistics, TRVector& rmAsymptotics, const int piSNPs, ETRVector& pmU, ETRMatrix& pmV, ETRVector& pmP, TRMatrix& pmCorrelation, TRMatrix& pmVBlitz, ETRVector& pmEvals, ETRMatrix& pmTemp) {
    TReal   lrP;
    string  lsP, lsN;
    TResult lrResult;
    //  UminP - Conneely and Boehnke
    pmP = pmU.cwiseAbs().cwiseQuotient(pmV.diagonal().cwiseSqrt());
    lrP = pmP.maxCoeff();
    if (1 == soProcess->miR.rinside && 1 == soProcess->miR.rcpp && 1 == soProcess->miR.mvtnorm) {
        lsN = boost::lexical_cast<std::string>(piSNPs);
        lsP = boost::lexical_cast<std::string>(lrP);
        pmVBlitz = CUtility::makeBlitz(pmV);
        pmCorrelation = CStatistics::correlation(pmVBlitz, true);
        CInterfaceR::setMatrix("mmCorrelation",pmCorrelation);
        rmStatistics(0)  = lrP;
        rmAsymptotics(0) = CInterfaceR::runCommandR("as.numeric(1 - pmvnorm(lower=c(rep(-" + lsP + "," + lsN + ")),upper=c(rep(" + lsP + "," + lsN + ")),mean=c(rep(0," + lsN + ")),sigma=mmCorrelation))");
        CInterfaceR::runCommand("remove(list=ls())");
    } else {
        rmStatistics(0)  = lrP;
        rmAsymptotics(0) = 1.0;
    }
    //  SSU
    moSolver1->compute(pmV,Eigen::EigenvaluesOnly);
    pmEvals          = moSolver1->eigenvalues();
    lrResult         = chiSquared(pmEvals, pmU.squaredNorm());
    rmStatistics(1)  = lrResult(0);
    rmAsymptotics(1) = lrResult(1);
    //  SSUw
    pmTemp.setConstant(0);
    //  Collapse these two lines.
    pmTemp.diagonal() = pmV.diagonal();
    pmTemp = pmTemp.inverse();
    moSolver2->compute(pmV * pmTemp,Eigen::EigenvaluesOnly);
    pmEvals          = real(moSolver2->eigenvalues().array());
    lrResult         = chiSquared(pmEvals, pmU.dot(pmU.cwiseQuotient(pmV.diagonal())));
    rmStatistics(2)  = lrResult(0);
    rmAsymptotics(2) = lrResult(1);
    //  Score
    rmStatistics(3)  = pmU.dot(pmV.inverse() * pmU);
    rmAsymptotics(3) = CStatistics::pvalueChiSquared(rmStatistics(3), piSNPs);
    //  Sum
    rmStatistics(4)  = pmU.sum() * pmU.sum() / pmV.sum();
    rmAsymptotics(4) = CStatistics::pvalueChiSquared(rmStatistics(4), 1);
//cout << rmStatistics << endl << rmAsymptotics << endl;
}

void CAnalysisScore::test() {
    int         i;
    TRVector    lmResult;
    ETRVector   lmY;
    TRMatrix    lmX;

    cout << setw(24) << "Test Score";

    lmX.resize(miSubjects,miSNPs);
    lmY.resize(miSubjects,1);
    lmResult.resize(5);

    lmX = moPopulation->getHaplotypes()->getDosage();
    lmY = CUtility::makeEigen(moPopulation->getPhenotypes());
    lmY.array() -= 1;

    if (1 == soProcess->miR.rinside && 1 == soProcess->miR.rcpp && 1 == soProcess->miR.mvtnorm) {
        CInterfaceR::setMatrix("lmX", lmX);
        CInterfaceR::setVector("lmY", lmY);
        CInterfaceR::runCommand("source('analysis_score.R')");
        lmResult = CInterfaceR::runCommandRV("analysis_score(lmY,lmX)");

        if (sum(abs(mmDataLine - lmResult)) < 0.0001) {
            cout << "Pass" << endl;
        } else {
            cout << "Fail" << endl;
            cerr << setprecision(4) << left << setw(24) << "     Calculated";
            for (i = 0; i < mmDataLine.rows(); i++) {
                cerr << setw(10) << mmDataLine(i);
            }
            cerr << endl << setw(24) << "     Reference";
            for (i = 0; i < lmResult.rows(); i++) {
                cerr << setw(10) << lmResult(i);
            }
            cerr << endl;
        }
    } else {
        cout << "Untested" << endl;
    }
}

TResult CAnalysisScore::chiSquared (const ETRVector& pmEvals, const TReal prScore) {
    TReal   lrSum1, lrSum2, lrSum3, lrA, lrB, lrD;
    TResult rrChiSquared;
    lrSum1 = pmEvals.sum();
    lrSum2 = pmEvals.array().pow(2).sum();
    lrSum3 = pmEvals.array().pow(3).sum();
    lrA = lrSum3 / lrSum2;
    lrB = lrSum1 - lrSum2*lrSum2 / lrSum3;
    lrD = lrSum2*lrSum2*lrSum2 / lrSum3/lrSum3;
    rrChiSquared(0) = fabs((prScore - lrB) / lrA);
    rrChiSquared(1) = CStatistics::pvalueChiSquared(fabs((prScore - lrB) / lrA), lrD);
    return rrChiSquared;
}

THetero CAnalysisScore::score (ESNPCode peSNP, ESNPCode peFather, ESNPCode peMother) {
    ESNPCombined    leCombined;
    THetero         raBC;
    raBC.miB = 0;
    raBC.miC = 0;
    if (ceSNPHetero == peFather || ceSNPHetero == peMother) {
        leCombined = (ESNPCombined) (((ceSNPHetero == peFather) ? (peMother << 2) : (peFather << 2)) + peSNP);
        switch (leCombined) {
            case ceSNPMajorMajor:
            case ceSNPMinorHetero:  raBC.miC = 1; break;
            case ceSNPMajorHetero:
            case ceSNPMinorMinor:   raBC.miB = 1; break;
            case ceSNPHeteroMajor:  raBC.miC = 2; break;
            case ceSNPHeteroMinor:  raBC.miB = 2; break;
            case ceSNPHeteroHetero: raBC.miB = 1; raBC.miC = 1; break;
            default: break;
        }
    }
    return raBC;
}
