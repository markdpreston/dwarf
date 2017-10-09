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
#include "process.h"
#include "command.h"
#include "statistics.h"
#include "command_load.h"
#include "command_power.h"
#include "command_simulation.h"

void CCommandPower::run(const vector<string>& paParameters) {
    int             i;
    unsigned int    liSNPs;
    string          lsMAFEmapFile;
    TReal           lrAlpha, lrBeta;
    TUnits          loUnit;
    blitz::Array<TPower,1> lmPower;
    CCommandLoad    loLoad;
    CPopulation     loPopulation;

    //  Get and check parameters.
    liSNPs = CUtility::getParameterInteger(paParameters, "snps", 0);
    lrAlpha = CUtility::getParameterReal(paParameters, "alpha", 0.05);
    lrBeta = CUtility::getParameterReal(paParameters, "beta", 0.8);
    lsMAFEmapFile = CUtility::getParameterString(paParameters, "emap", "");
    loUnit = CCommandSimulation::getUnit(paParameters);
    if (0 == liSNPs) {
        msError << "Missing snps parameter";
        error();
    }
    if ("" == lsMAFEmapFile) {
        msError << "Missing emap file";
        error();
    }
    if (0 == loUnit.miCaseSiblings + loUnit.miCaseParents + loUnit.miCaseUnrelateds) {
        msError << "No cases";
        error();
    }

    //  Set the population up.
    loPopulation.setType(cePopulationHaplotypes);
    loPopulation.resizeSNPs(liSNPs);

    //  Load SNP information.
    loLoad.setPopulation(&loPopulation);
    loLoad.run("emap", lsMAFEmapFile);

    //  Manually set-up the data.
    lmPower.resize(liSNPs);

    //  Analyse the snps.
    for (i = 0; i < (int) liSNPs; i++) {
        lmPower(i) = power(loPopulation.getSNPInfo(i), loUnit, lrAlpha, lrBeta);
    }

    //  Output the results.
    cout << setw(10) << "SNP" << setw(10) << "Count" << setw(10) << "P" << setw(10) << "Count" << setw(10) << "P" << endl;
    for (i = 0; i < (int) liSNPs; i++) {
        cout << setw(10) << loPopulation.getSNPInfo(i).getName() << setw(10) << loUnit.miCount << setw(10) << lmPower(i).mrPower << setw(10) << lmPower(i).miCount << setw(10) << lrBeta << endl;
    }
}

TPower CCommandPower::power (const CSnip& poSNP, const TUnits& poUnit, const TReal prAlpha, const TReal prBeta) {
    int             i, j;
    TReal           lrMu, lrSigma2, lrMAF, lrZalpha, lrZbeta, lrDeviation;
    TRVector        lmF(3), lmP(3), lmPi(3), lmPsi2(3), lmEpsilon(3), lmTau2(3);
    TRMatrix        lmM(3,3);
    TPower          roPower;

    //  Default return values.
    roPower.miCount = 0;
    roPower.mrPower = 0.0;

    //  Normal significance level.
    lrZalpha = sqrt(2.0) * boost::math::erfc_inv(2.0 * prAlpha);
    lrZbeta  = sqrt(2.0) * boost::math::erfc_inv(2.0 * prBeta);

    //  Calculation variables.
    lrMu        = 0.0;
    lrSigma2    = 0.0;
    lrDeviation = 0.0;
    lmF         = poSNP.getRelativeRisk(0.05);
    lmP         = poSNP.getHWE();
    lmPi        = pi(lmF);
    lmPsi2      = psi2(lmF);
    lmM         = m(lmF, lmP, poUnit.miCaseSiblings, poUnit.miCaseParents);

    //  Basic statistics of r affected siblings and 2 parents.
    if (0 == poUnit.miControlUnrelateds && 0 == poUnit.miCaseUnrelateds) {
        lrMu     = lmPi(0)*lmM(0,1) + lmPi(0)*lmM(1,0) + lmPi(1)*lmM(1,1) + lmPi(2)*lmM(1,2) + lmPi(2)*lmM(2,1);
        lrSigma2 = (     lmPsi2(0)*(lmM(0,1)+lmM(1,0)) +       lmPsi2(1)*lmM(1,1) +       lmPsi2(2)*(lmM(1,2)+lmM(2,1))) / poUnit.miCaseSiblings
                 + lmPi(0)*lmPi(0)*(lmM(0,1)+lmM(1,0)) + lmPi(1)*lmPi(1)*lmM(1,1) + lmPi(2)*lmPi(2)*(lmM(1,2)+lmM(2,1)) - lrMu*lrMu;
        lrMAF    = 0.0;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                lrMAF += 0.25 * (i + j) * lmM(i,j);
            }
        }
        if (0 == poUnit.miCaseParents && 2 == poUnit.miControlParents && 0 == poUnit.miControlSiblings) {
            lrDeviation = sqrt(lrMAF * (1 - lrMAF) / 4 / poUnit.miCaseSiblings);
        } else if (0 == poUnit.miCaseParents && 0 == poUnit.miControlParents && 0 < poUnit.miControlSiblings) {
            lrSigma2   += (lmM(1,0) + lmM(1,1) + lmM(1,2)) / 8 / poUnit.miControlSiblings;
            lrMAF      += lrMu * poUnit.miCaseSiblings / (poUnit.miCaseSiblings + poUnit.miControlSiblings);
            lrDeviation = sqrt(poUnit.miCaseSiblings + poUnit.miControlSiblings) * sqrt(lrMAF * (1 - lrMAF)) / sqrt(4 * poUnit.miCaseSiblings * poUnit.miControlSiblings);
        } else if (0 == poUnit.miCaseParents && 1 == poUnit.miControlParents && 0 == poUnit.miControlSiblings) {
            lrSigma2   += (lmM(1,0) + 2 * lmM(1,1) + lmM(1,2)) / 8 / (poUnit.miControlSiblings+1)
                        + (4 * lmM(2,0) - lmM(1,1)) / 8 / (poUnit.miControlSiblings+1) / (poUnit.miControlSiblings+1);
            lrMAF      += lrMu * poUnit.miCaseSiblings / (poUnit.miCaseSiblings + poUnit.miControlSiblings + 1);
            lrDeviation = sqrt(poUnit.miCaseSiblings + poUnit.miControlSiblings) * sqrt(lrMAF * (1 - lrMAF)) / sqrt(4 * poUnit.miCaseSiblings * poUnit.miControlSiblings);
        } else if (1 == poUnit.miCaseParents && 1 == poUnit.miControlParents && 0 == poUnit.miControlSiblings) {
            lrDeviation = sqrt(lrMAF * (1 - lrMAF) / 4 / poUnit.miCaseSiblings);
        } else {
            msLog << "No family power available";
            log();
        }
    } else {
        if (0 == poUnit.miCaseParents && 0 == poUnit.miControlParents && 0 == poUnit.miControlSiblings && 0 == poUnit.miCaseUnrelateds) {
            lmEpsilon(0) = 0.5 * lmF(1) / (lmF(1) + 1.0);                 // 01 10
            lmEpsilon(1) = (lmF(1) + lmF(2)) / (1.0 + 2*lmF(1) + lmF(2)); // 11
            lmEpsilon(2) = 0.5 * (2*lmF(2) + lmF(1)) / (lmF(2) + lmF(1)); // 12 21
            lmTau2(0)    = 0.25 * lmF(1) / pow(1.0 + lmF(1),2);                                // 01 10
            lmTau2(1)    = 0.5  * (lmF(2)*lmF(1) + 2*lmF(2) + lmF(1)) / pow(1.0 + 2*lmF(1) + lmF(2),2); // 11
            lmTau2(2)    = 0.25 * lmF(1) * lmF(2) / pow(lmF(1) + lmF(2),2);                             // 12 21
            lrMu         = lmEpsilon(0)*(lmM(0,1)+lmM(1,0)) + 0.5*(lmM(0,2)+lmM(2,0)) + lmEpsilon(1)*lmM(1,1) + lmEpsilon(2)*(lmM(1,2)+lmM(2,1)) + lmM(2,2) - poSNP.getMAF();
            lrSigma2     = (lmTau2(0)*(lmM(0,1)+lmM(1,0)) + lmTau2(1)*lmM(1,1) + lmTau2(2)*(lmM(1,2)+lmM(2,1))) / poUnit.miCaseSiblings
                           + lmEpsilon(0)*lmEpsilon(0)*(lmM(0,1)+lmM(1,0)) + 0.25*(lmM(0,2)+lmM(2,0)) + lmEpsilon(1)*lmEpsilon(1)*lmM(1,1) + lmEpsilon(2)*lmEpsilon(2)*(lmM(1,2)+lmM(2,1)) + lmM(2,2)
                           - pow(lrMu + poSNP.getMAF(),2) + 0.5*poSNP.getMAF()*(1-poSNP.getMAF())/poUnit.miControlUnrelateds;
            lrMAF        = poSNP.getMAF() + 2*poUnit.miCaseSiblings*lrMu / (poUnit.miCaseSiblings*poUnit.miControlUnrelateds + 2.0*poUnit.miCaseSiblings + poUnit.miControlUnrelateds);
            lrDeviation  = 0.5 * sqrt(poUnit.miCaseSiblings*poUnit.miControlUnrelateds + 2.0*poUnit.miCaseSiblings + poUnit.miControlUnrelateds) * sqrt(lrMAF * (1 - lrMAF))
                           / sqrt(poUnit.miCaseSiblings*poUnit.miControlUnrelateds);
        } else {
            msLog << "No case/control power available";
            log();
        }
    }

    roPower.miCount = (unsigned int) ceil(pow((lrZbeta * sqrt(lrSigma2) - lrZalpha * lrDeviation) / lrMu, (TReal) 2.0));
    roPower.mrPower = 1.0 - CStatistics::pvalueNormal((lrZalpha * lrDeviation - sqrt(poUnit.miCount) * lrMu) / sqrt(lrSigma2));

    return roPower;
}

TRVector CCommandPower::pi(const TRVector& pmF) {
    TRVector    rmPi(3);
    rmPi(0) = 0.25 * (pmF(1) - 1.0)    / (1.0 +   pmF(1)         ); // 10
    rmPi(1) = 0.5  * (pmF(2) - 1.0)    / (1.0 + 2*pmF(1) + pmF(2)); // 11
    rmPi(2) = 0.25 * (pmF(2) - pmF(1)) / (        pmF(1) + pmF(2)); // 12
    return rmPi;
}

TRVector CCommandPower::psi2(const TRVector& pmF) {
    TRVector    rmPsi2(3);
    rmPsi2(0) = 0.25 *  pmF(1)                             / pow(1.0 +   pmF(1),         2); // 10
    rmPsi2(1) = 0.5  * (pmF(1) + pmF(1)*pmF(2) + 2*pmF(2)) / pow(1.0 + 2*pmF(1) + pmF(2),2); // 11
    rmPsi2(2) = 0.25 *  pmF(1) * pmF(2)                    / pow(        pmF(1) + pmF(2),2); // 12
    return rmPsi2;
}

TRMatrix CCommandPower::m(const TRVector pmF, const TRVector pmP, const unsigned int piCaseSiblings, const unsigned int piCaseParents) {
    TReal       lrM;
    TRMatrix    rmM;

    rmM.resize(3,3);

    rmM(0,0) = pmP(0)*pmP(0);
    rmM(0,1) = pmP(0)*pmP(1) * pow(0.5 * (1.0 + pmF(1)), (TReal) piCaseSiblings);
    rmM(0,2) = pmP(0)*pmP(2) * pow(pmF(1), (TReal) piCaseSiblings);
    rmM(1,0) = rmM(0,1);
    rmM(1,1) = pmP(1)*pmP(1) * pow(0.25 * (1.0 + 2*pmF(1) + pmF(2)), (TReal) piCaseSiblings);
    rmM(1,2) = pmP(1)*pmP(2) * pow(0.5 * (pmF(1)+pmF(2)), (TReal) piCaseSiblings);
    rmM(2,0) = rmM(0,2);
    rmM(2,1) = rmM(1,2);
    rmM(2,2) = pmP(2)*pmP(2) * pow(pmF(2), (TReal) piCaseSiblings);

    if (1 == piCaseParents || 2 == piCaseParents) {
        rmM(1,0) = pmF(1) * rmM(1,0);
        rmM(1,1) = pmF(1) * rmM(1,1);
        rmM(1,2) = pmF(1) * rmM(1,2);
        rmM(2,0) = pmF(2) * rmM(2,0);
        rmM(2,1) = pmF(2) * rmM(2,1);
        rmM(2,2) = pmF(2) * rmM(2,2);
    }
    if (2 == piCaseParents) {
        rmM(0,1) = pmF(1) * rmM(1,1);
        rmM(0,2) = pmF(2) * rmM(1,2);
        rmM(1,1) = pmF(1) * rmM(1,1);
        rmM(1,2) = pmF(2) * rmM(1,2);
        rmM(2,1) = pmF(1) * rmM(2,1);
        rmM(2,2) = pmF(2) * rmM(2,2);
    }

    blitz::secondIndex j;
    lrM = sum(sum(rmM,j));
//cout << setw(10) << lrM << endl;
//cout << setw(10) << "preM" << rmM;
    rmM /= lrM;

//lrM = sum(sum(rmM,j));
//cout << setw(10) << lrM << endl;

    return rmM;
}
