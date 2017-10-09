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
#include "subject.h"
#include "statistics.h"

void CAnalysisTDT::statistics(TRVector& rmStatistics, const bool pbShuffle) {
    if (pbShuffle) { CRandom::shuffle(mmPhenotypes); }
    calculation(rmStatistics, true);
}

void CAnalysisTDT::asymptotics(TRVector& rmStatistics) {
    calculation(rmStatistics, false);
}

void CAnalysisTDT::calculation(TRVector& rmCalculated, const bool pbStatistic) {
    unsigned int    i;
    TSubjectVector *lmSubjects;
    TSNPVector      lmSNPs;
    TResult         lmChiSquared;

    // Perform tdt tests.
    lmSubjects = moPopulation->getSubjects();
    lmSNPs.resize(miSubjects);
    for (i = 0; i < miSNPs; i++) {
        switch (moPopulation->getType()) {
            case cePopulationGenotypes:  lmSNPs = moPopulation->getGenotypes()->getSNP(i); break;
            case cePopulationHaplotypes: lmSNPs = CHaplotype::transform(moPopulation->getHaplotypes()->getSNP(i)); break;
            default: msError << "Population type not implemented"; error();
        }
        lmChiSquared = pvalue(lmSubjects, lmSNPs, mmPhenotypes);
//cout << "TDT " << lmChiSquared << endl;
        if (-1 != lmChiSquared(0)) {
            rmCalculated(i) = lmChiSquared(pbStatistic ? 0 : 1);
        } else {
            rmCalculated(i) = 1.0;
        }
    }
}

TResult CAnalysisTDT::pvalue(const TSubjectVector *pmSubjects, const TSNPVector& pmSNPs, const TIVector& pmPhenotypes) {
    int             i, liB, liC;
    ESNPCode        leSNP, leFather, leMother;
    ESNPCombined    leCombined;
    CSubject        loSubject;
    TResult         rmChiSquared;

    rmChiSquared(0) = -1.0;
    rmChiSquared(1) =  1.0;
    liB = 0;
    liC = 0;
    for (i = 0; i < pmSubjects->rows(); i++) {
        loSubject = (*pmSubjects)(i);
        if (2 != loSubject.countParents() || 2 != pmPhenotypes(i)) { continue; }
        leSNP    = pmSNPs(loSubject.getIdNumber());
        leFather = pmSNPs(loSubject.getFather()->getIdNumber());
        leMother = pmSNPs(loSubject.getMother()->getIdNumber());
        if (ceHaplotype10 == leSNP)    { leSNP    = ceHaplotype01; }
        if (ceHaplotype10 == leFather) { leFather = ceHaplotype01; }
        if (ceHaplotype10 == leMother) { leMother = ceHaplotype01; }
        if (ceHaplotype01 != leFather && ceHaplotype01 != leMother) { continue; }
        leCombined = (ESNPCombined) (((ceHaplotype01 == leFather) ? (leMother << 2) : (leFather << 2)) + leSNP);
        switch (leCombined) {
            case ceSNPMajorMajor:
            case ceSNPMinorHetero:            liC += 1; break;
            case ceSNPMajorHetero:
            case ceSNPMinorMinor:   liB += 1;           break;
            case ceSNPHeteroMajor:            liC += 2; break;
            case ceSNPHeteroMinor:  liB += 2;           break;
            case ceSNPHeteroHetero: liB += 1; liC += 1; break;
            default: break;
        }
//cout << "lmScore(liSnp," << i << ") = " << liB - liC << ";" << endl;
    }


    if (0 != liB || 0 != liC) {
        rmChiSquared(0) = pow(liB - liC,2.0) / (liB + liC);
        rmChiSquared(1) = CStatistics::pvalueChiSquared(pow(liB - liC,2.0) / (liB + liC),1.0);
    }

    return rmChiSquared;
}

void CAnalysisTDT::test () {
    int             i;
    TRVector        lmDataLine;
    ifstream        fin;
    vector<string>  laParameters;
    CCommandSave    loSave;
    CCommandSystem  loSystem;
    lmDataLine.resize(miSNPs);
    //  Output dwarf results.
    cout << setw(24) << "Test TDT";
    //  Save the data out to a temp ped and map file.
    loSave.setPopulation(moPopulation);
    laParameters.push_back("ped");
    laParameters.push_back("testtdt.ped");
    loSave.run(laParameters);
    laParameters.clear();
    laParameters.push_back("map");
    laParameters.push_back("testtdt.map");
    loSave.run(laParameters);
    laParameters.clear();
    //  Run plink.
    laParameters.push_back("plink --tdt --file testtdt --out testtdt > /dev/null");
    loSystem.run(laParameters);
    laParameters.clear();
    //  Print the p-values.
    laParameters.push_back("awk -F' ' '{ print $10 }' < testtdt.tdt | tr '\n' ' ' | cut -b3- > testtdt.out");
    loSystem.run(laParameters);
    laParameters.clear();
    //  Read in the p-values.
    fin.open("testtdt.out");
    for (i = 0; i < mmDataLine.rows(); i++) {
        fin >> lmDataLine(i);
    }
    fin.close();
    if (sum(abs(lmDataLine - mmDataLine)) < 0.001) {
        cout << "Pass" << endl;
    } else {
        cout << "Fail" << endl;
        cerr << setprecision(4);
        cerr << setw(24) << "     Calculated";
        for (i = 0; i < mmDataLine.rows(); i++) {
            cerr << " " << mmDataLine(i);
        }
        cerr << endl << setw(24) << "     Reference";
        for (i = 0; i < mmDataLine.rows(); i++) {
            cerr << " " << lmDataLine(i);
        }
        cerr << endl << endl;
    }
    //  Delete the temp files.
    laParameters.push_back("rm -f testtdt.*");
    loSystem.run(laParameters);
    laParameters.clear();
}
