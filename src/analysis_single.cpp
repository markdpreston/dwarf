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

bool CAnalysisSingle::preRun(const vector<string>& paParameters) {
    msMethod = CUtility::getParameterString(paParameters, "method", "unrelated");
    return CAnalysis::preRun(paParameters);
}

void CAnalysisSingle::statistics(TRVector& rmStatistics, const bool pbShuffle) {
    if (pbShuffle) { CRandom::shuffle(mmPhenotypes); }
    calculation(rmStatistics, true);
}

void CAnalysisSingle::asymptotics(TRVector& rmAsymptotics) {
    calculation(rmAsymptotics, false);
}

void CAnalysisSingle::calculation(TRVector& rmCalculation, const bool pbStatistic) {
    unsigned int    i;
    string          lsMethod;
    TRVector        lmP;
    TSubjectVector *lmSubjects;
    TSNPVector      lmSNPs;
    TResult         loChiSquared;

    lmSNPs.resize(miSubjects);
    lmP.resize(miSNPs);
    lmSubjects = moPopulation->getSubjects();

    for (i = 0; i < miSNPs; i++) {
        switch (mePopulation) {
            case cePopulationGenotypes:  lmSNPs = moPopulation->getGenotypes()->getSNP(i); break;
            case cePopulationHaplotypes: lmSNPs = CHaplotype::transform(moPopulation->getHaplotypes()->getSNP(i)); break;
            default: msError << "Not supported population type"; error(); break;
        }
        if ("unrelated" == msMethod) {
            loChiSquared = CAnalysisAssociation::pvalue(lmSNPs, mmPhenotypes);
        } else if ("family" == msMethod) {
            loChiSquared = CAnalysisTDT::pvalue(lmSubjects, lmSNPs, mmPhenotypes);
        } else {
            msError << "Bad method";
            error();
        }
        lmP(i) = loChiSquared(pbStatistic ? 0 : 1);
    }

    mmDataLine(0) = min(lmP);
    mmDataLine(1) = mmDataLine(0) * miSNPs;
    mmDataLine(2) = 1 - pow(1 - mmDataLine(0), (TReal) miSNPs);
}

void CAnalysisSingle::test() {
    unsigned int    i;
    TReal           lrP;
    TRVector        lmP;
    ifstream        fin;
    vector<string>  laParameters;
    CCommandSave    loSave;
    CCommandSystem  loSystem;
    //  Output dwarf results.
    cout << setw(24) << "Test Single";
    //  Save the data out to a temp ped and map file.
    loSave.setPopulation(moPopulation);
    laParameters.push_back("ped");
    laParameters.push_back("testsingle.ped");
    loSave.run(laParameters);
    laParameters.clear();
    laParameters.push_back("map");
    laParameters.push_back("testsingle.map");
    loSave.run(laParameters);
    laParameters.clear();
    //  Run plink.
    if ("family" == msMethod) {
      laParameters.push_back("plink --tdt --file testsingle --out testsingle > /dev/null");
      loSystem.run(laParameters);
      laParameters.clear();
      //  Extract the p-values.
      laParameters.push_back("awk -F' ' '{ print $10 }' < testsingle.tdt | tr '\n' ' ' | cut -b3- > testsingle.out");
      loSystem.run(laParameters);
      laParameters.clear();
    } else {
      laParameters.push_back("plink --assoc --file testsingle --out testsingle > /dev/null");
      loSystem.run(laParameters);
      laParameters.clear();
      //  Extract the p-values.
      laParameters.push_back("awk -F' ' '{ print $9 }' < testsingle.assoc | tr '\n' ' ' | cut -b 3- > testsingle.out");
      loSystem.run(laParameters);
      laParameters.clear();
    }
    //  Read in the p-values.
    lmP.resize(miSNPs);
    fin.open("testsingle.out");
    for (i = 0; i < miSNPs; i++) {
      fin >> lmP(i);
    }
    fin.close();
    lrP = min(lmP);
    if (abs(lrP - mmDataLine(0)) < 0.001) {
        cout << "Pass" << endl;
    } else {
        cout << "Fail" << endl;
        cerr << setprecision(4) << left;
        cerr << setw(24) << "     Calculated";
        for (i = 0; i < (unsigned int) mmDataLine.rows(); i++) {
            cerr << setw(10) << mmDataLine(i);
        }
        cerr << endl << setw(24) << "     Reference" << lrP << endl;
    }

    //  Delete the temp files.
    laParameters.push_back("rm -f testsingle.*");
    loSystem.run(laParameters);
    laParameters.clear();
}
