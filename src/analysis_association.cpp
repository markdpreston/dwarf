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

void CAnalysisAssociation::statistics(TRVector& rmStatistics, const bool pbShuffle) {
    if (pbShuffle) { CRandom::shuffle(mmPhenotypes); }
    calculation(rmStatistics, true);
}

void CAnalysisAssociation::asymptotics(TRVector& rmStatistics) {
    calculation(rmStatistics, false);
}

void CAnalysisAssociation::calculation(TRVector& rmCalculated, const bool pbStatistic) {
    unsigned int    i;
    CSnip           loSNP;
    TSNPVector      lmSNPs;
    TResult         lmChiSquared;

    mmPhenotypeControls.resize(miSubjects);
    mmPhenotypeCases.resize(miSubjects);
    mmSNPMajor.resize(miSubjects);
    mmSNPHetero.resize(miSubjects);
    mmSNPMinor.resize(miSubjects);
    mmContingency.resize(2,2);

    mmPhenotypeCases    = (mmPhenotypes == 2);
    mmPhenotypeControls = (mmPhenotypes != 2);

    //  Header line.
    if (! pbStatistic && out) {
        (*out) << "CHR   SNP           BP          A1   F_A       F_U       A2        CHISQ              P        OR" << endl;
        (*out) << fixed << setprecision(4);
    }

    //  Perform simple association tests.
    lmSNPs.resize(miSubjects);
    for (i = 0; i < miSNPs; i++) {
        switch (moPopulation->getType()) {
            case cePopulationGenotypes:  lmSNPs = moPopulation->getGenotypes()->getSNP(i); break;
            case cePopulationHaplotypes: lmSNPs = CHaplotype::transform(moPopulation->getHaplotypes()->getSNP(i)); break;
            default: msError << "Population type not implemented"; error();
        }
        lmChiSquared = pvalue(lmSNPs);
        if (-1 != lmChiSquared(0)) {
            rmCalculated(i) = lmChiSquared(pbStatistic ? 0 : 1);
        } else {
            rmCalculated(i) = 1.0;
        }
        if (! pbStatistic && out) {
            loSNP = moPopulation->getSNPInfo(i);
            (*out) << setw(6)  << left  << loSNP.getChromosome()
                   << setw(14) << left  << loSNP.getName()
                   << setw(10) << right << loSNP.getPosition() << "  "
                   << setw(5)  << left  << CSnip::transform(loSNP.getAllele(1))
                   << setw(10) << left  << mmContingency(1,1) / miSubjects
                   << setw(10) << left  << mmContingency(1,0) / miSubjects
                   << setw(5)  << left  << CSnip::transform(loSNP.getAllele(2))
                   << setw(10) << right << lmChiSquared(0);
            out->unsetf(ios_base::floatfield);
            (*out) << setw(15) << right << lmChiSquared(1)
                   << setw(10) << right << fixed << mmContingency(0,0) * mmContingency(1,1) / mmContingency(0,1) / mmContingency(1,0)
                   << endl;
        }
    }
    if (! pbStatistic && out) {
    }

}

TResult CAnalysisAssociation::pvalue(const TSNPVector& pmSNPs) {
    TResult     rmChiSquared;
    mmSNPMinor  = (pmSNPs == ceSNPMinor);
    mmSNPHetero = (pmSNPs == ceSNPHetero);
    mmSNPMajor  = (pmSNPs == ceSNPMajor);
    mmContingency(0,0) = 2.0 * sum(mmSNPMinor && mmPhenotypeControls) + sum(mmSNPHetero && mmPhenotypeControls);
    mmContingency(0,1) = 2.0 * sum(mmSNPMinor && mmPhenotypeCases)    + sum(mmSNPHetero && mmPhenotypeCases);
    mmContingency(1,0) = 2.0 * sum(mmSNPMajor && mmPhenotypeControls) + sum(mmSNPHetero && mmPhenotypeControls);
    mmContingency(1,1) = 2.0 * sum(mmSNPMajor && mmPhenotypeCases)    + sum(mmSNPHetero && mmPhenotypeCases);
    rmChiSquared       = CStatistics::chiSquared(mmContingency);
    return rmChiSquared;
}

TResult CAnalysisAssociation::pvalue(const TSNPVector& pmSNPs, const TIVector& pmPhenotypes) {
    TRMatrix    lmContingency;
    TResult     rmChiSquared;
    lmContingency.resize(2,2);
    lmContingency = 0;
    lmContingency(0,0) = 2.0 * sum((pmSNPs == ceSNPMinor) && (pmPhenotypes != 2)) + sum((pmSNPs == ceSNPHetero) && (pmPhenotypes != 2));
    lmContingency(0,1) = 2.0 * sum((pmSNPs == ceSNPMinor) && (pmPhenotypes == 2)) + sum((pmSNPs == ceSNPHetero) && (pmPhenotypes == 2));
    lmContingency(1,0) = 2.0 * sum((pmSNPs == ceSNPMajor) && (pmPhenotypes != 2)) + sum((pmSNPs == ceSNPHetero) && (pmPhenotypes != 2));
    lmContingency(1,1) = 2.0 * sum((pmSNPs == ceSNPMajor) && (pmPhenotypes == 2)) + sum((pmSNPs == ceSNPHetero) && (pmPhenotypes == 2));
    rmChiSquared       = CStatistics::chiSquared(lmContingency);
    return rmChiSquared;
}

void CAnalysisAssociation::test () {
    unsigned int    i;
    TRVector        lmDataLine;
    ifstream        fin;
    vector<string>  laParameters;
    CCommandSave    loSave;
    CCommandSystem  loSystem;
    lmDataLine.resize(miSNPs);
    //  Output dwarf results.
    cout << setw(24) << "Test Association";
    //  Save the data out to a temp ped and map file.
    loSave.setPopulation(moPopulation);
    laParameters.push_back("ped");
    laParameters.push_back("testassoc.ped");
    loSave.run(laParameters);
    laParameters.clear();
    laParameters.push_back("map");
    laParameters.push_back("testassoc.map");
    loSave.run(laParameters);
    laParameters.clear();
    //  Run plink.
    laParameters.push_back("plink --assoc --file testassoc --out testassoc > /dev/null");
    loSystem.run(laParameters);
    laParameters.clear();
    //  Extract the p-values.
    laParameters.push_back("awk -F' ' '{ print $9 }' < testassoc.assoc | tr '\n' ' ' | cut -b 3- > testassoc.out");
    loSystem.run(laParameters);
    laParameters.clear();
    //  Read in the p-values.
    fin.open("testassoc.out");
    for (i = 0; i < miSNPs; i++) {
        fin >> lmDataLine(i);
    }
    fin.close();
    //  Test the p-values.
    if (sum(abs(lmDataLine - mmDataLine)) < 0.001*miSNPs) {
        cout << "Pass" << endl;
    } else {
        cout << "Fail" << endl;
        cerr << setprecision(4);
        cerr << setw(24) << "     Calculated";
        for (i = 0; i < miSNPs; i++) {
            cerr << " " << mmDataLine(i);
        }
        cerr << endl << setw(24) << "     Reference";
        for (i = 0; i < miSNPs; i++) {
            cerr << " " << lmDataLine(i);
        }
        cerr << endl << endl;
    }
    //  Delete the temp files.
    laParameters.push_back("rm -f testassoc.*");
    loSystem.run(laParameters);
    laParameters.clear();
}
