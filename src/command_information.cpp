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
 *  it under the terms of the GNU General PublmContingency(1,0) lmContingency(1,0)ense as published by
 *  the Free Software Foundation, either version 3 of the lmContingency(1,0)ense, or
 *  (at your option) any later version.
 *
 *  DWARF is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General PublmContingency(1,0) lmContingency(1,0)ense for more details.
 *
 *  You should have received a copy of the GNU General PublmContingency(1,0) lmContingency(1,0)ense
 *  along with DWARF.  If not, see <http://www.gnu.org/lmContingency(1,0)enses/>.
 *
 */
#include "command.h"
#include "command_information.h"
#include "population.h"
#include "utility.h"

enum {
    ceN = 0,
    ceA,
    ceC,
    ceG,
    ceT,
    ceI,
    ceILength,
    ceD,
    ceDLength,
    ceTotal
};

void CCommandInformation::run (const vector<string>& paParameters) {
    bool    lbPrint;

    preRun(paParameters);

    lbPrint = CUtility::getParameterString(paParameters,"print","false") != "false";

    mmData.resize(miSNPs,ceTotal);
    mmData = 0;

    (*out) << left << setw(10) << msPopulation << right << endl;
    (*out) << setw(10) << "Type"                  << setw(10) << "Subjects"                    << setw(10) << "SNPs" << endl;
    (*out) << setw(10) << moPopulation->getType() << setw(10) << moPopulation->countSubjects() << setw(10) << moPopulation->countSNPs() << endl;

    if (lbPrint) {
        switch (moPopulation->getType()) {
            case cePopulationNone:                     break;
            case cePopulationDosages:    dosages();    break;
            case cePopulationGenotypes:  genotypes();  break;
            case cePopulationHaplotypes: haplotypes(); break;
            case cePopulationPolyData:   polydata();   break;
        }
    }

    postRun();
}

void CCommandInformation::dosages() {
}

void CCommandInformation::haplotypes() {
    unsigned int    i, j;
    bool            lbPrint = false;
    ESNPCode        leSNP, leFather, leMother;
    ESNPCombined    leCombined;
    CSubject        loSubject;
    TRMatrix        lmContingency;
    TSNPVector      lmSNPs;

    lmSNPs.resize(miSubjects);
    lmContingency.resize(2,2);
    for (i = 0; i < miSNPs; i++) {
        switch (moPopulation->getType()) {
            case cePopulationGenotypes:  lmSNPs = moPopulation->getGenotypes()->getSNP(i); break;
            case cePopulationHaplotypes: lmSNPs = CHaplotype::transform(moPopulation->getHaplotypes()->getSNP(i)); break;
            default: msError << "Population type not implemented"; error();
        }
        lmContingency = 0;
        lmContingency(0,0) = 2.0 * sum((lmSNPs == ceSNPMinor) && (mmPhenotypes != 2)) + sum((lmSNPs == ceSNPHetero) && (mmPhenotypes != 2));
        lmContingency(0,1) = 2.0 * sum((lmSNPs == ceSNPMinor) && (mmPhenotypes == 2)) + sum((lmSNPs == ceSNPHetero) && (mmPhenotypes == 2));
        lmContingency(1,0) = 2.0 * sum((lmSNPs == ceSNPMajor) && (mmPhenotypes != 2)) + sum((lmSNPs == ceSNPHetero) && (mmPhenotypes != 2));
        lmContingency(1,1) = 2.0 * sum((lmSNPs == ceSNPMajor) && (mmPhenotypes == 2)) + sum((lmSNPs == ceSNPHetero) && (mmPhenotypes == 2));
        cout << lmContingency << endl;

        lmContingency = 0;
        for (j = 0; j < miSubjects; j++) {
            loSubject = *(moPopulation->getSubjectP(j));
            if (2 != loSubject.countParents() || 2 != mmPhenotypes(j)) { continue; }
            lbPrint = true;
            leSNP    = lmSNPs(loSubject.getIdNumber());
            leFather = lmSNPs(loSubject.getFather()->getIdNumber());
            leMother = lmSNPs(loSubject.getMother()->getIdNumber());
            if (ceHaplotype10 == leSNP)    { leSNP    = ceHaplotype01; }
            if (ceHaplotype10 == leFather) { leFather = ceHaplotype01; }
            if (ceHaplotype10 == leMother) { leMother = ceHaplotype01; }
            if (ceHaplotype01 != leFather && ceHaplotype01 != leMother) {
                if (ceHaplotype00 == leFather) { lmContingency(0,0) += 2; }
                if (ceHaplotype11 == leFather) { lmContingency(1,1) += 2; }
                if (ceHaplotype00 == leMother) { lmContingency(0,0) += 2; }
                if (ceHaplotype11 == leMother) { lmContingency(1,1) += 2; }
                continue;
            }
            leCombined = (ESNPCombined) (((ceHaplotype01 == leFather) ? (leMother << 2) : (leFather << 2)) + leSNP);
            switch (leCombined) {
                case ceSNPMajorMajor:
                case ceSNPMinorHetero:  lmContingency(1,0) += 1; lmContingency(0,0) += 1; break;
                case ceSNPMajorHetero:
                case ceSNPMinorMinor:   lmContingency(0,1) += 1; lmContingency(0,0) += 1; break;
                case ceSNPHeteroMajor:  lmContingency(1,0) += 2; break;
                case ceSNPHeteroMinor:  lmContingency(0,1) += 2; break;
                case ceSNPHeteroHetero: lmContingency(0,1) += 1; lmContingency(1,0) += 1; break;
                default: break;
            }
        }
        if (lbPrint) { cout << lmContingency << endl; }
    }
/*
    for (i = 0; i < miSubjects; i++) {
        cout << setw(5) << moPopulation->getSubjectP(i)->getPedigree() << " ";
        cout << setw(8) << moPopulation->getSubjectP(i)->getId() << " ";
        cout << setw(8) << moPopulation->getSubjectP(i)->getParent1() << " ";
        cout << setw(8) << moPopulation->getSubjectP(i)->getParent2() << " ";
        cout << moPopulation->getSubjectP(i)->getGender() << " ";
        cout << moPopulation->getPhenotype(i) << "   ";
        for (j = 0; j < miSNPs; j++) {
            cout << moPopulation->getHaplotypes()->get((int)i,(int)j) << " ";
        }
        cout << "   ";
        for (j = 0; j < miSNPs; j++) {
            cout << moPopulation->getHaplotypes()->get((int)i,(int)j,0) << " ";
        }
        cout << "   ";
        for (j = 0; j < miSNPs; j++) {
            cout << moPopulation->getHaplotypes()->get((int)i,(int)j,1) << " ";
        }
        cout << endl;
    }
*/
}

void CCommandInformation::genotypes() {
    unsigned int    i, j;
    TSNPVector      lmSNPs;
    TRMatrix        lmDosage;

    lmSNPs.resize(miSubjects);

    for (i = 0 ; i < miSNPs; i++) {
        lmSNPs = moPopulation->getGenotypes()->getSNP(i);
        cout << setw(3) << i << setw(8) << sum(lmSNPs == 0) << setw(8) << sum(lmSNPs == 1) << setw(8) << sum(lmSNPs == 2) << setw(8) << sum(lmSNPs == 3) << endl;
    }

    for (i = 0; i < miSNPs; i++) {
        for (j = 0; j < miSubjects; j++) {
            cout << moPopulation->getGenotypes()->get((int)j,(int)i) << " ";
        }
        cout << "   ";
        for (j = 0; j < miSubjects; j++) {
            cout << moPopulation->getGenotypes()->get(j,i) << " ";
        }
        cout << endl;
    }
}

void CCommandInformation::polydata() {
    unsigned int    i, j, liSNP;
    int             liIndex;
    TReal           lrTotal, lrMAF;
    TIVector        lmCount;
    TRVector        lmMAF;
    CSnip           loSNP;
    blitz::TinyVector<TReal,4> lmAF;

    lmMAF.resize(miSNPs);
    lmMAF = 0.0;
    lmAF = 0.0;
    lmCount.resize(miSNPs);
    lmCount = -1;

    for (i = 0 ; i < miSNPs; i++) {
        for (j = 0 ; j < miSubjects; j++) {
            liSNP = moPopulation->getPolyData()->get(j,i);
            switch (CPolyData::getPolyType(liSNP)) {
                case 'I': {
                    mmData(i,ceI)++;
                    mmData(i,ceILength) += CPolyData::getIndelLength(liSNP);
                    break;
                }
                case 'D': {
                    mmData(i,ceD)++;
                    mmData(i,ceDLength) += CPolyData::getIndelLength(liSNP);
                    break;
                }
                case 'S': {
                    switch (CPolyData::getQuadData(liSNP,1)) {
                        case 0: mmData(i,ceN)++; break;
                        case 1: mmData(i,ceT)++; break;
                        case 2: mmData(i,ceG)++; break;
                        case 4: mmData(i,ceC)++; break;
                        case 8: mmData(i,ceA)++; break;
                    }
                    switch (CPolyData::getQuadData(liSNP,2)) {
                        case 0: mmData(i,ceN)++; break;
                        case 1: mmData(i,ceT)++; break;
                        case 2: mmData(i,ceG)++; break;
                        case 4: mmData(i,ceC)++; break;
                        case 8: mmData(i,ceA)++; break;
                    }
                    break;
                }
            }
        }
    }

    (*out) << left << setw(5) << "Chr" << setw(20) << "Position" << right << setw(8) << "N" << setw(8) << "A" << setw(8) << "C" << setw(8) << "G" << setw(8) << "T" << setw(8) << "In" << setw(8) << "Del" << left << setw(8) << "AFs" << endl;
    for (i = 0; i < miSNPs; i++) {
        loSNP = moPopulation->getSNPInfo(i);
        (*out) << left << setw(5) << loSNP.getChromosome() << setw(20) << loSNP.getPosition()
               << right << setw(8) << mmData(i,ceN)
               << setw(8) << mmData(i,ceA) << setw(8) << mmData(i,ceC) << setw(8) << mmData(i,ceG) << setw(8) << mmData(i,ceT)
               << setw(8) << mmData(i,ceI) << setw(8) << mmData(i,ceD) << "  ";
        lmAF = 0.0;
        lrTotal = mmData(i,ceA) + mmData(i,ceC) + mmData(i,ceG) + mmData(i,ceT);
        lmAF(0) = mmData(i,ceA) / lrTotal;
        lmAF(1) = mmData(i,ceC) / lrTotal;
        lmAF(2) = mmData(i,ceG) / lrTotal;
        lmAF(3) = mmData(i,ceT) / lrTotal;
        out->precision(4);
        (*out) << left << fixed;
        lrMAF = 0.0;
        for (j = 0; j < 4; j++) {
            liIndex = blitz::maxIndex(lmAF);
            if (0.00001 < lmAF(liIndex)) {
                (*out) << setw(8) << lmAF(liIndex);
                lrMAF = lmAF(liIndex);
                lmAF(liIndex) = 0.0;
            } else {
                lmCount(i) = j;
                lmMAF(i) = lrMAF;
                break;
            }
        }
        if (-1 == lmCount(i)) {
            lmCount(i) = 4;
        }
        out->unsetf(ios_base::floatfield);
        (*out) << endl;
    }
    (*out) << "lmCount = [";
    for (i = 0; i < miSNPs; i++) {
        (*out) << lmCount(i) << ",";
    }
    (*out) << "];" << endl;
    (*out) << "lmMAF = [";
    for (i = 0; i < miSNPs; i++) {
        (*out) << lmMAF(i) << ",";
    }
    (*out) << "];" << endl;

}

