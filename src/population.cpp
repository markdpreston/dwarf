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
#include "population.h"
#include "snp.h"

CPopulation::CPopulation () {
    clear();
}

CPopulation::~CPopulation () {
    clear();
}

void CPopulation::clear() {
    meType = cePopulationNone;
    miPedigrees = 0;
    moDosages.clear();
    moGenotypes.clear();
    moHaplotypes.clear();
    moPolyData.clear();
    moSNPs.clear();
    msTrait = "binary";
    resizeSubjects(0);
    resizeSNPs(0);
}

void CPopulation::resizeSubjects (const unsigned int piSubjects) {
    mmSubjects.resize(piSubjects);
    mmPhenotypes.resize(piSubjects);
    mmPhenotypes = 0;
    mmTraits.resize(piSubjects);
    mmTraits = 0.0;
}

void CPopulation::resizeSNPs (const unsigned int piSNPs) {
    if (piSNPs != moSNPs.countSNPs()) {
        moSNPs.resize(piSNPs);
    }
}

void CPopulation::resizeData () {
    switch (meType) {
        case cePopulationNone:       break;
        case cePopulationDosages:                           moGenotypes.resize(0,0); moHaplotypes.resize(0,0); moPolyData.resize(0,0); moDosages.resize   (countSubjects(),countSNPs()); break;
        case cePopulationGenotypes:  moDosages.resize(0,0);                          moHaplotypes.resize(0,0); moPolyData.resize(0,0); moGenotypes.resize (countSubjects(),countSNPs()); break;
        case cePopulationHaplotypes: moDosages.resize(0,0); moGenotypes.resize(0,0);                           moPolyData.resize(0,0); moHaplotypes.resize(countSubjects(),countSNPs()); break;
        case cePopulationPolyData:   moDosages.resize(0,0); moGenotypes.resize(0,0); moHaplotypes.resize(0,0);                         moPolyData.resize  (countSubjects(),countSNPs()); break;
    }
}

/*
void CPopulation::setMAFs () {
    unsigned int    i;
    mmMAFs.resize(miSNPs);
    mmMAFs = 0;
    for (i = 0; i < miSNPs; i++) {
        mmMAFs(i) = sum(moHaplotypes.getHaplotype(i,0)) + sum(moHaplotypes.getHaplotype(i,1));
    }
    mmMAFs /= 2*miSubjects;
}
*/
CSubject* CPopulation::find (const string psFamily, const string psId) {
    unsigned int        i;
    size_t              liHash;
    boost::hash<string> loHash;
    if (psFamily.length() > 0 && psFamily != "0" && psId.length() > 0 && psId != "0") {
        liHash = loHash(psFamily + psId);
        for (i = 0; i < countSubjects(); i++) {
            if (liHash == mmSubjects(i).getHash()) {
                if (psFamily == mmSubjects(i).getPedigree() && psId == mmSubjects(i).getId()) {
                    return &(mmSubjects(i));
                }
            }
        }
    }
    return NULL;
}

CPedigree CPopulation::getPedigree (const int piPedigree) {
    unsigned int    i;
    CPedigree       roPedigree;
    roPedigree.setPedigree(piPedigree);
    for (i = 0; i < countSubjects(); i++) {
        if (piPedigree == mmSubjects(i).getPedigreeNumber()) {
            roPedigree.add(&mmSubjects(i),mmPhenotypes(i));
        }
    }
    return roPedigree;
}

int CPopulation::getPedigreeMaximum () {
    unsigned int    i;
    int             liPedigree, riPedigree = 0;
    for (i = 0; i < countSubjects(); i++) {
        liPedigree = mmSubjects(i).getPedigreeNumber();
        if (liPedigree > riPedigree) {
            riPedigree = liPedigree;
        }
    }
    return riPedigree;
}

TUVector CPopulation::getPhenotypes (const int piPhenotype) {
    unsigned int    i, liCount;
    TUVector        rmPhenotypes;
    rmPhenotypes.resize(countPhenotypes(piPhenotype));
    liCount = 0;
    for (i = 0; i < countSubjects(); i++) {
        if (piPhenotype == mmPhenotypes((int)i)) {
            rmPhenotypes(liCount) = i;
            liCount++;
        }
    }
    return rmPhenotypes;
}

TSubjectPVector CPopulation::getAffected () {
    unsigned int    i, liCount;
    TSubjectPVector rmAffected;
    rmAffected.resize(countSubjects());
    liCount = 0;
    for (i = 0; i < countSubjects(); i++) {
        if (mmPhenotypes(i) == 2) {
            rmAffected(liCount) = &mmSubjects(i);
            liCount++;
        }
    }
    rmAffected.resizeAndPreserve(liCount);
    return rmAffected;
}

void CPopulation::keepRemove (const bool pbKeep, const bool pbSubjects, const TUVector& pmIndices) {
//    cout << "    keepRemove: " << psType << " " << psFile << endl;
cout << "Fix miPedigress and miAffected" << endl;
    unsigned int    i, liNewLength, liCount;
    TUVector        lmIndices;
    if (pbKeep) {
        if (pbSubjects) {
            liNewLength = pmIndices.rows();
            for (i = 0; i < liNewLength; i++) {
                 if (pmIndices(i) != i) { mmSubjects(i) = mmSubjects(pmIndices(i)); }
            }
            mmSubjects.resizeAndPreserve(liNewLength);
            moDosages.keepSubjects(pmIndices);
            moGenotypes.keepSubjects(pmIndices);
            moHaplotypes.keepSubjects(pmIndices);
        } else {
            liNewLength = pmIndices.rows();
            for (i = 0; i < liNewLength; i++) {
                 if (pmIndices(i) != i) { moSNPs.setSNP(i, moSNPs.getSNP(pmIndices(i))); }
            }
            moSNPs.resizeAndPreserve(liNewLength);
            moDosages.keepSNPs(pmIndices);
            moGenotypes.keepSNPs(pmIndices);
            moHaplotypes.keepSNPs(pmIndices);
        }
    } else {
        if (pbSubjects) {
            liNewLength = countSubjects() - pmIndices.rows() + sum(pmIndices >= countSubjects());
            lmIndices.resize(liNewLength);
            liCount = 0;
            for (i = 0; i < countSubjects(); i++) {
                 if (0 == sum(pmIndices == i)) { lmIndices(liCount) = i; liCount++; }
            }
            for (i = 0; i < liNewLength; i++) {
                 if (lmIndices(i) != i) { mmSubjects(i) = mmSubjects(lmIndices(i)); }
            }
            mmSubjects.resizeAndPreserve(liNewLength);
            moDosages.removeSubjects(pmIndices);
            moGenotypes.removeSubjects(pmIndices);
            moHaplotypes.removeSubjects(pmIndices);
        } else {
            liNewLength = countSNPs() - pmIndices.rows() + sum(pmIndices >= countSNPs());
            lmIndices.resize(liNewLength);
            liCount = 0;
            for (i = 0; i < countSNPs(); i++) {
                 if (0 == sum(pmIndices == i)) { lmIndices(liCount) = i; liCount++; }
            }
            for (i = 0; i < liNewLength; i++) {
                 if (lmIndices(i) != i) { moSNPs.setSNP(i, moSNPs.getSNP(lmIndices(i))); }
            }
            moSNPs.resizeAndPreserve(liNewLength);
            moDosages.removeSNPs(pmIndices);
            moGenotypes.removeSNPs(pmIndices);
            moHaplotypes.removeSNPs(pmIndices);
        }
    }
}

void CPopulation::setGenotypes (const CHaplotype& poHaplotypes) {
    moGenotypes.resize(poHaplotypes.miSubjects,poHaplotypes.miSNPs);
    memcpy(moGenotypes.maData,poHaplotypes.maData,moGenotypes.miSize);
}

/*
{
    int             i, j, liCount;
    string          lsLineAll, lsPedigree, lsId;
    stringstream    lsLine;
    ifstream        fin;
    CSnip           loSNP;
    CSubject        loSubject;
    TIVector        lmIndex;
    blitz::firstIndex   loIndex;
    fin.open(psFile.c_str());
    if (! fin.good()) {
        cerr << "Data file " << psFile << " not opened" << endl;
        exit(0);
    }
    getline(fin, lsLineAll);
    lsLine.clear();
    lsLine.str(lsLineAll);
    lsLine >> liCount;
    lmIndex.resize(liCount);
    lmIndex = -1;
    if (psType == "snps") {
        for (i = 0; i < liCount; i++) {
            lsLine >> loSNP.msChromosome >> loSNP.msName;
            for (j = 0; j < miSNPs; j++) { if (mmSNPs(j) == loSNP) { lmIndex(i) = j; break; } }
            if (lmIndex(i) < 0) { cout << "Not found: " << loSNP.msChromosome << " " << loSNP.msName << endl; }
        }
        if (pbKeep) {
            maDosage.keepSNPs(lmIndex);
            moGenotypes.keepSNPs(lmIndex);
            maHaplotype.keepSNPs(lmIndex);
            maSimulated.keepSNPs(lmIndex);
        } else {
            maDosage.removeSNPs(lmIndex);
            moGenotypes.removeSNPs(lmIndex);
            maHaplotype.removeSNPs(lmIndex);
            maSimulated.removeSNPs(lmIndex);
        }
    }
   if (psType == "subjects") {
        for (i = 0; i < liCount; i++) {
            lsLine >> lsPedigree >> lsId;
            loSubject.setPedigree(lsPedigree);
            loSubject.setId(lsId);
            for (j = 0; j < miSubjects; j++) { if (mmSubjects(j) == loSubject) { lmIndex(i) = j; break; } }
            if (lmIndex(i) < 0) { cout << "Not found: " << loSubject.getPedigree() << " " << loSubject.getId() << endl; }
        }
        if (pbKeep) {
            maDosage.keepSubjects(lmIndex);
            moGenotypes.keepSubjects(lmIndex);
            maHaplotype.keepSubjects(lmIndex);
            maSimulated.keepSubjects(lmIndex);
        } else {
            maDosage.removeSubjects(lmIndex);
            moGenotypes.removeSubjects(lmIndex);
            maHaplotype.removeSubjects(lmIndex);
            maSimulated.removeSubjects(lmIndex);
        }
    }
    fin.close();
}
 */
