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
//#include <unordered_map>
#include "population.h"

void CPopulation::setRelations () {
    unsigned int    i, j;
    CSubject       *loParent1, *loParent2, *loFather, *loMother;
    vector<string>  laPedigreeNames;
    vector<string>::iterator    liPedigreeEnd;
    blitz::Array<size_t,1>      laPedigrees;
    std::map<string, int>       loPedigreeMap;
    std::map<string,int>::iterator  liPedigree;

    laPedigreeNames.reserve(countSubjects());
    for (i = 0; i < countSubjects(); i++) {
        laPedigreeNames.push_back(mmSubjects(i).getPedigree());
        mmSubjects(i).setIdNumber(i);
        loParent1 = find(mmSubjects(i).getPedigree(), mmSubjects(i).getParent1());
        loParent2 = find(mmSubjects(i).getPedigree(), mmSubjects(i).getParent2());
        mmSubjects(i).setParents(loParent1, loParent2);
        if (loParent1) { loParent1->addChild(&mmSubjects(i)); }
        if (loParent2) { loParent2->addChild(&mmSubjects(i)); }
    }
    sort(laPedigreeNames.begin(), laPedigreeNames.end());
    liPedigreeEnd = unique(laPedigreeNames.begin(), laPedigreeNames.end());
    laPedigreeNames.resize(liPedigreeEnd - laPedigreeNames.begin());
    miPedigrees = laPedigreeNames.size();
    laPedigrees.resize(miPedigrees);

    for (i = 0; i < miPedigrees; i++) {
        loPedigreeMap.insert(std::pair<string,int>(laPedigreeNames[i],i));
    }
    for (i = 0; i < countSubjects(); i++) {
        mmSubjects(i).setPedigreeNumber(-1);
        liPedigree = loPedigreeMap.find(laPedigreeNames[i]);
        if (liPedigree != loPedigreeMap.end()) {
            mmSubjects(i).setPedigreeNumber(liPedigree->second);
        }
    }
    for (i = 0; i < countSubjects(); i++) {
        if (2 == mmSubjects(i).countParents()) {
            loFather = mmSubjects(i).getFather();
            loMother = mmSubjects(i).getMother();
            for (j = i+1; j < countSubjects(); j++) {
                if (mmSubjects(j).getFather() == loFather && mmSubjects(j).getMother() == loMother) {
                    mmSubjects(i).addSibling(&mmSubjects(j));
                    mmSubjects(j).addSibling(&mmSubjects(i));
                }
            }
        }
    }
}

void CPopulation::clean (const bool pbRelations) {
    unsigned int    i;
    stringstream    loId;
    CSubject       *loSubject, *loFather, *loMother;
    map<string,int> lmFamilies;
    pair<map<string,int>::iterator,bool> loFamily;

    if (pbRelations) { setRelations(); }

    //  Change the family and individual ids.
    for (i = 0; i < countSubjects(); i++) {
        loSubject = getSubjectP(i);
        loFamily  = lmFamilies.insert(pair<string,int>(loSubject->getPedigree(),lmFamilies.size()));
        loId << "F" << setfill('0') << setw(4) << loFamily.first->second;
        loSubject->setPedigree(loId.str());
        loId.str("");
        loId << "I" << setfill('0') << setw(6) << i;
        loSubject->setId(loId.str());
        loId.str("");
    }
    //  Update the parent ids.
    for (i = 0; i < countSubjects(); i++) {
        loSubject = getSubjectP(i);
        loFather  = loSubject->getFather();
        loMother  = loSubject->getMother();
        loSubject->setParent1(loFather != NULL ? loFather->getId() : "0");
        loSubject->setParent2(loMother != NULL ? loMother->getId() : "0");
    }
    //  Check that there are not too many families or individuals.
    if (100000000 < lmFamilies.size() || 100000000 < countSubjects()) {
        cerr << "Too many families (" << lmFamilies.size() << ") or individuals (" << countSubjects() << ")" << endl;
        exit(0);
    }
    //  Clean the phenotype data.
    for (i = 0; i < countSubjects(); i++) {
        if (mmPhenotypes(i) == -9 || mmPhenotypes(i) == 0) {
            mmPhenotypes(i) = 1;
        }
    }
}
