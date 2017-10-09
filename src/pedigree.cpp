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
#include "pedigree.h"

CPedigree::CPedigree () {
    clear();
}

void CPedigree::clear() {
    miSize = 0;
    mmSubjects.resize(miSize);
    mmPhenotypes.resize(miSize);
}

void CPedigree::add(CSubject *poSubject, const int piPhenotype) {
    miSize++;
    mmSubjects.resizeAndPreserve(miSize);
    mmPhenotypes.resizeAndPreserve(miSize);
    mmSubjects(miSize - 1) = poSubject;
    mmPhenotypes(miSize - 1) = piPhenotype;
}

int CPedigree::countParents() {
    int i, riParents = 0;
    for (i = 0; i < miSize; i++) {
        if (0 < mmSubjects(i)->countChildren()) {
            riParents++;
        }
    }
    return riParents;
}

CSubject* CPedigree::getParent(const int piParent) {
    int i, liParents = 0;
    for (i = 0; i < miSize; i++) {
        if (0 < mmSubjects(i)->countChildren()) {
            if (piParent == liParents) { return mmSubjects(i); }
            liParents++;
        }
    }
    return NULL;
}

int CPedigree::countAffectedSiblings() {
    int i, riSiblings = 0;
    for (i = 0; i < miSize; i++) {
        if (2 == mmPhenotypes(i) && 0 < mmSubjects(i)->countParents()) {
            riSiblings++;
        }
    }
    return riSiblings;
}

CSubject* CPedigree::getAffectedSibling(const int piSibling) {
    int i, liSiblings = 0;
    for (i = 0; i < miSize; i++) {
        if (2 == mmPhenotypes(i) && 0 < mmSubjects(i)->countParents()) {
            if (piSibling == liSiblings) { return mmSubjects(i); }
            liSiblings++;
        }
    }
    return NULL;
}

CPedigree& CPedigree::operator=(const CPedigree& poPedigree) {
    miPedigree = poPedigree.miPedigree;
    miSize     = poPedigree.miSize;
    mmSubjects.resize(miSize);
    mmPhenotypes.resize(miSize);
    mmSubjects   = poPedigree.mmSubjects;
    mmPhenotypes = poPedigree.mmPhenotypes;
    return *this;
}
