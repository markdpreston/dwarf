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
#include "subject.h"
#include "random.h"

CSubject::CSubject () {
    clear();
}

CSubject::CSubject (string psPedigree, string psId, string psParent1, string psParent2, int piGender) {
    clear();
    msPedigree = psPedigree;
    msId       = psId;
    msParent1  = psParent1;
    msParent2  = psParent2;
    miGender   = piGender;
    setHash();
}

void CSubject::initialise (string psPedigree, string psId, string psParent1, string psParent2, int piGender) {
    clear();
    msPedigree = psPedigree;
    msId       = psId;
    msParent1  = psParent1;
    msParent2  = psParent2;
    miGender   = piGender;
    setHash();
}

void CSubject::clear () {
    miPedigree = -1;
    miId       = -1;
    msPedigree = "";
    msId       = "";
    msParent1  = "";
    msParent2  = "";
    miGender   = 1;
    moParent1  = NULL;
    moParent2  = NULL;
    maChildren.resize(0);
    setHash();
}

void CSubject::setHash () {
    boost::hash<string> loHash;
    miHash = loHash(msPedigree + msId);
    miHashPedigree = loHash(msPedigree);
}

void CSubject::setParents (CSubject *poParent1, CSubject *poParent2) {
    moParent1 = poParent1;
    moParent2 = poParent2;
}

CSubject* CSubject::getFather () {
     if (moParent1 && moParent1->miGender == 1) { return moParent1; }
     if (moParent2 && moParent2->miGender == 1) { return moParent2; }
     return NULL;
}

CSubject* CSubject::getMother () {
     if (moParent1 && moParent1->miGender == 2) { return moParent1; }
     if (moParent2 && moParent2->miGender == 2) { return moParent2; }
     return NULL;
}

int CSubject::countParents () {
     if ( moParent1 &&  moParent2) { return 2; }
     if (!moParent1 && !moParent2) { return 0; }
     return 1;
}

void CSubject::addChild(CSubject *poChild) {
    maChildren.resizeAndPreserve(maChildren.size() + 1);
    maChildren(maChildren.size() - 1) = poChild;
}

int CSubject::countChildren() {
    return maChildren.rows();
}

void CSubject::addSibling(CSubject *poSibling) {
    maSiblings.resizeAndPreserve(maSiblings.size() + 1);
    maSiblings(maSiblings.size() - 1) = poSibling;
}

bool CSubject::operator== (CSubject poSubject) {
    return (miHash == poSubject.miHash) && (msPedigree == poSubject.msPedigree) && (msId == poSubject.msId);
}

CSubject CSubject::operator= (CSubject poSubject) {
    miPedigree = poSubject.miPedigree;
    miId       = poSubject.miId;
    msPedigree = poSubject.msPedigree;
    msId       = poSubject.msId;
    msParent1  = poSubject.msParent1;
    msParent2  = poSubject.msParent2;
    miGender   = poSubject.miGender;
    moParent1  = poSubject.moParent1;
    moParent2  = poSubject.moParent2;
    maChildren.resize(poSubject.maChildren.rows());
    maChildren = poSubject.maChildren;
    maSiblings.resize(poSubject.maSiblings.rows());
    maSiblings = poSubject.maSiblings;
    setHash();
    return *this;
}

istream& operator>> (istream& in, CSubject& poSubject) {
    in >> poSubject.msPedigree >> poSubject.msId >> poSubject.msParent1 >> poSubject.msParent2 >> poSubject.miGender;
    poSubject.setHash();
    return in;
}

ostream& operator<< (ostream& out, const CSubject& poSubject) {
    out << left;
    out << setw(12) << poSubject.msPedigree;
    out << setw(10) << poSubject.msId;
    out << setw(10) << poSubject.msParent1;
    out << setw(10) << poSubject.msParent2;
    out << setw( 4) << poSubject.miGender;
    out << right;
    return out;
}
