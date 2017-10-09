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
#ifndef SUBJECT_H
#define SUBJECT_H

#include "types.h"

class CSubject;

typedef    blitz::Array<CSubject*,1>   TSubjectPVector;
typedef    blitz::Array<CSubject,1>    TSubjectVector;
typedef    blitz::Array<CSubject,2>    TSubjectMatrix;

class CSubject {
    public:
    //  Utility
                    CSubject            ();
                    CSubject            (string psFamily, string psId, string psParent1, string psParent2, int piGender);
    void            initialise          (string psFamily, string psId, string psParent1, string psParent2, int piGender);
    void            clear               ();
    bool            operator==          (CSubject poSubject);
    CSubject        operator=           (CSubject poSubject);
    friend istream& operator>>          (istream& in, CSubject& poSubject);
    friend ostream& operator<<          (ostream& out, const CSubject& poSubject);
    //  Getters and setters.
    void            setHash             ();
    void            setPedigreeNumber   (const int piPedigree)    { miPedigree = piPedigree; }
    void            setIdNumber         (const int piId)          { miId       = piId; }
    void            setPedigree         (const string psPedigree) { msPedigree = psPedigree; setHash(); }
    void            setId               (const string psId)       { msId       = psId;       setHash(); }
    void            setParent1          (const string psParent1)  { msParent1  = psParent1; }
    void            setParent2          (const string psParent2)  { msParent2  = psParent2; }
    void            setGender           (const int piGender)      { miGender   = piGender; }
    size_t          getHash             () { return miHash; }
    size_t          getHashPedigree     () { return miHashPedigree; }
    int             getPedigreeNumber   () { return miPedigree; }
    int             getIdNumber         () { return miId; }
    string          getPedigree         () { return msPedigree; }
    string          getId               () { return msId; }
    string          getParent1          () { return msParent1; }
    string          getParent2          () { return msParent2; }
    int             getGender           () { return miGender; }
    //  Relationship: parents
    void            setParents          (CSubject *poParent1, CSubject *poParent2);
    CSubject*       getFather           ();
    CSubject*       getMother           ();
    int             countParents        ();
    //  Relationships: children
    void            addChild            (CSubject *poChild);
    CSubject*       getChild            (int piChild)       { return maChildren(piChild); }
    int             countChildren       ();
    //  Relationships: sublings
    void            addSibling          (CSubject *poSibling);
    CSubject*       getSibling          (int piSibling)       { return maSiblings(piSibling); }
    int             countSiblings       ()                    { return maSiblings.rows(); }
    private:
    size_t          miHash;
    size_t          miHashPedigree;
    int             miPedigree;
    int             miId;
    string          msPedigree;
    string          msId;
    string          msParent1;
    string          msParent2;
    int             miGender;
    CSubject       *moParent1;
    CSubject       *moParent2;
    TSubjectPVector maChildren;
    TSubjectPVector maSiblings;
};

    int             countCommon         (const TSubjectPVector& pmA, const TSubjectPVector& pmB);
    int             countCommon         (const TUVector& pmA, const TUVector& pmB);

    TSubjectPVector getSubjects         (const TSubjectPVector& pmA, const TUVector& pmB);
    TUVector        getIndices          (const TSubjectPVector& pmA);

    TUVector        getParentCount      (const TSubjectPVector& pmA);
    TUVector        getSiblingCount     (const TSubjectPVector& pmA);
    TUVector        getChildrenCount    (const TSubjectPVector& pmA);

#endif
