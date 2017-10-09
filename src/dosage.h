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
#ifndef DOSAGE_H
#define DOSAGE_H

#include "types.h"
#include "pedigree.h"

typedef struct {
    TReal           mrMajor;
    TReal           mrHetero;
    TReal           mrMinor;
} TDosage;

typedef blitz::Array<TDosage,1> TDVector;
typedef blitz::Array<TDosage,2> TDMatrix;
const TDosage coDosageNull = {0, 0, 0};

class CDosage {
    public:
                    CDosage             ();
                   ~CDosage             ();
    void            clear               ();
    //  Correct for pedigrees.
    void            pedigreeDosage      (const CPedigree& poPedigree);
    //  Size.
    unsigned int    countSubjects       () { return miSubjects; }
    unsigned int    countSNPs           () { return miSNPs; }
    unsigned int    resize              (const unsigned int piSubjects, const unsigned int piSNPs);
    //  Setters.
    void            setSubject          (const unsigned int piSubject, const TDVector& pmValues);
    void            setSNP              (const unsigned int piSNP, const TDVector& pmValues);
    void            set                 (const unsigned int piSubject, const unsigned int piSNP, const TDosage& prValue);
    void            set                 (const unsigned int piSubject, const TUVector& pmSNPs, const TDVector& pmValues);
    void            set                 (const TUVector& pmSubjects, const unsigned int piSNP, const TDVector& pmValues);
    void            set                 (const TUVector& pmSubjects, const TUVector& pmSNPs, const TDMatrix& pmValues);
    //  Getters.
    TRMatrix        getDosage           () const;
    TDVector        getSubject          (const unsigned int piSubject) const;
    TDVector        getSNP              (const unsigned int piSNP) const;
    TReal           getDosage           (const unsigned int piSubject, const unsigned int piSNP) const;
    TDosage         get                 (const unsigned int piSubject, const unsigned int piSNP) const;
    TDVector        get                 (const unsigned int piSubject, const TUVector& pmSNPs) const;
    TDVector        get                 (const TUVector& pmSubjects, const unsigned int piSNP) const;
    TDMatrix        get                 (const TUVector& pmSubjects, const TUVector& pmSNPs) const;
    //  Data.
    void            keepSNPs            (const TUVector& pmSNPs);
    void            keepSubjects        (const TUVector& pmSubjects);
    void            removeSNPs          (const TUVector& pmSNPs);
    void            removeSubjects      (const TUVector& pmSubjects);
    //  IO.
    friend istream& operator>>          (istream& in, CDosage& poDosage);
    friend ostream& operator<<          (ostream& out, const CDosage& poDosage);
//    private:
    unsigned int    miSubjects;
    unsigned int    miSNPs;
    TUVector        mmAllSubjects;
    TUVector        mmAllSNPs;
    TRMatrix        mmMajorData;
    TRMatrix        mmMinorData;
    static TReal    convertToDosage     (const TDosage& poD);
    static TRVector convertToDosage     (const TDVector& pmD);
    static TRVector convertToDosage     (const TSNPVector& pmS, const bool pbError = false);
    static TRVector convertToDosage     (const THVector& pm0, const THVector& pm1);
    static TIVector convertToIntegerDosage (const TRVector& pmD);
};
#endif
