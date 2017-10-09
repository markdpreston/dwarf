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
#ifndef POLYDATA_H
#define POLYDATA_H

#include "types.h"
#include "snp.h"

class CPolyData {
    public:
                    CPolyData           ();
                   ~CPolyData           ();
    void            clear               ();
    //  Size.
    unsigned int    countSubjects       () { return miSubjects; }
    unsigned int    countSNPs           () { return miSNPs; }
    unsigned int    resize              (const unsigned int piSubjects, const unsigned int piSNPs);
    //  Setters.
    void            setSubject          (const unsigned int piSubject, const TUVector& pmValues);
    void            setSNP              (const unsigned int piSNP, const TUVector& pmValues);
    void            set                 (const unsigned int piSubject, const unsigned int piSNP, const unsigned char pcPolyType, unsigned char pcQuadData);
    void            set                 (const unsigned int piSubject, const unsigned int piSNP, const unsigned int prValue);
    void            set                 (const unsigned int piSubject, const TUVector& pmSNPs, const TUVector& pmValues);
    void            set                 (const TUVector& pmSubjects, const unsigned int piSNP, const TUVector& pmValues);
    void            set                 (const TUVector& pmSubjects, const TUVector& pmSNPs, const TUMatrix& pmValues);
    //  Getters.
    TUVector        getSubject          (const unsigned int piSubject) const;
    TUVector        getSNP              (const unsigned int piSNP) const;
    ESNPCode        getSNPCode          (const unsigned int piSubject, const unsigned int piSNP, const CSnip poSNP) const;
    unsigned int    get                 (const unsigned int piSubject, const unsigned int piSNP) const;
    TUVector        get                 (const unsigned int piSubject, const TUVector& pmSNPs) const;
    TUVector        get                 (const TUVector& pmSubjects, const unsigned int piSNP) const;
    TUMatrix        get                 (const TUVector& pmSubjects, const TUVector& pmSNPs) const;
    //  Data.
    void            keepSNPs            (const TUVector& pmSNPs);
    void            keepSubjects        (const TUVector& pmSubjects);
    void            removeSNPs          (const TUVector& pmSNPs);
    void            removeSubjects      (const TUVector& pmSubjects);
    //  Utility methods.
    static string           output      (const unsigned int pcA, bool pbSpace = false);
    static unsigned int     transform   (const string pcA, const string pcB);
//    static unsigned char    transform   (const unsigned char pcA, const unsigned char pcB);

    static const unsigned int      ciSNPOptions  = 12;
    static const unsigned int      ciQuadOptions = 10;
    static unsigned char    cmQuadOption    (unsigned int piIndex);
    static unsigned char    getPolyType     (const unsigned int piData);
    static unsigned int     getIndelLength  (const unsigned int piData);
    static unsigned char    getQuadData     (const unsigned int piData, const unsigned int piAllele = 0);
    static string           quadOutput      (const unsigned char pcA, bool pbSpace = false);
    static unsigned char    quadTransform   (const string pcA, const string pcB);
    static unsigned char    quadTransform   (const unsigned char pcA, const unsigned char pcB);

    private:
    unsigned int    miSubjects;
    unsigned int    miSNPs;
    TUVector        mmAllSubjects;
    TUVector        mmAllSNPs;
    TUMatrix        mmData;
};

#endif
