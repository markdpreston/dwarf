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
#ifndef BINARYDATA_H
#define BINARYDATA_H

#include "types.h"
#include "data.h"

class CBinaryData {
    public:
                    CBinaryData         ();
                   ~CBinaryData         ();
    void            clear               ();
    //  Size.
    unsigned int    countSubjects       () { return miSubjects; }
    unsigned int    countSNPs           () { return miSNPs; }
    unsigned int    resize              (const unsigned int piSubjects, const unsigned int piSNPs);
    unsigned int    resizeAndPreserve   (const unsigned int piSubjects, const unsigned int piSNPs);
    //  Setters.
    void            setSubject          (const unsigned int piSubject, const TSNPVector& pmCodes);
    void            setSNP              (const unsigned int piSNP, const TSNPVector& pmCodes);
    void            set                 (const unsigned int piSubject, const unsigned int piSNP, ESNPCode piCode);
    void            set                 (const unsigned int piSubject, const TUVector& pmSNPs, const TSNPVector& pmCodes);
    void            set                 (const TUVector& pmSubjects, const unsigned int piSNP, const TSNPVector& pmCodes);
    void            set                 (const TUVector& pmSubjects, const TUVector& pmSNPs, const TSNPMatrix& pmCodes);
    //  Getters.
    TSNPVector      getSubject          (const unsigned int piSubject) const;
    TSNPVector      getSNP              (const unsigned int piSNP) const;
    ESNPCode        get                 (const unsigned int piSubject, const unsigned int piSNP) const;
    TSNPVector      get                 (const unsigned int piSubject, const TUVector& pmSNPs) const;
    TSNPVector      get                 (const TUVector& pmSubjects, const unsigned int piSNP) const;
    TSNPMatrix      get                 (const TUVector& pmSubjects, const TUVector& pmSNPs) const;
    TRMatrix        getDosage           (const bool pbHaplotype = true) const;
    TReal           getDosage           (const unsigned int piSubject, const unsigned int piSNP, const bool pbHaplotype = true) const;
    static TReal    getDosage           (const ESNPCode peSNP, const bool pbHaplotype = true);
    //  Data.
    void            keepSubjects        (const TUVector& pmSubject);
    void            keepSNPs            (const TUVector& pmSNPs);
    void            removeSubjects      (const TUVector& pmSubject);
    void            removeSNPs          (const TUVector& pmSNPs);
    //  Assignment.
    CBinaryData&    operator=           (const CBinaryData& poBinaryData);
    //  Utility.
    void            transform           (const unsigned int piRecord, const unsigned char* paData, const unsigned int piBytes, const bool pbSNP);
    void            untransform         (const unsigned int piRecord, unsigned char* paData, const unsigned int piBytes, const bool pbSNP, const bool pbHaploid) const;



    unsigned int    miSize;

    friend class CPopulation;

    protected:
    void            internalError       (const int piError) const;
    unsigned char   maMasks[4];
    ESNPCode        maReverse[4];
    unsigned int    miSubjects;
    unsigned int    miSNPs;
    TUVector        mmAllSubjects;
    TUVector        mmAllSNPs;
    unsigned char  *maData;
};

inline ESNPCode CBinaryData::get(const unsigned int piSubject, const unsigned int piSNP) const {
    unsigned int    liByte, liOffset;
    unsigned char  *lcChar, rcValue;

    liByte   = (piSubject * miSNPs + piSNP) / 4;
    liOffset = (piSubject * miSNPs + piSNP) % 4;

    if (piSubject >= miSubjects) { internalError(1); }
    if (piSNP     >= miSNPs)     { internalError(2); }
    if (liByte    >= miSize)     { internalError(3); }

    lcChar = maData;
    rcValue = (lcChar[liByte] & maMasks[liOffset]) >> (6 - 2 * liOffset);

    return (ESNPCode) rcValue;
}

#endif
