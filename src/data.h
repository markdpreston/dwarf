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
#ifndef DATA_H
#define DATA_H

#include "types.h"

class CData {
    public:
                    CData         () { }
                   ~CData         () { }
    virtual void    clear         ();
    //  Size.
    unsigned int    countSubjects       () { return miSubjects; }
    unsigned int    countSNPs           () { return miSNPs; }
    virtual unsigned int    resize              (const unsigned int piSubjects, const unsigned int piSNPs) = 0;
    virtual unsigned int    resizeAndPreserve   (const unsigned int piSubjects, const unsigned int piSNPs) = 0;
    //  Setters.
    virtual void    setSubject          (const unsigned int piSubject, const TSNPVector& pmCodes) = 0;
    virtual void    setSNP              (const unsigned int piSNP, const TSNPVector& pmCodes) = 0;
    virtual void    set                 (const unsigned int piSubject, const unsigned int piSNP, ESNPCode piCode) = 0;
    virtual void    set                 (const unsigned int piSubject, const TUVector& pmSNPs, const TSNPVector& pmCodes) = 0;
    virtual void    set                 (const TUVector& pmSubjects, const unsigned int piSNP, const TSNPVector& pmCodes) = 0;
    virtual void    set                 (const TUVector& pmSubjects, const TUVector& pmSNPs, const TSNPMatrix& pmCodes) = 0;
    //  Getters.
    virtual TSNPVector      getSubject          (const unsigned int piSubject) const = 0;
    virtual TSNPVector      getSNP              (const unsigned int piSNP) const = 0;
    virtual ESNPCode        get                 (const unsigned int piSubject, const unsigned int piSNP) const = 0;
    virtual TSNPVector      get                 (const unsigned int piSubject, const TUVector& pmSNPs) const = 0;
    virtual TSNPVector      get                 (const TUVector& pmSubjects, const unsigned int piSNP) const = 0;
    virtual TSNPMatrix      get                 (const TUVector& pmSubjects, const TUVector& pmSNPs) const = 0;
    virtual TRMatrix        getDosage           (const bool pbHaplotype = true) const = 0;
    virtual TReal           getDosage           (const unsigned int piSubject, const unsigned int piSNP, const bool pbHaplotype = true) const = 0;
//    static virtual TReal    getDosage           (const ESNPCode peSNP, const bool pbHaplotype = true) = 0;
    //  Data.
    virtual void    keepSubjects        (const TUVector& pmSubject) = 0;
    virtual void    keepSNPs            (const TUVector& pmSNPs) = 0;
    virtual void    removeSubjects      (const TUVector& pmSubject) = 0;
    virtual void    removeSNPs          (const TUVector& pmSNPs) = 0;
    //  Assignment.
//    CData&    operator=           (const CData& poData) = 0;
    //  Utility.
    virtual void    transform           (const unsigned int piRecord, const unsigned char* paData, const unsigned int piBytes, const bool pbSNP) = 0;
    virtual void    untransform         (const unsigned int piRecord, unsigned char* paData, const unsigned int piBytes, const bool pbSNP, const bool pbHaploid) const = 0;

    unsigned int    miSize;

    protected:
    virtual void    internalError       (const int piError) const = 0;
    unsigned int    miSubjects;
    unsigned int    miSNPs;
    unsigned char  *maData;
};

#endif
