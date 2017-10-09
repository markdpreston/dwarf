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
#ifndef HAPLOTYPE_H
#define HAPLOTYPE_H

#include "types.h"
#include "binarydata.h"

class CHaplotype : public CBinaryData {
    public:
    ESNPCode        get                 (const unsigned int piSubject, const unsigned int piSNP);
    int             get                 (const unsigned int piSubject, const unsigned int piSNP, const unsigned int piStrand);
    void            set                 (const unsigned int piSubject, const unsigned int piSNP, const int piHaplotype0, const int piHaplotype1);
    using CBinaryData::getSubject;
    using CBinaryData::setSubject;
    void            setSubject          (const unsigned int piSubject, const THVector& pmCodes0, const THVector& pmCodes1);
    void            setSubject          (const unsigned int piSubject, const unsigned int piStrand, const THVector& pmCodes);
    void            setHaplotype        (const unsigned int piSNP, const unsigned int piStrand, const THVector& pmCodes);
    THVector        getSubject          (const unsigned int piSubject, const unsigned int piStrand) const;
    THVector        getHaplotype        (const unsigned int piSNP, const unsigned int piStrand) const;
    static ESNPCode   transform         (const ESNPCode piCode) { return (piCode != ceSNPError) ? piCode : ceSNPHetero; }
    static ESNPCode   transform         (const bool pbHaplotype0, const bool pbHaplotype1);
    static TSNPVector transform         (const TSNPVector& pmSNPs);
    static TSNPVector transform         (const THVector& pmHaplotype0, const THVector& pmHaplotype1);
    static THVector split               (const TSNPVector& pmSNPs, const int piStrand);
    static THVector concatenate         (const THVector& pmHaplotype0, const THVector& pmHaplotype1);
};

#endif
