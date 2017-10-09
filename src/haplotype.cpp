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
#include "haplotype.h"

ESNPCode CHaplotype::get(const unsigned int piSubject, const unsigned int piSNP) {
    return CBinaryData::get(piSubject, piSNP);
}

int CHaplotype::get(const unsigned int piSubject, const unsigned int piSNP, const unsigned int piStrand) {
    unsigned int    liStrand;
    ESNPCode        leSNP;
    leSNP = CBinaryData::get(piSubject,piSNP);
    if (0 == piStrand) { liStrand = 1; } else { liStrand = 2; }
    return (leSNP & liStrand) ? 1 : 0;
}

void CHaplotype::set (const unsigned int piSubject, const unsigned int piSNP, const int piAllele1, const int piAllele2) {
    CBinaryData::set(piSubject, piSNP, (ESNPCode) (2*piAllele1 + piAllele2));
}

void CHaplotype::setSubject (const unsigned int piSubject, const THVector& pmCodes0, const THVector& pmCodes1) {
    setSubject(piSubject, 0, pmCodes0);
    setSubject(piSubject, 1, pmCodes1);
}

void CHaplotype::setSubject (const unsigned int piSubject, const unsigned int piStrand, const THVector& pmCodes) {
    unsigned int    i, liStrand;
    TSNPVector      lmSNPs;
    if (miSNPs != (unsigned int) pmCodes.rows()) { cerr << "Bad haplotype length" << endl; }
    lmSNPs.resize(miSNPs);
    lmSNPs = CBinaryData::getSubject(piSubject);
    //  Mask the opposite strand.
    if (0 == piStrand) { liStrand = 0x2; } else { liStrand = 0x1; }
    if (piStrand > 1) { cerr << "Strand too large: " << piStrand << endl; exit(0); }
    for (i = 0; i < miSNPs; i++) {
        CBinaryData::set(piSubject, i, (ESNPCode) ((lmSNPs(i) & liStrand) + (pmCodes(i) << piStrand)));
    }
}

void CHaplotype::setHaplotype (const unsigned int piSNP, const unsigned int piStrand, const THVector& pmCodes) {
    unsigned int    i, liStrand;
    TSNPVector      lmSNPs;
    if (miSubjects != (unsigned int) pmCodes.rows()) { cerr << "Bad haplotype length" << endl; }
    lmSNPs.resize(miSubjects);
    lmSNPs = getSNP(piSNP);
    if (0 == piStrand) { liStrand = 1; } else { liStrand = 2; }
    for (i = 0; i < miSubjects; i++) {
        CBinaryData::set(i, piSNP, (ESNPCode) ((lmSNPs(i) & liStrand) + (pmCodes(i) << piStrand)));
    }
}

THVector CHaplotype::getSubject (const unsigned int piSubject, const unsigned int piStrand) const {
    unsigned int    i, liStrand;
    TSNPVector      lmSNPs;
    THVector        rmHaplotype;
    rmHaplotype.resize(miSNPs);
    lmSNPs.resize(miSNPs);
    lmSNPs = CBinaryData::getSubject(piSubject);
    if (0 == piStrand) { liStrand = 1; } else { liStrand = 2; }
    for (i = 0; i < miSNPs; i++) {
        rmHaplotype(i) = (lmSNPs(i) & liStrand) ? 1 : 0;
    }
    return rmHaplotype;
}

THVector CHaplotype::getHaplotype (const unsigned int piSNP, const unsigned int piStrand) const {
    unsigned int    i, liStrand;
    TSNPVector      lmSNP;
    THVector        rmHaplotype;
    rmHaplotype.resize(miSubjects);
    lmSNP.resize(miSubjects);
    lmSNP = getSNP(piSNP);
    if (0 == piStrand) { liStrand = 1; } else { liStrand = 2; }
    for (i = 0; i < miSubjects; i++) {
        rmHaplotype(i) = (lmSNP(i) & liStrand) ? 1 : 0;
    }
    return rmHaplotype;
}

ESNPCode CHaplotype::transform (const bool pbHaplotype1, const bool pbHaplotype2) {
    if (pbHaplotype1 != pbHaplotype2) {
        return ceSNPHetero;
    } else if (false == pbHaplotype1) {
        return ceSNPMajor;
    } else {
        return ceSNPMinor;
    }
}

TSNPVector CHaplotype::transform (const TSNPVector& pmSNPs) {
    unsigned int    i, liSize;
    TSNPVector      rmGenotype;
    liSize = pmSNPs.rows();
    rmGenotype.resize(liSize);
    rmGenotype = ceSNPError;
    for (i = 0; i < liSize; i++) {
        rmGenotype(i) = transform(pmSNPs(i));
    }
    return rmGenotype;
}

TSNPVector CHaplotype::transform (const THVector& pmHaplotype1, const THVector& pmHaplotype2) {
    unsigned int    i, liSize;
    TSNPVector      rmGenotype;
    liSize = pmHaplotype1.rows();
    rmGenotype.resize(liSize);
    rmGenotype = ceSNPError;
    for (i = 0; i < liSize; i++) {
        rmGenotype(i) = transform(pmHaplotype1(i), pmHaplotype2(i));
    }
    return rmGenotype;
}

THVector CHaplotype::split (const TSNPVector& pmSNPs, const int piStrand) {
    unsigned int    i, liSize;
    THVector        rmHaplotype;
    liSize = pmSNPs.rows();
    rmHaplotype.resize(liSize);
    for (i = 0; i < liSize; i++) {
        rmHaplotype(i) = (ceHaplotype11 == pmSNPs(i) || ceHaplotype10 == pmSNPs(i)) && 1 == piStrand;
    }
    return rmHaplotype;
}

THVector CHaplotype::concatenate (const THVector& pmHaplotype1, const THVector& pmHaplotype2) {
    unsigned int    i, liSize;
    THVector        rmHaplotype;
    liSize = pmHaplotype1.rows();
    rmHaplotype.resize(2*liSize);
    for (i = 0; i < liSize; i++) {
        rmHaplotype(i)        = pmHaplotype1(i);
        rmHaplotype(liSize+i) = pmHaplotype2(i);
    }
    return rmHaplotype;
}
