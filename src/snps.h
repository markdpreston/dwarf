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
#ifndef SNPS_H
#define SNPS_H

#include "types.h"
#include "gene.h"
#include "snp.h"
#include "snpdata.h"
#include "binarydata.h"
#include "haplotype.h"
#include "random.h"

class CSNPs {
    public:
                    CSNPs               () { clear(); }
                   ~CSNPs               () { clear(); }
    void            clear               ();
    unsigned int    countSNPs           () const { return miSNPs; }
    void            resize              (const unsigned int piSNPs);
    void            resizeAndPreserve   (const unsigned int piSNPs);
    void            createMAFs          (const string psMAFDistribution, const TReal prMAFBound);
    void            setAffecting        (const unsigned int piAffectedCount, const TReal prOddsRatio);
    void            setBaseLine         (const TReal prBaseLine);
    void            setBounds           (const TReal prControlMin, const TReal prControlMax, const TReal prCaseMin, const TReal prCaseMax);

    TReal           getBaseLine         () const { return mrBaseLine; }
    CSnip           getSNP              (const unsigned int piSNP) const { return mmSNPs((int) piSNP); }
    void            setSNP              (const unsigned int piSNP, const CSnip& poSNP) { mmSNPs((int) piSNP) = poSNP; }
    void            copySNPs            (const CSNPs* poSNPs, const int piSNPs = -1);
    void            addSNPs             (const CSNPs* poSNPs);

    TUVector        getGene             (const CGene& poGene) const;

    //  Creation
    THVector        createHaplotype     () const;
    void            createOffspring     (const THVector pmParents[], THVector pmOffspring[]) const;
    void            untransmitted       (const THVector pmParents[], const THVector pmOffspring[], THVector pmUntransmitted[]) const;
    void            untransmitted       (const TSNPVector& pmParent1, const TSNPVector& pmParent2, const TSNPVector& pmOffspring, TSNPVector& pmUntransmitted) const;
    TReal           trait               (const THVector& pmHaplotype0, const THVector& pmHaplotype1);
    TReal           trait               (const TSNPVector& pmSubject) const;
    bool            affected            (const TReal prTrait) const;

    private:
    unsigned int    miSNPs;
    TSnipVector     mmSNPs;
    TSNPVector      mmTempSNPs;
    TReal           mrCaseMin;
    TReal           mrCaseMax;
    TReal           mrControlMin;
    TReal           mrControlMax;
    TReal           mrBaseLine;
    TReal           mrLogBaseLine;
};

inline THVector CSNPs::createHaplotype() const {
    unsigned int    i;
    THVector        rmHaplotype;
    rmHaplotype.resize(miSNPs);
    rmHaplotype = false;
    for (i = 0; i < miSNPs; i++) {
        if (CRandom::uniform(mmSNPs(i).mrMAF)) { rmHaplotype(i) = true; }
    }
    return rmHaplotype;
}

inline void CSNPs::createOffspring (const THVector pmParents[], THVector pmOffspring[]) const {
    unsigned int    i, liStrand0, liStrand1, liMask;
    static boost::uniform_smallint<> soBoostSmallInt(0,pow(2,miSNPs)-1);
    static boost::variate_generator<boost::mt19937&,boost::uniform_smallint<> > soBoostOffspring(CRandom::soBoostRandom,soBoostSmallInt);
    liStrand0 = soBoostOffspring();
    liStrand1 = soBoostOffspring();
    liMask = 1;
    for (i = 0; i < miSNPs; i++) {
        pmOffspring[0](i) = (liStrand0 & liMask) ? pmParents[0](i) : pmParents[1](i);
        pmOffspring[1](i) = (liStrand1 & liMask) ? pmParents[2](i) : pmParents[3](i);
        liMask <<= 1;
    }
}

inline TReal CSNPs::trait (const THVector& pmHaplotype0, const THVector& pmHaplotype1) {
    mmTempSNPs = CHaplotype::transform(pmHaplotype0, pmHaplotype1);
    return trait(mmTempSNPs);
}

inline TReal CSNPs::trait (const TSNPVector& pmSubject) const {
    unsigned int    i;
    TReal           lrLogOdds = mrLogBaseLine;
    for (i = 0; i < miSNPs; i++) {
        lrLogOdds += mmSNPs(i).mrLogOddsRatio * CBinaryData::getDosage(pmSubject(i),true);
    }
//    for (i = 0; i < miSNPs; i++) {
//        switch (mmSNPs(i).msEffect[0]) {
//            case 'R': lrLogOdds += mmSNPs(i).mrLogOddsRatio * (CBinaryData::getDosage(pmSubject(i),true) >= 0.999); break;
//            case 'D': lrLogOdds += mmSNPs(i).mrLogOddsRatio * (CBinaryData::getDosage(pmSubject(i),true) >= 1.999); break;
//            case 'A': lrLogOdds += mmSNPs(i).mrLogOddsRatio * CBinaryData::getDosage(pmSubject(i),true); break;
//            case 'M': lrLogOdds += mmSNPs(i).mrLogOddsRatio * CBinaryData::getDosage(pmSubject(i),true); break;
//        }
//    }
    return CRandom::uniform(1 / (1 + exp(-lrLogOdds)));
}

inline bool CSNPs::affected (const TReal prTrait) const {
    return ((mrCaseMin + crEpsilon) <= prTrait) && (prTrait <= mrCaseMax);
}

#endif
