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
#include "snps.h"
#include "random.h"
#include "haplotype.h"

void CSNPs::clear() {
    miSNPs        = 0;
    mmSNPs.resize(0);
    mmTempSNPs.resize(0);
    mrBaseLine    = 0;
    mrLogBaseLine = 0;
    mrControlMin  = -crInfinity;
    mrControlMax  = 0.0;
    mrCaseMin     = 0.0;
    mrCaseMax     = crInfinity;
}

void CSNPs::resize(const unsigned int piSNPs) {
    unsigned int    i;
    string          lsNumber;
    miSNPs = piSNPs;
    mmSNPs.resize(piSNPs);
    mmTempSNPs.resize(piSNPs);
    mmTempSNPs = ceSNPMajor;
    for (i = 0; i < miSNPs; i++) {
        lsNumber = CUtility::intToString(i+1);
        mmSNPs(i).initialise("S", "rs"+lsNumber, "dis"+lsNumber, i+1);
    }
//    mmTempSNPs
//    mmTempDosage
}

void CSNPs::resizeAndPreserve(const unsigned int piSNPs) {
    miSNPs = piSNPs;
    mmSNPs.resizeAndPreserve(piSNPs);
}

void CSNPs::createMAFs (const string psMAFDistribution, const TReal prMAFBound) {
    unsigned int    i;
    TReal           lrMAF;
    if ("existing" == psMAFDistribution) {
    } else if ("uniform" == psMAFDistribution) {
        for (i = 0; i < miSNPs; i++) {
            mmSNPs(i).mrMAF = prMAFBound * CRandom::uniform();
        }
    } else if ("beta" == psMAFDistribution) {
        static TGammaD soGammaDX(0.345), soGammaDY(1.058);
        static TGammaG soGammaX(CRandom::soBoostRandom,soGammaDX), soGammaY(CRandom::soBoostRandom,soGammaDY);
        for (i = 0; i < miSNPs; i++) {
            do {
                lrMAF = 0.5 * CRandom::beta(soGammaX, soGammaY);
            } while (lrMAF > prMAFBound);
            mmSNPs(i).mrMAF = lrMAF;
        }
    } else {
        cerr << "Bad MAF distribution: " << psMAFDistribution << endl;
        exit(0);
    }
}

void CSNPs::setAffecting (const unsigned int piAffectingCount, const TReal prOddsRatio) {
    unsigned int    i, liAffecting;
    liAffecting = piAffectingCount;
    if (liAffecting > miSNPs) {
        cerr << "Too many affecting SNPs (" << liAffecting << "), limiting to number of SNPs (" << miSNPs << ")" << endl;
        liAffecting = miSNPs;
    }
    for (i = 0; i < liAffecting; i++) {
        mmSNPs(i).msEffect       = "M";
        mmSNPs(i).mrOddsRatio    = prOddsRatio;
        mmSNPs(i).mrLogOddsRatio = log(prOddsRatio);
    }
    for (i = piAffectingCount; i < miSNPs; i++) {
        mmSNPs(i).msEffect       = "0";
        mmSNPs(i).mrOddsRatio    = 1.0;
        mmSNPs(i).mrLogOddsRatio = 0.0;
    }
}

void CSNPs::setBounds(const TReal prControlMin, const TReal prControlMax, const TReal prCaseMin, const TReal prCaseMax) {
    mrControlMin = prControlMin;
    mrControlMax = prControlMax;
    mrCaseMin    = prCaseMin;
    mrCaseMax    = prCaseMax;
}

void CSNPs::setBaseLine(const TReal prBaseLine) {
    mrBaseLine    = prBaseLine;
    mrLogBaseLine = log(prBaseLine / (1 - prBaseLine));
}

void CSNPs::copySNPs(const CSNPs* poSNPs, const int piSNPs) {
    if (-1 == piSNPs || piSNPs > (int) poSNPs->miSNPs) {
        miSNPs = poSNPs->miSNPs;
        mmSNPs.resize(miSNPs);
        mmSNPs = poSNPs->mmSNPs;
    } else {
        miSNPs = piSNPs;
        mmSNPs.resize(miSNPs);
        mmSNPs = poSNPs->mmSNPs(blitz::Range(0,piSNPs-1));
    }
    mrBaseLine    = poSNPs->mrBaseLine;
    mrLogBaseLine = poSNPs->mrLogBaseLine;
    mrControlMin  = poSNPs->mrControlMin;
    mrControlMax  = poSNPs->mrControlMax;
    mrCaseMin     = poSNPs->mrCaseMin;
    mrCaseMax     = poSNPs->mrCaseMax;
}

void CSNPs::addSNPs(const CSNPs* poSNPs) {
    mmSNPs.resizeAndPreserve(miSNPs + poSNPs->miSNPs);
    mmSNPs(blitz::Range(miSNPs,blitz::toEnd)) = poSNPs->mmSNPs;
    miSNPs = miSNPs + poSNPs->miSNPs;
}

TUVector CSNPs::getGene(const CGene& poGene) const {
    unsigned int  i, liCount, liPositionStart, liPositionEnd;
    string        lsChromosome;
    TUVector      rmSNPs;
    lsChromosome    = poGene.getChromosome();
    liPositionStart = poGene.getPositionStart();
    liPositionEnd   = poGene.getPositionEnd();
    liCount = 0;
    rmSNPs.resize(0);
    for (i = 0; i < miSNPs; i++) {
        if (lsChromosome    == mmSNPs((int)i).msChromosome &&
            liPositionStart <= mmSNPs((int)i).miPosition &&
            liPositionEnd   >= mmSNPs((int)i).miPosition) {
            rmSNPs.resizeAndPreserve(liCount+1);
            rmSNPs(liCount) = i;
            liCount++;
        }
    }
    return rmSNPs;
}

void CSNPs::untransmitted (const THVector pmParents[], const THVector pmOffspring[], THVector pmUntransmitted[]) const {
    unsigned int    i;
    int             liMinors;
    pmUntransmitted[0] = false;
    pmUntransmitted[1] = false;
    for (i = 0; i < miSNPs; i++) {
        liMinors = (pmParents[0](i) ? 1 : 0) + (pmParents[1](i) ? 1 : 0) + (pmParents[2](i) ? 1 : 0) + (pmParents[3](i) ? 1 : 0) - (pmOffspring[0](i) ? 1 : 0) - (pmOffspring[1](i) ? 1 : 0);
        switch (liMinors) {
            case 0: break;
            case 1: pmUntransmitted[CRandom::binary() ? 1 : 0](i) = true; break;
            case 2: pmUntransmitted[0](i) = true; pmUntransmitted[1](i) = true; break;
            default: cerr << "Bad haplotypes: " << pmParents[0](i) << " " << pmParents[1](i) << " " << pmParents[2](i) << " " << pmParents[3](i) << " " << pmOffspring[0](i) << " " << pmOffspring[1](i) << " " << liMinors << " " << i << endl; exit(0);
        }
    }
}

void CSNPs::untransmitted (const TSNPVector& pmParent1, const TSNPVector& pmParent2, const TSNPVector& pmOffspring, TSNPVector& pmUntransmitted) const {
    unsigned int    i;
    int             liMinors;
    pmUntransmitted = ceSNPMajor;
    for (i = 0; i < miSNPs; i++) {
        liMinors = pmParent1(i) + pmParent2(i) - pmOffspring(i);
        switch (liMinors) {
            case 0: break;
            case 1: pmUntransmitted(i) = ceSNPHetero; break;
            case 2: pmUntransmitted(i) = ceSNPMinor; break;
            default: cerr << "Bad haplotypes: " << pmParent1(i) << " " << pmParent2(i) << " " << pmOffspring(i) << " " << liMinors << " " << i << endl; exit(0);
        }
    }
}
