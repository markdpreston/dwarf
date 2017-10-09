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
/*
#include "simulated.h"

void CSimulated::resizeMS (const unsigned int piHaplotypes, const unsigned int piChromosomes) {
    miSubjects = piHaplotypes / 2;
    mmChromosomes.resize(piChromosomes);
    mmChromosomes = 0;
}

void CSimulated::initialiseMS () {
    CHaplotype::resize(miSubjects, getChromosomes());
}

unsigned int CSimulated::getChromosomes () const {
    return sum(mmChromosomes);
}

unsigned int CSimulated::getChromosome (const unsigned int piChromosome) const {
    return mmChromosomes(piChromosome);
}

void CSimulated::setChromosome (const unsigned int piChromosome, const unsigned int piSNPs) {
    mmChromosomes(piChromosome) = piSNPs;
}

TReal CSimulated::getPosition (const unsigned int piChromosome, const unsigned int piSNP) const {
    return mmPositions(piChromosome * piSNP);
}

void CSimulated::setPosition (const unsigned int piChromosome, const unsigned int piSNP, const TReal prPosition) {
    int     liPosition = 0;
    if (0 != piChromosome) {
        liPosition = sum(mmChromosomes(blitz::Range(0,piChromosome-1)));
    }
    liPosition += piSNP;
    mmPositions.resizeAndPreserve(liPosition+1);
    mmPositions(liPosition) = prPosition;
}

void CSimulated::setHaplotype (const unsigned int piChromosome, const unsigned int piSubject, const unsigned int piStrand, const string psHaplotype) {
    int         i, liLength;
    THVector    lmHaplotype;
    liLength = psHaplotype.length();
    lmHaplotype.resize(liLength);
    for (i = 0; i < liLength; i++) {
        lmHaplotype(i) = (psHaplotype[i] == '0') ? 0 : 1;
    }
    CHaplotype::setSubject(piSubject, piStrand, lmHaplotype);
}
*/
