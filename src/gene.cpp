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
#include "gene.h"

CGene::CGene () {
    clear();
}

void CGene::clear () {
    msChromosome    = "0";
    msName          = "0";
    miPositionStart = 0;
    miPositionEnd   = 0;
    msNotes         = "";
}

void CGene::initialise (string psChromosome, string psName, int piPositionStart, int piPositionEnd, string psNotes) {
    msChromosome    = psChromosome;
    msName          = psName;
    miPositionStart = piPositionStart;
    miPositionEnd   = piPositionEnd;
    msNotes         = "";
}

void CGene::input (string psSNPLine) {
    stringstream lsSNPLine;
    lsSNPLine << psSNPLine;
    lsSNPLine >> msChromosome >> msName >> miPositionStart >> miPositionEnd >> msNotes;
}

CGene CGene::operator= (const CGene& poSNP) {
    msChromosome    = poSNP.msChromosome;
    msName          = poSNP.msName;
    miPositionStart = poSNP.miPositionStart;
    miPositionEnd   = poSNP.miPositionEnd;
    msNotes         = poSNP.msNotes;
    return *this;
}
