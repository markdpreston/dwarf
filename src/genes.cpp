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
#include "genes.h"

void CGenes::clear() {
    miGenes = 0;
    mmGenes.resize(0);
}

void CGenes::resize(const unsigned int piGenes) {
    miGenes = piGenes;
    mmGenes.resize(piGenes);
}

void CGenes::resizeAndPreserve(const unsigned int piGenes) {
    miGenes = piGenes;
    mmGenes.resizeAndPreserve(piGenes);
}

void CGenes::copyGenes (const CGenes* poGenes) {
    miGenes = poGenes->miGenes;
    mmGenes.resize(miGenes);
    mmGenes = poGenes->mmGenes;
}

void CGenes::addGenes (const CGenes* poGenes) {
    mmGenes.resizeAndPreserve(miGenes + poGenes->miGenes);
    mmGenes(blitz::Range(miGenes,blitz::toEnd)) = poGenes->mmGenes;
    miGenes = miGenes + poGenes->miGenes;
}
