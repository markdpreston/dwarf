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
#include "snp.h"

CSnip::CSnip () {
    clear();
}

CSnip::CSnip (string psSNPLine) {
    clear();
    input(psSNPLine);
}

CSnip::CSnip (string psChromosome, string psName, string psDistance, int piPosition, string psAllele1, string psAllele2) {
    clear();
    initialise(psChromosome, psName, psDistance, piPosition, psAllele1, psAllele2);
}

void CSnip::clear () {
    msChromosome   = "0";
    msName         = "0";
    msDistance     = "0";
    miPosition     = 0;
    miAllele1      = 1;
    miAllele2      = 2;
    mrMAF          = 0.0;
    msEffect       = "0";
    mrOddsRatio    = 1.0;
    mrLogOddsRatio = 0.0;
    msNotes        = "";
}

void CSnip::initialise (string psChromosome, string psName, string psDistance, int piPosition, string psAllele1, string psAllele2, string psMAF, string psEffect, string psOddsRatio) {
    msChromosome   = psChromosome;
    msName         = psName;
    msDistance     = psDistance;
    miPosition     = piPosition;
    miAllele1      = transform(psAllele1[0]);
    miAllele2      = transform(psAllele2[0]);
    mrMAF          = CUtility::stringToReal(psMAF);
    msEffect       = psEffect;
    mrOddsRatio    = CUtility::stringToReal(psOddsRatio);
    mrLogOddsRatio = log(mrOddsRatio);
    msNotes        = "";
}

void CSnip::input (string psSNPLine) {
    string       lsAllele1, lsAllele2;
    stringstream lsSNPLine;
    lsSNPLine << psSNPLine;
    lsSNPLine >> msChromosome >> msName >> msDistance >> miPosition >> lsAllele1 >> lsAllele2;
    miAllele1    = transform(lsAllele1[0], false);
    miAllele2    = transform(lsAllele2[0], false);
    msNotes      = "";
}

int CSnip::transform (const char pcCode, const bool pbExists) {
    int riCode = 0;
    switch (pcCode) {
    case '0':           riCode = 0; break;
    case '1': case 'A': riCode = 1; break;
    case '2': case 'C': riCode = 2; break;
    case '3': case 'G': riCode = 3; break;
    case '4': case 'T': riCode = 4; break;
    default:
        if (pbExists) {
           cerr << "Bad SNP code: " << pcCode << endl;
           exit(0);
           break;
        }
    }
    return riCode;
}

string CSnip::transform (const int piCode) {
    string rsCode;
    switch (piCode) {
    case 0: rsCode = "0"; break;
    case 1: rsCode = "1"; break;
    case 2: rsCode = "2"; break;
    case 3: rsCode = "3"; break;
    case 4: rsCode = "4"; break;
    }
    return rsCode;
}

int CSnip::getAllele(const unsigned int piAllele) const {
    switch (piAllele) {
        case 1:  return miAllele1;
        case 2:  return miAllele2;
        default: return -1;
    }
}

TRVector CSnip::getRelativeRisk(const TReal prBaseLine) const {
    TReal       lrRR;
    TRVector    rmRR(3);
    lrRR = mrOddsRatio / (1 - mrOddsRatio * (1 - prBaseLine));
    rmRR(0) = 1.0;
    rmRR(1) = 1.0;
    rmRR(2) = 1.0;
    switch (msEffect[0]) {
        case 'R': rmRR(1) = 1.0;  rmRR(2) = lrRR;       break;
        case 'D': rmRR(1) = lrRR; rmRR(2) = lrRR;       break;
        case 'A': rmRR(1) = lrRR; rmRR(2) = 2*lrRR - 1; break;
        case 'M': rmRR(1) = lrRR; rmRR(2) = lrRR*lrRR;  break;
    }
    return rmRR;
}

TRVector CSnip::getHWE() const {
    TRVector    rmHWE(3);
    rmHWE(0) = (1 - mrMAF) * (1 - mrMAF); // 00
    rmHWE(1) =  2 * mrMAF  * (1 - mrMAF); // 10 or 01
    rmHWE(2) =      mrMAF  *      mrMAF;  // 11
    return rmHWE;
}

bool CSnip::operator== (const CSnip& poSNP) const {
    return (msChromosome == msChromosome) && (msName == msName);
}

CSnip CSnip::operator= (const CSnip& poSNP) {
    msChromosome   = poSNP.getChromosome();
    msName         = poSNP.getName();
    msDistance     = poSNP.getDistance();
    miPosition     = poSNP.getPosition();
    miAllele1      = poSNP.getAllele(1);
    miAllele2      = poSNP.getAllele(2);
    mrMAF          = poSNP.getMAF();
    msEffect       = poSNP.getEffect();
    mrOddsRatio    = poSNP.getOddsRatio();
    mrLogOddsRatio = poSNP.getLogOddsRatio();
    msNotes        = poSNP.getNotes();
    return *this;
}

ostream& operator<< (ostream& out, const CSnip& poSNP) {
    out << left;
    out << setw( 4) << poSNP.getChromosome();
    out << setw(20) << poSNP.getName();
    out << setw(20) << poSNP.getDistance();
    out << setw(10) << poSNP.getPosition();
    out << setw( 4) << poSNP.getAllele(1);
    out << setw( 4) << poSNP.getAllele(2);
    out << right;
    return out;
}

