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
#include "polydata.h"

CPolyData::CPolyData () {
    clear();
}

CPolyData::~CPolyData () {
    clear();
}

void CPolyData::clear() {
    resize(0,0);
    miSubjects = 0;
    miSNPs = 0;
}

unsigned int CPolyData::resize (const unsigned int piSubjects, const unsigned int piSNPs) {
    unsigned int        riSize;
    blitz::firstIndex   loIndex1;

    miSubjects = piSubjects;
    miSNPs     = piSNPs;
    riSize     = miSubjects * miSNPs;

    mmData.resize(miSubjects,miSNPs);
    mmData = 0;
    mmAllSubjects.resize(miSubjects);
    mmAllSubjects = loIndex1;
    mmAllSNPs.resize(miSNPs);
    mmAllSNPs = loIndex1;

    return riSize;
}

ESNPCode CPolyData::getSNPCode (const unsigned int piSubject, const unsigned int piSNP, const CSnip poSNP) const {
    return ceSNPMajor;
}

unsigned char CPolyData::cmQuadOption (unsigned int piIndex) {
    blitz::TinyVector<unsigned char,10> rmHaplotypes;
    rmHaplotypes = 128 + 8, 128 + 4, 128 + 2, 128 + 1, 64 + 4, 64 + 2, 64 + 1, 32 + 2, 32 + 1, 16 + 1;
    return rmHaplotypes(piIndex);
}

unsigned char CPolyData::getPolyType (const unsigned int piData) {
    return piData >> 24;
}

unsigned int CPolyData::getIndelLength (const unsigned int piData) {
    if ('I' == getPolyType(piData) || 'D' == getPolyType(piData)) {
        return piData & 0x00FFFFFF;
    }
    return 0;
}

unsigned char CPolyData::getQuadData (const unsigned int piData, const unsigned int piAllele) {
    if ('S' == getPolyType(piData)) {
        switch (piAllele) {
            case 0: return piData & 0x000000FF;
            case 1: return (piData & 0x000000F0) >> 4;
            case 2: return piData & 0x0000000F;
        }
    }
    return 0;
}


string CPolyData::quadOutput(unsigned char pcA, const bool pbSpace) {
    stringstream rsOutput;
    switch ((unsigned char)(pcA & 0xF0)) {
        case 128: rsOutput << "A"; break;
        case 64:  rsOutput << "C"; break;
        case 32:  rsOutput << "G"; break;
        case 16:  rsOutput << "T"; break;
        default:  rsOutput << "0";
    }
    if (pbSpace) { rsOutput << " "; }
    switch ((unsigned char)(pcA & 0x0F)) {
        case 8:  rsOutput << "A"; break;
        case 4:  rsOutput << "C"; break;
        case 2:  rsOutput << "G"; break;
        case 1:  rsOutput << "T"; break;
        default: rsOutput << "0";
    }
    return rsOutput.str();
}

unsigned char CPolyData::quadTransform (const string pcA, const string pcB) {
    unsigned char   lcA = 0, lcB = 0, rcA = 0;
    if ("A" == pcA) { lcA = 8; } else if ("C" == pcA) { lcA = 4; } else if ("G" == pcA) { lcA = 2; } else if ("T" == pcA) { lcA = 1; } else { lcA = 0; }
    if ("A" == pcB) { lcB = 8; } else if ("C" == pcB) { lcB = 4; } else if ("G" == pcB) { lcB = 2; } else if ("T" == pcB) { lcB = 1; } else { lcB = 0; }
//    if (0 == lcA || 0 == lcB) { lcA = 0; lcB = 0; }
    if (lcA >= lcB) {
       rcA = (lcA << 4) + lcB;
    } else {
       rcA = (lcB << 4) + lcA;
    }
    return rcA;
}

unsigned char CPolyData::quadTransform (const unsigned char pcA, const unsigned char pcB) {
    unsigned char   lcA = 0, lcB = 0, rcA = 0;
    switch (pcA) {
        case 'A': lcA = 8; break;
        case 'C': lcA = 4; break;
        case 'G': lcA = 2; break;
        case 'T': lcA = 1; break;
    }
    switch (pcB) {
        case 'A': lcB = 8; break;
        case 'C': lcB = 4; break;
        case 'G': lcB = 2; break;
        case 'T': lcB = 1; break;
    }
//    if (0 == lcA || 0 == lcB) { lcA = 0; lcB = 0; }
    if (lcA >= lcB) {
       rcA = (lcA << 4) + lcB;
    } else {
       rcA = (lcB << 4) + lcA;
    }
    return rcA;
}

string CPolyData::output(unsigned int pcA, const bool pbSpace) {
    stringstream rsOutput;
    switch ((0xFF000000 & pcA) >> 24) {
        case 'D': rsOutput << "D" << (int) (0xFFFFFF & pcA); break;
        case 'I': rsOutput << "I" << (int) (0xFFFFFF & pcA); break;
        case 'S': rsOutput << CPolyData::quadOutput(0xFF & pcA, pbSpace); break;
        default:  cerr << "Bad polydata output" << (pcA >> 24) << endl; exit(0);
    }
    return rsOutput.str();
}

unsigned int CPolyData::transform (const string pcA, const string pcB) {
    unsigned int    riA;
    riA = 0;
    if ('I' == pcA[0] || 'D' == pcA[0]) {
        riA = (pcA[0] << 24) + CUtility::stringToInt(pcB);
    } else {
        riA = ('S' << 24) + CPolyData::quadTransform(pcA, pcB);
    }
    return riA;
}

void CPolyData::set (const unsigned int piSubject, const unsigned int piSNP, const unsigned char pcPolyType, const unsigned char pcQuadData) {
    mmData((int)piSubject,(int)piSNP) = (pcPolyType << 24) + pcQuadData;
}

void CPolyData::set (const unsigned int piSubject, const unsigned int piSNP, const unsigned int piData) {
    if (piSubject >= miSubjects || piSNP >= miSNPs) {
        cout << "polydata::set " << piSNP << "  "<< piSubject << " : " << miSNPs << "  "<< miSubjects << endl;
        exit(1);
    }
    mmData((int)piSubject,(int)piSNP) = piData;
}

unsigned int CPolyData::get (const unsigned int piSubject, const unsigned int piSNP) const {
    if (piSubject >= miSubjects || piSNP >= miSNPs) {
        cout << "polydata::get " << piSNP << "  "<< piSubject << " : " << miSNPs << "  "<< miSubjects << endl;
        exit(1);
    }
    return mmData((int)piSubject,(int)piSNP);
}

void CPolyData::setSubject (const unsigned int piSubject, const TUVector& pmSNPs) {
    set(piSubject, mmAllSNPs, pmSNPs);
}

void CPolyData::setSNP (const unsigned int piSNP, const TUVector& pmSNPs) {
    set(mmAllSubjects, piSNP, pmSNPs);
}

TUVector CPolyData::getSubject (const unsigned int piSubject) const {
    return get(piSubject, mmAllSNPs);
}

TUVector CPolyData::getSNP (const unsigned int piSNP) const {
    return get(mmAllSubjects, piSNP);
}

void CPolyData::set (const unsigned int piSubject, const TUVector& pmSNPs, const TUVector& pmValues) {
    int     i;
    for (i = 0; i < pmSNPs.rows(); i++) {
        set(piSubject,pmSNPs(i),pmValues(i));
    }
}

void CPolyData::set (const TUVector& pmSubjects, const unsigned int piSNP, const TUVector& pmValues) {
    int     i;
    for (i = 0; i < pmSubjects.rows(); i++) {
        set(pmSubjects(i),piSNP,pmValues(i));
    }
}

void CPolyData::set (const TUVector& pmSubjects, const TUVector& pmSNPs, const TUMatrix& pmValues) {
    int     i, j;
    for (i = 0; i < pmSubjects.rows(); i++) {
        for (j = 0; j < pmSNPs.rows(); j++) {
            set(pmSubjects(i),pmSNPs(j),pmValues(i,j));
        }
    }
}

TUVector CPolyData::get (const unsigned int piSubject, const TUVector& pmSNPs) const {
    int         i;
    TUVector    rmValues;
    rmValues.resize(pmSNPs.rows());
    rmValues = 0;
    for (i = 0; i < pmSNPs.rows(); i++) {
        rmValues(i) = get(piSubject,pmSNPs(i));
    }
    return rmValues;
}

TUVector CPolyData::get (const TUVector& pmSubjects, const unsigned int piSNP) const {
    int         i;
    TUVector    rmValues;
    rmValues.resize(pmSubjects.rows());
    rmValues = 0;
    for (i = 0; i < pmSubjects.rows(); i++) {
        rmValues(i) = get(pmSubjects(i),piSNP);
    }
    return rmValues;
}

TUMatrix CPolyData::get (const TUVector& pmSubjects, const TUVector& pmSNPs) const {
    int         i, j;
    TUMatrix    rmValues;
    rmValues.resize(pmSubjects.rows(),pmSubjects.rows());
    rmValues = 0;
    for (i = 0; i < pmSubjects.rows(); i++) {
        for (j = 0; j < pmSNPs.rows(); j++) {
            rmValues(i,j) = get(pmSubjects(i),pmSNPs(j));
        }
    }
    return rmValues;
}

void CPolyData::keepSNPs (const TUVector& pmSNPs) {
}

void CPolyData::keepSubjects (const TUVector& pmSubjects) {
}

void CPolyData::removeSNPs (const TUVector& pmSNPs) {
}

void CPolyData::removeSubjects (const TUVector& pmSubjects) {
}
