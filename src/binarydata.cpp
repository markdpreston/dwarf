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
#include "binarydata.h"

CBinaryData::CBinaryData () {
    maMasks[0] = 0xC0;
    maMasks[1] = 0x30;
    maMasks[2] = 0x0C;
    maMasks[3] = 0x03;
    maReverse[0] = ceSNPMajor;
    maReverse[1] = ceSNPError;
    maReverse[2] = ceSNPHetero;
    maReverse[3] = ceSNPMinor;
    maData = NULL;
    clear();
}

CBinaryData::~CBinaryData () {
    clear();
}

void CBinaryData::clear () {
    miSize     = 0;
    miSubjects = 0;
    miSNPs     = 0;
    if (maData != NULL) {
        free(maData);
        maData = NULL;
    }
    resize(0,0);
}

unsigned int CBinaryData::resize (const unsigned int piSubjects, const unsigned int piSNPs) {
    bool                lbClear;
    blitz::firstIndex   loIndex1;

    miSubjects = piSubjects;
    miSNPs     = piSNPs;
    miSize     = (miSubjects * miSNPs) / 4 + ((miSubjects * miSNPs) % 4 ? 1 : 0);
    if (miSize > 0) {
        if (0 != miSNPs && piSNPs != miSNPs) {
            cerr << "In subject-major data, illegal change of SNP count." << endl;
            exit(1);
        }
        lbClear = (NULL == maData);
        maData = (unsigned char*) realloc(maData,miSize);
        if (NULL == maData) {
            cerr << "Failed to reallocate memory with " << miSize << "bytes (" << miSubjects << " subjects and " << miSNPs << "SNPs)." << endl;
            exit(0);
        }
        if (lbClear) { memset(maData, 0, miSize); }
    } else {
        maData = NULL;
    }
    mmAllSubjects.resize(miSubjects);
    mmAllSubjects = loIndex1;
    mmAllSNPs.resize(miSNPs);
    mmAllSNPs = loIndex1;
    return miSize;
}

unsigned int CBinaryData::resizeAndPreserve (const unsigned int piSubjects, const unsigned int piSNPs) {
    //  Check that only the subjects are changing due to the storage order.
   if (piSNPs != miSNPs) {
       cerr << "Failed to resize and preserve: " << piSNPs << " != " << miSNPs << endl;
       exit(0);
    }
    //  Resize (realloc => preservation of data).
    return resize(piSubjects, piSNPs);
}

CBinaryData& CBinaryData::operator= (const CBinaryData& poBinaryData) {
    resize(poBinaryData.miSubjects,poBinaryData.miSNPs);
    memcpy(maData,poBinaryData.maData,miSize);
    return *this;
}

void CBinaryData::set (const unsigned int piSubject, const unsigned int piSNP, ESNPCode peValue) {
    unsigned int    liByte, liOffset;
    unsigned char  *lcChar, lcValue;

    liByte   = (piSubject * miSNPs + piSNP) / 4;
    liOffset = (piSubject * miSNPs + piSNP) % 4;

//    if (0 == miSubjects || piSubject >= miSubjects) { cerr << "Attempt to write subject too large: " << piSubject << "; max: " << miSubjects << endl; exit(0); }
//    if (0 == miSNPs || piSNP >= miSNPs) { cerr << "Attempt to write SNP too large: " << piSNP << "; max: " << miSNPs << endl; exit(0); }
//    if (liByte >= miSize) {cerr << "Attempt to write byte too large: " << liByte << "; max: " << miSize << endl; exit(0); }

    //  Convert ESPCode to enable it to be <<ed.
    lcValue = (unsigned char) peValue;
    lcValue <<= (6 - 2 * liOffset);

    lcChar = maData;
    lcChar[liByte] = (lcChar[liByte] & (255 - maMasks[liOffset])) + lcValue;
}

void CBinaryData::setSubject (const unsigned int piSubject, const TSNPVector& pmSNPs) {
    set(piSubject, mmAllSNPs, pmSNPs);
}

void CBinaryData::setSNP (const unsigned int piSNP, const TSNPVector& pmSNPs) {
    set(mmAllSubjects, piSNP, pmSNPs);
}

TSNPVector CBinaryData::getSubject (const unsigned int piSubject) const {
    return get(piSubject, mmAllSNPs);
}

TSNPVector CBinaryData::getSNP (const unsigned int piSNP) const {
    return get(mmAllSubjects, piSNP);
}

void CBinaryData::set (const unsigned int piSubject, const TUVector& pmSNPs, const TSNPVector& pmValues) {
    int     i;
    for (i = 0; i < pmSNPs.rows(); i++) {
        set(piSubject,pmSNPs(i),pmValues(i));
    }
}

void CBinaryData::set (const TUVector& pmSubjects, const unsigned int piSNP, const TSNPVector& pmValues) {
    int     i;
    for (i = 0; i < pmSubjects.rows(); i++) {
        set(pmSubjects(i),piSNP,pmValues(i));
    }
}

void CBinaryData::set (const TUVector& pmSubjects, const TUVector& pmSNPs, const TSNPMatrix& pmValues) {
    int     i, j;
    for (i = 0; i < pmSubjects.rows(); i++) {
        for (j = 0; j < pmSNPs.rows(); j++) {
            set(pmSubjects(i),pmSNPs(j),pmValues(i,j));
        }
    }
}

TSNPVector CBinaryData::get (const unsigned int piSubject, const TUVector& pmSNPs) const {
    int         i;
    TSNPVector  rmValues;
    rmValues.resize(pmSNPs.rows());
    rmValues = ceSNPError;
//cout << "BD " << piSubject << " " << miSNPs << " " << miSubjects << endl << flush;
    for (i = 0; i < pmSNPs.rows(); i++) {
        rmValues(i) = get(piSubject,pmSNPs(i));
    }
//cout << "BD" << endl << flush;
    return rmValues;
}

TSNPVector CBinaryData::get (const TUVector& pmSubjects, const unsigned int piSNP) const {
    int         i;
    TSNPVector  rmValues;
    rmValues.resize(pmSubjects.rows());
    rmValues = ceSNPError;
    for (i = 0; i < pmSubjects.rows(); i++) {
        rmValues(i) = get(pmSubjects(i),piSNP);
    }
    return rmValues;
}

TSNPMatrix CBinaryData::get (const TUVector& pmSubjects, const TUVector& pmSNPs) const {
    int         i, j;
    TSNPMatrix  rmValues;
    rmValues.resize(pmSubjects.rows(),pmSubjects.rows());
    rmValues = ceSNPError;
    for (i = 0; i < pmSubjects.rows(); i++) {
        for (j = 0; j < pmSNPs.rows(); j++) {
            rmValues(i,j) = get(pmSubjects(i),pmSNPs(j));
        }
    }
    return rmValues;
}

void CBinaryData::keepSNPs(const TUVector& pmSNPs) {
    unsigned int    i, liNewSNPs;
    CBinaryData     loData;
    liNewSNPs = sum(0 <= pmSNPs && pmSNPs < miSNPs);
    if (0 == miSubjects || 0 == miSNPs || 0 == liNewSNPs) { return; }
    loData.resize(miSubjects,liNewSNPs);
    for (i = 0; i < liNewSNPs; i++) {
        loData.setSNP(i,getSNP(pmSNPs(i)));
    }
    resize(miSubjects,liNewSNPs);
    memcpy(maData,loData.maData,miSize);
}

void CBinaryData::keepSubjects(const TUVector& pmSubjects) {
    unsigned int    i, liNewSubjects;
    CBinaryData     loData;
    liNewSubjects = sum(0 <= pmSubjects && pmSubjects < miSubjects);
    if (0 == miSubjects || 0 == miSubjects || 0 == liNewSubjects) { return; }
    loData.resize(liNewSubjects,miSNPs);
    for (i = 0; i < liNewSubjects; i++) {
        loData.setSubject(i,getSubject(pmSubjects(i)));
    }
    resize(liNewSubjects,miSNPs);
    memcpy(maData,loData.maData,miSize);
}

void CBinaryData::removeSNPs(const TUVector& pmSNPs) {
    unsigned int    i, liCount;
    TUVector        lmSNPs;
    if (0 == miSubjects || 0 == miSubjects) { return; }
    lmSNPs.resize(miSNPs - pmSNPs.rows() + sum(pmSNPs >= miSNPs));
    liCount = 0;
    for (i = 0; i < miSNPs; i++) {
        if (0 == sum(pmSNPs == i)) {
            lmSNPs(liCount) = i;
            liCount++;
        }
    }
    keepSNPs(lmSNPs);
}

void CBinaryData::removeSubjects(const TUVector& pmSubjects) {
    unsigned int    i, liCount;
    TUVector        lmSubjects;
    if (0 == miSubjects || 0 == miSubjects) { return; }
    lmSubjects.resize(miSubjects - pmSubjects.rows() + sum(pmSubjects >= miSubjects));
    liCount = 0;
    for (i = 0; i < miSubjects; i++) {
        if (0 == sum(pmSubjects == i)) {
            lmSubjects(liCount) = i;
            liCount++;
        }
    }
    keepSubjects(lmSubjects);
}

TRMatrix CBinaryData::getDosage (const bool pbHaplotypes) const {
    unsigned int    i, j;
    TRMatrix        rmDosage;
    rmDosage.resize(miSubjects,miSNPs);
    for (i = 0; i < miSubjects; i++) {
        for (j = 0; j < miSNPs; j++) {
            rmDosage((int)i,(int)j) = getDosage(i,j,pbHaplotypes);
        }
    }
    return rmDosage;
}

TReal CBinaryData::getDosage (const unsigned int piSubject, const unsigned int piSNP, const bool pbHaplotypes) const {
    ESNPCode    leSNP;
    leSNP = get(piSubject, piSNP);
    return getDosage(leSNP, pbHaplotypes);
}

TReal CBinaryData::getDosage (const ESNPCode peSNP, const bool pbHaplotypes) {
    TReal   rrDosage = 0.0;
    if (pbHaplotypes) {
        switch (peSNP) {
            case ceHaplotype00: rrDosage = 0.0; break;
            case ceHaplotype01: rrDosage = 1.0; break;
            case ceHaplotype10: rrDosage = 1.0; break;
            case ceHaplotype11: rrDosage = 2.0; break;
        }
    } else {
        switch (peSNP) {
            case ceSNPMajor:  rrDosage = 0.0; break;
            case ceSNPHetero: rrDosage = 1.0; break;
            case ceSNPError:  rrDosage = 0.0; break;
            case ceSNPMinor:  rrDosage = 2.0; break;
        }
    }
    return rrDosage;
}

void CBinaryData::internalError(const int piError) const {
    switch (piError) {
//        case 1: cerr << "Attempt to read subject too large: " << piSubject << "; max: " << miSubjects << endl; break;
//        case 2: cerr << "Attempt to read SNP too large: " << piSNP << "; max: " << miSNPs << endl; break;
//        case 3: cerr << "Attempt to write byte too large: " << liByte << "; max: " << miSize << endl; break;
    }
    exit(0);
}

void CBinaryData::transform (const unsigned int piRecord, const unsigned char* paBytes, const unsigned int piBytes, const bool pbSNP) {
    unsigned int    i, j;
    unsigned char   liSNP;
    for (i = 0; i < piBytes; i++) {
        for (j = 0; j < 4; j++) {
            liSNP = (paBytes[i] & maMasks[3-j]) >> (2*j);
            if (pbSNP) {
                if (4*i+j < miSubjects) {
                    set(4*i+j, piRecord, maReverse[liSNP]);
                }
            } else {
                if (4*i+j < miSNPs) {
                    set(piRecord, 4*i+j, maReverse[liSNP]);
                }
            }
        }
    }
}

void CBinaryData::untransform (const unsigned int piRecord, unsigned char* paBytes, const unsigned int piBytes, const bool pbSNP, const bool pbHaploid) const {
    unsigned int    i, j, leSNP;
    unsigned char   liSNP;
    for (i = 0; i < piBytes; i++) {
        paBytes[i] = 0;
        for (j = 0; j < 4; j++) {
            leSNP = 0;
            if (pbSNP) {
                if (4*i+j < miSubjects) {
                    leSNP = get(4*i+j, piRecord);
                }
            } else {
                if (4*i+j < miSNPs) {
                    leSNP = get(piRecord,4*i+j);
                }
            }
//            if (leSNP == 2) { leSNP = 1; }
            liSNP = maReverse[leSNP];
            paBytes[i] = paBytes[i] | liSNP << (2*j);
        }
    }
}
