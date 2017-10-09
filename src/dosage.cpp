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
#include "dosage.h"

CDosage::CDosage () {
    clear();
}

CDosage::~CDosage () {
    clear();
}

void CDosage::clear() {
    resize(0,0);
    miSubjects = 0;
    miSNPs = 0;
}

unsigned int CDosage::resize (const unsigned int piSubjects, const unsigned int piSNPs) {
    unsigned int        riSize;
    blitz::firstIndex   loIndex1;

    miSubjects = piSubjects;
    miSNPs     = piSNPs;
    riSize     = miSubjects * miSNPs;

    mmMajorData.resize(miSubjects,miSNPs);
    mmMajorData = 0;
    mmMinorData.resize(miSubjects,miSNPs);
    mmMinorData = 0;
    mmAllSubjects.resize(miSubjects);
    mmAllSubjects = loIndex1;
    mmAllSNPs.resize(miSNPs);
    mmAllSNPs = loIndex1;

    return riSize;
}

void CDosage::set (const unsigned int piSubject, const unsigned int piSNP, const TDosage& poDosage) {
    mmMajorData((int)piSubject,(int)piSNP) = poDosage.mrMajor;
    mmMinorData((int)piSubject,(int)piSNP) = poDosage.mrMinor;
}

TRMatrix CDosage::getDosage () const {
    unsigned int    i, j;
    TRMatrix        rmDosage;
    rmDosage.resize(miSubjects,miSNPs);
    for (i = 0; i < miSubjects; i++) {
        for (j = 0; j < miSNPs; j++) {
            rmDosage((int)i,(int)j) = getDosage(i,j);
        }
    }
    return rmDosage;
}

TDosage CDosage::get (const unsigned int piSubject, const unsigned int piSNP) const {
    TDosage roDosage;
    if (piSubject >= miSubjects || piSNP >= miSNPs) {
        cout << "dosage::get " << piSNP << "  "<< piSubject << " : " << miSNPs << "  "<< miSubjects << endl;
        exit(1);
    }
    roDosage.mrMajor = mmMajorData((int)piSubject,(int)piSNP);
    roDosage.mrMinor = mmMinorData((int)piSubject,(int)piSNP);
    roDosage.mrHetero = 1.0 - roDosage.mrMajor - roDosage.mrMinor;
//cout << piSubject << " " << piSNP << " " << roDosage.mrMajor << " " << roDosage.mrHetero << " " << roDosage.mrMinor << endl;
    return roDosage;
}

TReal CDosage::getDosage (const unsigned int piSubject, const unsigned int piSNP) const {
    if (piSubject >= miSubjects || piSNP >= miSNPs) {
        cout << "dosage::get " << piSNP << "  "<< piSubject << " : " << miSNPs << "  "<< miSubjects << endl;
        exit(1);
    }
    return 1.0 - mmMajorData((int)piSubject,(int)piSNP) + mmMinorData((int)piSubject,(int)piSNP);
}

istream& operator>> (istream& in, CDosage& poData) {
/*
    int             i, j, liSNP;
    double          lrPMajor, lrPHetero, lrPMinor;
    string          lsSNP, lsSample;
    blitz::Array<int,1>    lmSample;

    lmSample.resize (piSamples);
    mmDosage.resize (piSNPs, piSamples);

    //  Read the header line.
    in >> lsSNP;
    for (i = 0; i < piSamples; i++) {
        //  Store the sample id order.
        in >> lsSample;
        lmSample(i) = i;  /. ind the actual index.
        //  Ignore the two repeats.
        in >> lsSample;
        in >> lsSample;
    }

    //  Read the dosage data in for each SNP.
    for (i = 0; i < piSNPs; i++) {
        //  Read the SNP.
        in >> lsSNP;
        liSNP = i;  // Find the actual index.
        for (j = 0; j < piSamples; j++) {
            //  Read in the dosages.
            in >> lrPMajor >> lrPHetero >> lrPMinor;
            //  Store the major and minor homozygous probabilities.
            mmDosage(liSNP, lmSample(j)) = 2 * lrPMajor + lrPHetero;
        }
    }

    in.close();
*/
    return in;
/*
    unsigned int    i, liBytes;
    char            lcCode[3];
    unsigned char  *laBytes;
    in.read(lcCode,3);
    //  Individual order, i.e. all SNPs for one individual.
    if (0 == lcCode[2]) {
        liBytes = (poData.miSNPs / 4) + ((poData.miSNPs % 4) == 0 ? 0 : 1);
        laBytes = (unsigned char*) malloc (liBytes);
        for (i = 0; i < poData.miSubjects; i++) {
            in.read((char*) laBytes, liBytes);
            poData.transform(i, laBytes, liBytes, false);
        }
    //  SNP order, i.e. all individual for each SNP.
    } else {
        liBytes = (poData.miSubjects / 4) + ((poData.miSubjects % 4) == 0 ? 0 : 1);
        laBytes = (unsigned char*) malloc (liBytes);
        for (i = 0; i < poData.miSNPs; i++) {
            in.read((char*) laBytes, liBytes);
            poData.transform(i, laBytes, liBytes, true);
        }
    }
    free(laBytes);
*/
    return in;
}

ostream& operator<< (ostream& out, const CDosage& poData) {
/*
    unsigned int    i, j;
    cout << "Genotypic data: subjects " << poData.miSubjects << ", SNPs " << poData.miSNPs << endl;
    for (i = 0; i < poData.miSubjects; i++) {
        for (j = 0; j < poData.miSNPs; j++) {
            switch (poData.get(i,j)) {
            case 0: out << "22 "; break;
            case 1: out << "21 "; break;
            case 2: out << "00 "; break;
            case 3: out << "11 "; break;
            }
        }
        out << endl;
    }
    out << endl;
*/
    return out;
}

void CDosage::setSubject (const unsigned int piSubject, const TDVector& pmSNPs) {
    set(piSubject, mmAllSNPs, pmSNPs);
}

void CDosage::setSNP (const unsigned int piSNP, const TDVector& pmSNPs) {
    set(mmAllSubjects, piSNP, pmSNPs);
}

TDVector CDosage::getSubject (const unsigned int piSubject) const {
    return get(piSubject, mmAllSNPs);
}

TDVector CDosage::getSNP (const unsigned int piSNP) const {
    return get(mmAllSubjects, piSNP);
}

void CDosage::set (const unsigned int piSubject, const TUVector& pmSNPs, const TDVector& pmValues) {
    int     i;
    for (i = 0; i < pmSNPs.rows(); i++) {
        set(piSubject,pmSNPs(i),pmValues(i));
    }
}

void CDosage::set (const TUVector& pmSubjects, const unsigned int piSNP, const TDVector& pmValues) {
    int     i;
    for (i = 0; i < pmSubjects.rows(); i++) {
        set(pmSubjects(i),piSNP,pmValues(i));
    }
}

void CDosage::set (const TUVector& pmSubjects, const TUVector& pmSNPs, const TDMatrix& pmValues) {
    int     i, j;
    for (i = 0; i < pmSubjects.rows(); i++) {
        for (j = 0; j < pmSNPs.rows(); j++) {
            set(pmSubjects(i),pmSNPs(j),pmValues(i,j));
        }
    }
}

TDVector CDosage::get (const unsigned int piSubject, const TUVector& pmSNPs) const {
    int         i;
    TDVector    rmValues;
    rmValues.resize(pmSNPs.rows());
    rmValues = coDosageNull;
    for (i = 0; i < pmSNPs.rows(); i++) {
        rmValues(i) = get(piSubject,pmSNPs(i));
    }
    return rmValues;
}

TDVector CDosage::get (const TUVector& pmSubjects, const unsigned int piSNP) const {
    int         i;
    TDVector    rmValues;
    rmValues.resize(pmSubjects.rows());
    rmValues = coDosageNull;
    for (i = 0; i < pmSubjects.rows(); i++) {
        rmValues(i) = get(pmSubjects(i),piSNP);
    }
    return rmValues;
}

TDMatrix CDosage::get (const TUVector& pmSubjects, const TUVector& pmSNPs) const {
    int         i, j;
    TDMatrix    rmValues;
    rmValues.resize(pmSubjects.rows(),pmSubjects.rows());
    rmValues = coDosageNull;
    for (i = 0; i < pmSubjects.rows(); i++) {
        for (j = 0; j < pmSNPs.rows(); j++) {
            rmValues(i,j) = get(pmSubjects(i),pmSNPs(j));
        }
    }
    return rmValues;
}

void CDosage::keepSNPs (const TUVector& pmSNPs) {
}

void CDosage::keepSubjects (const TUVector& pmSubjects) {
}

void CDosage::removeSNPs (const TUVector& pmSNPs) {
}

void CDosage::removeSubjects (const TUVector& pmSubjects) {
}

TReal CDosage::convertToDosage (const TDosage& poD) {
    return poD.mrHetero + 2.0 * poD.mrMinor;
}

TRVector CDosage::convertToDosage (const TDVector& pmD) {
    int i, liR = pmD.rows();
    TRVector rmD;
    rmD.resize(liR);
    for (i = 0; i < liR; i++) {
        rmD(i) = convertToDosage(pmD(i));
    }
    return rmD;
}

TRVector CDosage::convertToDosage (const TSNPVector& pmS, const bool pbError) {
    TRVector rmD;
    rmD.resize(pmS.rows());
    for (int i = 0; i < pmS.rows(); i++) {
        switch (pmS(i)) {
            case ceSNPMajor:  rmD(i) = 0.0; break;
            case ceSNPHetero: rmD(i) = 1.0; break;
            case ceSNPError:  rmD(i) = pbError ? -1.0 : 1.0; break;
            case ceSNPMinor:  rmD(i) = 2.0; break;
        }
    }
    return rmD;
}

TRVector CDosage::convertToDosage (const THVector& pm0, const THVector& pm1) {
    TRVector rmD;
    rmD.resize(pm0.rows());
    for (int i = 0; i < pm0.rows(); i++) {
        rmD(i) = pm0(i) + pm1(i);
    }
    return rmD;
}

TIVector CDosage::convertToIntegerDosage (const TRVector& pmD) {
    TIVector rmD;
    rmD.resize(pmD.rows());
    for (int i = 0; i < pmD.rows(); i++) {
        rmD(i) = (int) round(pmD(i));
    }
    return rmD;
}
