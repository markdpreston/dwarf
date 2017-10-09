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
#include "random.h"
#include "utility.h"

//  Class static member declarations
stringstream     CUtility::ss;
blitz::Range     CUtility::soRangeAll = blitz::Range::all();
time_t           CUtility::soTime = 0;
map<string, EFileType> CUtility::smFileTypes;

string CUtility::getParameterString(const vector<string>& paParameters, const string psName, const string psDefault) {
    string         lsParameter;
    vector<string> lsNameValue;
    foreach(lsParameter, paParameters) {
        boost::split(lsNameValue, lsParameter, boost::is_any_of("="), boost::algorithm::token_compress_on);
        if (lsNameValue[0] == psName) {
            if (2 == lsNameValue.size()) {
                return lsNameValue[1];
            } else {
                return "";
            }
        }
    }
    return psDefault;
}

bool CUtility::getParameterBoolean(const vector<string>& paParameters, const string psName, const bool pbDefault) {
    string  lsValue;
    lsValue = getParameterString(paParameters, psName);
    if ("true" == lsValue) {
        return true;
    } else if ("false" == lsValue) {
        return false;
    } else {
        return pbDefault;
    }
}

int CUtility::getParameterInteger(const vector<string>& paParameters, const string psName, const int piDefault) {
    string  lsValue;
    lsValue = getParameterString(paParameters, psName);
    if ("" != lsValue) {
        return stringToInt(lsValue);
    } else {
        return piDefault;
    }
}

TReal CUtility::getParameterReal(const vector<string>& paParameters, const string psName, const TReal prDefault) {
    string  lsValue;
    lsValue = getParameterString(paParameters, psName);
    if ("" != lsValue) {
        return stringToReal(lsValue);
    } else {
        return prDefault;
    }
}

TRVector CUtility::makeReal(const TIVector &pmA)   {
    TRVector rmA;
    rmA.resize(pmA.rows());
    for (int i = 0; i < pmA.rows(); i++) {
        rmA(i) = pmA(i);
    }
    return rmA;
}

TRMatrix CUtility::makeReal(const TIMatrix &pmA)   {
    TRMatrix rmA;
    rmA.resize(pmA.shape());
    for (int i = 0; i < pmA.rows(); i++) {
        for (int j = 0; j < pmA.cols(); j++) {
            rmA(i, j) = pmA(i, j);
        }
    }
    return rmA;
}

ETRVector CUtility::makeEigen(const TIVector &pmA)   {
    ETRVector rmA;
    rmA.resize(pmA.rows());
    for (int i = 0; i < pmA.rows(); i++) {
        rmA(i) = pmA(i);
    }
    return rmA;
}

ETRMatrix CUtility::makeEigen(const TIMatrix &pmA)   {
    ETRMatrix rmA;
    rmA.resize(pmA.rows(), pmA.cols());
    for (int i = 0; i < pmA.rows(); i++) {
        for (int j = 0; j < pmA.cols(); j++) {
            rmA(i, j) = pmA(i, j);
        }
    }
    return rmA;
}

ETRVector CUtility::makeEigen(const TRVector &pmA)   {
    ETRVector rmA;
    rmA.resize(pmA.rows());
    for (int i = 0; i < pmA.rows(); i++) {
        rmA(i) = pmA(i);
    }
    return rmA;
}

ETRMatrix CUtility::makeEigen(const TRMatrix &pmA)   {
    ETRMatrix rmA;
    rmA.resize(pmA.rows(), pmA.cols());
    for (int i = 0; i < pmA.rows(); i++) {
        for (int j = 0; j < pmA.cols(); j++) {
            rmA(i, j) = pmA(i, j);
        }
    }
    return rmA;
}

TRVector  CUtility::makeBlitz(const ETRVector &pmA)  {
    TRVector  rmA;
    rmA.resize(pmA.rows());
    for (int i = 0; i < pmA.rows(); i++) {
        rmA(i) = pmA(i);
    }
    return rmA;
}

TRMatrix  CUtility::makeBlitz(const ETRMatrix &pmA)  {
    TRMatrix  rmA;
    rmA.resize(pmA.rows(), pmA.cols());
    for (int i = 0; i < pmA.rows(); i++) {
        for (int j = 0; j < pmA.cols(); j++) {
            rmA(i, j) = pmA(i, j);
        }
    }
    return rmA;
}

void CUtility::setEigen(ETRVector &pmVector, const TIVector &pmA)   {
    if (pmVector.rows() != pmA.rows()) {
        pmVector.resize(pmA.rows());
        pmVector = makeEigen(pmA);
    } else {
        for (int i = 0; i < pmA.rows(); i++) {
            pmVector(i) = pmA(i);
        }
    }
}

void CUtility::setEigen(ETRMatrix& pmMatrix, const TIMatrix& pmA) {
    if (pmMatrix.rows() != pmA.rows() || pmMatrix.cols() != pmA.cols()) {
        pmMatrix.resize(pmA.rows(),pmA.cols());
        pmMatrix = makeEigen(pmA);
    } else {
        for (int i = 0; i < pmA.rows(); i++) {
            for (int j = 0; j < pmA.cols(); j++) {
                pmMatrix(i, j) = pmA(i, j);
            }
        }
    }
}
void CUtility::setEigen(ETRVector& pmVector, const TRVector& pmA) {
    if (pmVector.rows() != pmA.rows()) {
        pmVector.resize(pmA.rows());
        pmVector = makeEigen(pmA);
    } else {
        for (int i = 0; i < pmA.rows(); i++) {
            pmVector(i) = pmA(i);
        }
    }
}

void CUtility::setEigen(ETRMatrix& pmMatrix, const TRMatrix& pmA) {
    if (pmMatrix.rows() != pmA.rows() || pmMatrix.cols() != pmA.cols()) {
        pmMatrix.resize(pmA.rows(),pmA.cols());
        pmMatrix = makeEigen(pmA);
    } else {
        for (int i = 0; i < pmA.rows(); i++) {
            for (int j = 0; j < pmA.cols(); j++) {
                pmMatrix(i, j) = pmA(i, j);
            }
        }
    }
}

void CUtility::initialise() {
    smFileTypes.insert(pair<string, EFileType>("bed",      EFileTypeBed));
    smFileTypes.insert(pair<string, EFileType>("bim",      EFileTypeBim));
    smFileTypes.insert(pair<string, EFileType>("emap",     EFileTypeEmap));
    smFileTypes.insert(pair<string, EFileType>("fam",      EFileTypeFam));
    smFileTypes.insert(pair<string, EFileType>("gens",     EFileTypeGens));
    smFileTypes.insert(pair<string, EFileType>("genes",    EFileTypeGenes));
    smFileTypes.insert(pair<string, EFileType>("haps",     EFileTypeHaps));
    smFileTypes.insert(pair<string, EFileType>("kbac",     EFileTypeKbac));
    smFileTypes.insert(pair<string, EFileType>("legend",   EFileTypeLegend));
    smFileTypes.insert(pair<string, EFileType>("map",      EFileTypeMap));
    smFileTypes.insert(pair<string, EFileType>("markers",  EFileTypeMarkers));
    smFileTypes.insert(pair<string, EFileType>("mito",     EFileTypeMito));
    smFileTypes.insert(pair<string, EFileType>("ms",       EFileTypeMS));
    smFileTypes.insert(pair<string, EFileType>("ped",      EFileTypePed));
    smFileTypes.insert(pair<string, EFileType>("poly",     EFileTypePoly));
    smFileTypes.insert(pair<string, EFileType>("samples",  EFileTypeSamples));
    smFileTypes.insert(pair<string, EFileType>("tab",      EFileTypeTab));
    smFileTypes.insert(pair<string, EFileType>("transmit", EFileTypeTransmit));
    smFileTypes.insert(pair<string, EFileType>("vcf",      EFileTypeVCF));
}

EFileType CUtility::getFileType(const string psType) {
    EFileType   reType = EFileTypeNone;
    map<string, EFileType>::const_iterator loType(smFileTypes.find(psType));
    if (smFileTypes.end() != loType) {
        reType = loType->second;
    }
    return reType;
}
