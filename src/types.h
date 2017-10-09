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
#ifndef TYPES_H
#define TYPES_H

#include <ctime>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <blitz/array.h>
#include <eigen3/Eigen/Dense>
#include "boost.h"

using namespace std;

//  0 is the norm, 1 is the variant.
typedef enum {
    ceSNPMajor    = 0,
    ceSNPHetero   = 1,
    ceSNPError    = 2,
    ceSNPMinor    = 3,
    ceHaplotype00 = 0,
    ceHaplotype01 = 1,
    ceHaplotype10 = 2,
    ceHaplotype11 = 3
} ESNPCode;
typedef enum {
    ceSNPMajorMajor   = (ceSNPMajor << 2)  + ceSNPMajor,
    ceSNPMajorHetero  = (ceSNPMajor << 2)  + ceSNPHetero,
    ceSNPMajorMinor   = (ceSNPMajor << 2)  + ceSNPMinor,
    ceSNPHeteroMajor  = (ceSNPHetero << 2) + ceSNPMajor,
    ceSNPHeteroHetero = (ceSNPHetero << 2) + ceSNPHetero,
    ceSNPHeteroMinor  = (ceSNPHetero << 2) + ceSNPMinor,
    ceSNPMinorMajor   = (ceSNPMinor << 2)  + ceSNPMajor,
    ceSNPMinorHetero  = (ceSNPMinor << 2)  + ceSNPHetero,
    ceSNPMinorMinor   = (ceSNPMinor << 2)  + ceSNPMinor
} ESNPCombined;

typedef long double                     TReal;
typedef blitz::Array<bool,1>            TBVector;
typedef blitz::Array<bool,2>            TBMatrix;
typedef blitz::Array<unsigned char,1>   TCVector;
typedef blitz::Array<unsigned char,2>   TCMatrix;
typedef blitz::Array<bool,1>            THVector;
typedef blitz::Array<bool,2>            THMatrix;
typedef blitz::Array<int,1>             TIVector;
typedef blitz::Array<int,2>             TIMatrix;
typedef blitz::Array<TReal,1>           TRVector;
typedef blitz::Array<TReal,2>           TRMatrix;
typedef blitz::Array<string,1>          TSVector;
typedef blitz::Array<string,2>          TSMatrix;
typedef blitz::Array<ESNPCode,1>        TSNPVector;
typedef blitz::Array<ESNPCode,2>        TSNPMatrix;
typedef blitz::Array<unsigned int,1>    TUVector;
typedef blitz::Array<unsigned int,2>    TUMatrix;
typedef blitz::TinyVector<TReal,2>      TResult;

typedef Eigen::Matrix<int,Eigen::Dynamic,1>                 ETIVector;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>    ETIMatrix;
typedef Eigen::Matrix<TReal,Eigen::Dynamic,1>               ETRVector;
typedef Eigen::Matrix<TReal,Eigen::Dynamic,Eigen::Dynamic>  ETRMatrix;

const TReal crEpsilon  = 0.000001;
const TReal crInfinity = numeric_limits<TReal>::infinity();

#include "utility.h"

#endif
