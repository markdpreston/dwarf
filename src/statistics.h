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
#ifndef STATISTICS_H
#define STATISTICS_H

#include "types.h"

extern"C" {
    void sadmvn_(int *n, double *lower, double *upper, int *inf, double *correl, int *max, double *abseps, double *releps, double *error, double *value, int *inform);
}

class CStatistics {
    public:
    //  Basic statistical methods.
    static TReal    mean                (const TRVector& pmA);
    static TRVector mean                (const TRMatrix& pmA);
    static TReal    variance            (const TRVector& pmA);
    static TRVector variance            (const TRMatrix& pmA);
    static TReal    covariance          (const TRVector& pmA, const TRVector& pmB);
    static TRMatrix covariance          (const TRMatrix& pmA);
    static TRMatrix correlation         (const TRMatrix& pmA, const bool pbCovariance);
    //  Statistical tests.
    static TResult  chiSquared          (const TRMatrix& pmA);
    static TResult  hotellingTSquared   (const TRMatrix& pmA, const TRMatrix& pmB);
    static TResult  normal              (const TReal prValue, const TRVector& pmA);
    //  Integer versions.
    static TReal    mean                (const TIVector& pmA) { return mean(CUtility::makeReal(pmA)); }
    static TRVector mean                (const TIMatrix& pmA) { return mean(CUtility::makeReal(pmA)); }
    static TReal    variance            (const TIVector& pmA) { return variance(CUtility::makeReal(pmA)); }
    static TReal    covariance          (const TIVector& pmA, const TIVector& pmB) { return covariance(CUtility::makeReal(pmA), CUtility::makeReal(pmB)); }
    static TRMatrix covariance          (const TIMatrix& pmA) { return covariance(CUtility::makeReal(pmA)); }
    static TResult  chiSquared          (const TIMatrix& pmA) { return chiSquared(CUtility::makeReal(pmA)); }
    static TResult  hotellingTSquared   (const TIMatrix& pmA, const TIMatrix& pmB) { return hotellingTSquared(CUtility::makeReal(pmA),CUtility::makeReal(pmB)); }
    static TResult  normal              (const int piValue,   const TIVector& pmA) { return normal(1.0*piValue, CUtility::makeReal(pmA)); }
    //  P-value calculations.
    static TReal    pvalueChiSquared    (const TReal prValue, const TReal prDoF);
    static TReal    pvalueFisherF       (const TReal prValue, const unsigned int piDoF1, const unsigned int piDoF2);
    static TReal    pvalueNormal        (const TReal prValue);
    static TReal    pvalueStudentsT     (const TReal prValue, const TReal prDoF);
    static TReal    pvalueMultiNormal   (const TRVector& pmLower, const TRMatrix& pmCorrelation);
    private:
};

#endif
