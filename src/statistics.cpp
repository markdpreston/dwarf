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
#include "statistics.h"
#include "maths.h"

TReal CStatistics::mean (const TRVector& pmA) {
    return blitz::mean(pmA);
}

TRVector CStatistics::mean (const TRMatrix& pmA) {
    int         i;
    TRVector    rmMean;
    rmMean.resize(pmA.cols());
    for (i = 0; i < pmA.cols(); i++) {
        rmMean(i) = blitz::mean(pmA(CUtility::getAll(),i));
    }
    return rmMean;
}

TReal CStatistics::variance (const TRVector& pmA) {
    TReal   lrMean = mean(pmA);
    return blitz::mean(sqr(pmA - lrMean));
}

TRVector CStatistics::variance (const TRMatrix& pmA) {
    int         i;
    TRVector    lmMean = mean(pmA), rmVariance;
    rmVariance.resize(pmA.cols());
    for (i = 0; i < pmA.cols(); i++) {
        rmVariance(i) = blitz::mean(sqr(pmA(CUtility::getAll(),i) - lmMean(i)));
    }
    return rmVariance;
}

TReal CStatistics::covariance (const TRVector& pmA, const TRVector& pmB) {
    TReal   lrMeanA = mean(pmA), lrMeanB = mean(pmB);
    return sum((pmA - lrMeanA) * (pmB - lrMeanB)) / pmA.rows();
}


TRMatrix CStatistics::covariance (const TRMatrix& pmA) {
    int         i, j, liRows, liCols;
    TRVector    lmA, lmB;
    TRMatrix    rmCovariance;
    liRows = pmA.rows();
    liCols = pmA.cols();
    lmA.resize(liRows);
    lmB.resize(liRows);
    rmCovariance.resize(liCols, liCols);
    for (i = 0; i < liCols; i++) {
        lmA = pmA(CUtility::getAll(),i);
        for (j = i; j < liCols; j++) {
            lmB = pmA(CUtility::getAll(),j);
            rmCovariance(i,j) = covariance(lmA, lmB);
            rmCovariance(j,i) = rmCovariance(i,j);
        }
    }
    return rmCovariance;
}

TRMatrix CStatistics::correlation (const TRMatrix& pmA, const bool pbCovariance) {
    int         i, j, liRows;
    TRVector    lmSD;
    TRMatrix    rmCorrelation;
    liRows = pmA.rows();
    lmSD.resize(liRows);
    rmCorrelation.resize(liRows, liRows);
    if (pbCovariance) {
       rmCorrelation = pmA;
    } else {
       rmCorrelation = covariance(pmA);
    }
    for (i = 0; i < liRows; i++) {
        lmSD(i) = sqrt(rmCorrelation(i,i));
    }
    for (i = 0; i < liRows; i++) {
        for (j = i; j < liRows; j++) {
            rmCorrelation(i,j) /= lmSD(i) * lmSD(j);
            rmCorrelation(j,i)  = rmCorrelation(i,j);
        }
    }
    return rmCorrelation;
}

TResult CStatistics::chiSquared (const TRMatrix& pmA) {
    TReal               lrTotal;
    TRVector            lmRowSum, lmColSum;
    TRMatrix            lmPredicted;
    blitz::firstIndex   loIndex1;
    blitz::secondIndex  loIndex2;
    TResult             rmChiSquared;
    //  Shape the vectors and matrices.
    lmRowSum.resize(pmA.rows());
    lmColSum.resize(pmA.cols());
    lmPredicted.resize(pmA.shape());
    //  Get the sums.
    lmRowSum = sum(pmA, loIndex2);
    lmColSum = sum(pmA(loIndex2,loIndex1), loIndex2);
    lrTotal  = sum(pmA);
    //  Create the frequencies table.
    lmPredicted = lmRowSum(loIndex1) * lmColSum(loIndex2) / lrTotal;
    if (any(lmPredicted == 0)) {
        rmChiSquared(0) = -1;
        rmChiSquared(1) = 1;
    } else {
        //  Work out the chi-squared value.
        rmChiSquared(0) = sum(sqr(pmA - lmPredicted) / lmPredicted);
        rmChiSquared(1) = pvalueChiSquared(rmChiSquared(0), (pmA.rows() - 1) * (pmA.cols() - 1));
    }
    return rmChiSquared;
}

TResult CStatistics::hotellingTSquared (const TRMatrix& pmA, const TRMatrix& pmB) {
    int                 liA, liB, liCols;
    TRVector            lmMean, lmTemp;
    TRMatrix            lmCovA, lmCovB, lmS;
    blitz::firstIndex   loIndex1;
    blitz::secondIndex  loIndex2;
    TResult             rmHotellingTSquared;
    //  Prepare variables.
    liA = pmA.rows();
    liB = pmB.rows();
    liCols = pmA.cols();
    lmMean.resize(liCols);
    lmTemp.resize(liCols);
    lmCovA.resize(liCols, liCols);
    lmCovB.resize(liCols, liCols);
    lmS.resize(liCols, liCols);
    //  Get the variant means and covariances.
    lmMean = CStatistics::mean(pmA) - CStatistics::mean(pmB);
    lmCovA = CStatistics::covariance(pmA);
    lmCovB = CStatistics::covariance(pmB);
    lmS = (liA * lmCovA + liB * lmCovB) / (liA + liB - 2);
    lmS = CMaths::matrixInversion(lmS);
    //  Calculate Hotelling's T^2 statistic and p value.
    lmTemp = sum(lmS(loIndex1, loIndex2) * lmMean(loIndex2), loIndex2);
    rmHotellingTSquared(0) = sum(lmTemp * lmMean) * liA * liB / (liA + liB) * (liA + liB - liCols - 1) / liCols / (liA + liB - 2);
    rmHotellingTSquared(1) = pvalueFisherF(rmHotellingTSquared(0), liCols, liA + liB - liCols - 1);
    return rmHotellingTSquared;
}

TResult CStatistics::normal (const TReal prValue, const TRVector& pmValues) {
    TResult rmNormal;
    rmNormal(0) = max(pmValues) > min(pmValues) ? (prValue - mean(pmValues)) / sqrt (variance(pmValues)) : -10000;
    rmNormal(1) = pvalueNormal(rmNormal(0));
    return rmNormal;
}

TReal CStatistics::pvalueChiSquared (const TReal prValue, const TReal prDF) {
    TReal   rrPvalue;
    if (! boost::math::isnan(prValue) && ! boost::math::isnan(prDF) && prValue > -1e-9) {
        boost::math::chi_squared loChiSquared(prDF);
        rrPvalue = 1 - boost::math::cdf(loChiSquared,prValue);
    } else {
        rrPvalue = 1;
        cerr << "Chi";
    }
    return rrPvalue;
}

TReal CStatistics::pvalueFisherF (const TReal prValue, const unsigned int piDoF1, const unsigned int piDoF2) {
    boost::math::fisher_f loFisher(piDoF1, piDoF2);
    return 1 - boost::math::cdf(loFisher,prValue);
}

TReal CStatistics::pvalueNormal (const TReal prValue) {
    boost::math::normal loNormal;
    TReal rrCDF;
    if (! boost::math::isnan(prValue)) {
        rrCDF = boost::math::cdf(loNormal,prValue);
    } else {
        rrCDF = 1.0;
    }
//    return lrCDF < 0.5 ? lrCDF : 1.0 - lrCDF;
    return rrCDF;
}

TReal CStatistics::pvalueStudentsT (const TReal prValue, const TReal piDoF) {
    boost::math::students_t loStudentT(piDoF);
    return 1 - boost::math::cdf(loStudentT,prValue);
}

//  The storage order for the covariance doesn't matter as it is symmetric anyway.  :)
TReal CStatistics::pvalueMultiNormal (const TRVector& pmLower, const TRMatrix& pmCovariance) {
/*
    int         i, j, liN, liMax, liReturn;
    double      lrAbsEps, lrRelEps, lrError, lrValue;
    TIVector    lmLimits;
    blitz::Array<double,1> lmLower, lmUpper;
    blitz::Array<double,2> lmCovariance;

    liN = pmLower.rows();
    liMax = 2000*liN;

    lrAbsEps = 1e-6;
    lrRelEps = 0.0;

    lmLower.resize(liN);
    lmUpper.resize(liN);
    lmLimits.resize(liN);
    lmCovariance.resize(liN,liN);

    for (i = 0; i < liN; i++) {
        lmLower(i)  = pmLower(i);
        lmUpper(i)  = numeric_limits<double>::infinity();
        lmLimits(i) = 1;
        for (j = 0; j < liN; j++) {
            lmCovariance(i,j) = pmCovariance(i,j);
        }
    }

//    sadmvn_(&liN, lmLower.data(), lmUpper.data(), lmLimits.data(), lmCovariance.data(), &liMax, &lrAbsEps, &lrRelEps, &lrError, &lrValue, &liReturn);
cout << setw(12) << lrError << setw(12) << liReturn << setw(12) << lrValue << endl;

    return lrValue;
*/
    return 1.0;
}
