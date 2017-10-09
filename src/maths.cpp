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
#include "maths.h"

int      CMaths::siBinomialCount = 0;
TIVector CMaths::smBinomialIndices;
vector<TRVector>  CMaths::smBinomial;

TRMatrix CMaths::matrixInversion (const TRMatrix& pmMatrix) {
    int         i, j, liRows;
    TRMatrix    rmInverse;
    liRows = pmMatrix.rows();
    rmInverse.resize(liRows,liRows);
    rmInverse = 0;

    ublas::matrix<TReal>              lmWork(liRows,liRows), lmInverse(liRows,liRows);
    ublas::permutation_matrix<size_t> lmPermutation(liRows);

    for (i = 0; i < liRows; i++) {
        for (j = 0; j < liRows; j++) {
            lmWork(i,j) = pmMatrix(i,j);
        }
    }

    if (0 != lu_factorize(lmWork,lmPermutation)) {
        cerr << "Matrix inversion failed." << endl;
        exit(0);
    }

    lmInverse.assign(ublas::identity_matrix<TReal>(liRows));
    lu_substitute(lmWork, lmPermutation, lmInverse);

    for (i = 0; i < liRows; i++) {
        for (j = 0; j < liRows; j++) {
            rmInverse(i,j) = lmInverse(i,j);
        }
    }

    return rmInverse;
}

TReal CMaths::binomial(const int piN, const int piK) {
    int         i, liIndex;
    TRVector    lmBinomial;
    TReal       rrBinomial;

    if (piK > piN) {
        return 0.0;
    }

    if (0 == siBinomialCount) {
        smBinomialIndices.resize(ciBinomialCount);
        smBinomialIndices = -1;
    }

    liIndex = -1;
    for (i = 0; i < siBinomialCount; i++) {
        if (piN == smBinomialIndices(i)) {
            liIndex = i;
            break;
        }
    }
    //liIndex = blitz::find(piN,smBinomialIndices);

    if (0 > liIndex) {
        if (siBinomialCount < ciBinomialCount) {
            smBinomialIndices(siBinomialCount) = piN;
            siBinomialCount++;
            lmBinomial.resize(piN+1);
            for (i = 0; i <= piN; i++) {
                lmBinomial(i) = binomialInternal(piN, i);
            }
            smBinomial.push_back(lmBinomial);
            rrBinomial = lmBinomial(piK);
        } else {
            rrBinomial = binomialInternal(piN,piK);
        }
    } else {
        rrBinomial = smBinomial[liIndex](piK);
    }
    return rrBinomial;
}

TReal CMaths::binomialInternal(const int piN, const int piK) {
    TReal   rrBinomial;
    try {
        rrBinomial = boost::math::binomial_coefficient<long double>(piN, piK);
    } catch (...) {
        rrBinomial = exp(lgamma(piN) - lgamma(piK) - lgamma(piN - piK));
    }
    return rrBinomial;
}
