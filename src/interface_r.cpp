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
#include "types.h"
#include "interface_r.h"

#ifdef INCLUDE_R

RInside CInterfaceR::soR(0,NULL);
SEXP    CInterfaceR::soAnswer;

void CInterfaceR::runCommand(const string psCommand) {
    try {
        soAnswer = soR.parseEval(psCommand);
    } catch (...) {
        cerr << "R error: " << psCommand << endl;
        exit(0);
    }
}
int CInterfaceR::runCommandI(const string psCommand) {
    runCommand(psCommand);
    return Rcpp::as<int>(soAnswer);
}
TReal CInterfaceR::runCommandR(const string psCommand) {
    runCommand(psCommand);
    return Rcpp::as<TReal>(soAnswer);
}

TRVector CInterfaceR::runCommandRV(const string psCommand) {
    runCommand(psCommand);
    return getRVector();
}
TRMatrix CInterfaceR::runCommandRM(const string psCommand) {
    runCommand(psCommand);
    return getRMatrix();
}
/*
TSVector CInterfaceR::runCommandSV(const string psCommand) {
    runCommand(psCommand);
    return getSVector();
}
*/

Rcpp::NumericVector CInterfaceR::createRVector(const TRVector& pmA) {
    int i;
    Rcpp::NumericVector rmVector(pmA.rows());
    for (i = 0; i < pmA.rows(); i++) {
        rmVector(i) = pmA(i);
    }
    return rmVector;
}

Rcpp::NumericMatrix CInterfaceR::createRMatrix(const TRMatrix& pmA) {
    int i, j;
    Rcpp::NumericMatrix rmMatrix(pmA.rows(),pmA.cols());
    for (i = 0; i < pmA.rows(); i++) {
        for (j = 0; j < pmA.cols(); j++) {
            rmMatrix(i,j) = pmA(i,j);
        }
    }
    return rmMatrix;
}

Rcpp::NumericVector CInterfaceR::createRVector(const ETRVector& pmA) {
    int i;
    Rcpp::NumericVector rmVector(pmA.rows());
    for (i = 0; i < pmA.rows(); i++) {
        rmVector(i) = pmA(i);
    }
    return rmVector;
}

Rcpp::NumericMatrix CInterfaceR::createRMatrix(const ETRMatrix& pmA) {
    int i, j;
    Rcpp::NumericMatrix rmMatrix(pmA.rows(),pmA.cols());
    for (i = 0; i < pmA.rows(); i++) {
        for (j = 0; j < pmA.cols(); j++) {
            rmMatrix(i,j) = pmA(i,j);
        }
    }
    return rmMatrix;
}

void CInterfaceR::setRVector(Rcpp::NumericVector& pmVector, const TRVector& pmA) {
    int i;
    if (pmVector.size() != pmA.rows()) {
        cerr << "No R resize";
        exit(0);
//        pmVector.resize(pmA.rows());
//        pmVector = createRVector(pmA);
    } else {
        for (i = 0; i < pmA.rows(); i++) {
            pmVector(i) = pmA(i);
        }
    }
}

void CInterfaceR::setRMatrix(Rcpp::NumericMatrix& pmMatrix, const TRMatrix& pmA) {
    int i, j;
    if (pmMatrix.nrow() != pmA.rows() || pmMatrix.ncol() != pmA.cols()) {
        cerr << "No R resize";
        exit(0);
    } else {
        for (i = 0; i < pmA.rows(); i++) {
            for (j = 0; j < pmA.cols(); j++) {
                pmMatrix(i,j) = pmA(i,j);
            }
        }
    }
}

void CInterfaceR::setRVector(Rcpp::NumericVector& pmVector, const ETRVector& pmA) {
    int i;
    if (pmVector.size() != pmA.rows()) {
        cerr << "No R resize";
        exit(0);
    } else {
        for (i = 0; i < pmA.rows(); i++) {
            pmVector(i) = pmA(i);
        }
    }
}

void CInterfaceR::setRMatrix(Rcpp::NumericMatrix& pmMatrix, const ETRMatrix& pmA) {
    int i, j;
    if (pmMatrix.nrow() != pmA.rows() || pmMatrix.ncol() != pmA.cols()) {
        cerr << "No R resize";
        exit(0);
    } else {
        for (i = 0; i < pmA.rows(); i++) {
            for (j = 0; j < pmA.cols(); j++) {
                pmMatrix(i,j) = pmA(i,j);
            }
        }
    }
}

TRVector CInterfaceR::extractRVector(const Rcpp::NumericVector& pmA) {
    int i;
    TRVector rmVector;
    rmVector.resize(pmA.size());
    for (i = 0; i < pmA.size(); i++) {
        rmVector(i) = pmA(i);
    }
    return rmVector;
}

TRMatrix CInterfaceR::extractRMatrix(const Rcpp::NumericMatrix& pmA) {
    int i, j;
    TRMatrix rmMatrix;
    rmMatrix.resize(pmA.nrow(),pmA.ncol());
    for (i = 0; i < pmA.nrow(); i++) {
        for (j = 0; j < pmA.ncol(); j++) {
            rmMatrix(i,j) = pmA(i,j);
        }
    }
    return rmMatrix;
}
/*
TSVector CInterfaceR::extractSVector(const Rcpp::CharacterVector& pmA) {
    int i;
    TSVector rmVector;
    rmVector.resize(pmA.size());
    for (i = 0; i < pmA.size(); i++) {
        rmVector(i) = (string) pmA(i);
    }
    return rmVector;
}
TSMatrix CInterfaceR::extractSMatrix(const Rcpp::CharacterMatrix& pmA) {
    int i, j;
    TSMatrix rmMatrix;
    rmMatrix.resize(pmA.nrow(),pmA.ncol());
    for (i = 0; i < pmA.nrow(); i++) {
        for (j = 0; j < pmA.ncol(); j++) {
            rmMatrix(i,j) = pmA(i,j);
        }
    }
    return rmMatrix;
}
*/

void CInterfaceR::getMatrix(const string psName, TRMatrix& pmA) {
    unsigned int    i, j, liCols, liRows;
    if (psName.length()) {
        runCommand(psName);
    }
    liRows = Rf_nrows(CInterfaceR::soAnswer);
    liCols = Rf_ncols(CInterfaceR::soAnswer);
    pmA.resize(liRows,liCols);
    if (1 == liRows && 1 == liCols) {
        pmA(0,0) = Rcpp::as<TReal>(soAnswer);
    } else if (1 == liRows || 1 == liCols) {
        Rcpp::NumericVector lmB(CInterfaceR::soAnswer);
        pmA.resize(max(liRows,liCols),1);
        for (i = 0; i < max(liRows,liCols); i++) {
            pmA((int)i,0) = lmB(i);
        }
    } else {
        Rcpp::NumericMatrix lmB(CInterfaceR::soAnswer);
        pmA.resize(liRows,liCols);
        for (i = 0; i < liRows; i++) {
            for (j = 0; j < liCols; j++) {
                pmA((int)i,(int)j) = lmB(i,j);
            }
        }
    }
}

#endif
