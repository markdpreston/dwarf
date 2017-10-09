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
#ifndef INTERFACE_R_H
#define INTERFACE_R_H

#include "types.h"

#define INCLUDE_R

#ifdef INCLUDE_R

#include "RInside.h"

class CInterfaceR {
    public:
//    static void     setReal             (const string psName, const TReal prA) { soR.assign(prA,psName); }
    static void     setVector           (const string psName, const ETRVector& pmA) { soR.assign(createRVector(pmA),psName); }
    static void     setMatrix           (const string psName, const ETRMatrix& pmA) { soR.assign(createRMatrix(pmA),psName); }
    static void     setVector           (const string psName, const TRVector& pmA) { soR.assign(createRVector(pmA),psName); }
    static void     setMatrix           (const string psName, const TRMatrix& pmA) { soR.assign(createRMatrix(pmA),psName); }
    static void     setVector           (const string psName, const Rcpp::NumericVector& pmA) { soR.assign(pmA,psName); }
    static void     setMatrix           (const string psName, const Rcpp::NumericMatrix& pmA) { soR.assign(pmA,psName); }
    static void     getMatrix           (const string psName, TRMatrix& pmA);
    static SEXP*    getAnswer           () { return &soAnswer; }
    static void     removeData          (const string psName) { Rcpp::Environment loR; loR.remove(psName); }
    static void     runCommand          (const string psCommand);
    static int      runCommandI         (const string psCommand);
    static TReal    runCommandR         (const string psCommand);
    static TRVector runCommandRV        (const string psCommand);
    static TRMatrix runCommandRM        (const string psCommand);
//    static TSVector runCommandSV        (const string psCommand);
    private:
    static RInside  soR;
    static SEXP     soAnswer;
    static Rcpp::NumericVector createRVector (const ETRVector& pmA);
    static Rcpp::NumericMatrix createRMatrix (const ETRMatrix& pmA);
    static Rcpp::NumericVector createRVector (const TRVector& pmA);
    static Rcpp::NumericMatrix createRMatrix (const TRMatrix& pmA);
    static TRVector getRVector          () { Rcpp::NumericVector lmA(soAnswer); return extractRVector(lmA); }
    static TRMatrix getRMatrix          () { Rcpp::NumericMatrix lmA(soAnswer); return extractRMatrix(lmA); }
//    static TSVector getSVector          () { Rcpp::CharacterVector lmA(soAnswer); return extractSVector(lmA); }
//    static TSMatrix getSMatrix          () { Rcpp::CharacterMatrix lmA(soAnswer); return extractSMatrix(lmA); }
    static void     setRVector          (Rcpp::NumericVector& pmVector, const TRVector& pmA);
    static void     setRMatrix          (Rcpp::NumericMatrix& pmMatrix, const TRMatrix& pmA);
    static void     setRVector          (Rcpp::NumericVector& pmVector, const ETRVector& pmA);
    static void     setRMatrix          (Rcpp::NumericMatrix& pmMatrix, const ETRMatrix& pmA);
    static TRVector extractRVector      (const Rcpp::NumericVector& pmA);
    static TRMatrix extractRMatrix      (const Rcpp::NumericMatrix& pmA);
//    static TSVector extractSVector      (const Rcpp::CharacterVector& pmA);
//    static TSMatrix extractSMatrix      (const Rcpp::CharacterMatrix& pmA);
};

#endif

#endif
