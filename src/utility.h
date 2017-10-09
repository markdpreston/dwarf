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
#ifndef UTILITY_H
#define UTILITY_H

#include "types.h"

typedef enum {
    EFileTypeNone = 0,
    EFileTypeBed,
    EFileTypeBim,
    EFileTypeEmap,
    EFileTypeFam,
    EFileTypeGens,
    EFileTypeGenes,
    EFileTypeHaps,
    EFileTypeKbac,
    EFileTypeLegend,
    EFileTypeMap,
    EFileTypeMarkers,
    EFileTypeMito,
    EFileTypeMS,
    EFileTypePed,
    EFileTypePoly,
    EFileTypeSamples,
    EFileTypeTab,
    EFileTypeTransmit,
    EFileTypeVCF
} EFileType;

class CUtility {
    public:
    static void     initialise      ();

    template<class T>
    static void     addToVector     (blitz::Array<T,1>& a, T b)  { int i = a.rows(); a.resizeAndPreserve(i + 1); a(i) = b; }

    static TRVector makeReal        (const TIVector& pmA);
    static TRMatrix makeReal        (const TIMatrix& pmA);

    static ETRVector makeEigen      (const TIVector& pmA);
    static ETRMatrix makeEigen      (const TIMatrix& pmA);
    static ETRVector makeEigen      (const TRVector& pmA);
    static ETRMatrix makeEigen      (const TRMatrix& pmA);
    static TRVector  makeBlitz      (const ETRVector& pmA);
    static TRMatrix  makeBlitz      (const ETRMatrix& pmA);

    static void      setEigen      (ETRVector& pmVector, const TIVector& pmA);
    static void      setEigen      (ETRMatrix& pmMatrix, const TIMatrix& pmA);
    static void      setEigen      (ETRVector& pmVector, const TRVector& pmA);
    static void      setEigen      (ETRMatrix& pmMatrix, const TRMatrix& pmA);

    static blitz::Range getAll      ()                      { return soRangeAll; }

    static void     startTimer      ()                      { soTime = time(NULL); }
    static time_t   getTimer        ()                      { if (0 == soTime) { startTimer(); return 0; } else { return time(NULL) - soTime; } }

    static int      stringToInt     (const string ps, const int pi = 0) { if (0 == ps.length()) { return pi; } int ri; ss.clear();  ss.str(ps); ss >> ri; return ri; }
    static TReal    stringToReal    (const string ps)       { TReal   rr; ss.clear(); ss.str(ps); ss >> rr; return rr; }
    static string   intToString     (const int pi, const int piSize = 0) { ss.clear(); ss.str(""); ss << setfill('0') << setw(piSize) << pi; return ss.str(); }
    static string   realToString    (const TReal pr)        { ss.clear(); ss.str(""); ss << pr; return ss.str(); }

    static string   getParameterString  (const vector<string>& paParameters, const string psName, const string psDefault = "");
    static bool     getParameterBoolean (const vector<string>& paParameters, const string psName, const bool pbDefault = false);
    static int      getParameterInteger (const vector<string>& paParameters, const string psName, const int piDefault = 0);
    static TReal    getParameterReal    (const vector<string>& paParameters, const string psName, const TReal prDefault = 0.0);

    static EFileType getFileType    (const string psType);

    private:
    static blitz::Range soRangeAll;
    static time_t       soTime;
    static stringstream ss;
    static map<string,EFileType>    smFileTypes;
};

//inline TPVector     operator&&      (const TBVector& pmA, const TBVector& pmB) { TPVector rmC; rmC.resize(pmA.rows()); for (int i = 0; i < pmA.rows(); i++) { rmC(i) = (pmA(i) && pmB(i)); } return rmC; }
//inline TPVector     operator&&      (const TPVector& pmA, const TPVector& pmB) { TPVector rmC; rmC.resize(pmA.rows()); for (int i = 0; i < pmA.rows(); i++) { rmC(i) = (pmA(i) && pmB(i)); } return rmC; }
//inline TPVector     operator||      (const TBVector& pmA, const TBVector& pmB) { TPVector rmC; rmC.resize(pmA.rows()); for (int i = 0; i < pmA.rows(); i++) { rmC(i) = (pmA(i) || pmB(i)); } return rmC; }
//inline TPVector     operator||      (const TPVector& pmA, const TPVector& pmB) { TPVector rmC; rmC.resize(pmA.rows()); for (int i = 0; i < pmA.rows(); i++) { rmC(i) = (pmA(i) || pmB(i)); } return rmC; }
//inline TPVector     cwiseEquals     (const TPVector& pmA, const int piB)        { TPVector rmC; rmC.resize(pmA.rows()); for (int i = 0; i < pmA.rows(); i++) { rmC(i) = (pmA(i) == piB); } return rmC; }
//inline TPVector     cwiseNotEquals  (const TSNPVector& pmA, const ESNPCode piB) { TPVector rmC; rmC.resize(pmA.rows()); for (int i = 0; i < pmA.rows(); i++) { rmC(i) = (pmA(i) != piB); } return rmC; }
//inline TPVector     cwiseNotEquals  (const TPVector& pmA, const int piB)        { TPVector rmC; rmC.resize(pmA.rows()); for (int i = 0; i < pmA.rows(); i++) { rmC(i) = (pmA(i) != piB); } return rmC; }
inline bool         operator==      (const THVector& pmA, const THVector& pmB)     { if (pmA.rows() != pmB.rows()) { return false; } for (int i = 0; i < pmA.rows(); i++) { if (pmA(i) != pmB(i)) { return false; } } return true; }
inline bool         operator==      (const TIVector& pmA, const TIVector& pmB)     { if (pmA.rows() != pmB.rows()) { return false; } for (int i = 0; i < pmA.rows(); i++) { if (pmA(i) != pmB(i)) { return false; } } return true; }
inline bool         operator==      (const TSNPVector& pmA, const TSNPVector& pmB) { if (pmA.rows() != pmB.rows()) { return false; } for (int i = 0; i < pmA.rows(); i++) { if (pmA(i) != pmB(i)) { return false; } } return true; }

#endif
