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
#ifndef ANALYSIS_SCORE_H
#define ANALYSIS_SCORE_H

#include "types.h"
#include "analysis.h"

typedef struct {
    int     miB;
    int     miC;
} THetero;

class CAnalysisScore : public CAnalysis {
    public:
                    CAnalysisScore      ();
                   ~CAnalysisScore      ();
    void            initialise          () { msCommand = "score"; miStatistics = 5; mbSNPs = false; }
    void            clean               ();
    virtual bool    preRun              (const vector<string>& paParameters);
    void            statistics          (TRVector& rmStatistics, const bool pbShuffle);
    void            asymptotics         (TRVector& rmAsymptotics);
    void            test                ();

    private:
    void            calculation         (TRVector& rmStatistics, TRVector& rmAsymptotics);
    void            pvalues             (TRVector& rmStatistics, TRVector& rmAsymptotics, const int piSNPs, ETRVector& pmU, ETRMatrix& pmV, ETRVector& pmP, TRMatrix& pmCorrelation, TRMatrix& pmVBlitz, ETRVector& pmEvals, ETRMatrix& pmTemp);
    TResult         chiSquared          (const ETRVector& pmEvals, const TReal prScore);
    THetero         score               (ESNPCode peSNP, ESNPCode peFather, ESNPCode peMother);
    //  General computational variables.
    string          msMethod;
    TRVector        mmZeros;
    TRMatrix        mmCorrelation;
    TRMatrix        mmVBlitz;
    ETRVector       mmU;
    ETRVector       mmEvals;
    ETRVector       mmP;
    ETRMatrix       mmV;
    ETRMatrix       mmTemp;
    Eigen::SelfAdjointEigenSolver<ETRMatrix> *moSolver1;
    Eigen::ComplexEigenSolver<ETRMatrix>     *moSolver2;
    //  Unrelated computational variables.
    TRVector        mmClusters;
    ETRVector       mmY;
    ETRMatrix       mmX;
    ETRMatrix       mmMeanX;
    //  Family computational variables.
    unsigned int    miAffected;
    TSubjectPVector mmAffected;
    ETRMatrix       mmBminusC;
    ETRMatrix       mmMean;
};

#endif
