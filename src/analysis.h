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
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "types.h"
#include "command.h"

class CAnalysis : public CCommand {
    public:
                    CAnalysis           ();
    virtual        ~CAnalysis           ();
    virtual void    clean               ();
    virtual void    run                 (const vector<string>& paParameters);
    virtual void    statistics          (TRVector& pmStatistics, const bool pbShuffle = false) = 0;
    virtual void    asymptotics         (TRVector& pmAsymptotics) = 0;
    protected:
    unsigned int    cluster             (TRVector& pmDosage, TIVector& pmPhenotypes);
    unsigned int    cluster             (TRMatrix& pmDosage, TIVector& pmPhenotypes);
    unsigned int    cluster             (ETRMatrix& pmDosage, TIVector& pmPhenotypes);
    string          msData;
    TRMatrix*       mmData;
    unsigned int    miPermutations;
    TRMatrix        mmPermutations;
    int             miDataEntry;
    unsigned int    miStatistics;
    bool            mbCluster;
    TIVector        mmClusters;
    bool            mbSNPs;
    TRVector        mmStatistics;
    int             miDataFields;
    int             miDataSize;
    TRVector        mmDataLine;
    TRVector        mmPermutationLine;
    TReal           mrAlpha;
    virtual bool    preRun              (const vector<string>& paParameters);
    virtual bool    postRun             (const bool pbSuppress = false);
};

#include "analysis_association.h"
#include "analysis_calpha.h"
#include "analysis_kbac.h"
#include "analysis_regression.h"
#include "analysis_score.h"
#include "analysis_single.h"
#include "analysis_skat.h"
#include "analysis_tdt.h"

#endif
