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
#ifndef COMMAND_SIMULATION_H
#define COMMAND_SIMULATION_H

#include "types.h"
#include "population.h"
#include "command.h"

typedef struct {
    unsigned int         miCount;
    unsigned int         miSiblings;
    unsigned int         miParents;
    unsigned int         miUnrelateds;
    unsigned int         miUnitSizeFalse;
    unsigned int         miUnitSizeTrue;

    unsigned int         miUnknownSiblings;
    unsigned int         miUnknownParents;
    unsigned int         miUnknownUnrelateds;
    unsigned int         miCaseSiblings;
    unsigned int         miCaseParents;
    unsigned int         miCaseUnrelateds;
    unsigned int         miControlSiblings;
    unsigned int         miControlParents;
    unsigned int         miControlUnrelateds;
} TUnits;

class CCommandSimulation  : public CCommand {
    public:
    void            initialise          ();
    void            run                 (const vector<string>& paParameters);
    static TUnits   getUnit             (const vector<string>& paParameters);

    private:
    CPopulation    *moPopulationSample;
    void            getSubjects         (const TIVector& pmPhenotypes, TIVector& pmSubjects, TIVector& pmPedigrees, const int piStatus = 0);
    void            snps                ();
    void            affection           ();
    void            haplotypes          ();
    void            subjects            (const TUnits poUnit);
    void            replace             ();
    bool            parentsOk           (const TUnits poUnit, const unsigned int piCases);
    bool            siblingsOk          (const TUnits poUnit, const unsigned int piCases);
    bool            unrelatedsOk        (const TUnits poUnit, const unsigned int piCases);
};

#endif
