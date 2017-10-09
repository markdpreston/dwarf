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
#ifndef PROCESS_H
#define PROCESS_H

#include "types.h"
#include "command.h"
#include "analysis.h"
#include "population.h"

typedef struct {
    int     rcpp;
    int     rinside;
    int     mvtnorm;
    int     skat;
    int     kbac;
} TRLibraries;

typedef struct {
    string          msVariable;
    string          msValue;
    char            mcFill;
    unsigned int    miFill;
} TVariable;

class CProcess {
    public:
                    CProcess            ();
                   ~CProcess            ();
    void            run                 (const string psScriptFile);
    //  Command control.
    unsigned int    countCommands       ();
    CCommand*       getCommand          (const unsigned int piCommand);
    CCommand*       getCommand          (const string psCommand);
    CPopulation*    getPopulation       (const string psPopulation);
    TRMatrix*       getData             (const string psData);
    vector<string>  getDataNames        ();
//    bool            isCommand           (const string psCommand)    { return 0 != maCommand.count(psCommand); }
    bool            isPopulation        (const string psPopulation) { return 0 != maPopulations.count(psPopulation); }
    bool            isData              (const string psData)       { return 0 != maData.count(psData); }
//    void            load                (CCommand* poCommand);
    //  File handling.
    static ostream* open                (const string psFile = "null", const bool pbOverwrite = false);
    static void     close               (ostream* out);
    TRLibraries     miR;
    private:
    void            initialise          ();
    static bool     mbQuiet;
    map<string,CPopulation*>    maPopulations;
    map<string,TRMatrix*>       maData;
    map<string,TVariable*>      maVariables;
    vector<CCommand*>           maCommands;
    map<string,string>          maConstants;
};

#endif
