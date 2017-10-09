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
#ifndef COMMAND_H
#define COMMAND_H

#include "types.h"
#include "population.h"
#include "interface_r.h"

class CProcess;

class CCommand {
    public:
    static CProcess *soProcess;
                    CCommand            ();
    virtual        ~CCommand            ();
    virtual void    initialise          () = 0;
    virtual void    run                 (const vector<string>& paParameters) = 0;
    virtual void    clean               ();
    void            setPopulation       (CPopulation* poPopulation) { moPopulation = poPopulation; }
    virtual void    test                () { cout << "Test " << setw(19) << msCommand << "NA - no test defined" << endl; }
    void            help                () { string lsParameter; cout << "Help: " << msCommand << ": " << msHelp << endl << msCommand; foreach(lsParameter, maParameters) { cout << " " << lsParameter; } cout << endl;}
    void            log                 () { cout << "Log: " << msCommand << endl << msLog.str() << endl; msLog.str(""); }
    void            warning             () { cerr << "Warning: " << msCommand << endl << msWarning.str() << endl; msWarning.str(""); }
    void            error               () { cerr << "Error: " << msCommand << endl << msError.str() << endl; msError.str(""); exit(0); }
    string          msCommand;
    string          msHelp;
    protected:
    vector<string>  maParameters;
    bool            mbTest;
    string          msPopulation;
    CPopulation*    moPopulation;
    EPopulationType mePopulation;
    unsigned int    miSubjects;
    unsigned int    miPedigrees;
    unsigned int    miSNPs;
    unsigned int    miCases;
    unsigned int    miControls;
    TIVector        mmPhenotypes;
    TIVector        mmPedigrees;
    string          msFile;
    stringstream    msLog;
    stringstream    msWarning;
    stringstream    msError;
    ostream*        out;
    virtual bool    preRun              (const vector<string>& paParameters);
    virtual bool    postRun             ();
};

class CCommandClean       : public CCommand { public: void initialise() { msCommand = "clean"; }     void run(const vector<string>& paParameters); };
class CCommandClear       : public CCommand { public: void initialise() { msCommand = "clear"; }     void run(const vector<string>& paParameters); };
class CCommandCluster     : public CCommand { public: void initialise() { msCommand = "cluster"; }   void run(const vector<string>& paParameters); };
class CCommandEcho        : public CCommand { public: void initialise() { msCommand = "echo"; }      void run(const vector<string>& paParameters); };
class CCommandEnrich      : public CCommand { public: void initialise() { msCommand = "enrich"; }    void run(const vector<string>& paParameters); };
class CCommandFilter      : public CCommand { public: void initialise() { msCommand = "filter"; }    void run(const vector<string>& paParameters); };
class CCommandHelp        : public CCommand { public: void initialise() { msCommand = "help"; }      void run(const vector<string>& paParameters); };
class CCommandKeep        : public CCommand { public: void initialise() { msCommand = "keep"; }      void run(const vector<string>& paParameters); };
class CCommandMerge       : public CCommand { public: void initialise() { msCommand = "merge"; }     void run(const vector<string>& paParameters); };
class CCommandPseudo      : public CCommand { public: void initialise() { msCommand = "pseudo"; }    void run(const vector<string>& paParameters); };
class CCommandR           : public CCommand { public: void initialise() { msCommand = "r"; }         void run(const vector<string>& paParameters); };
class CCommandRemove      : public CCommand { public: void initialise() { msCommand = "remove"; }    void run(const vector<string>& paParameters); };
class CCommandRun         : public CCommand { public: void initialise() { msCommand = "run"; }       void run(const vector<string>& paParameters); };
class CCommandSeed        : public CCommand { public: void initialise() { msCommand = "seed"; }      void run(const vector<string>& paParameters); };
class CCommandSystem      : public CCommand { public: void initialise() { msCommand = "system"; }    void run(const vector<string>& paParameters); };
class CCommandTick        : public CCommand { public: void initialise() { msCommand = "tick"; }      void run(const vector<string>& paParameters); };
class CCommandTime        : public CCommand { public: void initialise() { msCommand = "time"; }      void run(const vector<string>& paParameters); };
class CCommandTransform   : public CCommand { public: void initialise() { msCommand = "transform"; } void run(const vector<string>& paParameters); };

#include "command_data.h"
#include "command_load.h"
#include "command_information.h"
#include "command_population.h"
#include "command_power.h"
#include "command_save.h"
#include "command_simulation.h"

#endif
