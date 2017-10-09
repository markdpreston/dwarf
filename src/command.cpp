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
#include "command.h"
#include "process.h"
#include "random.h"
#include "statistics.h"

CProcess* CCommand::soProcess;

bool keepRemove (CPopulation* poPopulation, const vector<string>& paParameters, bool pbKeep);

CCommand::CCommand () {
    out       = NULL;
    msCommand = "";
    msHelp    = "";
    clean();
}

CCommand::~CCommand () {
    msCommand = "";
    msHelp    = "";
    clean();
}

void CCommand::clean () {
    miSubjects = 0;
    miSNPs     = 0;
    maParameters.clear();
    CProcess::close(out);
    out = NULL;
}

bool CCommand::preRun(const vector<string>& paParameters) {
    unsigned int    i;
    bool            lbOverwrite = false;
    //  Get the common parameters, if they exist, and open an output file.
    maParameters = paParameters;
    mbTest = CUtility::getParameterString(paParameters,"test","false") != "false";
    msPopulation = CUtility::getParameterString(paParameters,"population","default");
    moPopulation = soProcess->getPopulation(msPopulation);
    mePopulation = moPopulation->getType();
    if (msFile != CUtility::getParameterString(paParameters,"file","null")) {
        CProcess::close(out);
        msFile = CUtility::getParameterString(paParameters,"file","null");
        lbOverwrite = CUtility::getParameterString(paParameters,"overwrite","false") != "false";
        out = CProcess::open(msFile, lbOverwrite);
        if (! out) { cerr << msFile << " NOT OPENED!" << endl; exit(0); }
    }

    //  Initialise the size variables.
    miSubjects   = moPopulation->countSubjects();
    miPedigrees  = moPopulation->countPedigrees();
    miSNPs       = moPopulation->countSNPs();
    miControls   = moPopulation->countPhenotypes(1);
    miCases      = moPopulation->countPhenotypes(2);
    mmPhenotypes.resize(miSubjects);
    mmPhenotypes = moPopulation->getPhenotypes();
    mmPedigrees.resize(miSubjects);
    for (i = 0; i < miSubjects; i++) {
        mmPedigrees(i) = moPopulation->getSubjectP(i)->getPedigreeNumber();
    }

    return true;
}

//  Test and clean object if necessary.
bool CCommand::postRun () {
    if (mbTest) { test(); }
    return true;
}

void CCommandClean::run(const vector<string>& paParameters) {
    preRun(paParameters);
    moPopulation->clean(true);
    postRun();
}

void CCommandClear::run(const vector<string>& paParameters) {
    preRun(paParameters);
    moPopulation->clear();
    postRun();
}

void CCommandEcho::run(const vector<string>& paParameters) {
    unsigned int i;
    for (i = 0; i < paParameters.size(); i++) { cout << paParameters[i] << " "; } cout << endl;
}

void CCommandHelp::run(const vector<string>& paParameters) {
    unsigned int  i;
    CCommand     *loCommand;
    if (paParameters.size() == 0) {
        help();
        for (i = 0; i < soProcess->countCommands(); i++) {
            cout << "  " << soProcess->getCommand(i)->msCommand << endl;
        }
    } else {
        loCommand = soProcess->getCommand(paParameters[0]);
        if (loCommand) {
           cout << "Help: " << paParameters[0] << endl << loCommand->msHelp << endl;
        } else {
           cout << "Help: " << paParameters[0] << ": command not recognised" << endl;
        }
    }
}

void CCommandKeep::run(const vector<string>& paParameters) {
    preRun(paParameters);
    if (! keepRemove(moPopulation, paParameters, true)) { help(); }
    postRun();
}

//void CCommandLoad::run(const vector<string>& paParameters) {
//    preRun(paParameters);
//    paParameters.size() >= 2 ? moPopulation->load(paParameters[0], paParameters[1]) : help();
//    postRun();
//}

void CCommandRemove::run(const vector<string>& paParameters) {
    preRun(paParameters);
    if (! keepRemove(moPopulation, paParameters, false)) { help(); }
    postRun();
//    keepRemove(poPopulation, paParameters, false) == true ? 0 : help();
}

void CCommandRun::run(const vector<string>& paParameters) {
    paParameters.size() >= 1 ? soProcess->run(paParameters[0]) : help();
}

void CCommandSeed::run(const vector<string>& paParameters) {
    paParameters.size() >= 1 ? CRandom::seed(CUtility::stringToInt(paParameters[0])) : help();
}

void CCommandSystem::run(const vector<string>& psCommand) {
    unsigned int    i;
    int             liReturn;
    stringstream    lsCommand;
    for (i = 0; i < psCommand.size(); i++) {
        lsCommand << psCommand[i] << " ";
    }
    liReturn = std::system(lsCommand.str().c_str());
    if (0 != liReturn) {
        cerr << "system: " << lsCommand.str() << " : returned " << liReturn << endl;
    }
}

void CCommandTick::run(const vector<string>& paParameters) {
    int     liDivisor, liLoop, liTotal;
    liTotal   = CUtility::stringToInt(paParameters[0]);
    liLoop    = CUtility::stringToInt(paParameters[1]);
    liDivisor = CUtility::stringToInt(paParameters[2]);
    if (liLoop % liDivisor == 0) {
        cout << "tick      " << right << setw(10) << liTotal << setw(10) << liLoop << setw(10) << CUtility::getTimer() << left << endl;
    }
}

void CCommandTime::run(const vector<string>& paParameters) {
    cout << "time      " << right << setw(30) << CUtility::getTimer() << left << endl;
}

bool keepRemove (CPopulation* poPopulation, const vector<string>& paParameters, bool pbKeep) {
    bool            lbSubjects;
    unsigned int    i;
    TUVector        lmIndices;
    lmIndices.resize(paParameters.size() - 1);
    if (paParameters.size() >= 2) {
        lbSubjects = ("subjects" == paParameters[0]);
        for (i = 0; i < paParameters.size() - 1; i++) {
            lmIndices(i) = CUtility::stringToInt(paParameters[i+1]);
        }
        poPopulation->keepRemove(pbKeep, lbSubjects, lmIndices);
        return true;
    } else {
        return false;
    }
}

