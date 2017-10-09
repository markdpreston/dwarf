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
#include "process.h"
#include "subject.h"
#include "statistics.h"

bool      CProcess::mbQuiet = false;

CProcess::CProcess () {
    initialise();
    CCommand::soProcess = this;
}

CProcess::~CProcess () {
    while (! maPopulations.empty()) {
        delete maPopulations.begin()->second;
        maPopulations.erase(maPopulations.begin());
    }
    while (! maData.empty()) {
        delete maData.begin()->second;
        maData.erase(maData.begin());
    }
    while (! maCommands.empty()) {
        delete *(maCommands.begin());
        maCommands.erase(maCommands.begin());
    }
    while (! maVariables.empty()) {
        maVariables.erase(maVariables.begin());
    }
    maConstants.clear();
}

void CProcess::run (const string psScriptFile) {
    int              i, liLoop, liLoopStart[10], liLoopStop[10], liLoopWidth[10];
    bool             lbQuit;
    string           lsLine, lsLoopLine, lsCommand, lsLoopVariable[10];
    stringstream     lsLoop;
    size_t           liComment;
    vector<string>   lsFile, lsCommands[10], lsWords;
    map<string,string>::iterator liConstant;
    ifstream         fin;
    CCommand        *loCommand;
    //  Clear the loop variables.
    liLoop = 0;
    for (i = 0; i < 10; i++) {
        liLoopStart[i]    = 0;
        liLoopStop[i]     = 0;
        lsLoopVariable[i] = "";
    }
    //  Open the command file.
    if ("" != psScriptFile) {
        fin.open(psScriptFile.c_str());
        if (!fin) {
            cerr << "Bad file: " << psScriptFile << endl;
            exit(0);
        }
    }
    //  Parse and run the command lines.
    lbQuit = false;
    do {
        //  Read in line
        if ("" != psScriptFile) {
            getline(fin, lsLine);
        } else {
            cout << "> "; getline(cin, lsLine);
        }
        //  Substitute any constants in.
        for (liConstant = maConstants.begin(); liConstant != maConstants.end(); liConstant++) {
            boost::algorithm::replace_all(lsLine,liConstant->first,liConstant->second);
        }
        //  Remove comments, trim and split into tokens/words.
        liComment = lsLine.find_first_of("#");
        if (string::npos != liComment) { lsLine.erase(liComment); }
        boost::algorithm::trim(lsLine);
        boost::split(lsWords, lsLine, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
        if (0 == lsWords.size()) {
            continue;
        }
        //  Process language constructs: constant, for and next.
        if ("constant" == lsWords[0]) {
            maConstants["%" + lsWords[1] + "%"] = lsWords[2];
        } else if ("for" == lsWords[0]) {
            if (3 < lsWords.size()) {
                liLoop++;
                lsLoopVariable[liLoop] = lsWords[1];
                liLoopStart[liLoop]    = CUtility::stringToInt(lsWords[2]);
                liLoopStop[liLoop]     = CUtility::stringToInt(lsWords[3]);
                if (4 < lsWords.size()) {
                    liLoopWidth[liLoop] = CUtility::stringToInt(lsWords[4]);
                 } else {
                    liLoopWidth[liLoop] = lsWords[3].length();
                 }
            } else {
                cerr << "Bad for loop: " << lsLine << endl;
                exit(1);
            }
        } else if ("next" == lsWords[0]) {
            if (liLoopStart[liLoop] <= liLoopStop[liLoop]) {
                for (i = liLoopStart[liLoop]; i <= liLoopStop[liLoop]; i++) {
                    foreach(lsLoopLine, lsCommands[liLoop]) {
                        lsLoop.str("");
                        lsLoop << setfill('0') << setw(liLoopWidth[liLoop]) << i;
                        boost::algorithm::replace_all(lsLoopLine,"%" + lsLoopVariable[liLoop] + "%",lsLoop.str());
                        lsCommands[liLoop-1].push_back(lsLoopLine);
                    }
                }
            } else {
                for (i = liLoopStart[liLoop]; i >= liLoopStop[liLoop]; i--) {
                    foreach(lsLoopLine, lsCommands[liLoop]) {
                        lsLoop.str("");
                        lsLoop << setfill('0') << setw(liLoopWidth[liLoop]) << i;
                        boost::algorithm::replace_all(lsLoopLine,"%" + lsLoopVariable[liLoop] + "%",lsLoop.str());
                        lsCommands[liLoop-1].push_back(lsLoopLine);
                    }
                }
            }
            lsCommands[liLoop].clear();
            liLoopStart[liLoop]    = 0;
            liLoopStop[liLoop]     = 0;
            liLoopWidth[liLoop]    = 0;
            lsLoopVariable[liLoop] = "";
            liLoop--;
        } else {
            lsCommands[liLoop].push_back(lsLine);
        }
        //  Run the commands if all loops have been unrolled.
        if (liLoop == 0) {
            foreach(lsLine, lsCommands[0]) {
                //  Process line.
                boost::split(lsWords, lsLine, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                if (lsWords[0].size() > 0) {
                    //  Exit the script/console.
                    if ("quit" == lsWords[0] || "exit" == lsWords[0]) { lbQuit = true; break; }
                    //  Script output.
                    if ("quiet" == lsWords[0])   { mbQuiet = true;  continue; }
                    if (! mbQuiet)               { cout << lsLine << endl; }
                    if ("verbose" == lsWords[0]) { mbQuiet = false; continue; }
                    if ("script" == lsWords[0])  { run(lsWords[1]); continue; }
                    //  Process the command.
                    lsCommand = lsWords[0];
                    lsWords.erase(lsWords.begin());
                    loCommand = getCommand(lsCommand);
                    if (loCommand) {
                        loCommand->run(lsWords);
                    } else {
                        cerr << "Command unrecognised: " << lsCommand << endl;
                    }
                }
                lsWords.clear();
            }
            lsCommands[0].clear();
        }
        if (! lbQuit) {
            lbQuit = ("" != psScriptFile) ? fin.eof() : cin.eof();
        }
    } while (! lbQuit);
    if ("" != psScriptFile) {
        fin.close();
    }
}

CPopulation* CProcess::getPopulation (const string psPopulation) {
    CPopulation*    roPopulation;
    if (0 == maPopulations.count(psPopulation)) {
        roPopulation = new CPopulation();
        maPopulations.insert(pair<string,CPopulation*>(psPopulation,roPopulation));
    } else {
        roPopulation = maPopulations.find(psPopulation)->second;
    }
    return roPopulation;
}

TRMatrix* CProcess::getData (const string psData) {
    TRMatrix*    rmData;
    if (0 == maData.count(psData)) {
        rmData = new TRMatrix();
        maData.insert(pair<string,TRMatrix*>(psData,rmData));
    } else {
        rmData = maData.find(psData)->second;
    }
    return rmData;
}

vector<string> CProcess::getDataNames () {
    map<string,TRMatrix*>::iterator i;
    vector<string>   rmNames;
    for (i = maData.begin(); i != maData.end(); i++) {
        rmNames.push_back(i->first);
    }
    return rmNames;
}

unsigned int CProcess::countCommands() {
    return maCommands.size();
}

CCommand* CProcess::getCommand(const unsigned int piCommand) {
    return piCommand < maCommands.size() ? maCommands[piCommand] : NULL;
}


CCommand* CProcess::getCommand (const string psCommand) {
    unsigned int    i;
    for (i = 0; i < maCommands.size(); i++) {
        if (psCommand == maCommands[i]->msCommand) { return maCommands[i]; }
    }
    return NULL;
}

ostream* CProcess::open (const string psFile, const bool pbOverwrite) {
    if ("" == psFile || "cout" == psFile) {
        cout << setprecision(4);
        return &cout;
    } else if ("null" == psFile) {
        return new ofstream("/dev/null");
    } else {
        ostream *out = new ofstream(psFile.c_str(),pbOverwrite ? ios_base::out : ios_base::app);
        (*out) << setprecision(4);
        return out;
    }
}

void CProcess::close (ostream* out) {
    if (NULL != out && &cout != out) {
        ((ofstream*) out)->close();
        delete out;
    }
}

void CProcess::initialise () {
    unsigned int    i;

    miR.rcpp    = 0;
    miR.rinside = 0;
    miR.mvtnorm = 0;
    miR.skat    = 0;
    miR.kbac    = 0;

    maCommands.push_back(new CCommandClean());
    maCommands.push_back(new CCommandClear());
    maCommands.push_back(new CCommandData());
    maCommands.push_back(new CCommandEcho());
    maCommands.push_back(new CCommandEnrich());
    maCommands.push_back(new CCommandFilter());
    maCommands.push_back(new CCommandHelp());
    maCommands.push_back(new CCommandInformation());
    maCommands.push_back(new CCommandKeep());
    maCommands.push_back(new CCommandLoad());
    maCommands.push_back(new CCommandMerge());
    maCommands.push_back(new CCommandPopulation());
    maCommands.push_back(new CCommandPower());
    maCommands.push_back(new CCommandPseudo());
    maCommands.push_back(new CCommandR());
    maCommands.push_back(new CCommandRemove());
    maCommands.push_back(new CCommandRun());
    maCommands.push_back(new CCommandSave());
    maCommands.push_back(new CCommandSeed());
    maCommands.push_back(new CCommandSimulation());
    maCommands.push_back(new CCommandSystem());
    maCommands.push_back(new CCommandTick());
    maCommands.push_back(new CCommandTime());
    maCommands.push_back(new CCommandTransform());
    maCommands.push_back(new CAnalysisAssociation());
    maCommands.push_back(new CAnalysisCAlpha());
    maCommands.push_back(new CAnalysisKBAC());
    maCommands.push_back(new CAnalysisRegression());
    maCommands.push_back(new CAnalysisScore());
    maCommands.push_back(new CAnalysisSKAT());
    maCommands.push_back(new CAnalysisSingle());
    maCommands.push_back(new CAnalysisTDT());

    for (i = 0; i < maCommands.size(); i++) {
        maCommands[i]->initialise();
    }
}
