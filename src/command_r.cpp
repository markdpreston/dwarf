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
#include "interface_r.h"
#include "process.h"

void CCommandR::run(const vector<string>& paParameters) {
    unsigned int    i;
    TRMatrix        lmDosages;
    string          lsLine, lsName;
    stringstream    lsCommand;
    ifstream        fin;
#ifdef INCLUDE_R
    if ("command" == paParameters[0]) {
        for (i = 1; i < paParameters.size(); i++) {
            lsCommand << paParameters[i] << " ";
        }
        CInterfaceR::runCommand(lsCommand.str());
    }
    if ("script" == paParameters[0]) {
        lsCommand << "source(\"" << paParameters[1] << "\")";
        CInterfaceR::runCommand(lsCommand.str());
    }
    if ("load" == paParameters[0]) {
        if (miSubjects == 0) {
            cerr << "No dosage or phenotypes loaded into R" << endl;
            return;
        }
        CInterfaceR::setVector("phenotypes", CUtility::makeReal(mmPhenotypes));
        if (miSNPs == 0) {
            cerr << "No dosage loaded into R" << endl;
            return;
        }
        lmDosages.resize(miSubjects, miSNPs);
        switch (mePopulation) {
            case cePopulationGenotypes:  lmDosages = moPopulation->getGenotypes()->getDosage(false); break;
            case cePopulationHaplotypes: lmDosages = moPopulation->getHaplotypes()->getDosage();     break;
            case cePopulationDosages:    lmDosages = moPopulation->getDosages()->getDosage();        break;
            default: cerr << "No dosage loaded into R" << endl; return;
        }
        CInterfaceR::setMatrix("dosage", lmDosages);
        CInterfaceR::runCommand(lsCommand.str());
    }
    if ("store" == paParameters[0]) {
        TRMatrix *lmA;
        if (1 == paParameters.size()) {
            cerr << "No name for storage" << endl;
            return;
        }
        lmA = soProcess->getData(paParameters[1]);
        CInterfaceR::getMatrix(paParameters[1],*lmA);
    }
    if ("output" == paParameters[0]) {
        TRMatrix lmA;
        if (1 < paParameters.size()) {
            CInterfaceR::getMatrix(paParameters[1],lmA);
            cout << paParameters[1] << ": ";
        } else {
            CInterfaceR::getMatrix("",lmA);
            cout << "Last answer: ";
        }
        cout << lmA << endl;
    }
#endif
}
