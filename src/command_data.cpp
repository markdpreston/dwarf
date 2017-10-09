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
#include "command.h"

void CCommandData::run(const vector<string>& paParameters) {
    unsigned int    i;
    int             liCount, liEvery;
    string          lsCommand, lsVariables;
    vector<string>  laVariables;

    preRun(paParameters);

    lsCommand   = paParameters[0];
    lsVariables = CUtility::getParameterString(paParameters,"variables","all");
    liCount = CUtility::getParameterInteger(paParameters,"count",1);
    liEvery = CUtility::getParameterInteger(paParameters,"every",1);
    if (liCount % liEvery == 0) {
        if ("all" == lsVariables) {
            laVariables = soProcess->getDataNames();
        } else {
            boost::algorithm::split(laVariables, lsVariables, boost::is_any_of(","), boost::algorithm::token_compress_on);
        }

        if ("save" == lsCommand) {
            for (i = 0; i < laVariables.size(); i++) {
                if (soProcess->isData(laVariables[i])) {
                    save(out, laVariables[i], soProcess->getData(laVariables[i]));
                }
            }
        }
    }

    postRun();
}

void CCommandData::save (ostream* out, const string psName, const TRMatrix* pmData) {
    int     i, j;
    if (pmData->size()) {
        (*out) << left;
        (*out) << "# name: " << psName << endl;
        (*out) << "# type: matrix" << endl;
        (*out) << "# rows: " << pmData->rows() << endl;
        (*out) << "# columns: " << pmData->cols() << endl;
        for (i = 0; i < pmData->rows(); i++) {
            for (j = 0; j < pmData->cols(); j++) {
                (*out) << setw(15) << (*pmData)(i,j);
            }
            (*out) << endl;
        }
    }
}
