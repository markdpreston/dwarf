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
#ifndef COMMAND_LOAD_H
#define COMMAND_LOAD_H

#include "types.h"
#include "command.h"

class CCommandLoad : public CCommand {
    public:
    void            initialise          () { msCommand = "load"; }
    void            run                 (const vector<string>& paParameters);
    void            run                 (const string psType, const string psFile, const unsigned int piSubjects = 0, const unsigned int piSNPs = 0, const bool pbHaplotypes = false);

    private:
    void            loadMarkers         (ifstream& fin, const EFileType peType);
    void            loadPed             (ifstream& fin, const EFileType peType);
    void            loadTabular         (ifstream& fin);
    void            loadBed             (ifstream& fin, const unsigned int piSubjects = 0, const unsigned int piSNPs = 0, const bool pbHaplotypes = false);
    void            loadGenes           (ifstream& fin);
    void            loadVCF             (ifstream& fin);
    void            loadMito            (ifstream& fin);
    void            loadMS              (ifstream& fin);
    void            getFileSize         (ifstream& fin, int& piLines, int& piFields, const int piOffset = 0);
    //  Internal data structures
    map<string,EFileType>    mmFileTypes;
};

#endif
