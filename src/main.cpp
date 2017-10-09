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
#include "interface_r.h"
#include "process.h"
#include "random.h"
#include "data.h"

int main (int argc, char *argv[]) {
    int         i;
    CProcess    loProcess;

#ifdef INCLUDE_R
    loProcess.miR.rcpp    = CInterfaceR::runCommandI("require(Rcpp)");
    loProcess.miR.rinside = CInterfaceR::runCommandI("require(RInside)");
    loProcess.miR.mvtnorm = CInterfaceR::runCommandI("require(mvtnorm)");
    loProcess.miR.skat    = CInterfaceR::runCommandI("require(SKAT)");
    loProcess.miR.kbac    = CInterfaceR::runCommandI("require(KBAC)");
    if (0 == loProcess.miR.rcpp)    { cerr << "R: missing library \"Rcpp\".     Lots of functionality requires this library." << endl; }
    if (0 == loProcess.miR.rinside) { cerr << "R: missing library \"RInside\".  Lots of functionality requires this library." << endl; }
    if (0 == loProcess.miR.mvtnorm) { cerr << "R: missing library \"mvtnorm\".  Please install for full functionality." << endl; }
    if (0 == loProcess.miR.skat)    { cerr << "R: missing library \"SKAT\".     Please install for full functionality." << endl; }
    if (0 == loProcess.miR.kbac)    { cerr << "R: missing library \"KBAC\".     Please install for full functionality." << endl; }
#else
    cerr << "Non-R version.  Please install R and use appropriate version for full functionality." << endl;
#endif
    CRandom::seed();
    CUtility::startTimer();
    CUtility::initialise();

    if (1 == argc) {
        loProcess.run("");
    } else {
        for (i = 1; i < argc; i++) {
            loProcess.run(argv[i]);
        }
    }

    return 1;
}
