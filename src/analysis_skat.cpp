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
#include "analysis_skat.h"
#include "process.h"
#include "interface_r.h"

bool CAnalysisSKAT::preRun(const vector<string>& paParameters) {
    bool    rbSuccess = false;
    rbSuccess = CAnalysis::preRun(paParameters);
    mmY.resize(miSubjects);
    mmZ.resize(miSubjects,miSNPs);
    if (0 == soProcess->miR.rinside || 0 == soProcess->miR.rcpp || 0 == soProcess->miR.skat) {
        msError << "R libraries required: RInside, Rcpp, SKAT";
        error();
    }

    return rbSuccess;
}

void CAnalysisSKAT::run(const vector<string>& paParameters) {
    int         liCount;

    if (! preRun(paParameters)) { return; }

#ifdef INCLUDE_R
    //  Expose the genotypes and phenotypes to R.
    switch (moPopulation->getType()) {
        case cePopulationDosages:    mmZ = moPopulation->getDosages()->getDosage(); break;
        case cePopulationGenotypes:  mmZ = moPopulation->getGenotypes()->getDosage(false); break;
        case cePopulationHaplotypes: mmZ = moPopulation->getHaplotypes()->getDosage(); break;
        default: msError << "Population type not implemented"; error();
    }

    if (mbCluster) {
        liCount = CAnalysis::cluster(mmZ, mmPhenotypes);
        mmY.resize(liCount);
    }

    mmY = mmPhenotypes - 1.0;
    CInterfaceR::setVector(string("skatY"), mmY);
    CInterfaceR::setMatrix(string("skatZ"), mmZ);
    CInterfaceR::runCommand("o <- SKAT_Null_Model(skatY ~ 1, out_type=\"D\")");
    mmDataLine(0) = CInterfaceR::runCommandR("SKAT(skatZ,o)$p.value");
#else
    mmDataLine(0) = 1.0;
#endif

    postRun();
}
