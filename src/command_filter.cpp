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
#include "polydata.h"

void CCommandFilter::run(const vector<string>& paParameters) {
    unsigned int    i, j, liSNPs, liSubjects, liSNPCount, laCount[4];
    unsigned char   lcSNP;
    string          lsPopulation;
    TUVector        lmFilterSNPs, lmSNPs;
    EPopulationType lePopulation;
    CSnip           loSNP;
    CPopulation*    loPopulation;

    preRun(paParameters);

    lsPopulation = paParameters[0];
    if (! soProcess->isPopulation(lsPopulation)) {
        msError << "No population to filter from: " << lsPopulation;
        error();
    }
    loPopulation = soProcess->getPopulation(lsPopulation);
    lePopulation = loPopulation->getType();
    liSubjects = loPopulation->countSubjects();
    liSNPs = loPopulation->countSNPs();
    lmSNPs.resize(liSubjects);
    if (lePopulation != cePopulationPolyData) {
        msError << "Bad population type: " << lePopulation;
        error();
    }

    //  Count SNPs that pass the filter.
    lmFilterSNPs.resize(liSNPs);
    lmFilterSNPs = 0;
    liSNPCount = 0;
    for (i = 0; i < liSNPs; i++) {
        lmSNPs = loPopulation->getPolyData()->getSNP(i);
        laCount[0] = 0;
        laCount[1] = 0;
        laCount[2] = 0;
        laCount[3] = 0;
        for (j = 0; j < liSubjects; j++) {
            lcSNP = CPolyData::getQuadData(lmSNPs(j));
            laCount[0] += ((lcSNP & 128) ? 1 : 0);
            laCount[0] += ((lcSNP &   8) ? 1 : 0);
            laCount[1] += ((lcSNP &  64) ? 1 : 0);
            laCount[1] += ((lcSNP &   4) ? 1 : 0);
            laCount[2] += ((lcSNP &  32) ? 1 : 0);
            laCount[2] += ((lcSNP &   2) ? 1 : 0);
            laCount[3] += ((lcSNP &  16) ? 1 : 0);
            laCount[3] += ((lcSNP &   1) ? 1 : 0);
        }
        //  Check to see if there are more than 3 alleles?
        if ((laCount[0] > 0 && laCount[1] > 0 && laCount[2] > 0) ||
            (laCount[0] > 0 && laCount[1] > 0 && laCount[3] > 0) ||
            (laCount[0] > 0 && laCount[2] > 0 && laCount[3] > 0) ||
            (laCount[1] > 0 && laCount[2] > 0 && laCount[3] > 0)) {
            loSNP = loPopulation->getSNPs()->getSNP(i);
            cout << setw(4) << loSNP.getChromosome() << setw(20) << loSNP.getName() << setw(20) << loSNP.getDistance() << setw(14) << loSNP.getPosition() << setw(10) << laCount[0] << setw(10) << laCount[1] << setw(10) << laCount[2] << setw(10) << laCount[3] << endl;
        } else {
            lmFilterSNPs(liSNPCount) = i;
            liSNPCount++;
        }
    }
    lmFilterSNPs.resizeAndPreserve(liSNPCount);

    //
    moPopulation->resizeSubjects(liSubjects);
    moPopulation->resizeSNPs(liSNPCount);
    moPopulation->setType(cePopulationPolyData);
    moPopulation->resizeData();

    //  Copy only filtered SNPs.
    for (i = 0; i < liSNPCount; i++) {
        moPopulation->getSNPs()->setSNP(i,loPopulation->getSNPs()->getSNP(lmFilterSNPs(i)));
        lmSNPs = loPopulation->getPolyData()->getSNP(lmFilterSNPs(i));
        moPopulation->getPolyData()->setSNP(i,lmSNPs);
    }
    moPopulation->setPhenotypes(loPopulation->getPhenotypes());
    moPopulation->setTraits(loPopulation->getTraits());
    moPopulation->clean(false);

    postRun();
}
