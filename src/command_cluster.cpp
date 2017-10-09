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
#include "haplotype.h"
#include "random.h"
#include "utility.h"

void CCommandCluster::run(const vector<string>& paParameters) {
    unsigned int    i, liCount, liCluster;
    blitz::firstIndex   loIndex1;
    string          lsPopulation;
    TIVector        lmPhenotypes, lmClusters;
    TRMatrix        lmDosages;
    CDosage         loDosages;
    CSNPs           loSNPs;
    CPopulation*    loPopulation;

    preRun(paParameters);

    //  The new population to create.
    lsPopulation = CUtility::getParameterString(paParameters, "to", "default");
    loPopulation = soProcess->getPopulation(lsPopulation);

    //  Prepare.
    lmDosages.resize(miSubjects, miSNPs);
    switch (mePopulation) {
        case cePopulationGenotypes:  lmDosages = moPopulation->getGenotypes()->getDosage(false); break;
        case cePopulationHaplotypes: lmDosages = moPopulation->getHaplotypes()->getDosage();     break;
        case cePopulationDosages:    lmDosages = moPopulation->getDosages()->getDosage();        break;
        default: msError << "Not supported population type"; error(); break;
    }
    lmClusters = loIndex1;

    //  Sum the dosage.
    for (i = 0; i < miSubjects; i++) {
        liCluster = first(mmPedigrees == mmPedigrees(i) && mmPhenotypes == mmPhenotypes(i));
        if (liCluster != i) {
            lmClusters(i) = -1;
            lmDosages(liCluster,CUtility::getAll()) += lmDosages(i,CUtility::getAll());
            lmDosages(i,CUtility::getAll()) = 0;
        }
    }

    //  Shrink the dosage and phenotypes.
    liCount = 0;
    for (i = 0; i < miSubjects; i++) {
        if (lmClusters(i) != -1) {
            lmDosages(liCount,CUtility::getAll()) = lmDosages(i,CUtility::getAll());
            lmPhenotypes(liCount) = lmPhenotypes(i);
            liCount++;
        }
    }
    lmDosages.resizeAndPreserve(liCount, miSNPs);
    lmPhenotypes.resizeAndPreserve(liCount);

    //  Set up the population object and create dummy names.
    loSNPs.copySNPs(moPopulation->getSNPs());
    loPopulation->clear();
    loPopulation->setType(cePopulationDosages);
    loPopulation->resizeSubjects(liCount);
    loPopulation->resizeSNPs(miSNPs);
    loPopulation->getSNPs()->copySNPs(&loSNPs);
    loPopulation->resizeData();
    loPopulation->setTraits(CUtility::makeReal(lmPhenotypes));
    loPopulation->setPhenotypes(lmPhenotypes);
    for (i = 0; i < liCount; i++) {
        loPopulation->getSubjectP(i)->setPedigree(CUtility::intToString(i,5));
        loPopulation->getSubjectP(i)->setPedigreeNumber(i);
    }


    loPopulation->setDosages(loDosages);


    loPopulation->clean(false);

    postRun();
}
