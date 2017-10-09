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
#include "pedigree.h"
#include "utility.h"

void CCommandEnrich::run(const vector<string>& paParameters) {
    unsigned int    i, j, k, liCount;
    string          lsPopulation;
    TIVector        lmPhenotypes, lmPedigree;
    THVector        lmHaplotype;
    TSNPVector      lmGenotype;
    TSubjectVector  lmSubjects;
    CHaplotype      loHaplotypes;
    CPedigree       loPedigree;
    CSubject       *loSubject;
    CSNPs           loSNPs;
    CPopulation*    loPopulation;

    preRun(paParameters);

    //  The new population to create.
    lsPopulation = CUtility::getParameterString(paParameters, "to", "default");
    loPopulation = soProcess->getPopulation(lsPopulation);
    lmHaplotype.resize(miSNPs);

    //  Count the ASP families and create a new haplotype.
    liCount = 0;
    for (i = 0; i < miPedigrees; i++) {
        loPedigree = moPopulation->getPedigree(i);
        if (2 == loPedigree.countParents() && 2 == loPedigree.countAffectedSiblings()) {
            liCount++;
        }
    }

    //  Set up the population object and create dummy names.
    loHaplotypes.resize(3*liCount, miSNPs);
    lmPhenotypes.resize(3*liCount);
    lmPedigree.resize(3*liCount);
    lmSubjects.resize(3*liCount);
    lmPhenotypes = 0;
    lmPedigree = 0;

    //  Loop over ASP families.
    liCount = 0;

    for (i = 0; i < miPedigrees; i++) {
        loPedigree = moPopulation->getPedigree(i);
        if (2 == loPedigree.countParents() && 2 == loPedigree.countAffectedSiblings()) {
            for (j = 0; j < 3; j++) {
                switch (j) {
                    case 0: lmPhenotypes(liCount) = 1; loSubject = loPedigree.getParent(0); break;
                    case 1: lmPhenotypes(liCount) = 1; loSubject = loPedigree.getParent(1); break;
                    case 2: lmPhenotypes(liCount) = 2; loSubject = loPedigree.getAffectedSibling(0); break;
                }
                switch (moPopulation->getType()) {
                    case cePopulationGenotypes: {
                        lmGenotype = moPopulation->getGenotypes()->getSubject(loSubject->getIdNumber());
                        loHaplotypes.setSubject(liCount, lmGenotype);
                        break;
                    }
                    case cePopulationHaplotypes: {
                        for (k = 0; k < 2; k++) {
                            lmHaplotype = moPopulation->getHaplotypes()->getSubject(loSubject->getIdNumber(), k);
                            loHaplotypes.setSubject(liCount, k, lmHaplotype);
                        }
                        break;
                    }
                    default: cerr << "Danger Will Robinson" << endl; exit(0);
                }
                lmSubjects(liCount) = *loSubject;
                liCount++;
            }
        }
    }

    loSNPs.copySNPs(moPopulation->getSNPs());
    loPopulation->clear();
    loPopulation->setType(cePopulationHaplotypes);
    loPopulation->resizeSubjects(liCount);
    loPopulation->resizeSNPs(miSNPs);
    loPopulation->getSNPs()->copySNPs(&loSNPs);
    loPopulation->resizeData();
    loPopulation->setTraits(CUtility::makeReal(lmPhenotypes));
    loPopulation->setPhenotypes(lmPhenotypes);
    for (i = 0; i < liCount; i++) {
        loPopulation->setSubject(i, lmSubjects(i));
    }
    loPopulation->setHaplotypes(loHaplotypes);
    loPopulation->clean(true);

    postRun();
}
