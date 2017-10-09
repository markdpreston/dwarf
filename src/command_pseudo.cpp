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

void CCommandPseudo::run(const vector<string>& paParameters) {
    unsigned int    i, liCount;
    string          lsPopulation;
    TIVector        lmPhenotypes, lmPedigree;
    THVector        lmCase[2], lmControl[2], lmParents[4];
    CHaplotype      loHaplotypes;
    CSubject        loSubject;
    CSNPs           loSNPs;
    CPopulation*    loPopulation;

    preRun(paParameters);

    //  The new population to create.
    lsPopulation = CUtility::getParameterString(paParameters, "to", "default");
    loPopulation = soProcess->getPopulation(lsPopulation);

    //  Prepare.
    lmCase[0].resize(miSNPs);
    lmCase[1].resize(miSNPs);
    lmControl[0].resize(miSNPs);
    lmControl[1].resize(miSNPs);
    lmParents[0].resize(miSNPs);
    lmParents[1].resize(miSNPs);
    lmParents[2].resize(miSNPs);
    lmParents[3].resize(miSNPs);

    //  Count the affected subjects with two parents?
    liCount = 0;
    for (i = 0; i < miSubjects; i++) {
//        if (moPopulation->getSubjectInfo(i).countParents() == 2 && moPopulation->getPhenotype(i) == 2) {
        if (moPopulation->getSubjectInfo(i).countParents() == 2) {
            liCount++;
        }
    }

    //  Create a new haplotype.
    loHaplotypes.resize(2*liCount, miSNPs);
    lmPhenotypes.resize(2*liCount);
    lmPhenotypes = 0;
    lmPedigree.resize(2*liCount);
    lmPedigree = 0;

    //  Loop over cases with two parents.
    liCount = 0;
    for (i = 0; i < miSubjects; i++) {
        loSubject = moPopulation->getSubjectInfo(i);
//        if (loSubject.countParents() == 2 && moPopulation->getPhenotype(i) == 2) {
        if (loSubject.countParents() == 2) {
            //  Get the case haplotypes, assume that there is no missing data in the genotype data.
            switch (mePopulation) {
//                case cePopulationGenotypes: {
//                    lmCase = moPopulation->getGenotypes()->getSubject(loSubject.getIdNumber());
//                    break;
//                }
                case cePopulationHaplotypes: {
                    lmCase[0] = moPopulation->getHaplotypes()->getSubject(loSubject.getIdNumber(),0);
                    lmCase[1] = moPopulation->getHaplotypes()->getSubject(loSubject.getIdNumber(),1);
                    //  Store case.
                    loHaplotypes.setSubject(2*liCount,0,lmCase[0]);
                    loHaplotypes.setSubject(2*liCount,1,lmCase[1]);
                    //  Create new subject based on untransmitted haplotype and store as pseudo control.
                    lmParents[0] = moPopulation->getHaplotypes()->getSubject(loSubject.getFather()->getIdNumber(),0);
                    lmParents[1] = moPopulation->getHaplotypes()->getSubject(loSubject.getFather()->getIdNumber(),1);
                    lmParents[2] = moPopulation->getHaplotypes()->getSubject(loSubject.getMother()->getIdNumber(),0);
                    lmParents[3] = moPopulation->getHaplotypes()->getSubject(loSubject.getMother()->getIdNumber(),1);
                    moPopulation->getSNPs()->untransmitted(lmParents, lmCase, lmControl);
                    loHaplotypes.setSubject(2*liCount+1,0,lmControl[0]);
                    loHaplotypes.setSubject(2*liCount+1,1,lmControl[1]);
                    break;
                }
                default: msError << "Not supported population type"; error(); break;
            }
            lmPhenotypes(2*liCount)   = 2;
            lmPhenotypes(2*liCount+1) = 1;
            lmPedigree(2*liCount)     = loSubject.getPedigreeNumber();
            lmPedigree(2*liCount+1)   = loSubject.getPedigreeNumber();
            liCount++;
        }
    }

    //  Set up the population object and create dummy names.
    loSNPs.copySNPs(moPopulation->getSNPs());
    loPopulation->clear();
    loPopulation->setType(cePopulationHaplotypes);
    loPopulation->resizeSubjects(2*liCount);
    loPopulation->resizeSNPs(miSNPs);
    loPopulation->getSNPs()->copySNPs(&loSNPs);
    loPopulation->resizeData();
    loPopulation->setTraits(CUtility::makeReal(lmPhenotypes));
    loPopulation->setPhenotypes(lmPhenotypes);
    for (i = 0; i < 2*liCount; i++) {
        loPopulation->getSubjectP(i)->setPedigree(CUtility::intToString(lmPedigree(i),5));
        loPopulation->getSubjectP(i)->setPedigreeNumber(lmPedigree(i));
    }
    loPopulation->setHaplotypes(loHaplotypes);
    loPopulation->clean(false);

    postRun();
}
