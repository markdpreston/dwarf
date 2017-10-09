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
#include "command_population.h"

void CCommandPopulation::run (const vector<string>& paParameters) {
    unsigned int    liSubjects, liSNPs;
    string          lsCommand, lsPopulationFrom, lsPopulationTo, lsJoin;

    preRun(paParameters);

    //  The command to run.
    lsCommand = paParameters[0];

    //  The new population (can be the same as the from population).
    lsPopulationFrom = CUtility::getParameterString(paParameters, "from", "default");
    lsPopulationTo   = CUtility::getParameterString(paParameters, "to", "default");
    lsJoin           = CUtility::getParameterString(paParameters, "by", "subjects");
    moPopulationFrom = soProcess->getPopulation(lsPopulationFrom);
    moPopulationTo   = soProcess->getPopulation(lsPopulationTo);

    //  Copy a population.
    if ("copy" == lsCommand) {
        moPopulation = soProcess->getPopulation("internalEmpty1985");
        join(false);
    }

    //  Join to the subjects or SNPs.
    if ("join" == lsCommand) {
        join(0 == lsJoin.compare("subjects"));
    }

    //  Keep the 'top-left' data.
    if ("keep" == lsCommand) {
        liSubjects   = CUtility::getParameterInteger(paParameters, "subjects", moPopulationFrom->countSubjects());
        liSNPs       = CUtility::getParameterInteger(paParameters, "snps", moPopulationFrom->countSNPs());
        keep(liSubjects, liSNPs);
    }

    postRun();
}

void CCommandPopulation::join (const bool pbSubjects) {
    unsigned int    i, j, k, liSubject;
    THVector        lmHaplotype, lmHaplotypeExtra;
    TSNPVector      lmGenotype, lmGenotypeExtra;
    CPopulation*    loPopulation;

    //  Calculate the new sizes and type.
    if (pbSubjects) {
        miNewSubjects = moPopulationFrom->countSubjects() + moPopulationTo->countSubjects();
        miNewSNPs = moPopulationFrom->countSNPs();
    } else {
        miNewSubjects = moPopulationFrom->countSubjects();
        miNewSNPs = moPopulationFrom->countSNPs() + moPopulationTo->countSNPs();
    }
    if (cePopulationHaplotypes == moPopulationFrom->getType() && cePopulationHaplotypes == moPopulationTo->getType()) {
        meNewPopulation = cePopulationHaplotypes;
    } else {
        meNewPopulation = cePopulationGenotypes;
    }

    //  Create a new haplotype.
    moNewHaplotypes.resize(miNewSubjects, miNewSNPs);
    mmNewSubjects.resize(miNewSubjects);
    mmNewPhenotypes.resize(miNewSubjects);
    mmNewPhenotypes = 0;

    //  Copy the haplotype, phenotype and pedigree data.
    for (i = 0; i < miNewSubjects; i++) {
        if (i < moPopulationTo->countSubjects()) {
            loPopulation = moPopulationTo;
            liSubject = i;
        } else {
            loPopulation = moPopulationFrom;
            liSubject = i - moPopulationTo->countSubjects();
        }
        mmNewSubjects(i) = loPopulation->getSubjectInfo(liSubject);
        mmNewPhenotypes(i) = loPopulation->getPhenotype(liSubject);
        if (i < moPopulationTo->countSubjects()) {
            mmNewSubjects(i).setPedigree(mmNewSubjects(i).getPedigree() + "NealXX");
        }
        switch (meNewPopulation) {
            case cePopulationHaplotypes: {
                for (j = 0; j < 2; j++) {
                    lmHaplotype.resize(loPopulation->countSNPs());
                    lmHaplotype = loPopulation->getHaplotypes()->getSubject(mmNewSubjects(i).getIdNumber(),j);
                    if (! pbSubjects) {
                        lmHaplotypeExtra.resize(moPopulationFrom->countSNPs());
                        lmHaplotypeExtra = moPopulationFrom->getHaplotypes()->getSubject(mmNewSubjects(i).getIdNumber(),j);
                        lmHaplotype.resizeAndPreserve(miNewSNPs);
                        for (k = 0; k < moPopulationFrom->countSNPs(); k++) {
                            lmHaplotype(moPopulationTo->countSNPs()+k) = lmHaplotypeExtra(k);
                        }
                    }
                    moNewHaplotypes.setSubject(i,j,lmHaplotype);
                }
                break;
            }
            case cePopulationGenotypes: {
                lmGenotype.resize(loPopulation->countSNPs());
                lmGenotype = loPopulation->getGenotypes()->getSubject(mmNewSubjects(i).getIdNumber());
                if (! pbSubjects) {
                    lmGenotypeExtra.resize(moPopulationFrom->countSNPs());
                    lmGenotypeExtra = moPopulationFrom->getGenotypes()->getSubject(mmNewSubjects(i).getIdNumber());
                    lmGenotype.resizeAndPreserve(miNewSNPs);
                    for (k = 0; k < moPopulationFrom->countSNPs(); k++) {
                        lmGenotype(moPopulationTo->countSNPs()+k) = lmGenotypeExtra(k);
                    }
                }
                moNewHaplotypes.setSubject(i,lmGenotype);
                break;
            }
            default:
                cerr << "Unable to copy/join - non-supported population type." << endl;
        }
    }

    //  Copy the SNP data.
    moNewSNPs.copySNPs(moPopulationTo->getSNPs());
    if (! pbSubjects) {
        moNewSNPs.addSNPs(moPopulationFrom->getSNPs());
    }

    //  Store this combined population.
    store();
}

void CCommandPopulation::keep (const unsigned int piSubjects, const unsigned int piSNPs) {
    unsigned int    i;
    TIVector        lmPhenotypes, lmPedigree;
    TSNPVector      lmGenotype;
    CBinaryData    *loBinaryData;
    CHaplotype      loHaplotypes;
    TSubjectVector  lmSubjects;
    CSNPs           loSNPs;

    //  Calculate the new sizes.
    meNewPopulation = moPopulationFrom->getType();
    miNewSubjects = piSubjects > moPopulationFrom->countSubjects() ? moPopulationFrom->countSubjects() : piSubjects;
    miNewSNPs = piSNPs > moPopulationFrom->countSNPs() ? moPopulationFrom->countSNPs() : piSNPs;

    //  Create a new haplotype.
    moNewHaplotypes.resize(miNewSubjects, miNewSNPs);
    mmNewSubjects.resize(miNewSubjects);
    mmNewPhenotypes.resize(miNewSubjects);
    mmNewPhenotypes = 0;

    //  Copy the haplotype, phenotype and pedigree data.
    switch (meNewPopulation) {
        case cePopulationGenotypes: {
            loBinaryData = moPopulationFrom->getGenotypes();
            break;
        }
        case cePopulationHaplotypes: {
            loBinaryData = moPopulationFrom->getHaplotypes();
            break;
        }
        default: {
           cerr << "Population not modified - non-supported type: " << meNewPopulation << endl;
           return;
        }
    }
    lmGenotype.resize(moPopulationFrom->countSNPs());
    for (i = 0; i < miNewSubjects; i++) {
        mmNewSubjects(i) = moPopulationFrom->getSubjectInfo(i);
        mmNewPhenotypes(i) = moPopulationFrom->getPhenotype(i);
        lmGenotype = loBinaryData->getSubject(mmNewSubjects(i).getIdNumber());
        lmGenotype.resizeAndPreserve(miNewSNPs);
        moNewHaplotypes.setSubject(i,lmGenotype);
    }

    //  Copy the SNP data.
    moNewSNPs.copySNPs(moPopulationFrom->getSNPs(),miNewSNPs);

    //  Store this combined population.
    store();
}

void CCommandPopulation::store() {
    //  Set up the population object and create dummy names.
    moPopulationTo->clear();
    moPopulationTo->setType(meNewPopulation);
    moPopulationTo->resizeSubjects(miNewSubjects);
    moPopulationTo->resizeSNPs(miNewSNPs);
    moPopulationTo->getSNPs()->copySNPs(&moNewSNPs);
    moPopulationTo->resizeData();
    moPopulationTo->setSubjects(mmNewSubjects);
    moPopulationTo->setPhenotypes(mmNewPhenotypes);
    moPopulationTo->setTraits(CUtility::makeReal(mmNewPhenotypes));
    switch (meNewPopulation) {
        case cePopulationGenotypes:  moPopulationTo->setGenotypes(moNewHaplotypes);  break;
        case cePopulationHaplotypes: moPopulationTo->setHaplotypes(moNewHaplotypes); break;
        default: cerr << "No new population genetic data copied." << endl;
    }
    moPopulationTo->clean(true);
}
