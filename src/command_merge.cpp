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

void CCommandMerge::run(const vector<string>& paParameters) {
    unsigned int    i, liSNPs;
    string          lsPopulation;
    EPopulationType lePopulation;
    CPopulation*    loPopulation;

    preRun(paParameters);

    lsPopulation = paParameters[0];
    if (! soProcess->isPopulation(lsPopulation)) {
        msError << "No population to merge from: " << lsPopulation;
        error();
    }
    loPopulation = soProcess->getPopulation(lsPopulation);
    lePopulation = loPopulation->getType();

    if (mePopulation != lePopulation) {
        msError << "Mismatch population types: " << mePopulation << " " << lePopulation;
        error();
    }

    liSNPs     = loPopulation->countSNPs();
    if (miSubjects != liSNPs) {
        msError << "Unequal subjects: " << miSubjects << " " << loPopulation->countSubjects();
        error();
    }

    loPopulation->resizeSubjects(miSubjects);
    loPopulation->resizeSNPs(miSNPs+liSNPs);
    loPopulation->setType(cePopulationHaplotypes);
    loPopulation->resizeData();

    //  Merge SNPs
    loPopulation->resizeSNPs(miSNPs+liSNPs);
    loPopulation->getSNPs()->copySNPs(moPopulation->getSNPs());
    for (i = 0; i < liSNPs; i++) {
        loPopulation->getSNPs()->setSNP(miSNPs+i,loPopulation->getSNPs()->getSNP(i));
    }

    //  Merge data
    switch (mePopulation) {
        case cePopulationNone: {
//            moPopulation.copy(loPopulation);
            break;
        }
        case cePopulationHaplotypes: {
            CHaplotype      loHaplotypes;
            THVector        lmHaplotype;
            lmHaplotype.resize(miSNPs+liSNPs);
            loHaplotypes.resize(miSubjects,miSNPs+liSNPs);
            for (i = 0; i < miSubjects; i++) {
                lmHaplotype(blitz::Range(0,miSNPs-1))   = moPopulation->getHaplotypes()->getSubject(i,0);
                lmHaplotype(blitz::Range(miSNPs,blitz::toEnd)) = loPopulation->getHaplotypes()->getSubject(i,0);
                loHaplotypes.setHaplotype(i,0,lmHaplotype);
                lmHaplotype(blitz::Range(0,miSNPs-1))   = moPopulation->getHaplotypes()->getSubject(i,1);
                lmHaplotype(blitz::Range(miSNPs,blitz::toEnd)) = loPopulation->getHaplotypes()->getSubject(i,1);
                loHaplotypes.setHaplotype(i,1,lmHaplotype);
            }
            loPopulation->setHaplotypes(loHaplotypes);
            break;
        }
        case cePopulationGenotypes: {
            CGenotype       loGenotypes;
            TSNPVector      lmGenotype;
            lmGenotype.resize(miSNPs+liSNPs);
            loGenotypes.resize(miSubjects,miSNPs+liSNPs);
            for (i = 0; i < miSubjects; i++) {
                lmGenotype(blitz::Range(0,miSNPs-1))   = moPopulation->getGenotypes()->getSubject(i);
                lmGenotype(blitz::Range(miSNPs,blitz::toEnd)) = loPopulation->getGenotypes()->getSubject(i);
                loGenotypes.setSubject(i,lmGenotype);
            }
            loPopulation->setGenotypes(loGenotypes);
            break;
        }
        case cePopulationDosages: {
            CDosage         loDosages;
            TDVector        lmDosage;
            lmDosage.resize(miSNPs+liSNPs);
            loDosages.resize(miSubjects,miSNPs+liSNPs);
            for (i = 0; i < miSubjects; i++) {
                lmDosage(blitz::Range(0,miSNPs-1))   = moPopulation->getDosages()->getSubject(i);
                lmDosage(blitz::Range(miSNPs,blitz::toEnd)) = loPopulation->getDosages()->getSubject(i);
                loDosages.setSubject(i,lmDosage);
            }
            loPopulation->setDosages(loDosages);
            break;
        }
        case cePopulationPolyData: {
            CPolyData       loPolyData;
            TUVector        lmPolyData;
            lmPolyData.resize(miSNPs+liSNPs);
            loPolyData.resize(miSubjects,miSNPs+liSNPs);
            for (i = 0; i < miSubjects; i++) {
                lmPolyData(blitz::Range(0,miSNPs-1))   = moPopulation->getPolyData()->getSubject(i);
                lmPolyData(blitz::Range(miSNPs,blitz::toEnd)) = loPopulation->getPolyData()->getSubject(i);
                loPolyData.setSubject(i,lmPolyData);
            }
            loPopulation->setPolyData(loPolyData);
            break;
        }
    }

    postRun();
}
