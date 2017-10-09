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
#include "command.h"
#include "population.h"

void CCommandTransform::run (const vector<string>& paParameters) {
    unsigned int    i, j;
    EPopulationType lePopulation;
    TDosage         loDosage;

    preRun(paParameters);

    lePopulation = cePopulationNone;
    if ("haplotypes" == paParameters[0]) { lePopulation = cePopulationHaplotypes; }
    if ("genotypes"  == paParameters[0]) { lePopulation = cePopulationGenotypes;  }
    if ("dosages"    == paParameters[0]) { lePopulation = cePopulationDosages;    }
    if ("polydata"   == paParameters[0]) { lePopulation = cePopulationPolyData;   }

    if (cePopulationNone == lePopulation || cePopulationDosages == lePopulation) {
        msError << "Bad population type: " << paParameters[0];
        error();
    }

    if (mePopulation == lePopulation) {
        msLog << "No transformation necessary: " << paParameters[0];
        postRun();
        return;
    }

    switch (mePopulation) {
        case cePopulationHaplotypes: {
            if (cePopulationGenotypes == lePopulation) {
                moPopulation->setType(cePopulationGenotypes);
                moPopulation->getGenotypes()->resize(miSubjects,miSNPs);
                for (i = 0; i < miSubjects; i++) {
                    for (j = 0; j < miSNPs; j++) {
                        switch (moPopulation->getHaplotypes()->get(i,j)) {
                            case ceHaplotype00: moPopulation->getGenotypes()->set(i,j,ceSNPMajor);  break;
                            case ceHaplotype01: moPopulation->getGenotypes()->set(i,j,ceSNPHetero); break;
                            case ceHaplotype10: moPopulation->getGenotypes()->set(i,j,ceSNPHetero); break;
                            case ceHaplotype11: moPopulation->getGenotypes()->set(i,j,ceSNPMinor);  break;
                        }
                    }
                }
                moPopulation->getHaplotypes()->clear();
            }
            if (cePopulationDosages == lePopulation) {
                moPopulation->setType(cePopulationDosages);
                moPopulation->getDosages()->resize(miSubjects,miSNPs);
                for (i = 0; i < miSubjects; i++) {
                    for (j = 0; j < miSNPs; j++) {
                        loDosage.mrMajor  = 0;
                        loDosage.mrHetero = 0;
                        loDosage.mrMinor  = 0;
                        switch (moPopulation->getHaplotypes()->get(i,j)) {
                            case ceHaplotype00: loDosage.mrMajor  = 1.0; moPopulation->getDosages()->set(i,j,loDosage);  break;
                            case ceHaplotype01:
                            case ceHaplotype10: loDosage.mrHetero = 1.0; moPopulation->getDosages()->set(i,j,loDosage); break;
                            case ceHaplotype11: loDosage.mrMinor  = 1.0; moPopulation->getDosages()->set(i,j,loDosage);  break;
                        }
                    }
                }
                moPopulation->getHaplotypes()->clear();
            }
            break;
        }
        case cePopulationGenotypes: {
            if (cePopulationHaplotypes == lePopulation) {
                moPopulation->setType(cePopulationHaplotypes);
                moPopulation->getHaplotypes()->resize(miSubjects,miSNPs);
                for (i = 0; i < miSubjects; i++) {
                    for (j = 0; j < miSNPs; j++) {
                        switch (moPopulation->getGenotypes()->get(i,j)) {
                            case ceSNPMajor:  moPopulation->getHaplotypes()->set(i,j,0,0); break;
                            case ceSNPHetero: moPopulation->getHaplotypes()->set(i,j,0,1); break;
                            case ceSNPError:  moPopulation->getHaplotypes()->set(i,j,1,0); break;
                            case ceSNPMinor:  moPopulation->getHaplotypes()->set(i,j,1,1); break;
                        }
                    }
                }
                moPopulation->getGenotypes()->clear();
            }
            if (cePopulationDosages == lePopulation) {
                moPopulation->setType(cePopulationDosages);
                moPopulation->getDosages()->resize(miSubjects,miSNPs);
                for (i = 0; i < miSubjects; i++) {
                    for (j = 0; j < miSNPs; j++) {
                        loDosage.mrMajor  = 0;
                        loDosage.mrHetero = 0;
                        loDosage.mrMinor  = 0;
                        switch (moPopulation->getGenotypes()->get(i,j)) {
                            case ceSNPMajor:  loDosage.mrMajor  = 1.0; moPopulation->getDosages()->set(i,j,loDosage); break;
                            case ceSNPHetero: loDosage.mrHetero = 1.0; moPopulation->getDosages()->set(i,j,loDosage); break;
                            case ceSNPMinor:  loDosage.mrMinor  = 1.0; moPopulation->getDosages()->set(i,j,loDosage); break;
                            case ceSNPError:  msError << "Bad transformation"; error();
                        }
                    }
                }
                moPopulation->getGenotypes()->clear();
            }
            break;
        }
        default: {}
    }

    if (lePopulation != moPopulation->getType()) {
        msError << "Transformation not taken place";
        error();
    }

    postRun();
}
