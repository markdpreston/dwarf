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
#include "command_save.h"
#include "utility.h"

void CCommandSave::run (const vector<string>& paParameters) {
    string      lsType, lsFile;
    EFileType   leType;
    ofstream    fout;

    if (! preRun(paParameters)) { return; }
    if (2 > paParameters.size()) { help(); return; }

    lsType = paParameters[0];
    lsFile = paParameters[1];

    fout.open(lsFile.c_str());
    if (! fout.good()) {
        msError << "Data file " << lsFile << " not opened";
        error();
    }
    leType = CUtility::getFileType(lsType);
    switch (leType) {
        case (EFileTypeNone):     fout.close(); return;
        case (EFileTypeBim):
        case (EFileTypeGens):
        case (EFileTypeLegend):
        case (EFileTypeMap):
        case (EFileTypeEmap):
        case (EFileTypeMarkers):  saveMarkers(fout, leType); break;
        case (EFileTypeFam):
        case (EFileTypeHaps):
        case (EFileTypeKbac):
        case (EFileTypePed):
        case (EFileTypePoly):
        case (EFileTypeSamples):
        case (EFileTypeTransmit): savePed(fout, leType);     break;
        case (EFileTypeBed):      saveBed(fout);             break;
        case (EFileTypeGenes):    saveGenes(fout);           break;
        default: msError << "Unrecognised file type: " << lsType; error();
    }
    fout.close();

    postRun();
}

void CCommandSave::saveMarkers (ofstream& fout, EFileType peType) {
//    cout << "    saveMarkers: " << miSNPs << endl;
    unsigned int    i, j;
    CSnip           loSNP;
    TDosage         loDosage;
    for (i = 0; i < miSNPs; i++) {
        loSNP = moPopulation->getSNPInfo(i);
        switch (peType) {
            case EFileTypeBim:     fout << loSNP.getChromosome() << " " << loSNP.getName() << " " << loSNP.getDistance() << " " << loSNP.getPosition() << " " << CSnip::transform(loSNP.getAllele(1)) << " " << CSnip::transform(loSNP.getAllele(2)) << endl; break;
            case EFileTypeEmap:    fout << loSNP.getChromosome() << " " << loSNP.getName() << " " << loSNP.getDistance() << " " << loSNP.getPosition() << " " << loSNP.getMAF() << " " << loSNP.getEffect() << " " << loSNP.getEffect() << "" << loSNP.getNotes() << endl; break;
            case EFileTypeGens:    fout << loSNP.getChromosome() << " " << loSNP.getName() << " " << loSNP.getPosition() << " " << CSnip::transform(loSNP.getAllele(1)) << " " << CSnip::transform(loSNP.getAllele(2)); break;
            case EFileTypeLegend:  fout << "- " << loSNP.getName() << " " << loSNP.getPosition() << " " << CSnip::transform(loSNP.getAllele(1)) << " " << CSnip::transform(loSNP.getAllele(2)) << endl; break;
            case EFileTypeMap:     fout << loSNP.getChromosome() << " " << loSNP.getName() << " " << loSNP.getDistance() << " " << loSNP.getPosition() << endl; break;
            case EFileTypeMarkers: fout << loSNP.getName() << " " << loSNP.getDistance() << " " << loSNP.getPosition() << endl; break;
            default: msError << "saveMarkers error" << miSNPs; error();
        }
        if (EFileTypeGens == peType) {
            for (j = 0; j < miSubjects; j++) {
                loDosage = moPopulation->getDosages()->get(j,i);
                fout << " " << loDosage.mrMajor << " " << loDosage.mrHetero << " " << loDosage.mrMinor;
            }
            fout << endl;
        }
    }
}

void CCommandSave::savePed (ofstream& fout, EFileType peType) {
    unsigned int    i, j;
    int             liPhenotype;
    CSNPData        loSNPData;
    CSubject        loSubject;
    switch (peType) {
        case EFileTypePed:      loSNPData.meEncode = ceSNPEncodePlink;    break;
        case EFileTypeTransmit: loSNPData.meEncode = ceSNPEncodeTransmit; break;
        default: break;
    }
    for (i = 0; i < miSubjects; i++) {
        loSubject   = moPopulation->getSubjectInfo(i);
        liPhenotype = moPopulation->getPhenotype(i);
        switch (peType) {
            case EFileTypeFam:
            case EFileTypePed:
            case EFileTypePoly:
            case EFileTypeTransmit: fout << setw(5) << loSubject.getPedigree() << setw(8) << loSubject.getId() << setw(8) << loSubject.getParent1() << setw(8) << loSubject.getParent2() << setw(2) << loSubject.getGender() << setw(2) << liPhenotype; break;
            case EFileTypeSamples:  fout << loSubject.getPedigree() << " " << loSubject.getId() << " - " << loSubject.getGender() << " " << liPhenotype; break;
            case EFileTypeHaps:
            case EFileTypeKbac:     fout << (liPhenotype == 2 ? "1 " : "0 "); break;
            default: msError << "savePed error"; error();
        }
        if (EFileTypePed == peType || EFileTypeTransmit == peType) {
            for (j = 0; j < miSNPs; j++) {
                switch (mePopulation) {
                    case (cePopulationGenotypes):  loSNPData.meCode = moPopulation->getGenotypes()->get(i,j); fout << loSNPData; break;
                    case (cePopulationHaplotypes): loSNPData.meCode = CHaplotype::transform(moPopulation->getHaplotypes()->get(i,j)); fout << loSNPData; break;
                    case (cePopulationPolyData):   fout << " " << CPolyData::quadOutput(CPolyData::getQuadData(moPopulation->getPolyData()->get(i,j)),true); break;
                    default: msError << "save ped: No genotype/haplotype to output"; error();
                }
            }
        }
        if (EFileTypePoly == peType) {
            for (j = 0; j < miSNPs; j++) {
                fout << " " << CPolyData::output(moPopulation->getPolyData()->get(i,j),true);
            }
        }
        if (EFileTypeKbac == peType) {
            if (0 == moPopulation->getHaplotypes()->miSize) {
                msError << "save kbac: No haplotype to output";
            }
            for (j = 0; j < miSNPs; j++) {
                 switch (moPopulation->getHaplotypes()->get(i,j)) {
                     case ceHaplotype00: fout << "0 "; break;
                     case ceHaplotype11: fout << "2 "; break;
                     default:            fout << "1 ";
                 }
             }
        }
        if (EFileTypeHaps == peType) {
            if (0 == moPopulation->getHaplotypes()->miSize) {
                msError << "save haps: No haplotype to output";
                error();
            }
            for (j = 0; j < miSNPs; j++) {
                 switch (moPopulation->getHaplotypes()->get(i,j)) {
                     case ceHaplotype00: fout << "00 "; break;
                     case ceHaplotype01: fout << "01 "; break;
                     case ceHaplotype10: fout << "10 "; break;
                     case ceHaplotype11: fout << "11 "; break;
                 }
             }
        }
        fout << endl;
    }
}

void CCommandSave::saveBed (ofstream& fout) {
    unsigned int    i, liBytes, liSNPs, liSubjects;
    char            lcCode[3];
    unsigned char  *laBytes;
    lcCode[0] = 0x6c;
    lcCode[1] = 0x1b;
    lcCode[2] = 0x01;
    fout.write(lcCode, 3);

    liSubjects = moPopulation->countSubjects();
    liSNPs     = moPopulation->countSNPs();
    liBytes    = (liSubjects / 4) + ((liSubjects % 4) == 0 ? 0 : 1);
    laBytes    = (unsigned char*) malloc (liBytes);
//    liBytes = (liSNPs / 4) + ((liSNPs % 4) == 0 ? 0 : 1);
//    laBytes = (unsigned char*) malloc (liBytes);
    for (i = 0; i < liSNPs; i++) {
        switch (moPopulation->getType()) {
            case cePopulationGenotypes:  moPopulation->getGenotypes()->untransform(i, laBytes, liBytes, true, false); break;
            case cePopulationHaplotypes: moPopulation->getHaplotypes()->untransform(i, laBytes, liBytes, true, true); break;
            default: msError << "Population type not implemented"; error();
        }
        fout.write((char*) laBytes, liBytes);
    }
    free(laBytes);
}

void CCommandSave::saveGenes (ofstream& fout) const {
}
