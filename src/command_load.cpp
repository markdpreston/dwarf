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
#include "command_load.h"
#include "utility.h"

void CCommandLoad::run (const vector<string>& paParameters) {
    bool            lbHaplotypes;
    unsigned int    liSubjects, liSNPs;

    if (! preRun(paParameters)) { return; }

    if (2 > paParameters.size()) { help(); return; }

    //  Parameters for BED files.
    liSubjects   = CUtility::getParameterInteger(paParameters, "subjects",   0);
    liSNPs       = CUtility::getParameterInteger(paParameters, "snps",       0);
    lbHaplotypes = CUtility::getParameterString(paParameters,  "haplotypes", "false") != "false";
    run (paParameters[0], paParameters[1], liSubjects, liSNPs, lbHaplotypes);

    postRun();
}

void CCommandLoad::run (const string psType, const string psFile, const unsigned int piSubjects, const unsigned int piSNPs, const bool pbHaplotypes) {
    EFileType   leType;
    ifstream    fin;

    fin.open(psFile.c_str());
    if (! fin.good()) {
        msError << "Data file " << psFile << " not opened";
        error();
    }
    leType = CUtility::getFileType(psType);
    switch (leType) {
        case (EFileTypeNone):     fin.close(); return;
        case (EFileTypeBim):
        case (EFileTypeEmap):
        case (EFileTypeGens):
        case (EFileTypeLegend):
        case (EFileTypeMap):
        case (EFileTypeMarkers):  loadMarkers(fin, leType); break;
        case (EFileTypeFam):
        case (EFileTypeHaps):
        case (EFileTypeKbac):
        case (EFileTypePed):
        case (EFileTypePoly):
        case (EFileTypeSamples):
        case (EFileTypeTransmit): loadPed(fin, leType); break;
        case (EFileTypeTab):      loadTabular(fin); break;
        case (EFileTypeGenes):    loadGenes(fin); break;
        case (EFileTypeBed):      loadBed(fin, piSubjects, piSNPs, pbHaplotypes); break;
        case (EFileTypeMS):       loadMS(fin); break;
        case (EFileTypeMito):     loadMito(fin); break;
        case (EFileTypeVCF):      loadVCF(fin); break;
        default: msError << "Unrecognised file type: " << psType; error();
    }
    fin.close();
}

void CCommandLoad::loadMarkers (ifstream& fin, EFileType peType) {
    unsigned int    i, j, liPosition = 0;
    int             liLines, liFields;
    char            lsNotes[100];
    string          lsLineAll, lsMarker, lsChromosome = "0", lsDistance = "0", lsAllele1 = "1", lsAllele2 = "2", lsMAF = "0.0", lsEffect = "0", lsOddsRatio = "1.0", lsCrap;
    stringstream    lsLine;
    CSnip           loSNP;
    TDosage         loDosage;
    //  Get the file parameters.
    getFileSize(fin, liLines, liFields, 1);
    //  Eat any header line.
    if (EFileTypeMap != peType && EFileTypeHaps != peType && EFileTypeBim != peType && EFileTypeEmap != peType) {
        getline(fin, lsLineAll);
        liLines--;
    }
    //  Set SNP count.
    moPopulation->resizeSNPs(liLines);
    //  Set dosage/genotype sizes.
    if (EFileTypeGens == peType) {
        if (0 == moPopulation->countSubjects()) {
            moPopulation->resizeSubjects((liFields - 3) / 2);
        }
        moPopulation->setType(cePopulationDosages);
        moPopulation->resizeData();
    }
    if (EFileTypeHaps == peType) {
        moPopulation->setType(cePopulationHaplotypes);
        moPopulation->resizeData();
    }
    //  Read in each marker.
    for (i = 0; i < moPopulation->countSNPs(); i++) {
        getline(fin, lsLineAll);
        lsLine.clear();
        lsLine.str(lsLineAll);
        switch (peType) {
            case EFileTypeBim:     lsLine >> lsChromosome >> lsMarker >> lsDistance >> liPosition >> lsAllele1 >> lsAllele2; break;
            case EFileTypeGens:    lsLine >> lsMarker >> lsAllele1 >> lsAllele2; break;
            case EFileTypeLegend:  lsLine >> lsCrap >> lsMarker >> liPosition >> lsAllele1 >> lsAllele2; break;
            case EFileTypeMap:
            case EFileTypeEmap:    lsLine >> lsChromosome >> lsMarker >> lsDistance >> liPosition >> lsMAF >> lsEffect >> lsOddsRatio; break;
            case EFileTypeMarkers: lsLine >> lsMarker >> lsDistance >> liPosition; break;
            default: msError << "loadMarkers error"; error(); break;
        }
        loSNP.initialise(lsChromosome, lsMarker, lsDistance, liPosition, lsAllele1, lsAllele2, lsMAF, lsEffect, lsOddsRatio);
        if (EFileTypeGens == peType) {
            for (j = 0; j < moPopulation->countSubjects(); j++) {
                lsLine >> loDosage.mrMajor >> loDosage.mrHetero;
                loDosage.mrMinor = 1 - loDosage.mrMajor - loDosage.mrHetero;
                moPopulation->getDosages()->set(j, i, loDosage);
            }
        }
        if (EFileTypeEmap == peType) {
            lsLine.get(lsNotes,100);
            loSNP.setNotes(lsNotes);
        }
        moPopulation->getSNPs()->setSNP(i,loSNP);
    }
}

void CCommandLoad::loadPed (ifstream& fin, EFileType peType) {
    unsigned int    i, j, liGender, liDosage, liSNPs;
    int             liLines, liFields, liPhenotype;
    string          lsLineAll, lsPedigree, lsId, lsParent1 = "0", lsParent2 = "0", lsMissing, lsPhenotype, lsCrap, lsAlleles, lsAllele1, lsAllele2;
    stringstream    lsLine;
    CSNPData        loSNPData;
    CSubject        loSubject;
    //  Get the file parameters.
    getFileSize(fin, liLines, liFields, 1);
    //  Set subject count.
    moPopulation->resizeSubjects(liLines);
    //  Set SNP count.
    switch (peType) {
        case EFileTypeHaps:     liSNPs = liFields;           break;
        case EFileTypeKbac:     liSNPs = liFields;           break;
        case EFileTypePed:      liSNPs = (liFields - 5) / 2; break;
        case EFileTypePoly:     liSNPs = (liFields - 5) / 2; break;
        case EFileTypeSamples:  liSNPs = (liFields - 4);     break;
        case EFileTypeTransmit: liSNPs = (liFields - 5);     break;
        default: liSNPs = moPopulation->countSNPs();
    }
    if (moPopulation->countSNPs() != liSNPs) {
        moPopulation->resizeSNPs(liSNPs);
    }
    //  Resize the data storage.
    switch (peType) {
        case EFileTypeHaps:     moPopulation->setType(cePopulationHaplotypes); moPopulation->resizeData(); break;
        case EFileTypeKbac:
        case EFileTypePed:
        case EFileTypeTransmit: moPopulation->setType(cePopulationGenotypes);  moPopulation->resizeData(); break;
        case EFileTypePoly:     moPopulation->setType(cePopulationPolyData);   moPopulation->resizeData(); break;
        default: break;
    }
    for (i = 0; i < moPopulation->countSubjects(); i++) {
        getline(fin, lsLineAll);
        lsLine.clear();
        lsLine.str(lsLineAll);
        switch (peType) {
            case EFileTypeHaps:
            case EFileTypeKbac:     lsPedigree = "0"; lsId = "0"; liGender = 0; break;
            case EFileTypeFam:
            case EFileTypePed:
            case EFileTypePoly:
            case EFileTypeTransmit: lsLine >> lsPedigree >> lsId >> lsParent1 >> lsParent2 >> liGender; break;
            case EFileTypeSamples:  lsLine >> lsPedigree >> lsId >> lsCrap >> liGender; break;
            default: msError << "loadPed error"; error(); break;
        }
        loSubject.initialise(lsPedigree, lsId, lsParent1, lsParent2, liGender);
        lsLine >> liPhenotype;
        if (liPhenotype == -9) { liPhenotype = 0; }
        if (EFileTypeKbac == peType) { liPhenotype++; }
        moPopulation->setSubject(i, loSubject);
        moPopulation->setTrait(i, liPhenotype);
        moPopulation->setPhenotype(i, liPhenotype);
        if (EFileTypeHaps == peType) {
            for (j = 0; j < moPopulation->countSNPs(); j++) {
                lsLine >> lsAlleles;
                if ("00" == lsAlleles) { moPopulation->getHaplotypes()->set(i, j, 0, 0); }
                if ("01" == lsAlleles) { moPopulation->getHaplotypes()->set(i, j, 0, 1); }
                if ("10" == lsAlleles) { moPopulation->getHaplotypes()->set(i, j, 1, 0); }
                if ("11" == lsAlleles) { moPopulation->getHaplotypes()->set(i, j, 1, 1); }
            }
        }
        if (EFileTypeKbac == peType) {
            for (j = 0; j < moPopulation->countSNPs(); j++) {
                lsLine >> liDosage;
                switch (liDosage) {
                    case 0: moPopulation->getGenotypes()->set(i, j, ceSNPMajor); break;
                    case 1: moPopulation->getGenotypes()->set(i, j, ceSNPHetero); break;
                    case 2: moPopulation->getGenotypes()->set(i, j, ceSNPMinor); break;
                }
            }
        }
        if (EFileTypePed == peType || EFileTypeTransmit == peType) {
            for (j = 0; j < moPopulation->countSNPs(); j++) {
                lsLine >> loSNPData;
                moPopulation->getGenotypes()->set(i, j, loSNPData.meCode);
            }
        }
        if (EFileTypePoly == peType) {
            for (j = 0; j < moPopulation->countSNPs(); j++) {
                lsLine >> lsAllele1 >> lsAllele2;
                moPopulation->getPolyData()->set(i, j, CPolyData::transform(lsAllele1, lsAllele2));
            }
        }
    }
    //  Set the internal relationships up.
    moPopulation->setRelations();
}

void CCommandLoad::loadTabular (ifstream& fin) {
    unsigned int    i, j, k, liLength, liPosition, liSubjects, liSNPs;
    string          lsLineAll, lsChromosome, lsReference, lsId, lsCrap, lsAlleles, lsAllele1, lsAllele2;
    stringstream    lsLine, lsPosition, lsDistance;
    CSnip           loSNP;
    CSubject        loSubject;
    vector<string>  lsSubjects;

    //  Calculate the number of subjects and SNPs in the tabular file.
    //  Read in all of the subject ids from the first line.
    getline(fin, lsLineAll);
    boost::split(lsSubjects, lsLineAll, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
    liSubjects = lsSubjects.size();
    if (3 > liSubjects) {
        msError << "Not enough fields in tabular file";
        error();
    }
    liSubjects -= 3;
    //  Loop over the remaining lines to count the references.
    liSNPs = 0;
    while (! fin.eof()) {
        getline(fin, lsLineAll);
        if (0 != lsLineAll.length()) {
            lsLine.clear();
            lsLine.str(lsLineAll);
            lsLine >> lsChromosome >> liPosition >> lsReference;
            liSNPs += lsReference.length();
        }
    }

    //  Initialise the data.
    moPopulation->setType(cePopulationPolyData);
    moPopulation->resizeSubjects(liSubjects);
    moPopulation->resizeSNPs(liSNPs);
    moPopulation->resizeData();

    //  Store the subjects.
    for (i = 0; i < liSubjects; i++) {
        loSubject.initialise(string("0"), lsSubjects[i+3], string("0"), string("0"), 0);
        moPopulation->setSubject(i, loSubject);
    }

    //  Move the file pointer back to the start of the file.
    fin.clear();
    fin.seekg(0, ios::beg);
    //  Ignore header/subject line.
    getline(fin, lsLineAll);
    //  Read in all of the SNPs and SNP data from the second line onwards.
    for (i = 0; i < liSNPs; i += liLength) {
        getline(fin, lsLineAll);
        lsLine.clear();
        lsLine.str(lsLineAll);
        lsLine >> lsChromosome >> liPosition >> lsReference;
        liLength = lsReference.length();
        //  Store the SNPs.
        for (j = 0; j < liLength; j++) {
            lsPosition.str(""); lsPosition << "pos"  << (liPosition+j);
            lsDistance.str(""); lsDistance << "dist" << (liPosition+j);
            loSNP.initialise(lsChromosome, lsPosition.str(), lsDistance.str(), liPosition + j);
            moPopulation->getSNPs()->setSNP(i + j,loSNP);
        }
        //  Store the data.
        for (j = 0; j < liSubjects; j++) {
            for (k = 0; k < liLength; k++) {
                lsLine >> lsAlleles;
                lsAllele1 = lsAlleles[k];
                lsAllele2 = lsAlleles[k + liLength + 1];
                moPopulation->getPolyData()->set(j, i + k, CPolyData::transform(lsAllele1, lsAllele2));
            }
        }
    }
}

void CCommandLoad::loadBed (ifstream& fin, const unsigned int piSubjects, const unsigned int piSNPs, const bool pbHaplotypes) {
    unsigned int    i, liBytes, liSNPs, liSubjects;
    char            lcCode[3];
    unsigned char  *laBytes;
    CBinaryData    *loData;
    fin.read(lcCode,3);
    if (0x6c != lcCode[0] || 0x1b != lcCode[1]) {
        msError << "Bad BED code";
        error();
    }
    liSubjects = piSubjects != 0 ? piSubjects : miSubjects;
    liSNPs     = piSNPs     != 0 ? piSNPs     : miSNPs;
    if (moPopulation->countSubjects() != liSubjects) {
       moPopulation->resizeSubjects(liSubjects);
    }
    if (moPopulation->countSNPs() != liSNPs) {
       moPopulation->resizeSNPs(liSNPs);
    }
    moPopulation->setType(pbHaplotypes ? cePopulationHaplotypes : cePopulationGenotypes);
    moPopulation->resizeData();
    loData = pbHaplotypes ? (CBinaryData*) moPopulation->getHaplotypes() : (CBinaryData*) moPopulation->getGenotypes();
    //  Individual order, i.e. all SNPs for one individual.
    if (0 == lcCode[2]) {
        liBytes = (liSNPs / 4) + ((liSNPs % 4) == 0 ? 0 : 1);
        laBytes = (unsigned char*) malloc (liBytes);
        for (i = 0; i < liSubjects; i++) {
            fin.read((char*) laBytes, liBytes);
            loData->transform(i, laBytes, liBytes, false);
        }
    //  SNP order, i.e. all individual for each SNP.
    } else {
        liBytes = (liSubjects / 4) + ((liSubjects % 4) == 0 ? 0 : 1);
        laBytes = (unsigned char*) malloc (liBytes);
        for (i = 0; i < liSNPs; i++) {
            fin.read((char*) laBytes, liBytes);
            loData->transform(i, laBytes, liBytes, true);
        }
    }
    free(laBytes);
}

void CCommandLoad::loadGenes (ifstream& fin) {
    unsigned int    i;
    int             liPositionStart, liPositionEnd, liLines, liFields;
    char            lsNotes[100];
    string          lsLineAll, lsChromosome, lsName;
    stringstream    lsLine;
    CGene           loGene;
    //  Get the file parameters.
    getFileSize(fin, liLines, liFields, 1);
    //  Resize the number of genes.
    moPopulation->getGenes()->resize(liLines);
    for (i = 0; i < moPopulation->countGenes(); i++) {
        getline(fin, lsLineAll);
        lsLine.clear();
        lsLine.str(lsLineAll);
        lsLine >> lsChromosome >> lsName >> liPositionStart >> liPositionEnd;
        lsLine.get(lsNotes, 100);
        loGene.initialise(lsChromosome, lsName, liPositionStart, liPositionEnd, lsNotes);
        moPopulation->getGenes()->setGene(i, loGene);
    }
}

void CCommandLoad::loadVCF (ifstream& fin) {
    int             i, j, liFields, liInfo, liLineCount;
    string          lsLineAll, lsChromosome, lsPosition, lsReference, lsName, lsGenotype, lsInfo, lsNull;
    stringstream    lsLine;
    vector<string>  laNames, laInfo, laGenotype;
    //  Skip the header.
    liLineCount = 0;
    do {
        getline(fin, lsLineAll);
        liLineCount++;
        lsLine.clear();
        lsLine.str(lsLineAll);
        lsLine >> lsNull;
    } while ("#CHROM" != lsNull);
    //  Ignore titles and read in line of subjects.
    for (i = 1; i < 9; i++) {
        lsLine >> lsNull;
    }
    liFields = 0;
    while (! lsLine.eof()) {
        liFields++;
        lsLine >> lsName;
        laNames.push_back(lsName);
    }
    //  Read in the genetic data.
    while (! fin.eof()) {
        getline(fin, lsLineAll);
        lsLine.clear();
        lsLine.str(lsLineAll);
        lsLine >> lsChromosome >> lsPosition >> lsReference >> lsNull >> lsNull >> lsNull >> lsNull >> lsNull >> lsInfo;
        boost::algorithm::split(laInfo, lsInfo, boost::is_any_of(":"));
        liInfo = -1;
        for (j = 0; j < laInfo.size(); j++) {
            if ("GT" == laInfo[j]) {
                liInfo = j;
                break;
            }
        }
        if (-1 == liInfo) {
            cerr << "No DP data in VCF file" << endl;
            exit(1);
        }
        for (j = 0; j < liFields; j++) {
            lsLine >> lsGenotype;
            boost::algorithm::split(laGenotype, lsGenotype, boost::is_any_of(":"));
            lsGenotype = laGenotype[liInfo];
        }
    }
}

void CCommandLoad::loadMito (ifstream& fin) {
    int                 i, j, liLines, liFields;
    string              lsLineAll, lsChromosome, lsPosition, lsReference, lsName, lsGenotype, lsNull;
    stringstream        lsLine;
    vector<string>      laNames, laChromosomes;
    vector<string>::iterator liChromosome;
    map<string,string>  laReferences;
    //  Get the file parameters.
    getFileSize(fin, liLines, liFields, 0);
cout << liLines << " " << liFields << endl;
    //  Read in the subject names.
    getline(fin, lsLineAll);
    lsLine.clear();
    lsLine.str(lsLineAll);
    for (i = 0; i < liFields; i++) {
        lsLine >> lsName;
        if (i > 9) {
            laNames.push_back(lsName);
        }
    }
cout << laNames.size() << endl << laNames[0] << endl;
    //  Resize the number of genes.
    for (i = 1; i < liLines; i++) {
        getline(fin, lsLineAll);
        lsLine.clear();
        lsLine.str(lsLineAll);
        lsLine >> lsChromosome >> lsPosition >> lsReference >> lsNull >> lsNull >> lsNull >> lsNull >> lsNull;
        laChromosomes.push_back(lsChromosome);
        laReferences[lsChromosome].append(lsReference);
        for (j = 10; j < liFields; j++) {
            lsLine >> lsGenotype;
            if ("-" == lsGenotype) { cout << "missing" << endl; }
            if (lsReference == lsGenotype) { cout << "major" << endl; }
        }
    }
    sort(laChromosomes.begin(), laChromosomes.end());
    liChromosome = unique(laChromosomes.begin(), laChromosomes.end());
    laChromosomes.resize(liChromosome - laChromosomes.begin());
cout << laChromosomes.size() << endl << laChromosomes[0] << endl;
cout << laReferences.size() << endl << laReferences[laChromosomes[0]] << endl;
}

void CCommandLoad::loadMS (ifstream& fin) {
/*
    unsigned int    i, j, liSamples, liChromosomes, liSegsites, liCount;
    TReal           lrPosition;
    string          lsLineAll, lsCrap;
    stringstream    lsLine, lsChromosome;
    CSnip           loSNP;
    //  Part 1: parse the file for matrix sizes and initialise the simulated data class.
    //  Command line.
    getline(fin, lsLineAll);
    lsLine.clear();
    lsLine.str(lsLineAll);
    lsLine >> lsCrap >> liSamples >> liChromosomes;
    moSimulated.resizeMS(liSamples, liChromosomes);
    //  Ignore random number line.
    getline(fin, lsLineAll);
    for (i = 0; i < liChromosomes; i++) {
        //  Ignore blank and // lines.
        getline(fin, lsLineAll);
        getline(fin, lsLineAll);
        //  Segsites line.
        getline(fin, lsLineAll);
        lsLine.clear();
        lsLine.str(lsLineAll);
        lsLine >> lsCrap >> liSegsites;
        moSimulated.setChromosome(i,liSegsites);
        //  Ignore the positions line and the hap lines.
        getline(fin, lsLineAll);
        for (j = 0; j < liSamples; j++) { getline(fin, lsLineAll); }
    }
    moSimulated.initialiseMS();
    //  Part 2: Read in the SNP and haplotype data.
    //  Store all of the SNPs.
    mmSNPs.resize(moSimulated.getChromosomes());
    liCount = 0;
    //  Ignore command and random lines.
    fin.seekg(ios_base::beg);
    getline(fin, lsLineAll);
    getline(fin, lsLineAll);
    for (i = 0; i < liChromosomes; i++) {
        //  Ignore blank, // and segsites lines.
        getline(fin, lsLineAll);
        getline(fin, lsLineAll);
        getline(fin, lsLineAll);
        //  Read positions lines.
        getline(fin, lsLineAll);
        lsLine.clear();
        lsLine.str(lsLineAll);
        lsLine >> lsCrap;
        lsChromosome.str("");
        lsChromosome << i;
        for (j = 0; j < moSimulated.getChromosome(i); j++) {
            lsLine >> lrPosition;
            loSNP.initialise(lsChromosome.str(), lexical_cast<string>(j), lexical_cast<string>(lrPosition), j, "0", "1");
            mmSNPs(liCount) = loSNP;
            liCount++;
        }
        //  Read the hap lines.
        for (j = 0; j < liSamples; j++) {
            getline(fin, lsLineAll);
            moSimulated.setHaplotype(i, j / 2, j % 2, lsLineAll);
        }
    }
*/
}

void CCommandLoad::getFileSize (ifstream& fin, int& piLines, int& piFields, const int piOffset) {
    string          lsLine;
    vector<string>  laWords;
    piLines  = 0;
    piFields = 0;
    //  Count the number of lines in the file.
    while (! fin.eof()) {
        getline(fin, lsLine);
        //  At the specific offset line count the number of fields.
        if (piOffset == piLines) {
            boost::split(laWords, lsLine, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
            piFields = laWords.size();
        }
        if (0 != lsLine.length()) {
            piLines++;
        }
    }
    //  Move the file pointer back to the start of the file.
    fin.clear();
    fin.seekg(0, ios::beg);
}


