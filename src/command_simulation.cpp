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
#include "command_load.h"
#include "haplotype.h"
#include "random.h"

void CCommandSimulation::initialise() {
    msCommand = "simulation";
}

void CCommandSimulation::run(const vector<string>& paParameters) {
    string  lsParameter = paParameters[0], lsPopulation;
    TUnits  loUnit;

    preRun(paParameters);

    if ("snps" == lsParameter) {
        lsPopulation = CUtility::getParameterString(paParameters,"population","sample");
        moPopulation = soProcess->getPopulation(lsPopulation);
        snps();
    } else if ("affection" == lsParameter) {
        lsPopulation = CUtility::getParameterString(paParameters,"population","sample");
        moPopulation = soProcess->getPopulation(lsPopulation);
        affection();
    } else if ("haplotypes" == lsParameter) {
        lsPopulation = CUtility::getParameterString(paParameters,"population","sample");
        moPopulation = soProcess->getPopulation(lsPopulation);
        haplotypes();
    } else if ("subjects" == lsParameter) {
        loUnit = CCommandSimulation::getUnit(maParameters);
        lsPopulation = CUtility::getParameterString(paParameters,"sample","sample");
        moPopulationSample = soProcess->getPopulation(lsPopulation);
        subjects(loUnit);
    } else if ("replace" == lsParameter) {
        lsPopulation = CUtility::getParameterString(paParameters,"sample","sample");
        moPopulationSample = soProcess->getPopulation(lsPopulation);
        replace();
    } else {
        help();
    }

    postRun();
}

void CCommandSimulation::getSubjects(const TIVector& pmPhenotypes, TIVector& pmSubjects, TIVector& pmPedigrees, const int piStatus) {
    int     i, liCount, liPedigree;
    pmSubjects.resize(sum(piStatus == pmPhenotypes));
    pmPedigrees.resize(pmSubjects.rows());
    liCount = 0;
    for (i = 0; i < pmPhenotypes.rows(); i++) {
        if (piStatus == pmPhenotypes(i)) {
            liPedigree = moPopulation->getSubjectP(i)->getPedigreeNumber();
            pmSubjects(liCount)  = i;
            pmPedigrees(liCount) = liPedigree > -1 ? liPedigree : i;
            liCount++;
        }
    }
}

void CCommandSimulation::snps() {
    int             liSNPs, liAffectingCount;
    TReal           lrMAFBound, lrOddsRatio;
    string          lsMAFDistribution, lsMAFEmapFile;
    CCommandLoad    loLoad;
    lsMAFEmapFile = CUtility::getParameterString(maParameters,"emap","");
    if ("" != lsMAFEmapFile) {
        //  Create the SNP profile from an EMAP file.
        loLoad.setPopulation(moPopulation);
        loLoad.run("emap",lsMAFEmapFile);
    } else {
        //  Extract the simulation parameters.
        liSNPs            = CUtility::getParameterInteger(maParameters,"count", 0);
        lsMAFDistribution = CUtility::getParameterString(maParameters,"mafdistribution", "uniform");
        lrMAFBound        = CUtility::getParameterReal(maParameters,"mafbound", 0.5);
        liAffectingCount  = CUtility::getParameterInteger(maParameters,"affectingcount", 0);
        lrOddsRatio       = CUtility::getParameterReal(maParameters,"oddsratio", 1);
        //  Create SNP information.
        moPopulation->resizeSNPs(liSNPs);
        moPopulation->getSNPs()->createMAFs(lsMAFDistribution, lrMAFBound);
        moPopulation->getSNPs()->setAffecting(liAffectingCount, lrOddsRatio);
    }
}

void CCommandSimulation::affection() {
    TReal   lrBaseLine, lrCaseMin, lrCaseMax, lrControlMin, lrControlMax;
    string  lsTrait;
    //  Extract the affected parameters.
    lsTrait      = CUtility::getParameterString(maParameters,"trait","binary");
    lrBaseLine   = CUtility::getParameterReal(maParameters,"baseline",    0.1);
    lrControlMin = CUtility::getParameterReal(maParameters,"controlMin", -crInfinity);
    lrControlMax = CUtility::getParameterReal(maParameters,"controlMax",  0.0);
    lrCaseMin    = CUtility::getParameterReal(maParameters,"caseMin",     0.0);
    lrCaseMax    = CUtility::getParameterReal(maParameters,"caseMax",     crInfinity);
    moPopulation->setTraitType(lsTrait);
    moPopulation->getSNPs()->setBaseLine(lrBaseLine);
    moPopulation->getSNPs()->setBounds(lrControlMin,lrControlMax,lrCaseMin,lrCaseMax);
}

void CCommandSimulation::haplotypes () {
    unsigned int    i, liSubjects;
    THVector        lmHaplotype;
    //  Set up the population.
    liSubjects = CUtility::stringToInt(maParameters[1]) / 2;
    moPopulation->setType(cePopulationHaplotypes);
    moPopulation->resizeSubjects(liSubjects);
    moPopulation->resizeData();
    //  Create random sample of haplotypes.
    lmHaplotype.resize(moPopulation->countSNPs());
    for (i = 0; i < moPopulation->countSubjects(); i++) {
        lmHaplotype = moPopulation->getSNPs()->createHaplotype();
        moPopulation->getHaplotypes()->setSubject(i, 0, lmHaplotype);
        lmHaplotype = moPopulation->getSNPs()->createHaplotype();
        moPopulation->getHaplotypes()->setSubject(i, 1, lmHaplotype);
    }
}

void CCommandSimulation::subjects (const TUnits poUnit) {
    unsigned int        i, j, liSNPs, liHaplotype, liOffset, liCount, liCaseParents, liCaseSiblings, liCaseUnrelateds;
    TBVector            lmCases;
    THVector           *lmHaplotype;
    TIVector            lmSubjects, lmPedigrees, lmIndices;
    TRVector            lmTraits;
    CSubject           *loSubject, *loParent0, *loParent1;
    blitz::firstIndex   loIndex1;

    liSNPs = moPopulationSample->countSNPs();
    lmHaplotype = new THVector[2*poUnit.miUnitSizeFalse];
    for (i = 0; i < 2*poUnit.miUnitSizeFalse; i++) {
        lmHaplotype[i].resize(liSNPs);
    }
    lmTraits.resize(poUnit.miUnitSizeFalse);
    lmCases.resize(poUnit.miUnitSizeFalse);
    lmIndices.resize(poUnit.miUnitSizeTrue);
    lmIndices = loIndex1 + 2 - poUnit.miParents;
    TUniformD loDistributionN(0,2*moPopulationSample->countSubjects()-1);
    TUniformG loRandomN(CRandom::soBoostRandom,loDistributionN);

    //  moPopulation must have enough room to create the cases/controls, SNP counts must match and moPopulationSample must contain haplotype data.
    if (cePopulationNone == moPopulation->getType()) {
        moPopulation->setType(cePopulationHaplotypes);
        moPopulation->resizeSubjects(poUnit.miUnitSizeTrue * poUnit.miCount);
        moPopulation->resizeSNPs(liSNPs);
        moPopulation->resizeData();
    }
    getSubjects(moPopulation->getPhenotypes(), lmSubjects, lmPedigrees);
    if (lmSubjects.rows() < (int) (poUnit.miUnitSizeTrue*poUnit.miCount)) {
        msError << "Not enough space to create subjects: " << moPopulation->countSubjects() << " " << lmSubjects.rows() << " " << poUnit.miUnitSizeTrue << " " << poUnit.miCount << endl;
        error();
    }
    if (0 == liSNPs || moPopulationSample->countSNPs() != liSNPs) {
        msError << "Wrong SNP count for simulation: " << moPopulationSample->countSNPs() << " " << liSNPs << endl;
        error();
    }
    if (cePopulationHaplotypes != moPopulationSample->getType()) {
        msError << "Wrong population type to sample from: " << moPopulationSample->getType() << endl;
        error();
    }
    if (cePopulationHaplotypes != moPopulation->getType()) {
        msError << "Wrong population type to create to: " << moPopulation->getType() << endl;
        error();
    }

    //  Create
    moPopulation->getSNPs()->copySNPs(moPopulationSample->getSNPs());
    for (i = 0; i < poUnit.miCount;) {
        //  Randomly select the parental haplotypes from the pool.
        liCount = 0;
        lmTraits = 0.0;
        lmCases = false;
        liCaseParents    = 0;
        liCaseSiblings   = 0;
        liCaseUnrelateds = 0;
        for (j = 0; j < 2; j++) {
            liHaplotype = loRandomN();
            lmHaplotype[2*liCount]   = moPopulationSample->getHaplotypes()->getSubject(liHaplotype/2,liHaplotype%2);
            liHaplotype = loRandomN();
            lmHaplotype[2*liCount+1] = moPopulationSample->getHaplotypes()->getSubject(liHaplotype/2,liHaplotype%2);
            lmTraits(liCount) = moPopulation->getSNPs()->trait(lmHaplotype[2*liCount], lmHaplotype[2*liCount+1]);
            lmCases(liCount) = moPopulation->getSNPs()->affected(lmTraits(liCount));
            liCaseParents += lmCases(liCount);
            liCount++;
        }
        //  Randomly create offspring from the parents.
        for (j = 0; j < poUnit.miSiblings; j++) {
            moPopulation->getSNPs()->createOffspring(lmHaplotype, &(lmHaplotype[2*liCount]));
            lmTraits(liCount) = moPopulation->getSNPs()->trait(lmHaplotype[2*liCount], lmHaplotype[2*liCount+1]);
            lmCases(liCount) = moPopulation->getSNPs()->affected(lmTraits(liCount));
            liCaseSiblings += lmCases(liCount);
            liCount++;
        }
        //  Randomly create unrelateds from the pool.
        for (j = 0; j < poUnit.miUnrelateds; j++) {
            liHaplotype = loRandomN();
            lmHaplotype[2*liCount]   = moPopulationSample->getHaplotypes()->getSubject(liHaplotype/2,liHaplotype%2);
            liHaplotype = loRandomN();
            lmHaplotype[2*liCount+1] = moPopulationSample->getHaplotypes()->getSubject(liHaplotype/2,liHaplotype%2);
            lmTraits(liCount) = moPopulation->getSNPs()->trait(lmHaplotype[2*liCount], lmHaplotype[2*liCount+1]);
            lmCases(liCount) = moPopulation->getSNPs()->affected(lmTraits(liCount));
            liCaseUnrelateds += lmCases(liCount);
            liCount++;
        }
        if (parentsOk(poUnit, liCaseParents) && siblingsOk(poUnit, liCaseSiblings) && unrelatedsOk(poUnit, liCaseUnrelateds)) {
            liOffset = poUnit.miUnitSizeTrue * i;
            loParent0 = NULL;
            loParent1 = NULL;
            switch (poUnit.miParents) {
            case 1:
                loParent0 = moPopulation->getSubjectP(lmSubjects(liOffset+0));
                loParent1 = NULL;
                lmIndices(0) = 0;
                if (1 == poUnit.miCaseParents    && ! lmCases(0)) { lmIndices(0) = 1; }
                if (1 == poUnit.miControlParents &&   lmCases(0)) { lmIndices(0) = 1; }
                break;
            case 2:
                loParent0 = moPopulation->getSubjectP(lmSubjects(liOffset+0));
                loParent1 = moPopulation->getSubjectP(lmSubjects(liOffset+1));
                break;
            }
            //  Store these new families in the population object.
            for (j = 0; j < poUnit.miUnitSizeTrue; j++) {
                loSubject = moPopulation->getSubjectP(lmSubjects(liOffset+j));
                loSubject->setId(CUtility::intToString(liOffset+j,5));
                loSubject->setIdNumber(lmSubjects(liOffset+j));
                loSubject->setPedigree(CUtility::intToString(lmPedigrees(liOffset),5));
                loSubject->setPedigreeNumber(lmPedigrees(liOffset));
                loSubject->setGender(j % 2 + 1);
                loSubject->setParents(loParent0,loParent1);
                if (loParent0) { loSubject->setParent1(loParent0->getId()); }
                if (loParent1) { loSubject->setParent2(loParent1->getId()); }
                moPopulation->setTrait(lmSubjects(liOffset+j), lmTraits(lmIndices(j)));
                moPopulation->setPhenotype(lmSubjects(liOffset+j), lmCases(lmIndices(j)) ? 2 : 1);
                moPopulation->getHaplotypes()->setSubject(lmSubjects(liOffset+j), lmHaplotype[2*lmIndices(j)],lmHaplotype[2*lmIndices(j)+1]);
            }
            //  Correct for overzealous loop above.
            if (loParent0) { loParent0->setParents(NULL,NULL); loParent0->setParent1(""); loParent0->setParent2(""); }
            if (loParent1) { loParent1->setParents(NULL,NULL); loParent1->setParent1(""); loParent1->setParent2(""); }
            i++;
        }
    }
    moPopulation->clean(true);
    delete[] lmHaplotype;
}

void CCommandSimulation::replace () {
    int         i;
    TIVector    lmSubjects, lmPedigrees;
    TUnits      loUnit;

    getSubjects(moPopulation->getPhenotypes(), lmSubjects, lmPedigrees, 1);
    for (i = 0; i < lmSubjects.rows(); i++) {
        moPopulation->setTrait(lmSubjects(i), 0.0);
        moPopulation->setPhenotype(lmSubjects(i), 0);
    }

    loUnit.miCount             = lmSubjects.rows();
    loUnit.miCaseSiblings      = 0;
    loUnit.miCaseParents       = 0;
    loUnit.miCaseUnrelateds    = 0;
    loUnit.miControlSiblings   = 0;
    loUnit.miControlParents    = 0;
    loUnit.miControlUnrelateds = 1;
    loUnit.miUnknownSiblings   = 0;
    loUnit.miUnknownParents    = 0;
    loUnit.miUnknownUnrelateds = 0;
    loUnit.miSiblings          = 0;
    loUnit.miParents           = 0;
    loUnit.miUnrelateds        = 1;
    loUnit.miUnitSizeFalse     = 3;
    loUnit.miUnitSizeTrue      = 1;
    subjects(loUnit);
}

TUnits CCommandSimulation::getUnit(const vector<string>& paParameters) {
    TUnits  roUnit;
    roUnit.miCount             = CUtility::getParameterInteger(paParameters, "units",             0);
    roUnit.miCaseSiblings      = CUtility::getParameterInteger(paParameters, "caseSiblings",      0);
    roUnit.miCaseParents       = CUtility::getParameterInteger(paParameters, "caseParents",       0);
    roUnit.miCaseUnrelateds    = CUtility::getParameterInteger(paParameters, "caseUnrelateds",    0);
    roUnit.miControlSiblings   = CUtility::getParameterInteger(paParameters, "controlSiblings",   0);
    roUnit.miControlParents    = CUtility::getParameterInteger(paParameters, "controlParents",    0);
    roUnit.miControlUnrelateds = CUtility::getParameterInteger(paParameters, "controlUnrelateds", 0);
    roUnit.miUnknownSiblings   = CUtility::getParameterInteger(paParameters, "unknownSiblings",   0);
    roUnit.miUnknownParents    = CUtility::getParameterInteger(paParameters, "unknownParents",    0);
    roUnit.miUnknownUnrelateds = CUtility::getParameterInteger(paParameters, "unknownUnrelateds", 0);

    roUnit.miSiblings      = roUnit.miCaseSiblings   + roUnit.miControlSiblings   + roUnit.miUnknownSiblings;
    roUnit.miParents       = roUnit.miCaseParents    + roUnit.miControlParents    + roUnit.miUnknownParents;
    roUnit.miUnrelateds    = roUnit.miCaseUnrelateds + roUnit.miControlUnrelateds + roUnit.miUnknownUnrelateds;
    roUnit.miUnitSizeFalse = 2 + roUnit.miSiblings + roUnit.miUnrelateds;
    roUnit.miUnitSizeTrue  = roUnit.miUnitSizeFalse + roUnit.miParents - 2;

    return roUnit;
}

bool CCommandSimulation::parentsOk(const TUnits poUnit, const unsigned int piCases) {
    switch (piCases) {
        case 0: return 0 == poUnit.miCaseParents;
        case 1: return 2 != poUnit.miCaseParents && 2 != poUnit.miControlParents;
        case 2: return 0 == poUnit.miControlParents;
    }
    return false;
}

bool CCommandSimulation::siblingsOk(const TUnits poUnit, const unsigned int piCases) {
    return piCases >= poUnit.miCaseSiblings && (poUnit.miSiblings - piCases) >= poUnit.miControlSiblings;
}

bool CCommandSimulation::unrelatedsOk(const TUnits poUnit, const unsigned int piCases) {
    return piCases >= poUnit.miCaseUnrelateds && (poUnit.miUnrelateds - piCases) >= poUnit.miControlUnrelateds;
}

