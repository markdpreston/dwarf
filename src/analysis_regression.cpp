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
#include "analysis.h"
#include "interface_r.h"

void CAnalysisRegression::run(const vector<string>& paParameters) {
    bool            lbBinary;
    unsigned int    i, liRows = 0;
    string          lsType, lsLink;
    TRVector        lmY;
    TRMatrix        lmRegressionData, lmBetas, lmData;
    TSVector        lmColumns;

    preRun(paParameters);
    lsType = CUtility::getParameterString(paParameters,"method","individual");
    lbBinary = "binary" == moPopulation->getTraitType();

    //  Expose the phenotypes to R.
    lmY.resize(miSubjects);
    if (lbBinary) {
        lmY = moPopulation->getPhenotypes() - 1.0;
    } else {
        lmY = moPopulation->getTraits();
    }
#ifdef INCLUDE_R
    CInterfaceR::setVector("y", lmY);
#endif

    //  Link function.
    lsLink = lbBinary ? "family=binomial(link=\"logit\"), " : "";

    //  Do individual regressions, or combined regression by default.
    if ("individual" == lsType) {
        //  Correct mmData size.
        lmRegressionData.resize(miSubjects,1);
        lmBetas.resize(2,4);
        for (i = 0; i < miSNPs; i++) {
#ifdef INCLUDE_R
            //  Extract and expose the dosage data for each SNP
            switch (moPopulation->getType()) {
                case cePopulationDosages:    lmRegressionData(CUtility::getAll(),0) = CDosage::convertToDosage(moPopulation->getDosages()->getSNP(i)); break;
                case cePopulationGenotypes:  lmRegressionData(CUtility::getAll(),0) = CDosage::convertToDosage(moPopulation->getGenotypes()->getSNP(i)); break;
                case cePopulationHaplotypes: lmRegressionData(CUtility::getAll(),0) = CDosage::convertToDosage(CHaplotype::transform(moPopulation->getHaplotypes()->getSNP(i))); break;
                default: msError << "Population type not implemented"; error();
            }
            CInterfaceR::setMatrix("x1", lmRegressionData);
            //  Run regression and extract data from R.
            CInterfaceR::runCommand("x2 <- as.data.frame(x1); remove(x1); r <- summary(glm(y ~ . , data = x2," + lsLink + " na.action=na.pass)); b <- coef(r);");
            liRows = CInterfaceR::runCommandI("nrow(b);");
            if (liRows == 2) {
                lmBetas = CInterfaceR::runCommandRM("b;");
                mmDataLine((int)i) = lmBetas(1,3);
            } else {
                mmDataLine((int)i) = 1.0;
            }
            CInterfaceR::runCommand("remove(x2); remove(r); remove(b)");
            liRows = miSNPs;
#endif
        }
    } else {
        lmRegressionData.resize(miSubjects,miSNPs);
        switch (moPopulation->getType()) {
            case cePopulationDosages:    lmRegressionData = moPopulation->getDosages()->getDosage(); break;
            case cePopulationGenotypes:  lmRegressionData = moPopulation->getGenotypes()->getDosage(false); break;
            case cePopulationHaplotypes: lmRegressionData = moPopulation->getHaplotypes()->getDosage(); break;
            default: msError << "Population type not implemented"; error();
        }
        //  Expose data to R.
#ifdef INCLUDE_R
        CInterfaceR::setMatrix("x1", lmRegressionData);
        //  Run regression and extract data from R.
        CInterfaceR::runCommand("x2 <- as.data.frame(x1); remove(x1); r <- summary(glm(y ~ . , data = x2, " + lsLink + ", na.action=na.pass)); b <- coef(r);"); // print (b);");
        lmRegressionData.resize(CInterfaceR::runCommandI("nrow(b);"), CInterfaceR::runCommandI("ncol(b);"));
        lmRegressionData = CInterfaceR::runCommandRM("b;");
        for (i = 0; i < miSNPs; i++) {
            mmDataLine((int)i) = lmRegressionData((int)i+1, 3);
        }
        //  Clean up R.
        lmRegressionData.resize(0,0);
        CInterfaceR::runCommand("remove(y); remove(x2); remove(r); remove(b)");
#endif
    }

    postRun();
}
