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
#ifndef POPULATION_H
#define POPULATION_H

#include "types.h"
#include "subject.h"
#include "snps.h"
#include "genes.h"
#include "pedigree.h"
#include "dosage.h"
#include "genotype.h"
#include "haplotype.h"
#include "polydata.h"

typedef enum {
    cePopulationNone,
    cePopulationDosages,
    cePopulationGenotypes,
    cePopulationHaplotypes,
    cePopulationPolyData
} EPopulationType;

class CPopulation {
    public:
    //  Population construction/destruction.
                    CPopulation         ();
                   ~CPopulation         ();
    //  Type.
    EPopulationType getType             () const { return meType; }
    void            setType             (const EPopulationType peType) { meType = peType; }
    //  Data.
    void            clear               ();
    void            clean               (const bool pbRelations);
    void            setRelations        ();
    void            keepRemove          (const bool pbKeep, const bool pbSubjects, const TUVector& pmIndices);
    void            resizeData          ();
    void            resizeSubjects      (const unsigned int piSubjects);
    void            resizeSNPs          (const unsigned int piSNPs);
    //  Counts.
    unsigned int    countSubjects       () const { return mmSubjects.rows(); }
    unsigned int    countSNPs           () const { return moSNPs.countSNPs(); }
    unsigned int    countGenes          () const { return moGenes.countGenes(); }
    unsigned int    countPedigrees      () const { return miPedigrees; }
    unsigned int    countAffected       () const { return sum(mmPhenotypes == 2); }
    unsigned int    countPhenotypes     (const int piPhenotype) const { return sum(mmPhenotypes == piPhenotype); }
    //  Access.
    CSnip           getSNPInfo          (const int piSNP)     const { return moSNPs.getSNP(piSNP); }
    CSubject        getSubjectInfo      (const int piSubject) const { return mmSubjects(piSubject); }
    CSubject*       getSubjectP         (const int piSubject) { return &mmSubjects(piSubject); }
    CSubject*       find                (const string psPedigree, const string psId);
    CGene           getGene             (const int piGene) { return moGenes.getGene(piGene); }
    int             getPhenotype        (const int piSubject) const { return mmPhenotypes(piSubject); }
    TIVector&       getPhenotypes       () { return mmPhenotypes; }
    TUVector        getPhenotypes       (const int piPhenotype);
    string          getTraitType        () { return msTrait; }
    TReal           getTrait            (const int piSubject) { return mmTraits(piSubject); }
    TRVector&       getTraits           () { return mmTraits; }
    TSubjectPVector getAffected         ();
    CPedigree       getPedigree         (const int piPedigree);
    int             getPedigreeMaximum  ();

    TSubjectVector* getSubjects         () { return &mmSubjects; }
    CSNPs*          getSNPs             () { return &moSNPs; }
    CGenes*         getGenes            () { return &moGenes; }
    CDosage*        getDosages          () { return &moDosages; }
    CGenotype*      getGenotypes        () { return &moGenotypes; }
    CHaplotype*     getHaplotypes       () { return &moHaplotypes; }
    CPolyData*      getPolyData         () { return &moPolyData; }

    void            setSubject          (const int piSubject, const CSubject poSubject) { mmSubjects(piSubject)   = poSubject; }
    void            setSubjects         (const TSubjectVector& pmSubjects) { mmSubjects = pmSubjects; }
    void            setPhenotype        (const int piSubject, const int piPhenotype)    { mmPhenotypes(piSubject) = piPhenotype; }
    void            setPhenotypes       (const TIVector& pmPhenotypes) { mmPhenotypes = pmPhenotypes; }
    void            setTraitType        (const string psTrait) { msTrait = psTrait; }
    void            setTrait            (const int piSubject, const TReal piTrait)    { mmTraits(piSubject) = piTrait; }
    void            setTraits           (const TRVector& pmTraits) { mmTraits = pmTraits; }
    void            setDosages          (const CDosage& poDosages)       { moDosages    = poDosages;    }
    void            setGenotypes        (const CGenotype& poGenotypes)   { moGenotypes  = poGenotypes;  }
    void            setGenotypes        (const CHaplotype& poHaplotypes);
    void            setHaplotypes       (const CHaplotype& poHaplotypes) { moHaplotypes = poHaplotypes; }
    void            setPolyData         (const CPolyData& poPolyData)    { moPolyData   = poPolyData;   }

    TSNPVector      getSNP              (const int piSNP) { return CHaplotype::transform(moHaplotypes.getHaplotype(piSNP,0),moHaplotypes.getHaplotype(piSNP,1)); }
    TSNPVector      getSubject          (const int piSubject) { return moGenotypes.getSubject(piSubject); }

    private:
    //  Field information.
    EPopulationType meType;
    unsigned int    miPedigrees;

    TSubjectVector  mmSubjects;
    string          msTrait;
    TIVector        mmPhenotypes;
    TRVector        mmTraits;

    CSNPs           moSNPs;
    CGenes          moGenes;
    CDosage         moDosages;
    CGenotype       moGenotypes;
    CHaplotype      moHaplotypes;
    CPolyData       moPolyData;
};

#endif
