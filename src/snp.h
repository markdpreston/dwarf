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
#ifndef SNP_H
#define SNP_H

#include "types.h"

class CSNPs;

class CSnip {
    public:
                    CSnip               ();
                    CSnip               (string psSNPLine);
                    CSnip               (string psChromosome, string psName, string msDistance, int piPosition, string psAllele1, string psAllele2);
    void            clear               ();
    void            initialise          (string psChromosome, string psName, string msDistance = "0", int piPosition = 0, string psAllele1 = "1", string psAllele2 = "2", string psMAF = "0.0", string psEffect = "0", string psOddsRatio = "1.0");
    void            input               (string psSNPLine);
    static int      transform           (const char pcCode, const bool pbExists = true);
    static string   transform           (const int piCode);
    //  Getters
    string          getChromosome       () const { return msChromosome; }
    string          getName             () const { return msName; }
    string          getDistance         () const { return msDistance; }
    unsigned int    getPosition         () const { return miPosition; }
    int             getAllele           (const unsigned int piAllele) const;
    string          getEffect           () const { return msEffect; }
    TReal           getOddsRatio        () const { return mrOddsRatio; }
    TReal           getLogOddsRatio     () const { return mrLogOddsRatio; }
    TReal           getMAF              () const { return mrMAF; }
    string          getNotes            () const { return msNotes; }
    TRVector        getHWE              () const ;
    TRVector        getRelativeRisk     (const TReal prBaseLine) const;
    //  Setters
    void            setNotes            (const string psNotes) { msNotes = psNotes; }
    //  Operators
    bool            operator==          (const CSnip& poSubject) const;
    CSnip           operator=           (const CSnip& poSubject);
    friend istream& operator>>          (istream& in, CSnip& poSubject);
    friend ostream& operator<<          (ostream& out, const CSnip& poSubject);
    friend class    CSNPs;

    protected:
    string          msChromosome;
    string          msName;
    string          msDistance;
    unsigned int    miPosition;
    int             miAllele1;
    int             miAllele2;
    TReal           mrMAF;
    string          msEffect;
    TReal           mrOddsRatio;
    TReal           mrLogOddsRatio;
    string          msNotes;
};

typedef blitz::Array<CSnip,1> TSnipVector;

#endif
