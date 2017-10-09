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
#ifndef PEDIGREE_H
#define PEDIGREE_H

#include "types.h"
#include "subject.h"

class CPedigree {
     public:
                        CPedigree               ();
     void               clear                   ();
     void               setPedigree             (const int piPedigree) { miPedigree = piPedigree; }
     void               add                     (CSubject *poSubject, const int piPhenotype);
     int                countParents            ();
     int                countAffectedSiblings   ();
     CSubject*          getParent               (const int piParent);
     CSubject*          getAffectedSibling      (const int piSibling);
     CPedigree&         operator=               (const CPedigree& poPedigree);
//     private:
     int                miSize;
     int                miPedigree;
     TSubjectPVector    mmSubjects;
     TIVector           mmPhenotypes;
};

#endif
