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
#ifndef SNPDATA_H
#define SNPDATA_H

#include "types.h"
#include "snp.h"

typedef enum {
    ceSNPEncodeBinary,
    ceSNPEncodeDwarf,
    ceSNPEncodePlink,
    ceSNPEncodeMendel,
    ceSNPEncodeBeagle,
    ceSNPEncodeTransmit
} CSNPEncode;

class CSNPData {
    public:
    CSNPEncode      meEncode;
    ESNPCode        meCode;
    CSnip           moSNP;
    friend istream& operator>>          (istream& in, CSNPData& poSNP) {
                                            int liCode1, liCode2;
                                            string lsCode;
                                            //  Get the characters.
                                            in >> lsCode;
                                            //  Check the length.
                                            switch (lsCode.size()) {
                                                case 1:
                                                    liCode1 = CSnip::transform(lsCode[0]);
                                                    in >> lsCode;
                                                    liCode2 = CSnip::transform(lsCode[0]);
                                                    break;
                                                case 2:
                                                    liCode1 = CSnip::transform(lsCode[0]);
                                                    liCode2 = CSnip::transform(lsCode[1]);
                                                    break;
                                                case 3:
                                                    liCode1 = CSnip::transform(lsCode[0]);
                                                    liCode2 = CSnip::transform(lsCode[2]);
                                                    break;
                                                default:
                                                    cerr << "Bad PED format (length): " << in.tellg() << " " << lsCode.size() << " " << lsCode << endl;
                                                    exit(0);
                                            }
                                            //  Test the codes.
                                            if (0 == liCode1 || 0 == liCode2) {
                                                //  Bad data, return 2.
                                                poSNP.meCode = ceSNPError;
                                            } else if (liCode1 != liCode2) {
                                                //  Heterozygote.
                                                poSNP.meCode = ceSNPHetero;
                                            } else if (liCode1 == poSNP.moSNP.getAllele(1) && liCode1 == liCode2) {
                                                //  First homozygote.
                                                poSNP.meCode = ceSNPMajor;
                                            } else if (liCode1 == poSNP.moSNP.getAllele(2) && liCode1 == liCode2) {
                                                //  Second homozygote.
                                                poSNP.meCode = ceSNPMinor;
                                            } else {
                                                cerr << "Bad SNP codes (impossible to reach error): " << lsCode << " " << liCode1 << " " << liCode2 << " " << poSNP.moSNP.getAllele(1) << " " << poSNP.moSNP.getAllele(2) << endl;
                                                exit(0);
                                            }
                                            return in;
                                        }
    friend ostream& operator<<          (ostream& out, const CSNPData& poSNP) {
                                            string  lsSeparator;
                                            switch (poSNP.meEncode) {
                                                case ceSNPEncodeBinary: return out;
                                                case ceSNPEncodeDwarf:    lsSeparator = "";  break;
                                                case ceSNPEncodePlink:    lsSeparator = " "; break;
                                                case ceSNPEncodeMendel:   lsSeparator = "|"; break;
                                                case ceSNPEncodeBeagle:   lsSeparator = " "; break;
                                                case ceSNPEncodeTransmit: lsSeparator = "/"; break;
                                            }
                                            out << " ";
                                            switch (poSNP.meCode) {
                                                case ceSNPMinor:  out << CSnip::transform(poSNP.moSNP.getAllele(2)) << lsSeparator << CSnip::transform(poSNP.moSNP.getAllele(2)); break;
                                                case ceSNPHetero: out << CSnip::transform(poSNP.moSNP.getAllele(1)) << lsSeparator << CSnip::transform(poSNP.moSNP.getAllele(2)); break;
                                                case ceSNPError:  out << "0" << lsSeparator << "0"; break;
                                                case ceSNPMajor:  out << CSnip::transform(poSNP.moSNP.getAllele(1)) << lsSeparator << CSnip::transform(poSNP.moSNP.getAllele(1)); break;
                                            }
                                            return out;
                                        }
};

#endif
