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
#ifndef RANDOM_H
#define RANDOM_H

#include "types.h"

typedef boost::binomial_distribution<int,TReal> TBinomialD;
typedef boost::variate_generator<boost::mt19937&,boost::binomial_distribution<int,TReal> > TBinomialG;

typedef boost::gamma_distribution<TReal> TGammaD;
typedef boost::variate_generator<boost::mt19937&,boost::gamma_distribution<TReal> > TGammaG;

typedef boost::normal_distribution<TReal> TNormalD;
typedef boost::variate_generator<boost::mt19937&,boost::normal_distribution<TReal> > TNormalG;

typedef boost::poisson_distribution<int,TReal> TPoissonD;
typedef boost::variate_generator<boost::mt19937&,boost::poisson_distribution<int,TReal> > TPoissonG;

typedef boost::uniform_int<> TUniformD;
typedef boost::variate_generator<boost::mt19937&,boost::uniform_int<> > TUniformG;

struct shuffleRNG : std::unary_function<unsigned, unsigned> {
//    shuffleRNG(boost::mt19937 &state) : _state(state) {}
//    boost::mt19937 &_state;
//    unsigned operator()(unsigned i) {
//        boost::uniform_int<> rng(0, i - 1);
//        return rng(_state);
//    }
};


class CRandom {
    public:
    static boost::mt19937               soBoostRandom;

    static void     seed                () {
                                            soBoostRandom.seed(static_cast<unsigned int>(time(0)));
                                        }
    static void     seed                (const int piSeed) {
                                            soBoostRandom.seed(piSeed);
                                        }
    static bool     binary              () {
                                            static boost::uniform_smallint<> soBoostSmallInt(0,1);
                                            static boost::variate_generator<boost::mt19937&,boost::uniform_smallint<> > soBoostBinary(soBoostRandom,soBoostSmallInt);
                                            return soBoostBinary() == 1;
                                        }
    static TReal    uniform             () {
                                            static boost::uniform_01<boost::mt19937&>  soBoostUniform(soBoostRandom);
                                            return soBoostUniform();
                                        }
    static bool     uniform             (TReal prBound) {
                                            return uniform() < prBound;
                                        }
    static TReal    beta                (TGammaG poGammaX, TGammaG poGammaY) {
                                            TReal lrX = poGammaX(), lrY = poGammaY();
                                            return lrX / (lrX + lrY);
                                        }
    static int      uniform             (int i) {
                                            static boost::uniform_int<> soBoostUniform(0, i-1);
                                            static boost::variate_generator<boost::mt19937&, boost::uniform_int<> > soBoostUniformGen(soBoostRandom, soBoostUniform);
                                            return soBoostUniformGen();
                                        }
    static ptrdiff_t shuffleRandom      (int i) { return uniform(i); }
    static void     shuffle             (TIVector& pmVector) {
                                            static ptrdiff_t (*shuffleRandomP)(int) = shuffleRandom;
                                            random_shuffle(pmVector.data(), pmVector.data() + pmVector.size(), shuffleRandomP);
                                        }
};

#endif
