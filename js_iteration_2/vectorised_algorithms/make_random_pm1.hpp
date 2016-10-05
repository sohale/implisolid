#pragma once

#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace mp5_implicit {

namespace vectorised_algorithms {
    int global_rng_seed = 12;
};

namespace vectorised_algorithms {

vectorized_vect make_random_pm1(vindex_t n, int dims, REAL amplitude) {
    //vectorized_vect result {n, dims};
    vectorized_vect result {boost::extents[n][dims]};
    //mt11213b r = boost::mt11213b();
    int seed = vectorised_algorithms::global_rng_seed;
    boost::random::mt11213b rngen(seed);
    //boost::random::uniform_
    boost::random::uniform_01<REAL> distr;
    for (vindex_t i = 0; i < n; ++i) {
        result[i][0] = (distr(rngen) * 2.0 - 1.0 ) * amplitude;
        result[i][1] = (distr(rngen) * 2.0 - 1.0 ) * amplitude;
        result[i][2] = (distr(rngen) * 2.0 - 1.0 ) * amplitude;
    }
    return result;
}

}
}
