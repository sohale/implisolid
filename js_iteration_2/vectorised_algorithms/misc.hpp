

#include <random>

namespace mp5_implicit {
namespace vectorised_algorithms {


// void indices_of_zero_normals(verts_t& normals)

void replace_zero_normals_with_gaussian_random(verts_t& normals) {
    assert(normals.shape()[1] == 3);

    /*
    int seed = vectorised_algorithms::global_rng_seed;
    boost::random::mt11213b rngen(seed);
    boost::random::uniform_01<REAL> distr;
    */
    std::default_random_engine generator;
    std::normal_distribution<REAL> distribution(0.0, 1.0);

    for (int i = 0, e = normals.shape()[0]; i < e; i++) {
        REAL norm2 = norm_squared(normals[i][0], normals[i][1], normals[i][2]);
        if (norm2 < 0.9999) {
            // mi++;

            double x = distribution(generator);
            double y = distribution(generator);
            double z = distribution(generator);

            normals[i][0] = x;
            normals[i][1] = y;
            normals[i][2] = z;

            // Note that it is intentionally not normalised
        }
    }
}


// Will not work: template <REAL x, REAL y, REAL z>
//inline void set_vector(verts_t& vect, REAL x, REAL y, REAL z) {
//    satic_assert(vect.shape()[1] == 3);
//                set_vector<1, 0, 0>(dxc, 1, 0, 0);

inline void fill_vector(verts_t& vect, REAL x, REAL y, REAL z) {
    assert(vect.shape()[1] == 3);
    for (auto i = vect.begin(), e = vect.end(); i < e; i++) {
        (*i)[0] = x;
        (*i)[1] = y;
        (*i)[2] = z;
    }

}
/*
inline void set_vector_from(verts_t& target, const verts_t& source) {
    assert(target.shape()[1] == 3);
    auto i = target.begin();
    auto e = target.end();
    auto si = source.begin();
    for (; i < e; ++i, ++si) {
        (*i)[0] = (*si)[0];
        (*i)[1] = (*si)[1];
        (*i)[2] = (*si)[2];
    }
}
*/


inline bool sizes_are_equal(const verts_t & A, const verts_t & B) {
    return A.shape()[0] == B.shape()[0]  && A.shape()[1] == B.shape()[1];
}

inline bool sizes_are_equal(const  & A, const faces_t & B) {
    return A.shape()[0] == B.shape()[0]  && A.shape()[1] == B.shape()[1];
}


}
}
