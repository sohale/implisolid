

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

inline bool sizes_are_equal(const verts_t & A, const faces_t & B) {
    return A.shape()[0] == B.shape()[0]  && A.shape()[1] == B.shape()[1];
}

// A[C] = B[C]
inline void set_a_b_if_c(vectorized_vect & A, const vectorized_vect & B, const vectorized_bool &C){
    assert(sizes_are_equal(A, B));
    assert(B.shape()[0] == C.shape()[0]);

    for (int i=0, n = B.shape()[0]; i < n; i++) {
        if (C[i]) {
            A[i][0] = B[i][0];
            A[i][1] = B[i][1];
            A[i][2] = B[i][2];
        } else {
        }
    }
}


// assign_using_fancy_indexing
 // A = B[CI];
inline void assign_vects_chosen_by_fancy_indexing(
    vectorized_vect & A,
    const vectorized_vect & B,
    const array_of_indices & CI,
    vectorized_vect::index  other_size_to_assert
) {
    vectorized_vect::index n = CI.size();
    assert( other_size_to_assert == n);
    for (int i=0; i < other_size_to_assert; i++) {
        auto k = CI[i];
        A[k][0] = B[k][0];
        A[k][1] = B[k][1];
        A[k][2] = B[k][2];
    }
}


/*
// template <typename T>
inline void assign_selected(T& a, const T& b, bool[] bool_arr, int n) {
    for (int i=0; i < n; i++) {
        if (bool_arr[i]) {
            a[i][0] = b[i][0];
            a[i][1] = b[i][1];
            a[i][2] = b[i][2];
        }
    }
}

assign_fancy_indexes(centroids, x_bisect, relevants_bool_indices)
inline void assign_fancy_indexes(T& a, const T& b, int[] indices_arr, int m) {
    for (int i=0; i < m; i++) {
        if (bool_arr[i]) {
            a[f[i]][0] = b[i][0];
            a[f[i]][1] = b[i][1];
            a[f[i]][2] = b[i][2];
        }
    }
}

*/



inline void bool_find_zero_scalars(vectorized_bool & zeros2_bool, const vectorized_scalar & f2, REAL ROOT_TOLERANCE) {
    auto n = f2.shape()[0];
    assert( n == zeros2_bool.shape()[0]);
    for (int i=0; i < n; ++i) {
        zeros2_bool[i] = std::abs(f2[i])<= ROOT_TOLERANCE;
        /*

        if (std::abs(f2[i])<= ROOT_TOLERANCE) {
            zeros2_bool[i] = b_true;
        } else {
            zeros2_bool[i] = b_false;
        }
        */
    }
}


inline array_of_indices build_range_array(array_of_indices::value_type n) {
    array_of_indices_shape shape = {n};
    array_of_indices  A(shape);
    int i=0;
    for (auto iter = A.begin(), e = A.end(); iter < e; ++iter) {
        //A[i] = i;
        *iter = i;
        ++i;
    }
    assert(i == n);
    if (n > 0) {
        assert(A[n/2] == n/2);
        assert(A[0] == 0);
        assert(A[n-1] == n-1);
    }
    return A; // std::move(A);
}

/*
template <typename T, T v>
inline void set_array_to_value(A, v) {
    for (auto iter = A.begin(), e = A.end(); iter < e; ++iter) {
        *iter = v;
    }
}
*/

// todo: write a TDD test for this. // not checked
inline void set_all_array_elements_to_a_boolean_value(vectorized_bool & A, vectorized_bool::value_type v) {
    for (auto iter = A.begin(), e = A.end(); iter < e; ++iter) {
        *iter = v;
    }
}

}
}

inline ostream& operator<<(ostream& os, const vectorized_vect::value_type& single_vector)
{
    os << single_vector[0] << ',' << single_vector[1] << ',' << single_vector[2];
    return os;
}
