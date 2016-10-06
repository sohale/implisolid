#pragma once

#include "vectorised_algorithms/normalise_inplace.hpp"
//only for compute_centroid_gradient  that uses normalize111
#include "vectorised_algorithms/assert_are_normalised.hpp"

using mp5_implicit::vectorised_algorithms::assert_are_normalised;

// #include "configs.hpp"

// imlpicit_vectorized_loops.hpp

namespace mp5_implicit {



vectorized_scalar get_signs(const vectorized_scalar& scalars, REAL ROOT_TOLERANCE) {
    /*
    auto scalar_shape = scalars.shape();  // boost::extents[n]
    vectorized_scalar signs_c(scalar_shape);
    */
    auto scalar_shape = boost::extents[scalars.shape()[0]];  // n]
    vectorized_scalar signs_c(scalar_shape);

    int e = scalars.shape()[0];
    for (int i=0; i < e; i++) {
        if (scalars[i] > ROOT_TOLERANCE) {
            signs_c[i] = +1.0;
        } else if (scalars[i] < -ROOT_TOLERANCE) {
            signs_c[i] = -1.0;
        } else {
            signs_c[i] = 0.0;
        }
    }
    // todo: make sure it uses move constructor
    return signs_c;
}


/*
    Normalisation policty:
        random vector with normal distribution, amplitude.
        constant: sqrt(1/3)

        leave 0:
        norm = (norm < mp5_implicit::CONFIG_C::center_projection::min_gradient_len)? norm: 1.0;

    normalise1_(v, min_len)
*/

vectorized_vect
produce_facet_normals(
    vectorized_faces faces,
    vectorized_vect  verts,
    bool force_normalisation)
{
    const REAL min_norm_sq = mp5_implicit::CONFIG_C::MIN_AREA * mp5_implicit::CONFIG_C::MIN_AREA;   // why squared?
    const REAL osqrt3 = 1.0 / std::sqrt((REAL)3.0);
    auto numfaces = faces.shape()[0];
    vectorized_vect facet_normals = vectorized_vect(boost::extents[numfaces][3]);
    for (int fi = 0; fi < numfaces; ++fi) {
        auto f0 = faces[fi][0];
        auto f1 = faces[fi][1];
        auto f2 = faces[fi][2];

        auto x1 = verts[f1][0] - verts[f0][0];
        auto y1 = verts[f1][1] - verts[f0][1];
        auto z1 = verts[f1][2] - verts[f0][2];

        auto x2 = verts[f2][0] - verts[f0][0];
        auto y2 = verts[f2][1] - verts[f0][1];
        auto z2 = verts[f2][2] - verts[f0][2];

        REAL x = y1 * z2 - z1 * y2;
        REAL y = z1 * x2 - x1 * z2;
        REAL z = x1 * y2 - y1 * x2;

        if (force_normalisation) {
            REAL n2 = x*x + y*y + z*z;
            if ( n2 < min_norm_sq) {
                x = osqrt3;
                y = osqrt3;
                z = osqrt3;

                cout << "norm3 = " << x*x + y*y + z*z << std::endl;
            }
            else {
                REAL n = std::sqrt(n2);
                x = x / n;
                y = y / n;
                z = z / n;
            }
        }

        auto& nx = facet_normals[fi][0];
        auto& ny = facet_normals[fi][1];
        auto& nz = facet_normals[fi][2];

        /*
        cout <<
        facet_normals[fi][0] << " " <<
        facet_normals[fi][1] << " " <<
        facet_normals[fi][2] << " " <<
        " -> ";
        */

        nx = x;
        ny = y;
        nz = z;

        /*
        cout <<
        facet_normals[fi][0] << " " <<
        facet_normals[fi][1] << " " <<
        facet_normals[fi][2] << " " <<
        std::endl;
        */

    }
    if (force_normalisation) {
        /*
            bads = np.logical_or(np.isnan(facet_areas), np.abs(facet_areas - 0.) < MIN_AREA)
            if not np.allclose(np.linalg.norm(facet_normals[np.logical_not(bads), :], axis=1), 1.):
                set_trace()
            assert np.allclose(np.linalg.norm(facet_normals[np.logical_not(bads), :], axis=1), 1.)  #is not zero
            facet_normals[bads, :] = 1./np.sqrt(3.)  # not tested
        */

        assert(assert_are_normalised(facet_normals));
    }

    return facet_normals;
}



void print_vector(const std::string heading, const vectorized_vect& v, int count) {
    clog << heading << ": ";
    int e = v.shape()[0];
    if (count < e) {
        e = count;
    }
    for (int i=0; i < e; i++) {
        clog << v[i][0] << "," << v[i][1] << "," << v[i][2] << "  ";
    }
    clog << std::endl;
}




void compute_centroid_gradient(const vectorized_vect& X, vectorized_vect& N, implicit_function* gradou) {
    std::clog << "Using this function is discouraged. Can cause bugs." << std::endl;

    gradou->eval_gradient(X, &N);
    vectorised_algorithms::normalize_1111(N);
}


// The only compute_centroids() function
void compute_centroids(const vectorized_faces& faces, const vectorized_vect& verts, vectorized_vect& centroids) {
    int nt = faces.shape()[0];
    for (int j = 0; j < nt; j++) {
        const int f0 = faces[j][0];
        const int f1 = faces[j][1];
        const int f2 = faces[j][2];
        for (int di = 0; di < 3; di++) {
            centroids[j][di] = (verts[f0][di] + verts[f1][di] + verts[f2][di]) / (REAL)(3.0);
        }
    }
}


vectorized_scalar  eval_implicit_on_selected_points_indexed(
        mp5_implicit::implicit_function* object,
        const vectorized_vect & X,
        const array_of_indices & active_indices,
        int count) {
    assert(count == active_indices.shape()[0]);
    vectorized_scalar_shape scalar_shape { count    };
    vectorized_scalar f_a(scalar_shape);
    // assign_vects_chosen_by_fancy_indexing(xa, X, active_indices, ???);
    vectorized_vect_shape shape = {count, 3};
    vectorized_vect   xa(shape);
    for (int j=0; j < count; j++) {
        auto k = active_indices[j];
        xa[j][0] = X[k][0];
        xa[j][1] = X[k][1];
        xa[j][2] = X[k][2];
    }
    //object->eval_implicit(xa.begin(), xa.begin()+count, f_a.begin());

    ///////////// TODO: WE TO DEFINE THIS?
    // vectorized_scalar f_a(scalar_shape);
    object->eval_implicit(xa, &f_a);
    // clog << "Evaluating(a) " << active_indices.shape()[0] << " points" << std::endl;

    return f_a;
}

vectorized_scalar  eval_implicit_on_selected_points_bool(mp5_implicit::implicit_function* object, const vectorized_vect & X, const vectorized_bool & selected_bool,  int count) {
    assert(count == selected_bool.shape()[0]);
    vectorized_scalar_shape scalar_shape { count    };
    vectorized_scalar f_a(scalar_shape);
    // assign_vects_chosen_by_fancy_indexing(xa, X, selected_bool, ???);
    vectorized_vect_shape shape = {count, 3};
    vectorized_vect   xa(shape);
    for (int j=0; j < count; j++) {
        auto k = selected_bool[j];
        xa[j][0] = X[k][0];
        xa[j][1] = X[k][1];
        xa[j][2] = X[k][2];
    }
    //object->eval_implicit(xa.begin(), xa.begin()+count, f_a.begin());

    ///////////// TODO: WE TO DEFINE THIS?
    // vectorized_scalar f_a(scalar_shape);
    object->eval_implicit(xa, &f_a);
    clog << "Evaluating(b) " << selected_bool.shape()[0] << " points" << std::endl;

    return f_a;
}

bool test_if_conjugate_opposite_signs_indexed(
    mp5_implicit::implicit_function* object,
    const vectorized_vect & A, const vectorized_vect & B,
    const vectorized_bool & already_success,
    REAL ROOT_TOLERANCE
    ) {
    // object->eval_implicit(xa4, &f_a);
    vectorized_scalar  f1 = eval_implicit_on_selected_points_bool(object, A, already_success,  already_success.shape()[0]);
    vectorized_scalar  f2 = eval_implicit_on_selected_points_bool(object, B, already_success,  already_success.shape()[0]);

    vectorized_scalar s1 = get_signs(f1, ROOT_TOLERANCE);
    vectorized_scalar s2 = get_signs(f2, ROOT_TOLERANCE);

    assert( s1.shape()[0] == s2.shape()[0]);
    assert( s1.shape()[0] == already_success.shape()[0]);
    int n = s1.shape()[0];
    for(int i = 0; i < n; ++i) {
        if ( s1[i] * s2[i] > 0.0) {
            int i0 = already_success[i];
            clog << " Test failed at point k=" << i << " i=" << i0 <<
                " point:" << A[i0][0] << "," << A[i0][1] << "," << A[i0][2] <<
                " versus:" << B[i0][0] << "," << B[i0][1] << "," << B[i0][2] <<
                std::endl;
            return false;
        }
    }

    return true;
}


bool check_all_are_root(mp5_implicit::implicit_function* object, const vectorized_vect & x_bisect, int m, REAL ROOT_TOLERANCE) {
    // bool evaluate_and_assert_sign<0>(x_bisect, *object);
    /*
    boost::array<int, 2> x1_relevant_shape = {m, 3};
    boost::array<int, 1> f1_relevant_shape = {m};
    */
    const vectorized_scalar::size_type nn = x_bisect.shape()[0];
    vectorized_scalar f(boost::array<vectorized_scalar::size_type, 1>{nn});
    object->eval_implicit(x_bisect, &f);

    bool ok = true;
    for (int i = 0; i < m; i++) {
        // ok = ok && std::abs(f1_relevants[i]) < ROOT_TOLERANCE;
        ok = ok && std::abs(f[i]) < ROOT_TOLERANCE;
        if (!ok) {
            clog << f[i] << " [" << i << "]" << std::endl;
            break;
        }
    }
    return ok;
}

}