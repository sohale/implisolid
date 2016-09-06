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

verts_t
produce_facet_normals(
    boost::multi_array<int, 2> faces,
    vectorized_vect  verts,
    bool force_normalisation)
{
    const REAL min_norm_sq = mp5_implicit::CONFIG_C::MIN_AREA * mp5_implicit::CONFIG_C::MIN_AREA;
    const REAL osqrt3 = 1.0 / std::sqrt((REAL)3.0);
    auto numfaces = faces.shape()[0];
    verts_t facet_normals = verts_t(boost::extents[numfaces][3]);
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



void print_vector(const std::string heading, const verts_t& v, int count) {
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




void compute_centroid_gradient(const verts_t& X, verts_t& N, implicit_function* gradou) {
    std::clog << "Using this function is discouraged. This cause a bug." << std::endl;

    gradou->eval_gradient(X, &N);
    vectorised_algorithms::normalize_1111(N);
}



void compute_centroids(faces_t& faces, verts_t& verts, verts_t& centroids){
  int nt = faces.shape()[0];
  for (int j=0; j<nt; j++){
    int f0 = faces[j][0];
    int f1 = faces[j][1];
    int f2 = faces[j][2];
    for (int di=0; di<3; di++){
        centroids[j][di] = (verts[f0][di] + verts[f1][di] + verts[f2][di])/3.;

    }
  }
}


}

