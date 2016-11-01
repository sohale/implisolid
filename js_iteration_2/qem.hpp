#pragma once

#include <cstddef>   // for std::nullptr only

#include "../js_iteration_2/basic_functions.hpp"

#include "../js_iteration_2/matrix_functions.hpp"

namespace mp5_implicit {

bool super_quiet =true;

// get the matrix A and b used in vertex_apply_qem
void get_A_b__old(const std::vector<faceindex_type> & nai, const vectorized_vect& centroids, const vectorized_vect& centroid_gradients, vectorized_vect* A, vectorized_scalar* b) {

    int m = nai.size();

    boost::array<int, 2> center_array_shape = {m, 3};
    vectorized_vect  center_array(center_array_shape);
    vectorized_vect  normals(center_array_shape);

    for (int i=0; i < m; i++) {
        vindex_t cn = nai[i];
        normals[i][0] = centroid_gradients[cn][0];
        normals[i][1] = centroid_gradients[cn][1];
        normals[i][2] = centroid_gradients[cn][2];

        center_array[i][0] = centroids[cn][0];
        center_array[i][1] = centroids[cn][1];
        center_array[i][2] = centroids[cn][2];
    }

    REAL nn00;
    REAL nn01;
    REAL nn02;
    REAL nn11;
    REAL nn12;
    REAL nn22;

    (*A)[0][0] = 0;
    (*A)[0][1] = 0;
    (*A)[0][2] = 0;
    (*A)[1][1] = 0;
    (*A)[1][2] = 0;
    (*A)[2][2] = 0;
    (*A)[1][0] = 0;
    (*A)[2][0] = 0;
    (*A)[2][1] = 0;

    (*b)[0] = 0;
    (*b)[1] = 0;
    (*b)[2] = 0;
    for (int j=0; j < m; j++) {
        /* the matrix is symmetric so we dont need to compute all 9 elements*/
        /* A = dot(normals, normals.T) */
        nn00 = normals[j][0] * normals[j][0];
        nn01 = normals[j][0] * normals[j][1];
        nn02 = normals[j][0] * normals[j][2];
        nn11 = normals[j][1] * normals[j][1];
        nn12 = normals[j][1] * normals[j][2];
        nn22 = normals[j][2] * normals[j][2];

        (*A)[0][0] += nn00;
        (*A)[0][1] += nn01;
        (*A)[0][2] += nn02;
        (*A)[1][0] += nn01;
        (*A)[1][1] += nn11;
        (*A)[1][2] += nn12;
        (*A)[2][2] += nn22;
        (*A)[2][0] += nn02;
        (*A)[2][1] += nn12;

        REAL Px = center_array[j][0];
        REAL Py = center_array[j][1];
        REAL Pz = center_array[j][2];

        (*b)[0] -= nn00 * Px + nn01 * Py + nn02 * Pz;
        (*b)[1] -= nn01 * Px + nn11 * Py + nn12 * Pz;
        (*b)[2] -= nn02 * Px + nn12 * Py + nn22 * Pz;

    }

}


/*  Function: vertex_apply_qem(verts, faces, centroids, vertex_neighbours_facelist, centroid_gradients, treated)

    Description:
        This function follows the projection finding and applies the QEM process, to an array of verts,
        in order to update them and get new ones, more accurate.

    Implementation Details:
        So, the main parameter we are interested in
        changing is verts, others are declared as const.

    Parameters:
      verts:
          mesh vertices of the current shape, mutable and updated inside this function

      faces:
          mesh faces (or triangles) of the current shape, immutable for this function

      centroids:
          the centroids of the mesh, immutable

      vertex_neighbours_facelist:
          Holds neighbouring faces/triangles for all the vertices.


      centroids_gradients:
          The gradient of the implicit function of the current shape for
          all the centroids.

      treated:
          A list of flags that indicates the process is complete for every index
 */

void vertex_apply_qem__old(
    vectorized_vect* verts, const vectorized_faces faces,
    const vectorized_vect centroids,
    const std::vector< std::vector<faceindex_type>> vertex_neighbours_facelist,
    const vectorized_vect centroid_gradients
    //const vectorized_bool& treated
    )
{
    assert (verts != nullptr);
    /* The next asserts have no significance in C++, only Python */

    /* assert centroids is not None */
    /* assert vertex_neighbours_facelist is not None */
    /* assert centroids_gradients is not None */

    int nverts = verts->shape()[0];
    /* assert nvert = len(vertex_neighbours_facelist) */

    boost::array<int, 2> A_shape = { 3 , 3 };

    boost::array<int, 2> b_shape = { 3 , 1 };
    vectorized_scalar b(b_shape);
    vectorized_scalar y(b_shape);
    vectorized_scalar utb(b_shape);
    vectorized_scalar new_x(b_shape);
    vectorized_vect u(A_shape);
    vectorized_vect s(A_shape);
    vectorized_vect v(A_shape);
    vectorized_vect A(A_shape);
    for (int vi=0; vi < nverts; vi++) {
        std::vector<faceindex_type> nlist;
        for (int i=0; i < vertex_neighbours_facelist[vi].size(); i++) {
            nlist.push_back(vertex_neighbours_facelist[vi][i]);
        }

        /*
        if (true) {
            std::clog << " remove 'treated' " << std::endl;
            bool skip = false;
            for (int g = 0; g < vertex_neighbours_facelist[vi].size()-1; g++) {
                if (!treated[vertex_neighbours_facelist[vi][g]]) {
                    skip = true;
                }
                for (int g1 = g+1 ; g1 < vertex_neighbours_facelist[vi].size()-1; g1++) {
                    if ((abs(abs(centroid_gradients[vertex_neighbours_facelist[vi][g]][0]) - abs(centroid_gradients[vertex_neighbours_facelist[vi][g1]][0])) < 0.000001)
                    &&(abs(abs(centroid_gradients[vertex_neighbours_facelist[vi][g]][1]) - abs(centroid_gradients[vertex_neighbours_facelist[vi][g1]][1])) < 0.000001)
                    &&(abs(abs(centroid_gradients[vertex_neighbours_facelist[vi][g]][2]) - abs(centroid_gradients[vertex_neighbours_facelist[vi][g1]][2])) < 0.000001))
                    skip = true;
                }
            }



            if (skip) {
                clog << vi << endl;
                continue;
            }
        }
        */

        get_A_b__old(nlist, centroids, centroid_gradients, &A, &b);
        SVD__old(A, u, s, v); // the SVD
        // assert(test_svd(A, u, s, v));

        // in python the SVD values of s are sorted by the svd function, this is a possible workaround
        // (we may need to keep the A=u*s*v equality, which is done this way)
        // assert(np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))) validation assert,
        // also note that u and v are supposed to be unitary
        // so we can add: assert(u.T = u^-1) and assert(v.T == v^-1)

        REAL max_s = max(s[0][0], max(s[1][1], s[2][2]));

        REAL tau = 680;
        int rank = 0;
        if (s[0][0]/max_s < 1./tau) {
            s[0][0] = 0.;
        } else {
            ++rank;
        }
        if (s[1][1]/max_s < 1./tau) {
            s[1][1] = 0.;
        } else {
            ++rank;
        }

        if (s[2][2]/max_s < 1./tau) {
            s[2][2] = 0.;
        } else {
            ++rank;
        }

        // assert s[0] == np.max(s)  asserts that SVD produces descending order eigenvalues
        const REAL x0 = (*verts)[vi][0];
        const REAL y0 = (*verts)[vi][1];
        const REAL z0 = (*verts)[vi][2];
        y[0] = v[0][0] * x0 + v[1][0] * y0 + v[2][0] * z0;
        y[1] = v[0][1] * x0 + v[1][1] * y0 + v[2][1] * z0;
        y[2] = v[0][2] * x0 + v[1][2] * y0 + v[2][2] * z0;

        utb[0] = - u[0][0] * b[0] - u[1][0] * b[1] - u[2][0] * b[2];
        utb[1] = - u[0][1] * b[0] - u[1][1] * b[1] - u[2][1] * b[2];
        utb[2] = - u[0][2] * b[0] - u[1][2] * b[1] - u[2][2] * b[2];

        for (int i=0; i < rank; i++) {
          if (s[i][i] != 0) {
              y[i] = utb[i] / s[i][i];
          } else {
              rank++;
          }
        }

        new_x[0] = v[0][0] * y[0] + v[0][1] * y[1] + v[0][2] * y[2];
        new_x[1] = v[1][0] * y[0] + v[1][1] * y[1] + v[1][2] * y[2];
        new_x[2] = v[2][0] * y[0] + v[2][1] * y[1] + v[2][2] * y[2];


        (*verts)[vi][0] = new_x[0];
        (*verts)[vi][1] = new_x[1];
        (*verts)[vi][2] = new_x[2];
    }
}

/*
*/
inline bool check_normality_for_qem(REAL nx, REAL ny, REAL nz) {
    // constexpr REAL N2 = std::pow(CONFIG_C::MIN_NORMAL_LEN, 2);
    bool ok = true;
    if (norm_squared(nx, ny, nz) < std::pow(CONFIG_C::MIN_NORMAL_LEN, 2)) {
        std::clog << "Error: bad normal: " << nx <<","<< ny <<","<< nz << std::endl;
    }
    ok = ok && !is_bad_number(nx);
    ok = ok && !is_bad_number(ny);
    ok = ok && !is_bad_number(nz);
    return ok;
}

// Calculates A nd b for a selection of centroids only (an umbrella).
// Todo: vectorized version.
void get_A_b(
    const std::vector<faceindex_type> & neighbours_faces,
    const vectorized_vect& centroids,
    const vectorized_vect& centroid_normals,
    const Matrix<REAL, 3, 1> & qem_origin,
    Matrix<REAL, 3, 3> *A, Matrix<REAL, 3, 1> * b) {

    // centroid_normals must be normalised already


    int m = neighbours_faces.size();
    // centroids[neighbours_faces]

    *A << 0,0,0, 0,0,0, 0,0,0;
    *b << 0,0,0;

    for (int i=0; i < m; i++) {

        vindex_t ni = neighbours_faces[i];

        REAL Ni_x = centroid_normals[ni][0];
        REAL Ni_y = centroid_normals[ni][1];
        REAL Ni_z = centroid_normals[ni][2];


        assert(check_normality_for_qem(Ni_x, Ni_y, Ni_z));


        // Pi - origin
        REAL Px = centroids[ni][0] - qem_origin(0);
        REAL Py = centroids[ni][1] - qem_origin(1);
        REAL Pz = centroids[ni][2] - qem_origin(2);

        // Matrix<REAL, 3, 1> P;
        // P << centroids[ni][0], centroids[ni][1], centroids[ni][2];

        /* the matrix is symmetric so we dont need to compute all 9 elements*/
        REAL nn00 = Ni_x * Ni_x;
        REAL nn01 = Ni_x * Ni_y;
        REAL nn02 = Ni_x * Ni_z;
        REAL nn11 = Ni_y * Ni_y;
        REAL nn12 = Ni_y * Ni_z;
        REAL nn22 = Ni_z * Ni_z;

        // A += dot(normals, normals.T)
        (*A)(0,0) += nn00;
        (*A)(0,1) += nn01;
        (*A)(0,2) += nn02;
        (*A)(1,0) += nn01;  // (*,R) is an Abelian group.
        (*A)(1,1) += nn11;
        (*A)(1,2) += nn12;
        (*A)(2,2) += nn22;
        (*A)(2,0) += nn02;
        (*A)(2,1) += nn12;

        // nnt * p_i
        (*b)[0] -= nn00 * Px + nn01 * Py + nn02 * Pz;
        (*b)[1] -= nn01 * Px + nn11 * Py + nn12 * Pz;
        (*b)[2] -= nn02 * Px + nn12 * Py + nn22 * Pz;
    }
}




void vertex_apply_qem(
    vectorized_vect* verts
    , const vectorized_faces faces
    , const vectorized_vect centroids
    , const std::vector< std::vector<faceindex_type>> vertex_neighbours_facelist
    , const vectorized_vect centroid_gradients
    , array_of_indices *ranks_output
    , REAL maximum_displacement_distance
    )
    //const vectorized_bool& treated)
{
    REAL tau = 680.0;  // 1200.0; // 680;
    REAL svd_threshold = 1.0 / tau;

    // note: vertex_neighbours_facelist contains face indices

    // assert (verts != nullptr);

    int nverts = verts->shape()[0];
    assert(nverts = vertex_neighbours_facelist.size());

    bool store_ranks = (ranks_output != nullptr);
    if (store_ranks) {
        assert(ranks_output->shape()[0] == centroids.shape()[0]);
    }
    assert(centroid_gradients.shape()[0] == centroids.shape()[0]);

    // print_vertex_neighbourhood(vertex_neighbours_facelist);
    /*
    // Seems incorrect.
        0: 0, 1, 2, 3.
        1: 0, 1, 4, 5.
        2: 0, 2, 4, 6.
        3: 1, 3, 5, 7.
        4: 2, 3, 6, 7.
        5: 4, 5, 6, 7.
    */

    /*
    for( int j = 0; j <vertex_neighbours_facelist.size(); j++) {
        // clog << vertex_neighbours_facelist[j] << " "
        for( auto v : vertex_neighbours_facelist[j]) {
            clog << v << " ";
        }
        clog << std::endl;
    }
    */

    // Matrix<REAL, 3, 1> y;
    // Matrix<REAL, 3, 1> utb;
    // Matrix<REAL, 3, 1> new_x;

    Matrix<REAL, 3, 3> A;
    Matrix<REAL, 3, 1> b;

    Matrix<REAL, 3, 3> U;
    // Matrix<REAL, 3, 3> S;  // const SingularValuesType' (aka 'const Eigen::Matrix<float, 3, 1, 0, 3, 1>')
    Matrix<REAL, 3, 3> V;

    for (int vi=0; vi < nverts; vi++) {

        // array of faces
        const std::vector<faceindex_type> & nlist = vertex_neighbours_facelist[vi];

        Matrix<REAL, 3, 1> qem_origin;
        qem_origin << (*verts)[vi][0], (*verts)[vi][1], (*verts)[vi][2];  // old verts

        get_A_b(nlist, centroids, centroid_gradients, qem_origin, &A, &b);

        /*
        int rank = SVD(A, U, S, V, svd_threshold); // the SVD
        */

        // Whether calculate the rank based on Eigen3's internal mechanism or manually check the S in SVD.
        #define RANK_MANUALLY  false

        Eigen::JacobiSVD< Matrix<REAL, 3, 3> > svd(A,  Eigen::ComputeFullU | Eigen::ComputeFullV);
        #if (!RANK_MANUALLY)
            svd.setThreshold( svd_threshold );
        #endif
        auto S = svd.singularValues();
        U = svd.matrixU();
        V = svd.matrixV().transpose();
        #if (!RANK_MANUALLY)
            int rank = svd.rank();
        #endif


        assert(S(0) >= S(1));
        assert(S(1) >= S(2));


        #if RANK_MANUALLY
        int rank = 0;
        // if (S(0)/S(0) < svd_threshold) {
        //     S(0) = 0.;
        // } else {
        //     ++rank;
        // }
        ++rank;

        if (S(1)/S(0) < svd_threshold) {
            S(1) = 0.;
        } else {
            ++rank;
        }

        if (S(2)/S(0) < svd_threshold) {
            S(2) = 0.;
        } else {
            ++rank;
        }
        // std::clog << "rank" << rank << std::endl;
        #endif



        // std::clog << "!1" << std::endl;

        // assert s[0] == np.max(s)  asserts that SVD produces descending order eigenvalues
        /*
        MIN_EIGENVALUE = 0.000001
        if not s[0] > MIN_EIGENVALUE:
            print("Warning! sigma_1 == 0" )
            print(s)
            print("A", A)

            #not tested
            result_verts_ranks[vi] = 0
            new_verts[vi, 0:3] = new_x[:, 0]

            rank = 0

        assert np.all(s[:rank]/s[0] >= svd_threshold)
        */


        Matrix<REAL, 3, 1> v0;
        v0 << (*verts)[vi][0] - qem_origin(0), (*verts)[vi][1] - qem_origin(1), (*verts)[vi][2] - qem_origin(2);
        /*
        const REAL x0 = (*verts)[vi][0];
        const REAL y0 = (*verts)[vi][1];
        const REAL z0 = (*verts)[vi][2];
        y[0] = v[0][0] * x0 + v[1][0] * y0 + v[2][0] * z0;
        y[1] = v[0][1] * x0 + v[1][1] * y0 + v[2][1] * z0;
        y[2] = v[0][2] * x0 + v[1][2] * y0 + v[2][2] * z0;
        */
        Matrix<REAL, 3, 1> y  // deafault value. It will be changed later.
            = V * v0;
        /*
        utb[0] = - u[0][0] * b[0] - u[1][0] * b[1] - u[2][0] * b[2];
        utb[1] = - u[0][1] * b[0] - u[1][1] * b[1] - u[2][1] * b[2];
        utb[2] = - u[0][2] * b[0] - u[1][2] * b[1] - u[2][2] * b[2];
        */
        Matrix<REAL, 3, 1> utb
            = -U.transpose() * b;

        // Solves Ax = b,   i.e. minimzes |Ax-b|
        // ?? svd.solve(utb);

        // std::clog << "rank=" << rank << " S( " << S(0) << ", " << S(1) << ", " << S(2) << " ) " << std::endl;

        if (rank < 3) {
            if (VERBOSE_QEM)
                std::clog << " > S(rank+1-1)" << S(rank+1-1) << std::endl;
            for (int i = rank; i < 3; ++i) {
                // S(i) = 0.0;  // not necessary really
            }
            // assert( S(rank+1-1) == 0.0 );
        }

        #if DEBUG_VERBOSE
        clog << " rank=" << rank << " ";
        #endif

        for (int i = rank; i < 3; ++i) {
            S(i) = 0.0;  // not necessary really
        }
        const REAL MIN_EIGENVALUE = 0.000001;
        if (S(0) <= MIN_EIGENVALUE) {
            //rank = 0
            //new_verts[vi, 0:3] = new_x[:, 0]
        }

        // std::clog << "!4" << std::endl;
        // std::clog << " > S(0)" << S(0) << std::endl;
        // std::clog << "!4.5" << std::endl;

        y = V * v0; //default value. Is this needed?
        for (int i=0; i < rank; i++) {
            //if (S(i) != 0.0)
            y(i) = utb(i) / S(i);
            /*
            Never happens:
            if (S(i) == 0.0) {
                // y(i) = 0.0;
                clog "utb(i) " << utb(i) << " ";
            }
            */
        }

        for (int i=rank; i < 3; i++) {
            // y(i) = utb(i) / S(i);
            //if (S(i) == 0.0) {  // always
            if (!super_quiet)
                clog << "utb(i) " << utb(i) << " ";
            // There are two solutions for this:
            // // 1- y(i) = 0.0
            // // 2- y(i) = default
            //if (!DEFAUT_SOLUTIION)
            // y(i) = 0.0;

            //}
        }

        // std::clog << "!5" << std::endl;

        /*
        new_x[0] = v[0][0] * y[0] + v[0][1] * y[1] + v[0][2] * y[2];
        new_x[1] = v[1][0] * y[0] + v[1][1] * y[1] + v[1][2] * y[2];
        new_x[2] = v[2][0] * y[0] + v[2][1] * y[1] + v[2][2] * y[2];

        (*verts)[vi][0] = new_x[0];
        (*verts)[vi][1] = new_x[1];
        (*verts)[vi][2] = new_x[2];
        */

        Matrix<REAL, 3, 1> new_x = V.transpose() * y + qem_origin;

        if (maximum_displacement_distance > 0) {
            REAL dx = new_x(0) - (*verts)[vi][0];
            REAL dy = new_x(1) - (*verts)[vi][1];
            REAL dz = new_x(2) - (*verts)[vi][2];

            REAL dist2 = norm_squared(dx,dy,dz);

            if (dist2 <= maximum_displacement_distance*maximum_displacement_distance) {
                (*verts)[vi][0] = new_x(0);
                (*verts)[vi][1] = new_x(1);
                (*verts)[vi][2] = new_x(2);
            } else {
                // move "a bit"

                REAL dist = std::sqrt(dist2);

                // factor 1.5: works better on large objects (large MC step)
                //  factor 0.9: problems on high-resolution objection (small MC step)
                // togo: put the factor in maximum_displacement_distance (now it has discontinuity)

                REAL displacement_len = maximum_displacement_distance * 1.5; //* 0 .9; // * 0.5;
                if (displacement_len > dist) {
                    displacement_len = dist;
                }


                (*verts)[vi][0] += dx / dist * displacement_len;
                (*verts)[vi][1] += dy / dist * displacement_len;
                (*verts)[vi][2] += dz / dist * displacement_len;
            }

        } else {
            (*verts)[vi][0] = new_x(0);
            (*verts)[vi][1] = new_x(1);
            (*verts)[vi][2] = new_x(2);
        }

        if (store_ranks) {
            (*ranks_output)[vi] = rank;
        }

    }

    #if DEBUG_VERBOSE
    clog << std::endl;
    #endif
    if (!super_quiet)
        clog << std::endl;

}

}  // namespace
