/*
Copyright 2016 MyMiniFactory Ltd.
author: Marc, Solene, Sohail
*/
#pragma once


#include <iostream>
#include <math.h>
#include <cassert>
#include <map>

// #include <vector>
#include <string>
#include <tuple>
#include <fstream>

#include "boost/multi_array.hpp"
#include "boost/array.hpp"

// #include "../js_iteration_2/basic_data_structures.hpp"

#include "../js_iteration_2/vectorised_algorithms/make_random_pm1.hpp"
#include "../js_iteration_2/vectorised_algorithms/add_inplace.hpp"
#include "../js_iteration_2/vectorised_algorithms/cross_product.hpp"
#include "../js_iteration_2/vectorised_algorithms/normalise_inplace.hpp"
#include "../js_iteration_2/vectorised_algorithms/assert_are_normalised.hpp"
using mp5_implicit::vectorised_algorithms::assert_are_normalised;


#include "../js_iteration_2/vectorised_algorithms/misc.hpp"
using mp5_implicit::vectorised_algorithms::replace_zero_normals_with_gaussian_random;
using mp5_implicit::vectorised_algorithms::fill_vector;
using mp5_implicit::vectorised_algorithms::set_a_b_if_c;
using mp5_implicit::vectorised_algorithms::assign_vects_chosen_by_fancy_indexing;
using mp5_implicit::vectorised_algorithms::build_range_array;
using mp5_implicit::vectorised_algorithms::set_array_to_boolean_value;


#include "../js_iteration_2/implicit_vectorised_algorithms.hpp"
using mp5_implicit::get_signs;
using mp5_implicit::produce_facet_normals;
using mp5_implicit::compute_centroid_gradient;

using mp5_implicit::eval_implicit_on_selected_points_indexed;
using mp5_implicit::test_if_conjugate_opposite_signs_indexed;
using mp5_implicit::check_all_are_root;

#include "../js_iteration_2/matrix_functions.hpp"

#include "../js_iteration_2/faces_verts_algorithms.hpp"

using namespace std;

namespace mp5_implicit {

//#undef ASSERT_USED
//#define ASSERT_USED 0
//#define assert(x) {}

REAL compute_average_edge_length(const faces_t& faces, const verts_t& verts) {
    int nfaces = faces.shape()[0];
    REAL edge_length;
    for (int j=0; j < nfaces; j++) {
        auto f0 = faces[j][0];
        auto f1 = faces[j][1];
        auto f2 = faces[j][2];
        edge_length += norm_2(verts[f0][0] - verts[f1][0], verts[f0][1] - verts[f1][1], verts[f0][2] - verts[f1][2]);
        edge_length += norm_2(verts[f0][0] - verts[f2][0], verts[f0][1] - verts[f2][1], verts[f0][2] - verts[f2][2]);
        edge_length += norm_2(verts[f2][0] - verts[f1][0], verts[f2][1] - verts[f1][1], verts[f2][2] - verts[f1][2]);
    }
    return edge_length/(3.*nfaces);
}


/*void compute_centroids(const faces_t& faces, const verts_t& verts, verts_t& centroids) {
    int nt = faces.shape()[0];
    for (int j=0; j < nt; j++) {
        auto f0 = faces[j][0];
        auto f1 = faces[j][1];
        auto f2 = faces[j][2];
        for (int di=0; di < 3; di++) {
                centroids[j][di] = (verts[f0][di] + verts[f1][di] + verts[f2][di])/3.;

        }
    }
}
*/
inline REAL my_sign(REAL v, REAL ROOT_TOLERANCE) {
//     return np.sign(v) * (np.abs(v) > ROOT_TOLERANCE)
  /*
    return
      (v > 0) ?
      (+1) * (v > ROOT_TOLERANCE) :
      (v < 0) ?
      (-1) * (-v > ROOT_TOLERANCE) :
      //(v==0)
      0.0; //(0) * (v > ROOT_TOLERANCE);
  */
      /*
      REAL sgn = (v > 0) ? +1.0 : (v < 0)? -1.0 : 0.0;
      REAL vabs = (v > 0) ? v : (-v);
      // bool (v >= 0) ? (v > ROOT_TOLERANCE) : (-v > ROOT_TOLERANCE);
      REAL (vabs > ROOT_TOLERANCE) ? +1 : -1;
      return (v > 0) ? sgn
      */
    /*
    REAL r;
    if ( v > +ROOT_TOLERANCE ) {
        r = +1;
    } else if ( v < -ROOT_TOLERANCE ) {
        r = -1;
    } else {
        r = 0.0;
    }
    return r;
    */
    return
        (v > +ROOT_TOLERANCE)?
            (+1) :
        (v < -ROOT_TOLERANCE)?
            (-1)
        :
            (0.0);
}


#include "../js_iteration_2/polygoniser/bisection.hpp"

/*
namespace mp5 {
    int global_rng_seed = 12;
};
verts_t make_random_pm1(vindex_t n, int dims, REAL amplitude) {
    //verts_t result {n, dims};
    verts_t result {boost::extents[n][dims]};
    //mt11213b r = boost::mt11213b();
    int seed = mp5::global_rng_seed;
    boost::random::mt11213b rngen(seed);
    //boost::random::uniform_
    boost::random::uniform_01<REAL> distr;
    for (vindex_t i = 0; i < n; ++i) {
        result[i][0] = distr(rngen) * amplitude;
        result[i][1] = distr(rngen) * amplitude;
        result[i][2] = distr(rngen) * amplitude;
    }

    for (vindex_t i = 0; i < 10; ++i) {
        std::clog
            << result[i][0]
            << result[i][1]
            << result[i][2]
            << std::endl;
    }
    return result;
}
*/


std::vector<REAL> make_alpha_list(REAL initial_step_size, REAL min_step_size, REAL max_dist, int max_iter, bool EXTREME_ALPHA) {
    // Prepare a list of step sizes.
    // vectorized_scalar alpha_list(scalar_shape);
    // max_dist is not the absolute strict maximum length. It is the 1- A measure of full search, 2- a unit of length.
    REAL unit_of_length = max_dist;
    std::vector<REAL> alpha_list;
    REAL step_size = initial_step_size;
    // int iter = 0;
    assert(step_size > min_step_size);
    while (step_size > min_step_size) {

        step_size = step_size * 0.5;
        // maximum number of steps to pave max_dist with the stepsize, if not limited by max_iter.
        int total_steps = (int)(std::floor(max_dist / std::abs(step_size)+0.001));
        int max_steps = min(max_iter, total_steps);

        for (int i=1; i < max_steps + 1; i += 2) {
            REAL alpha = ((REAL)i) * step_size;
            alpha_list.push_back(alpha / unit_of_length);
            alpha_list.push_back(-alpha / unit_of_length);
            /*
            alpha_list[i + iter] = alpha/unit_of_length;
            alpha_list[i + iter +1] = -alpha/unit_of_length;
            */
        }
        /*
        iter += max_steps;
        */
    }
    // alpha_list.resize(boost::extents[iter]);
    // However, the following violates the max_dim condision, which is necessary. Caused problems.
    if (EXTREME_ALPHA) {
        std::vector<REAL> more_alphas {+1, -1, +1.5, -1.5, +2, -2};
        alpha_list.insert(std::end(alpha_list), std::begin(more_alphas), std::end(more_alphas));
        // alpha_list.append(+1);
    }

    #if ASSERT_USED
        // todo: D.R.Y.
        // if (VERBOSE)
        std::clog << "Alphas: ";
        for (std::vector<REAL>::iterator i = std::begin(alpha_list), e=std::end(alpha_list); i < e; ++i) {
            std::clog << *i << " ";
        }
        std::clog << std::endl;
    #endif

    return alpha_list;
}

/* creates a bundle of directions on surface, one for each centroid. */
void create_directions_bundle(
    //inputs
    const int iter_type, const std::vector<REAL>& alpha_list_full,
    const vectorized_vect   & directions_basedon_gradient, const vectorized_vect   & facet_normals_directions,
    //outputs
    verts_t & dxc_output, std::string & name_output, std::vector<REAL> & alpha_list_output
) {
    /*
        soetimes alpha_list_output is not changed. It assumes that it has been changes in previous iterations and  has a correct value. Hence, does nto need update.
        The size of dxc_output is automatically checked (because of boost/STL) in statement: dxc_output = directions_basedon_gradient.
    */


    std::vector<REAL> alpha_list1_10 = std::vector<REAL>(alpha_list_full.begin(), alpha_list_full.begin()+10);
    std::vector<REAL> alpha_list0 = std::vector<REAL>(alpha_list_full.begin(), alpha_list_full.begin()+1);

    int n = directions_basedon_gradient.shape()[0];
    boost::array<int, 2> vector_shape = {n, 3};
    assert( directions_basedon_gradient.shape()[0] == n);
    assert( directions_basedon_gradient.shape()[0] == n);
    assert( facet_normals_directions.shape()[0] == n);
    assert( dxc_output.shape()[0] == n);

    if (iter_type == 0) {
        /*
            Examine the directions based on gradients (gradient-normals).
        */
        name_output = "implicit normals toward surface";
        dxc_output = directions_basedon_gradient;
        alpha_list_output = alpha_list_full;  // copy

    } else if (iter_type == 1) {
        /*
            Examine directions based on mesh (mesh-normals).
        */
        name_output = "mesh normals";
        dxc_output = facet_normals_directions; // copy
        alpha_list_output = alpha_list1_10;

    } else if (iter_type == 2) {
        /*
            Examine a direction (III) (almost-) perpendicular to both directions in part I and II.
            Find any direction perpendicular to the gradient-normals.
            This is created by adding some random direction to the gradient-normals
            and cross-product-ing it with the mesh normals. (Why mesh-normals?).
            The result will be perpendicular to both direction I and direction II.
        */
        name_output = "...2";
        int count = directions_basedon_gradient.shape()[0];
        REAL R = 0.000001; // too small
        //REAL R = 0.001;
        verts_t perturb = vectorised_algorithms::make_random_pm1(count, 3, R);
        vectorised_algorithms::add_inplace(perturb, directions_basedon_gradient);
        verts_t z = verts_t(vector_shape);
        assert(assert_are_normalised(facet_normals_directions) && "*mesh normals*");

        vectorised_algorithms::cross_product(facet_normals_directions, perturb, z);
        // print_vector("**meshnormals", facet_normals_directions, 10);
        // print_vector("**perturb", perturb, 10);
        // print_vector("*z", z, 10);
        vectorised_algorithms::normalize_1111(z);
        // print_vector("*z normalized", z, 10);
        /*
        #if ASSERT_USED
        vectorised_algorithms::assert_are_normalised(z);
        #endif
        */
        assert(assert_are_normalised(z) && "1");

        clog <<
            dxc_output.shape()[0] << " x " << dxc_output.shape()[1] << dxc_output.shape()[2] <<
            "  =?=  " <<
            facet_normals_directions.shape()[0] << " x " << facet_normals_directions.shape()[1] << facet_normals_directions.shape()[2] <<
            std::endl;

        assert(dxc_output.shape()[0] == facet_normals_directions.shape()[0]); // remove this and make it inplace directly in dxc_output. or return using move constructor.
        // assert(dxc_output.shape() == facet_normals_directions.shape());  // FAILS?!!! yes. Becasue these are pointers.
        assert( vectorised_algorithms::sizes_are_equal(dxc_output, facet_normals_directions) );
        dxc_output = facet_normals_directions; // copy. checks size.
        alpha_list_output = alpha_list1_10;

        // static auto last_dxc = dxc_output;

    } else if (iter_type == 3) {
        /*
        Direction IV:
        Find the normals based on III & II. The result will not be direction I.
        Reason: (?).
        */
        int count = directions_basedon_gradient.shape()[0];
        //facet_normals_directions: should be (already) normalised, but allowing some to have zero length
        //missing = indices_of_zero_normals(facet_normals_directions);
        verts_t mesh_normals_modifiable = facet_normals_directions;
        replace_zero_normals_with_gaussian_random(mesh_normals_modifiable);
        verts_t z2 = verts_t(vector_shape);
        assert(z2.shape()[0] == mesh_normals_modifiable.shape()[0]);
        assert(z2.shape()[1] == mesh_normals_modifiable.shape()[1]);
        //z = ???????????????????
        //set_vector_from(z, dxc_output); // last dxc_output
        //verts_t z = dxc_output;  // copy
        const verts_t& z = dxc_output; // last dxc_output
        // print_vector("mesh_normals_modifiable", mesh_normals_modifiable, 100);
        // print_vector("z", z, 10);
        vectorised_algorithms::cross_product(mesh_normals_modifiable, z, z2);
        // print_vector("z2", z2, 10);
        vectorised_algorithms::normalise_inplace(z2, mp5_implicit::CONFIG_C::center_projection::min_gradient_len);
        // print_vector("z2 after normalisation", z2, 10);
        assert(assert_are_normalised(z2));
        dxc_output = z2;  // copy
        //alpha_list_output = same as before

        //****


    } else if (iter_type == 4 || iter_type == 5 || iter_type == 6) {
        /*
        Directions V, VI, VII: simply parallel to X,Y,Z axes.
        */
        int dimi = iter_type - 4;
        if (dimi == 0) {
            fill_vector(dxc_output, 1, 0, 0);
        } else if (dimi == 1) {
            fill_vector(dxc_output, 0, 1, 0);
        } else if (dimi == 2) {
            fill_vector(dxc_output, 0, 0, 1);
        }
        //alpha_list_output = same as before

    } else {
        cerr << "Error. incorrect direction type." << std::endl;
        assert(0);
    }
    // dxc_output = direction

    // (alpha_list_output, dxc_output) = directions(iter_type, directions_basedon_gradient, facet_normals_directions, ***);

    assert(alpha_list_output.size() > 0);
}

// main function version v3s002
void  set_centers_on_surface(
      mp5_implicit::implicit_function* object,
      const verts_t & centroids,
      const REAL average_edge,
      //nones_map
      const verts_t & facet_normals_directions,
      vectorized_bool& treated,
      verts_t & centroids_output) {

    /*
    Note that it's allowed to use centroids for centroids_output
    */

    // intilization, objects creation

    constexpr REAL min_gradient_len = 0.000001;  // Gradients smaller than this are considered zero.
    constexpr int max_iter = 20;
    constexpr bool USE_MESH_NORMALS = true;
    constexpr bool EXTREME_ALPHA = false;  // keep false

    const verts_t& X = centroids;

    REAL max_dist = average_edge;

    const int n = X.shape()[0];
    vectorized_scalar_shape scalar_shape = {n};
    // auto scalar_shape = boost::extents[n];

    vectorized_scalar fc_a(scalar_shape);
    object->eval_implicit(centroids, &fc_a);


    /*Compute the gradients, and normalise them.*/
    // g_direction_a:
    // g_a:
    boost::array<int, 2> g_a_shape = { n , 3 };
    verts_t  g_a(g_a_shape);
    // compute_centroid_gradient(X, g_a, object);
    object->eval_gradient(X, &g_a);
    mp5_implicit::vectorised_algorithms::normalise_inplace(g_a, mp5_implicit::CONFIG_C::center_projection::min_gradient_len);
    // print_vector("g_a", g_a, 100);
    assert(assert_are_normalised(g_a) && "0");
    // now g_a --> g_direction
    // Note: Elements will be modified
    auto& g_direction_a = g_a;
    // how to undefine g_a here?


    // todo: remove: negative_f_c

    vectorized_scalar signs_c = get_signs(fc_a, ROOT_TOLERANCE);


    assert(g_direction_a.shape()[0] == fc_a.shape()[0]);

    for (int i = 0, e = g_direction_a.shape()[0]; i < e; i++) {
        if (signs_c[i] < 0.0) {
            // i.e., if fc_a[i] < -ROOT_TOLERANCE
            g_direction_a[i][0] = - g_direction_a[i][0];
            g_direction_a[i][1] = - g_direction_a[i][1];
            g_direction_a[i][2] = - g_direction_a[i][2];
        }
    }

    /*
        dx0_c_grad: the DIRECTION that moves toward the surface, based on GRADIENTs. With the hope that the point and point+direction will be on two different sides of the solid.
    */
    // Move toward surface: If inside (the f value is positive) move outwards (oppoosite the -gradient), and if outside, move inwards (+gradient).
    boost::array<int, 2> vector_shape = {n, 3};
    vectorized_vect  dx0_c_grad(vector_shape);
    for (int i=0; i < fc_a.shape()[0]; i++) {
        // Bug fixed: negative sign missing.
        dx0_c_grad[i][0] = - g_direction_a[i][0]*signs_c[i];
        dx0_c_grad[i][1] = - g_direction_a[i][1]*signs_c[i];
        dx0_c_grad[i][2] = - g_direction_a[i][2]*signs_c[i];
    }
    // Problem: directions_basedon_gradient is basically EQUAL to "g_a"
    const vectorized_vect   & directions_basedon_gradient = dx0_c_grad;



    // stepsize happens to be equal to the max_dist.
    const REAL length_factor = max_dist;  // average_edge; // max_dist
    REAL initial_step_size = max_dist * 1.0;
    std::vector<REAL> alpha_list_full = make_alpha_list(initial_step_size, 0.001, max_dist, max_iter, EXTREME_ALPHA);
    /*
        Interesting observation: alpha_list[4] always finds many points. We can bring it forward, after [1] or [2].
        If we use mesh normals first, the results seem better.
    */

    /* ******************************************************************************************** */
    // THE algorithm
    /* ******************************************************************************************** */

    // array definition
    vectorized_vect   best_result_x (vector_shape);
    assert(best_result_x.shape()[0] == n);

    array_of_indices  active_indices = build_range_array(n);

    /*
    boost::multi_array<vindex_t, 1> active_indices(scalar_shape);
    for (int i=0; i < n; i++) {
        active_indices[i] = i;
    }
    */

    int active_count = n;

    array_of_indices  still_nonsuccess_indices(scalar_shape);
    int still_nonsuccess_indices___s_n_s = active_indices.shape()[0];
    //TODO RENAME  still_nonsuccess_indices___s_n_s --->  still_nonsuccess_indices_size
    still_nonsuccess_indices = active_indices;  // copy
    assert(std::begin(still_nonsuccess_indices) != std::begin(active_indices));  // assure it's a correct copy

    /*
    vectorized_bool  already_success(scalar_shape);
    for (int i=0; i < n; i++) {
        already_success[i] = b_false;
    }


    vectorized_bool  success(scalar_shape);
    //vectorized_bool  already_success(scalar_shape);
    for (int i=0; i < n; i++) {
        success[i] = b_false;  // = already_success[i]
    }

    // (vectorized_bool::value_type)(b_false)
    //set_array_to_value<vectorized_bool::value_type, b_false>(already_success);
    //set_array_to_value<vectorized_bool::value_type, b_false>(success);
    */

    vectorized_bool  already_success(scalar_shape);
    vectorized_bool  success(scalar_shape);

    set_array_to_boolean_value(already_success, b_false);
    set_array_to_boolean_value(success, b_false);

    //already_success must be all false here.

    assert(USE_MESH_NORMALS);
    //if (USE_MESH_NORMALS) {  // true
        //const verts_t & facet_normals_directions = facet_normals_directions;

        assert(assert_are_normalised(facet_normals_directions));
        /*
        #ifdef ASSERT_USED
        assert_are_normalised(facet_normals_directions);
        #endif
        */

        /*
        // refactor: assert_normalised(facet_normals_directions)

        #ifdef ASSERT_USED
        for (vindex_t i = 0, e = facet_normals_directions.shape()[0]; i < e; i++) {
            REAL norm2 = norm_squared(facet_normals_directions[i][0], facet_normals_directions[i][1], facet_normals_directions[i][2]);
            //assert(norm2 != 0.0);
            assert(std::abs(norm2 - 1.0) < mp5_implicit::vectorised_algorithms::ALL_CLOSE_EPS);
        }
        #endif
        */
    //  }


    /*
    *******************
    todo: these variables.
    make a loop
    refactor into functions
    *******************
    */
//    vectorized_vect   xa4(vector_shape);
/*
    vectorized_vect x1_half(vector_shape);

    vectorized_scalar f_a(scalar_shape);
    vectorized_scalar signs_a(scalar_shape);

    // boolean
    vectorized_bool  success0(scalar_shape);

    // indices arrays
*/

    std::vector<int> cases;
    if (USE_MESH_NORMALS) {
        cases = std::vector<int>{0, 1, 2, 3, 4, 5, 6};
    } else {
        cases = std::vector<int>{0};
    }

    // todo: move alpha_list here
    for (auto iter_type : cases) {

        std::vector<REAL> alpha_list1; // = std::vector<int>(alpha_list.begin(), alpha_list.begin());
        // alpha_list0, alpha_list1_10, alpha_list_full

        // always copy
        verts_t dxc = verts_t(vector_shape);
        std::string name;

        /**********************
        ***********************/
        create_directions_bundle(
            //inputs
            iter_type, alpha_list_full,
            directions_basedon_gradient, facet_normals_directions,
            // outputs
            dxc, name, alpha_list1
        );

        assert(alpha_list1.size() > 0);
        /**********************
        ***********************/


        // ********************
        //  ALPHA LOOP
        // ********************
        // input: dxc = direction at centroids

        // ?
        //int still_nonsuccess_indices___s_n_s = 0;

        for (int alpha_i =0, counter = -1; alpha_i < alpha_list1.size(); alpha_i++) {
            counter += 1;
            assert(counter == alpha_i);

            REAL alpha = alpha_list1[alpha_i];

            ///////////// TODO: WE TO DEFINE THIS?
            // moved outside
            vectorized_vect  x1_half(vector_shape);
            REAL c = (length_factor * alpha);
            assert(X.shape()[0] == dxc.shape()[0]);

            // todo: ONLY ACTIVE INDICES

            // for (int j=0; j < n; j++) {
            for (int j=0; j < X.shape()[0]; j++) {
                x1_half[j][0] = X[j][0] + c * dxc[j][0];
                x1_half[j][1] = X[j][1] + c * dxc[j][1];
                x1_half[j][2] = X[j][2] + c * dxc[j][2];
            }

            // simply: active_indices = still_nonsuccess_indices
            clog << still_nonsuccess_indices___s_n_s << "==" << still_nonsuccess_indices.shape()[0] << std::endl;
            assert(still_nonsuccess_indices___s_n_s == still_nonsuccess_indices.shape()[0]);
            active_indices.resize(boost::extents[still_nonsuccess_indices.shape()[0]]);
            active_indices = still_nonsuccess_indices;

/*
            ////////////// REPEATED CODE.
            active_indices = still_nonsuccess_indices;
            int still_nonsuccess_indices__effective_size = 0;
            int active_indices_size = still_nonsuccess_indices__effective_size;
*/
            vectorized_scalar f_a = eval_implicit_on_selected_points_indexed(object, x1_half, active_indices, active_indices.shape()[0]);
            /*
            vectorized_scalar f_a(scalar_shape);
            {
                // assign_vects_chosen_by_fancy_indexing(xa4, x1_half, active_indices, ???);
                int ac = active_indices.shape()[0];
                vectorized_vect_shape shape = {ac, 3};
                vectorized_vect   xa4(shape);
                for (int j=0; j < ac; j++) {
                    auto k = active_indices[j];
                    xa4[j][0] = x1_half[k][0];
                    xa4[j][1] = x1_half[k][1];
                    xa4[j][2] = x1_half[k][2];
                }
                //object->eval_implicit(xa4.begin(), xa4.begin()+ac, f_a.begin());

                ///////////// TODO: WE TO DEFINE THIS?
                // vectorized_scalar f_a(scalar_shape);
                object->eval_implicit(xa4, &f_a);
                clog << "Evaluating " << active_indices.shape()[0] << " points" << std::endl;
            }
            */

            vectorized_scalar
                signs_a = get_signs(f_a, ROOT_TOLERANCE);



            // not checked
            /*
            for (int j=0; j < active_indices.shape()[0]; j++) {
                success[j] = b_false;

                if (signs_a[j] * signs_c[active_indices[j]] <=0) {
                    success0[j] = b_true;
                } else {
                    success0[j] = b_false;
                }

            }

            for (int j=0; j < active_indices.shape()[0]; j++) {
                success[active_indices[j]] = success0[j];
            }
            */
            /*
            // not checked
            for (int j=0; j < success.shape()[0]; j++) {
                success[j] = b_false;
            }
            */
            set_array_to_boolean_value(success, b_false);

            for (int j=0; j < active_indices.shape()[0]; j++) {
                success[active_indices[j]] = (signs_a[j] * signs_c[active_indices[j]]) <= 0;
            }

            // NOT CHECKED

            array_of_indices  new_success_indices(scalar_shape);
            int new_success_indices___n_s = 0;
            //array_of_indices_struct  new_success_indices(scalar_shape);

            // NOT CHECKED


            still_nonsuccess_indices___s_n_s = 0;
            for (int j=0; j < active_indices.shape()[0]; j++) {
                if (success[j] && !already_success[j]) {
                    new_success_indices[new_success_indices___n_s] = j;
                    new_success_indices___n_s ++;
                }
                if (!success[j] && !already_success[j]) {
                    still_nonsuccess_indices[still_nonsuccess_indices___s_n_s] = j;
                    still_nonsuccess_indices___s_n_s ++;
                }
            }
            // not really necessary
            still_nonsuccess_indices.resize(boost::extents[still_nonsuccess_indices___s_n_s]);


            for (int j=0; j < new_success_indices___n_s; j++) {
                auto k = new_success_indices[j];
                // todo: re-factor
                best_result_x[k][0] = x1_half[k][0];
                best_result_x[k][1] = x1_half[k][1];
                best_result_x[k][2] = x1_half[k][2];
            }

            // already_success = already_success OR success
            // Monotonically increasing. Accumulating success.
            assert(already_success.shape()[0] == n);
            assert(success.shape()[0] == n);
            for (int j=0; j < n; j++) {
                if (success[j] == b_true) {
                    already_success[j] = b_true;
                }
            }

            clog << "[" << counter << "](+" << new_success_indices___n_s << ")" << still_nonsuccess_indices___s_n_s << " ";

            #if ASSERT_USED
                assert(test_if_conjugate_opposite_signs_indexed(object, X, best_result_x, already_success, ROOT_TOLERANCE ));
            #endif

            // Code that is Different to bisection_vectorized5_()


            // MOVED BELOW
            //
            // if (still_nonsuccess_indices___s_n_s == 0) {
            //     break;
            // }

            //??
            //active_indices.resize(boost::extents[still_nonsuccess_indices___s_n_s]);
            //still_nonsuccess_indices.resize(boost::extents[still_nonsuccess_indices___s_n_s]);

            if (still_nonsuccess_indices___s_n_s == 0) {
                break;
            }
        }  // alphas loop

        if (still_nonsuccess_indices___s_n_s == 0) {
            break;
        }

    }  //

        //*********




        // main part of the algor

    /*
        for (.....) {
...



            for (int j=0; j < active_indices.shape()[0]; j++) {
                if (f_a[j] > ROOT_TOLERANCE) {
                    signs_a[j] = +1.;
                } else if (f_a[j] < -ROOT_TOLERANCE) {
                    signs_a[j] = -1.;
                } else {
                    signs_a[j] = 0.;
                }

                success[j] = b_false;

                if (signs_a[j] * signs_c[active_indices[j]] <=0) {
                    success0[j] = b_true;
                } else {
                    success0[j] = b_false;
                }

            }

            for (int j=0; j < active_indices.shape()[0]; j++) {
                success[active_indices[j]] = success0[j];
            }


            int new_success_indices___n_s = 0;
            still_nonsuccess_indices___s_n_s = 0;
            for (int j=0; j < active_indices.shape()[0]; j++) {
                if (success[j] == b_true && already_success[j] == b_false) {
                    new_success_indices[new_success_indices___n_s] = j;
                    new_success_indices___n_s ++;
                }
                //} else {
                //    still_nonsuccess_indices[still_nonsuccess_indices___s_n_s] = j;
                //    still_nonsuccess_indices___s_n_s ++;
                //}

                if (success[j] == b_false && already_success[j] == b_false) {
                    still_nonsuccess_indices[still_nonsuccess_indices___s_n_s] = j;
                    still_nonsuccess_indices___s_n_s ++;
                }
            }

            for (int j=0; j < new_success_indices___n_s; j++) {
                // todo: re-factor
                auto k = new_success_indices[j];
                best_result_x[k][0] = x1_half[k][0];
                best_result_x[k][1] = x1_half[k][1];
                best_result_x[k][2] = x1_half[k][2];
            }

            for (int j=0; j < n; j++) {
                if (success[j] == b_true) {
                    already_success[j] = b_true;
                }
            }

            / Code that is Different to bisection_vectorized5_()


            // MOVED BELOW
            //
            // if (still_nonsuccess_indices___s_n_s == 0) {
            //     break;
            // }

            active_indices.resize(boost::extents[still_nonsuccess_indices___s_n_s]);
            still_nonsuccess_indices.resize(boost::extents[still_nonsuccess_indices___s_n_s]);

            if (still_nonsuccess_indices___s_n_s == 0) {
                break;
            }

        }
    */

    /////////////////////////////////////////////////////////////////////////////////////

    // todo: still_nonsuccess_indices__effective_size
    assert(still_nonsuccess_indices___s_n_s == still_nonsuccess_indices.size());

    assign_vects_chosen_by_fancy_indexing(best_result_x, centroids, still_nonsuccess_indices, still_nonsuccess_indices___s_n_s);  // A = B[C];
    /*
    for (int i=0; i < still_nonsuccess_indices___s_n_s; i++) {
        auto k = still_nonsuccess_indices[i];
        best_result_x[k][0] = centroids[k][0];
        best_result_x[k][1] = centroids[k][1];
        best_result_x[k][2] = centroids[k][2];
    }
    */

    // *************************************

    // vectorized_scalar f1(scalar_shape);
    const auto & f1 = fc_a;

    // vectorized_vect  xa1(vector_shape);
    // vectorized_vect  xa2(vector_shape);

    // todo: use X
    // const auto & xa1 = centroids;  // = x0_v3
    const auto & xa2 = best_result_x;

    // todo: assert f1 == F(xa1)

    vectorized_scalar f2(scalar_shape);
    object->eval_implicit(xa2, &f2);

    vectorized_bool  zeros2_bool(scalar_shape);
    vectorized_bool  zeros1_bool(scalar_shape);
    vectorized_bool  zeros1or2(scalar_shape);
    array_of_indices  relevants_bool_indices(scalar_shape);

    // clog << "3" << endl;



    //****************************
    // UP TO HERE

    // bool_find_zero_scalars(zeros2_bool, f2, ROOT_TOLERANCE);
    for (int i=0; i < n; i++) {
        zeros2_bool[i] = std::abs(f2[i])<= ROOT_TOLERANCE;
        /*

        if (std::abs(f2[i])<= ROOT_TOLERANCE) {
            zeros2_bool[i] = b_true;
        } else {
            zeros2_bool[i] = b_false;
        }
        */
    }

    // todo: turn this into a canonical function
    // bool_find_zero_scalars(zeros1_bool, f1, ROOT_TOLERANCE);
    for (int i=0; i < n; i++) {
        zeros1_bool[i] = std::abs(f1[i]) <= ROOT_TOLERANCE;
        /*
        if (std::abs(f1[i])<= ROOT_TOLERANCE) {
            zeros1_bool[i] = b_true;
        } else {
            zeros1_bool[i] = b_false;
        }
        */
    }

    // todo: turn this into a canonical function
    // best_result_x = x0_v3[zeros1_bool]
    set_a_b_if_c(best_result_x, centroids, zeros1_bool);


    // todo: turn this into a canonical function
    for (int i=0; i < n; i++) {
        zeros1or2[i] = zeros2_bool[i] || zeros1_bool[i];
    /*
    if (zeros2_bool[i] || zeros1_bool[i]) {
        zeros1or2[i] = b_true;
    } else {
        zeros1or2[i] = b_false;
    }
    */
    }

    // bug fixed. Initialisation was missing.
    int r_b = 0;
    for (int i=0; i < n; i++) {
        if (already_success[i] && !zeros1or2[i]) {   // another bug detected here. "!" was missing
            relevants_bool_indices[r_b] = i;
            r_b ++;
        }
    }

    #if ASSERT_USED
    {
        assert(scalar_shape[0] == n);
        vectorized_scalar f_(scalar_shape);
        object->eval_implicit(best_result_x, &f_);
        bool everything_alright = true;
        assert(zeros1or2.size() == n);
        for (int i = 0, e = zeros1or2.size(); i < e; ++i) {
            if (zeros1or2[i]) {
                bool ok = f_[i] <= ROOT_TOLERANCE; // f(best_result_x[i]);
                everything_alright = everything_alright && ok;
            }
        }
    }
    #endif

    int relevants_bool_indices_size = r_b;
    relevants_bool_indices.resize(boost::extents[r_b]);
    int m = relevants_bool_indices.shape()[0];
    assert(m == r_b);

    // template<>
    // void lookup_<func>(f1, zeros1or2)
    #if ASSERT_USED
    {
        bool everything_alright = true;
        for (int i = 0, e = zeros1or2.size(); i < e; ++i) {
            int j = zeros1or2[i];
            bool ok = std::abs(f1[j]) <= ROOT_TOLERANCE;
            everything_alright = everything_alright && ok;
        }
    }
    #endif

    // check all are non-zero
    #if ASSERT_USED
    {
      bool everything_alright = true;
      for (int i = 0, e = relevants_bool_indices.size(); i < e; ++i) {
          int j = relevants_bool_indices[i];
          bool ok = std::abs(f1[j]) > ROOT_TOLERANCE;
          everything_alright = everything_alright && ok;
      }
    }
    #endif



    // Check if all x1 & x2 points have srinctly opposite signs (i.e. non zero)
    #if ASSERT_USED
    {

        /*
        int n_ = x_vectorized.shape()[0];
        boost::array<int, 1> v1_shape = {n_};
        vectorized_scalar v_arr(v1_shape);  // x_vectorized.shape());
        */

        // object.eval_implicit(x1111, &f1);

        bool everything_alright = true;
        for (int i = 0, e = relevants_bool_indices.size(); i < e; ++i) {
            int j = relevants_bool_indices[i];
            auto s1 = my_sign(f1[j], ROOT_TOLERANCE);
            auto s2 = my_sign(f2[j], ROOT_TOLERANCE);
            bool ok = (s1 * s2 < 0);
            everything_alright = everything_alright && ok;
        }
    }
    #endif

    // define and assign x2_relevant and x1_relevant.
    // input: centroids, best_result_x

    boost::array<int, 2> x1_relevant_shape = {m, 3};
    boost::array<int, 1> f1_relevant_shape = {m};
    vectorized_vect  x1_relevant(x1_relevant_shape);
    vectorized_vect  x2_relevant(x1_relevant_shape);

    vectorized_scalar f2_relevants(f1_relevant_shape);

    assert(m == x1_relevant_shape[0]);
    assert(x1_relevant.shape()[0] == m);

    for (int i=0; i < m; i++) {
        x1_relevant[i][0] = centroids[relevants_bool_indices[i]][0];
        x1_relevant[i][1] = centroids[relevants_bool_indices[i]][1];
        x1_relevant[i][2] = centroids[relevants_bool_indices[i]][2];

        x2_relevant[i][0] = best_result_x[relevants_bool_indices[i]][0];
        x2_relevant[i][1] = best_result_x[relevants_bool_indices[i]][1];
        x2_relevant[i][2] = best_result_x[relevants_bool_indices[i]][2];

        /*
        assert("problem");
        // if (0)
        assert(  // "problem" &&
          !(
              x1_relevant[i][0] == x2_relevant[i][0]
              &&
              x1_relevant[i][1] == x2_relevant[i][1]
              &&
              x1_relevant[i][2] == x2_relevant[i][2]
          ));
        */
    }
    // clog << "x1x2" << endl;

    // Now we have x1_relevant and x2_relevant which have the same size.

    object->eval_implicit(x2_relevant, &f2_relevants);

    // TODO: REFACTOR
    // Check the signs are opposite
    #if ASSERT_USED
        vectorized_scalar f1_relevants(f1_relevant_shape);
        // object->eval_implicit(x2_relevant, &f1_relevants_);
        object->eval_implicit(x1_relevant, &f1_relevants);

        clog << x1_relevant.shape()[0] << " " << x1_relevant.shape()[1] << " " << x1_relevant.shape() << "  ,  " <<
            " " << x2_relevant.shape()[0] << " " << x2_relevant.shape()[1] << " " <<
            x2_relevant.shape()[0]<<"/"<<x2_relevant.shape()[1]<<"/"<<x2_relevant.shape()[2]<<"/"<<x2_relevant.shape()[3]<<"/"<<x2_relevant.shape()[4]<<"/"<<x2_relevant.shape()[5]<<"/"<<x2_relevant.shape()[6]<<"/"<<x2_relevant.shape()[7]
            << "  " <<
            x1_relevant.shape()[0]<<"/"<<x1_relevant.shape()[1]<<"/"<<x1_relevant.shape()[2]<<"/"<<x1_relevant.shape()[3]<<"/"<<x1_relevant.shape()[4]<<"/"<<x1_relevant.shape()[5]<<"/"<<x1_relevant.shape()[6]<<"/"<<x1_relevant.shape()[7]
            << endl;
        // clog << "fff" << endl;
        // fails: assert(x1_relevant.shape() == x2_relevant.shape());
        // assert(x1_relevant.shape() == x2_relevant.shape());
        assert(x1_relevant.shape()[0] == x2_relevant.shape()[0]);
        assert(x1_relevant.shape()[1] == x2_relevant.shape()[1]);

        // clog << "ggg" << endl;
        // clog << f1_relevants.size() << " " << f2_relevants.size() << endl;
        assert(f1_relevants.size() == f2_relevants.size());
        // clog << "hhh" << endl;

        // assert np.all(f1_relevants*f2_relevants <= +THRESHOLD_zero_interval)

        for (int i=0; i < m; i++) {
            REAL mult = f2_relevants[i] * f1_relevants[i];
            if (!  (mult <= - ROOT_TOLERANCE*ROOT_TOLERANCE)) {
                clog << mult <<" = " << f2_relevants[i] << " * " << f1_relevants[i] << " tol=" << ROOT_TOLERANCE << "[" << i << "]"<< endl;
            }
            // if (0)
            assert(mult <= + ROOT_TOLERANCE*ROOT_TOLERANCE);
        }
    #endif

    // Swap
    /*
    REAL temp0;
    REAL temp1;
    REAL temp2;
    */

    // clog << "m=" << m << endl;
    int ctr = 0;
    for (int i=0; i < m; i++) {
        // If x2 is inside, swap it. => x1 has to be outside.
        if (f2_relevants[i] < -ROOT_TOLERANCE) {
            // ****************************
            // problem: sometimes they are exactly equal!!
            // if (ctr<10) clog << x2_relevant[i][0] << ", " << x1_relevant[i][0] << " <-> ";
            std::swap(x2_relevant[i][0], x1_relevant[i][0]);
            std::swap(x2_relevant[i][1], x1_relevant[i][1]);
            std::swap(x2_relevant[i][2], x1_relevant[i][2]);
            // if (ctr<10) clog << x2_relevant[i][0] << ", " << x1_relevant[i][0] << endl;

            ctr++;

            /*
            temp0 = x2_relevant[i][0];
            temp1 = x2_relevant[i][1];
            temp2 = x2_relevant[i][2];
            x2_relevant[i][0] = x1_relevant[i][0];
            x2_relevant[i][1] = x1_relevant[i][1];
            x2_relevant[i][2] = x1_relevant[i][2];
            x1_relevant[i][0] = temp0;
            x1_relevant[i][1] = temp1;
            x1_relevant[i][2] = temp2;
            */
        }
    }
    clog << "swapped: " << ctr << endl;

    #if ASSERT_USED
        assert(test_if_points_are_inside(x2_relevant, *object, ROOT_TOLERANCE, true));
        assert(test_if_points_are_outside(x1_relevant, *object, ROOT_TOLERANCE, true));
    #endif

    vectorized_vect  x_bisect(x1_relevant_shape);
    // calling the vectorized bisection
    bisection(object, x_bisect, x1_relevant, x2_relevant, ROOT_TOLERANCE, treated);
    // x1_relevant: outside, x2_relevant: inside


    assert(check_all_are_root(object, x_bisect, m, ROOT_TOLERANCE));
    /*
    #if ASSERT_USED
        // bool evaluate_and_assert_sign<0>(x_bisect, *object);

        // boost::array<int, 2> x1_relevant_shape = {m, 3};
        // boost::array<int, 1> f1_relevant_shape = {m};

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
        assert(ok);
    #endif
        */




    assert(relevants_bool_indices.size() == x_bisect.shape()[0]);

    // clog << "1-centroids_output.begin() " << static_cast<void*>(&centroids_output) << " != " << " centroids.begin():" << static_cast<const void*>(&centroids) << std::endl;
    bool same_address1 = centroids_output.begin() != centroids.begin();

    centroids_output = centroids;  // copy
    // clog << "centroids_output.begin() " << static_cast<void*>(centroids_output.begin()) << " != " << " centroids.begin():" << static_cast<void*>(centroids.begin()) << std::endl;
    // clog << "2-centroids_output.begin() " << static_cast<void*>(&centroids_output) << " != " << " centroids.begin():" << static_cast<const void*>(&centroids) << std::endl;
    // assert(centroids_output.begin() != centroids.begin() );
    bool same_address2 = centroids_output.begin() != centroids.begin();
    assert(same_address1 && same_address2  || !same_address1 && !same_address2);

    // changing the values of the centroids
    assert(m == relevants_bool_indices.size());
    for (int i=0; i < m; i++) {
        centroids_output[relevants_bool_indices[i]][0] = x_bisect[i][0];
        centroids_output[relevants_bool_indices[i]][1] = x_bisect[i][1];
        centroids_output[relevants_bool_indices[i]][2] = x_bisect[i][2];
    }

    assert(n == zeros1or2.size());
    for (int i=0; i < n; i++) {
        if (zeros1or2[i]) {
            centroids_output[i][0] = best_result_x[i][0];
            centroids_output[i][1] = best_result_x[i][1];
            centroids_output[i][2] = best_result_x[i][2];
        }
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



// get the matrix A and b used in vertex_apply_qem
void get_A_b__old(const std::vector<int> & nai, const verts_t& centroids, const verts_t& centroid_gradients, verts_t* A, vectorized_scalar* b) {

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


/*  Function: vertex_apply_qem(verts, faces, centroids, vertex_neighbours_list, centroid_gradients, treated)

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

      vertex_neighbours_list:
          Holds neighbouring faces/triangles for all the vertices.


      centroids_gradients:
          The gradient of the implicit function of the current shape for
          all the centroids.

      treated:
          A list of flags that indicates the process is complete for every index
 */

void vertex_apply_qem__old(
    verts_t* verts, const faces_t faces,
    const verts_t centroids,
    const std::vector< std::vector<int>> vertex_neighbours_list,
    const verts_t centroid_gradients,
    const vectorized_bool& treated)
{
    assert (verts != nullptr);
    /* The next asserts have no significance in C++, only Python */

    /* assert centroids is not None */
    /* assert vertex_neighbours_list is not None */
    /* assert centroids_gradients is not None */

    int nverts = verts->shape()[0];
    /* assert nvert = len(vertex_neighbours_list) */

    boost::array<int, 2> A_shape = { 3 , 3 };

    boost::array<int, 2> b_shape = { 3 , 1 };
    vectorized_scalar b(b_shape);
    vectorized_scalar y(b_shape);
    vectorized_scalar utb(b_shape);
    vectorized_scalar new_x(b_shape);
    verts_t u(A_shape);
    verts_t s(A_shape);
    verts_t v(A_shape);
    verts_t A(A_shape);
    for (int vi=0; vi < nverts; vi++) {
        std::vector<int> nlist;
        for (int i=0; i < vertex_neighbours_list[vi].size(); i++) {
            nlist.push_back(vertex_neighbours_list[vi][i]);
        }

        if (true) {
            std::clog << " remove 'treated' " << std::endl;
            bool skip = false;
            for (int g = 0; g < vertex_neighbours_list[vi].size()-1; g++) {
                if (!treated[vertex_neighbours_list[vi][g]]) {
                    skip = true;
                }
                for (int g1 = g+1 ; g1 < vertex_neighbours_list[vi].size()-1; g1++) {
                    if ((abs(abs(centroid_gradients[vertex_neighbours_list[vi][g]][0]) - abs(centroid_gradients[vertex_neighbours_list[vi][g1]][0])) < 0.000001)
                    &&(abs(abs(centroid_gradients[vertex_neighbours_list[vi][g]][1]) - abs(centroid_gradients[vertex_neighbours_list[vi][g1]][1])) < 0.000001)
                    &&(abs(abs(centroid_gradients[vertex_neighbours_list[vi][g]][2]) - abs(centroid_gradients[vertex_neighbours_list[vi][g1]][2])) < 0.000001))
                    skip = true;
                }
            }



            if (skip) {
                clog << vi << endl;
                continue;
            }
        }

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




// Calculates A nd b for a selection of centroids only (an umbrella).
// Todo: vectorized version.
void get_A_b(
    const std::vector<int> & neighbours_faces,
    const verts_t& centroids,
    const verts_t& centroid_normals,
    const Matrix<REAL, 3, 1> & qem_origin,
    Matrix<REAL, 3, 3> *A, Matrix<REAL, 3, 1> * b) {

    // centroid_normals must be normalised

    int m = neighbours_faces.size();
    // centroids[neighbours_faces]

    *A << 0,0,0, 0,0,0, 0,0,0;
    *b << 0,0,0;

    for (int i=0; i < m; i++) {

        vindex_t ni = neighbours_faces[i];

        REAL Ni_x = centroid_normals[ni][0];
        REAL Ni_y = centroid_normals[ni][1];
        REAL Ni_z = centroid_normals[ni][2];

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
        (*A)(1,0) += nn01;
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
    verts_t* verts, const faces_t faces,
    const verts_t centroids,
    const std::vector< std::vector<int>> vertex_neighbours_list,
    const verts_t centroid_gradients,
    const vectorized_bool& treated)
{
    REAL tau = 680;
    REAL svd_threshold = 1.0 / tau;

    assert (verts != nullptr);

    int nverts = verts->shape()[0];
    assert(nverts = vertex_neighbours_list.size());


    /*
    for( int j = 0; j <vertex_neighbours_list.size(); j++) {
        // clog << vertex_neighbours_list[j] << " "
        for( auto v : vertex_neighbours_list[j]) {
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

        const std::vector<int> & nlist = vertex_neighbours_list[vi];

        Matrix<REAL, 3, 1> qem_origin;
        qem_origin << (*verts)[vi][0], (*verts)[vi][1], (*verts)[vi][2];  // old verts

        get_A_b(nlist, centroids, centroid_gradients, qem_origin, &A, &b);

        /*
        int rank = SVD(A, U, S, V, svd_threshold); // the SVD
        */

        Eigen::JacobiSVD< Matrix<REAL, 3, 3> > svd(A,  Eigen::ComputeFullU | Eigen::ComputeFullV);
        svd.setThreshold( svd_threshold );
        auto S = svd.singularValues();
        U = svd.matrixU();
        V = svd.matrixV();
        int rank = svd.rank();

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

        assert np.all(s[:rank]/s[0] >= 1.0/tau)
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


        // std::clog << "!4" << std::endl;
        // std::clog << " > S(0)" << S(0) << std::endl;
        // std::clog << "!4.5" << std::endl;

        for (int i=0; i < rank; i++) {
            //if (S(i) != 0.0)
            y(i) = utb(i) / S(i);
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

        (*verts)[vi][0] = new_x(0);
        (*verts)[vi][1] = new_x(1);
        (*verts)[vi][2] = new_x(2);

    }

}


void centroids_projection(mp5_implicit::implicit_function* object, std::vector<REAL>& result_verts, const std::vector<int>& result_faces) {
    boost::array<unsigned int, 2> verts_shape = { (unsigned int)result_verts.size()/3 , 3 };
    vectorized_vect  verts(verts_shape);
    boost::array<unsigned int, 2> faces_shape = { (unsigned int)result_faces.size()/3 , 3 };
    boost::multi_array<int, 2> faces(faces_shape);

    int output_verts = 0;
    auto i = result_verts.begin();
    auto e = result_verts.end();
    for ( ; i != e; i++, output_verts++) {
        verts[output_verts][0] = (*i);
        i++;
        verts[output_verts][1] = (*i);
        i++;
        verts[output_verts][2] = (*i);
    }


    int output_faces = 0;
    auto i_f = result_faces.begin();
    auto e_f = result_faces.end();
    for ( ; i_f != e_f; i_f++, output_faces++) {
        faces[output_faces][0] = (*i_f);
        i_f++;
        faces[output_faces][1] = (*i_f);
        i_f++;
        faces[output_faces][2] = (*i_f);
    }

    REAL average_edge;
    average_edge = compute_average_edge_length(faces,  verts);

    vectorized_vect_shape  centroids_shape = { static_cast<vectorized_vect::index>(result_faces.size()/3) , 3 };
    vectorized_vect centroids(centroids_shape);

    if (VERBOSE_QEM)
        clog << "1.centroids: " << centroids.shape()[0] << "x" << centroids.shape()[1] << " : "  << centroids[0][0] << std::endl;

    compute_centroids(faces, verts, centroids);

    // clog << "2.centroids: " << centroids.shape()[0] << "x" << centroids.shape()[1] << " : "  << centroids[0][0] << std::endl;

    if (STORE_POINTSETS)
    {
    verts_t ps1 = centroids;
    if (VERBOSE_QEM)
        clog << "2.centroids: " << ps1.shape()[0] << "x" << ps1.shape()[1] << " : " << ps1[0][0] << std::endl;
    point_set_set.emplace(std::make_pair(std::string("pre_p_centroids"), ps1));
    }


    vectorized_bool_shape  treated_shape = {static_cast<vectorized_vect::index>(result_faces.size())};
    vectorized_bool  treated(treated_shape);
    // IS this necessary? It was missing.
    for (auto it = treated.begin(), e=treated.end(); it != e; ++it) {
        *it = b_false;
    }
    /*
    verts_t mesh_normals(centroids_shape);
    std::clog << "Error: mesh_normals is not initialised" << endl;
    //todo: initialise mesh_normals

    assert(assert_are_normalised(facet_normals));
    mp5_implicit::set_centers_on_surface(object, centroids, average_edge, mesh_normals, treated);
    */
    verts_t facet_normals = produce_facet_normals(faces, verts, true);
    assert(assert_are_normalised(facet_normals));
    mp5_implicit::set_centers_on_surface(object, centroids, average_edge, facet_normals, treated, centroids);


    if (STORE_POINTSETS) {
    verts_t ps2 = centroids;
    if (VERBOSE_QEM)
        clog << "3.centroids: " << ps2.shape()[0] << "x" << ps2.shape()[1] << " : " << ps2[0][0] << std::endl;
    point_set_set.emplace(std::make_pair(std::string("post_p_centroids"), ps2));
    }



    /*
        // temporary
    for (int i=0; i < verts.shape()[0]; i++) {
        int ti = i*2;
        result_verts[ti*3+0] = centroids[i][0];
        result_verts[ti*3+1] = centroids[i][1];
        result_verts[ti*3+2] = centroids[i][2];
    }
    clog << "demarkate" << std::endl;
    return ;
    */


    std::vector< std::vector<int>> vertex_neighbours_list;
    vertex_neighbours_list = make_neighbour_faces_of_vertex(verts, faces);

    vectorized_vect  centroid_gradients(centroids_shape);

    compute_centroid_gradient(centroids, centroid_gradients, object);


    if (STORE_POINTSETS)
    {
    verts_t ps1 = verts;
    if (VERBOSE_QEM)
        clog << "1.verts: " << ps1.shape()[0] << "x" << ps1.shape()[1] << " : " << ps1[0][0] << std::endl;
    point_set_set.emplace(std::make_pair(std::string("pre_qem_verts"), ps1));
    }

    vertex_apply_qem(&verts, faces, centroids, vertex_neighbours_list, centroid_gradients, treated);

    if (STORE_POINTSETS)
    {
    verts_t ps1 = verts;
    if (VERBOSE_QEM)
        clog << "2.verts: " << ps1.shape()[0] << "x" << ps1.shape()[1] << " : " << ps1[0][0] << std::endl;
    point_set_set.emplace(std::make_pair(std::string("post_qem_verts"), ps1));
    }

    // if APPLY_QEM_RESULT
    for (int i=0; i < verts.shape()[0]; i++) {
        result_verts[i*3+0] = verts[i][0];
        result_verts[i*3+1] = verts[i][1];
        result_verts[i*3+2] = verts[i][2];
    }


}

}  // namespace
