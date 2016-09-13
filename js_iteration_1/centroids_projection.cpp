/*
Copyright 2016 MyMiniFactory Ltd.
author: Marc, Solene, Sohail
*/
#pragma once

#include "../js_iteration_2/qem.hpp"


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
using mp5_implicit::vectorised_algorithms::set_all_array_elements_to_a_boolean_value;


#include "../js_iteration_2/implicit_vectorised_algorithms.hpp"
using mp5_implicit::get_signs;
using mp5_implicit::produce_facet_normals;
using mp5_implicit::compute_centroid_gradient;

using mp5_implicit::eval_implicit_on_selected_points_indexed;
using mp5_implicit::test_if_conjugate_opposite_signs_indexed;
using mp5_implicit::check_all_are_root;

// #include "../js_iteration_2/matrix_functions.hpp"

#include "../js_iteration_2/faces_verts_algorithms.hpp"

// #include "../js_iteration_2/qem.hpp"

using namespace std;

namespace mp5_implicit {


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

void print_alpha_list(const std::vector<REAL> & alpha_list, const std::string & title) {
        // todo: D.R.Y.
        // if (VERBOSE)
        std::clog << title << ": (count=" << alpha_list.size() << ") ";
        // std::vector<REAL>::iterator
        for (auto i = std::begin(alpha_list), e=std::end(alpha_list); i < e; ++i) {
            std::clog << *i << " ";
        }
        std::clog << std::endl;
}


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
    /*
        // todo: D.R.Y.
        // if (VERBOSE)
        std::clog << "Alphas: (count=" << alpha_list.size() << ") ";
        for (std::vector<REAL>::iterator i = std::begin(alpha_list), e=std::end(alpha_list); i < e; ++i) {
            std::clog << *i << " ";
        }
        std::clog << std::endl;
    */
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

    clog << "iter_type:" << iter_type << std::endl;
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

        // fixme: This can *FAIL*. The vector in z2 can be zero (it happaned).
        assert(assert_are_normalised(z2));
        dxc_output = z2;  // copy
        //alpha_list_output = same as before

        // BUG FIXED. SETTING  alpha_list_output was missing.
        // was missing
        alpha_list_output = alpha_list1_10;

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
        alpha_list_output = alpha_list1_10;

    } else {
        cerr << "Error. incorrect direction type." << std::endl;
        assert(0);
    }
    // dxc_output = direction

    // (alpha_list_output, dxc_output) = directions(iter_type, directions_basedon_gradient, facet_normals_directions, ***);

    assert(alpha_list_output.size() > 0);
}


void print_out_array_of_indices(const array_of_indices& A, int count_elements, const std::string & title) {

    clog << title << ": [";
    for (int j=0; j < count_elements; j++) {
        auto k = A[j];
        std::clog << k << " ";
    }
    std::clog << "]" << std::endl;
}


void print_out_array_of_bools(const vectorized_bool& A, int count_elements, const std::string & title) {
    clog << title << ": [";
    for (int j=0; j < count_elements; j++) {
        auto k = A[j]?"1":"0";
        std::clog << k << "";
    }
    std::clog << "]" << std::endl;
}

void print_out_array_of_signs(const vectorized_scalar& A, int count_elements, const std::string & title) {
    clog << title << ": [";
    for (int j=0; j < count_elements; j++) {
        auto k = (A[j] > 0)? "+":  (A[j] < 0)? "-" : "0";
        std::clog << k << "";
    }
    std::clog << "]" << std::endl;
}


void print_out_array_of_scalars(const vectorized_scalar& A, int count_elements, const std::string & title) {
    clog << title << ": [";
    for (int j=0; j < count_elements; j++) {
        auto k = A[j];
        k = ((REAL)((long)(k * 100.0)))/100.0;
        std::clog << k << " ";
    }
    std::clog << "]" << std::endl;
}

void print_out_array_of_vectors(const vectorized_vect& A, int count_elements, const std::string & title) {
    clog << title << ": [";
    for (int j=0; j < count_elements; j++) {
        auto x = A[j][0];
        auto y = A[j][1];
        auto z = A[j][2];
        x = ((REAL)((long)(x * 10000.0)))/10000.0;
        y = ((REAL)((long)(y * 10000.0)))/10000.0;
        z = ((REAL)((long)(z * 10000.0)))/10000.0;
        std::clog << x << "," << y << "," << z << " ";
    }
    std::clog << "]" << std::endl;
}


// main function version v3s002
void  set_centers_on_surface(
      mp5_implicit::implicit_function* object,
      const verts_t & centroids,
      const REAL average_edge,
      //nones_map
      const verts_t & facet_normals_directions,
      //vectorized_bool& treated,
      verts_t & centroids_output) {

    /*
    Note that it's allowed to use centroids for centroids_output
    */
    auto NN1 = centroids.shape()[0];
    auto NN2 = centroids_output.shape()[0];

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
    clog << " initial_step_size: " << initial_step_size
    << " max_dist:" << max_dist << std::endl;
    // Goog value equivalent to 0.001, scale
    REAL scale_min_len = 1.0; // initial_step_size / 1.27158;
    std::vector<REAL> alpha_list_full = make_alpha_list(initial_step_size, 0.001*scale_min_len, max_dist, max_iter, EXTREME_ALPHA);
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
    #if ASSERT_USED
        constexpr REAL  MAGIC_VALUE_FOR_DEBUG = -10000.0;
        for (int i = 0; i < best_result_x.shape()[0]; ++i) {
            best_result_x[i][0] = MAGIC_VALUE_FOR_DEBUG;
            best_result_x[i][1] = MAGIC_VALUE_FOR_DEBUG;
            best_result_x[i][2] = MAGIC_VALUE_FOR_DEBUG;
        }
    #endif
    array_of_indices  active_indices = build_range_array(n);
    array_of_indices::index active_indices_size = n;


    /*
    boost::multi_array<vindex_t, 1> active_indices(scalar_shape);
    for (int i=0; i < n; i++) {
        active_indices[i] = i;
    }
    */

    int active_count = n;
    assert( active_indices_size == active_count);

    array_of_indices  still_nonsuccess_indices(scalar_shape);
    assert(active_indices.shape()[0] == active_indices_size);
    int still_nonsuccess_indices___s_n_s = active_indices.shape()[0];
    //TODO RENAME  still_nonsuccess_indices___s_n_s --->  still_nonsuccess_indices_size
    still_nonsuccess_indices = active_indices;  // copy
    assert(std::begin(still_nonsuccess_indices) != std::begin(active_indices));  // assure it's a correct copy

    /*
    vectorized_bool  already_success_bool(scalar_shape);
    for (int i=0; i < n; i++) {
        already_success_bool[i] = b_false;
    }


    vectorized_bool  success(scalar_shape);
    //vectorized_bool  already_success_bool(scalar_shape);
    for (int i=0; i < n; i++) {
        success[i] = b_false;  // = already_success_bool[i]
    }

    // (vectorized_bool::value_type)(b_false)
    //set_array_to_value<vectorized_bool::value_type, b_false>(already_success_bool);
    //set_array_to_value<vectorized_bool::value_type, b_false>(success);
    */

    vectorized_bool  already_success_bool(scalar_shape);
    vectorized_bool  success(scalar_shape);

    set_all_array_elements_to_a_boolean_value(already_success_bool, b_false);
    set_all_array_elements_to_a_boolean_value(success, b_false);

    //already_success_bool must be all false here.

    assert(USE_MESH_NORMALS);
    //if (USE_MESH_NORMALS) {  // true
        //const verts_t & facet_normals_directions = facet_normals_directions;

        assert(assert_are_normalised(facet_normals_directions));
        /*
        #if ASSERT_USED
        assert_are_normalised(facet_normals_directions);
        #endif
        */

        /*
        // refactor: assert_normalised(facet_normals_directions)

        #if ASSERT_USED
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

        #if DEBUG_VERBOSE
        print_out_array_of_vectors(directions_basedon_gradient, directions_basedon_gradient.shape()[0], "directions_basedon_gradient");
        print_out_array_of_vectors(facet_normals_directions, facet_normals_directions.shape()[0], "facet_normals_directions");
        print_alpha_list(alpha_list_full, "Alphas FULL");
        #endif

        /**********************
        ***********************/
        create_directions_bundle(
            //inputs
            iter_type, alpha_list_full,
            directions_basedon_gradient, facet_normals_directions,
            // outputs
            dxc, name, alpha_list1
        );

        #if DEBUG_VERBOSE
        print_out_array_of_vectors(dxc, dxc.shape()[0], "dx");
        print_alpha_list(alpha_list1, "Alphas");
        #endif

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
            // std::clog << "alpha: " << alpha << std::endl;


            // TODO: FIXME:
            // WHEN THE OBJECT IS VERY SMALL, THE "length_factor" IS TOO SMALL.


            ///////////// TODO: WE TO DEFINE THIS?
            // moved outside
            vectorized_vect  x1_half(vector_shape);
            REAL c = (length_factor * alpha) * 4.0;
            assert(X.shape()[0] == dxc.shape()[0]);



            std::clog << "alpha: " << alpha << " c:" << c << " length_factor:" << length_factor << std::endl;

            // todo: ONLY ACTIVE INDICES

            // for (int j=0; j < n; j++) {
            for (int j=0; j < X.shape()[0]; j++) {
                // But we need this only for active points
                // todo: only on chosen ones. ( fixme: )  set_h1_half --> h1_half_of_selected_points
                x1_half[j][0] = X[j][0] + c * dxc[j][0];
                x1_half[j][1] = X[j][1] + c * dxc[j][1];
                x1_half[j][2] = X[j][2] + c * dxc[j][2];
            }

            #if DEBUG_VERBOSE
            print_out_array_of_vectors(x1_half, x1_half.shape()[0], "x1_half");
            #endif

            // simply: active_indices = still_nonsuccess_indices
            clog << still_nonsuccess_indices___s_n_s << "==" << still_nonsuccess_indices.shape()[0] << std::endl;
            assert(still_nonsuccess_indices___s_n_s == still_nonsuccess_indices.shape()[0]);
            active_indices.resize(boost::extents[still_nonsuccess_indices.shape()[0]]);
            active_indices_size = still_nonsuccess_indices___s_n_s;

            active_indices = still_nonsuccess_indices;  // just copy all contents. Tobe simplified later.

/*
            ////////////// REPEATED CODE.
            active_indices = still_nonsuccess_indices;
            int still_nonsuccess_indices__effective_size = 0;
            int active_indices_size = still_nonsuccess_indices__effective_size;
*/
            vectorized_scalar f_a = eval_implicit_on_selected_points_indexed(object, x1_half, active_indices, active_indices.shape()[0]);

            #if DEBUG_VERBOSE
            print_out_array_of_scalars(f_a, f_a.shape()[0], "f(active_indices)");
            // print_out_array_of_vectors(dxc, dxc.shape()[0], "dx");
            #endif

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
                // f_a[] of only active_indices[:]



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
            set_all_array_elements_to_a_boolean_value(success, b_false);

            for (int j=0; j < active_indices.shape()[0]; j++) {
                assert(signs_a.shape()[0] == active_indices_size);
                assert(signs_a.shape()[0] == active_indices.shape()[0]);
                success[active_indices[j]] = (signs_a[j] * signs_c[active_indices[j]]) <= 0;
                // signs_a[j] == object.f(x1_half[active_indices[j]])
                // signs_c: sign(f(x)) at centroid (i.e. before projection)

                    /*
                    // DEBUG THINGY
                    {
                        auto jj = active_indices[j];
                        if (best_result_x[jj][0] == MAGIC_VALUE_FOR_DEBUG || (best_result_x[jj][0] == 0 && best_result_x[jj][1] == 0 && best_result_x[j][2] == 0)) {
                            clog << "[d]is absolute zero. jj=" << jj << std::endl;
                            abort();
                        }
                    }
                    */
            }

            // NOT CHECKED

            array_of_indices  new_success_indices(scalar_shape);
            int new_success_indices___n_s = 0;
            //array_of_indices_struct  new_success_indices(scalar_shape);

            // NOT CHECKED


            still_nonsuccess_indices___s_n_s = 0;
            for (int j=0; j < active_indices.shape()[0]; j++) {
                if (success[j] && !already_success_bool[j]) {
                    new_success_indices[new_success_indices___n_s] = j;
                    new_success_indices___n_s ++;
                }
                if (!success[j] && !already_success_bool[j]) {
                    still_nonsuccess_indices[still_nonsuccess_indices___s_n_s] = j;
                    still_nonsuccess_indices___s_n_s ++;
                }
                    /*
                    // does not fail. why? Because the success[j] is never true
                    // DEBUG THINGY
                    {
                        if(success[j])
                        if (best_result_x[j][0] == MAGIC_VALUE_FOR_DEBUG || (best_result_x[j][0] == 0 && best_result_x[j][1] == 0 && best_result_x[j][2] == 0)) {
                            // happens
                            clog << "[c1]is absolute zero. j=" << j << " " << best_result_x[j][0] << std::endl;
                            abort();
                        }
                        if(already_success_bool[j])
                        if (best_result_x[j][0] == MAGIC_VALUE_FOR_DEBUG || (best_result_x[j][0] == 0 && best_result_x[j][1] == 0 && best_result_x[j][2] == 0)) {
                            clog << "[c2]is absolute zero. j=" << j << " " << best_result_x[j][0] << std::endl;
                            abort();
                        }
                    }
                    */

            }
            // not really necessary
            still_nonsuccess_indices.resize(boost::extents[still_nonsuccess_indices___s_n_s]);

            #if DEBUG_VERBOSE
            print_out_array_of_indices(still_nonsuccess_indices, still_nonsuccess_indices___s_n_s, "still_nonsuccess_indices");

            print_out_array_of_indices(new_success_indices, new_success_indices___n_s, "new_success_indices");

            print_out_array_of_bools(already_success_bool, already_success_bool.shape()[0], "already_success_bool");
            print_out_array_of_bools(success, already_success_bool.shape()[0], "success");

            print_out_array_of_signs(signs_a, signs_a.shape()[0], "signs_a");
            print_out_array_of_signs(signs_c, signs_c.shape()[0], "signs_c");

            print_out_array_of_indices(active_indices, active_indices_size, "active_indices");
            #endif




            // note: new_success_indices[] can be empty

            for (int j=0; j < new_success_indices___n_s; j++) {
                auto k = new_success_indices[j];
                // todo: re-factor
                best_result_x[k][0] = x1_half[k][0];
                best_result_x[k][1] = x1_half[k][1];
                best_result_x[k][2] = x1_half[k][2];
            }

            // already_success_bool = already_success_bool OR success
            // Monotonically increasing. Accumulating success.
            assert(already_success_bool.shape()[0] == n);
            assert(success.shape()[0] == n);
            for (int j=0; j < n; j++) {
                if (success[j]) {  // it should be new_success? ****
                    already_success_bool[j] = b_true;

                    // DEBUG THINGY
                    {
                        if (best_result_x[j][0] == MAGIC_VALUE_FOR_DEBUG || (best_result_x[j][0] == 0 && best_result_x[j][1] == 0 && best_result_x[j][2] == 0)) {
                            clog << "[2]is absolute zero. j=" << j << " " << best_result_x[j][0] << std::endl;
                            abort();
                        }
                    }
                }
            }

            clog << "[" << counter << "](+" << new_success_indices___n_s << ")" << still_nonsuccess_indices___s_n_s << " ";

            // DEBUG THINGY
            {
                std::clog << "already_success_bool.shape()[0] = " << already_success_bool.shape()[0] << std::endl;
                if (already_success_bool.shape()[0] > 0) {
                    for (int j = 0; j < already_success_bool.shape()[0]; ++j) {
                        if(already_success_bool[j])
                          if (best_result_x[j][0] == MAGIC_VALUE_FOR_DEBUG || (best_result_x[j][0] == 0 && best_result_x[j][1] == 0 && best_result_x[j][2] == 0)) {
                            clog << "[-] is absolute zero. j=" << j << " " << best_result_x[j][0] << std::endl;
                            abort();
                        }
                    }
                }
            }

            #if ASSERT_USED
                // fails
                assert(test_if_conjugate_opposite_signs_indexed(
                    object, X, best_result_x, already_success_bool, ROOT_TOLERANCE ));
                // X === centroids
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
                if (success[j] == b_true && already_success_bool[j] == b_false) {
                    new_success_indices[new_success_indices___n_s] = j;
                    new_success_indices___n_s ++;
                }
                //} else {
                //    still_nonsuccess_indices[still_nonsuccess_indices___s_n_s] = j;
                //    still_nonsuccess_indices___s_n_s ++;
                //}

                if (success[j] == b_false && already_success_bool[j] == b_false) {
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
                    already_success_bool[j] = b_true;
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

    /* The points that failed to find their conjugates, just copy the original point (centroids).*/
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




    // clog << "3" << endl;



    //****************************
    // UP TO HERE

    vectorized_bool  zeros2_bool(scalar_shape);
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

    vectorized_bool  zeros1_bool(scalar_shape);
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


    vectorized_bool  zeros1or2(scalar_shape);
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

    array_of_indices  relevants_bool_indices(scalar_shape);
    // bug fixed. Initialisation was missing.
    int r_b = 0;
    for (int i=0; i < n; i++) {
        if (already_success_bool[i] && !zeros1or2[i]) {   // another bug detected here. "!" was missing
            relevants_bool_indices[r_b] = i;
            r_b ++;
        }
    }

    int relevants_bool_indices_size = r_b;
    relevants_bool_indices.resize(boost::extents[relevants_bool_indices_size]);
    int m = relevants_bool_indices_size; //relevants_bool_indices.shape()[0];
    assert(m == r_b);

    // Invariant at this point:
    // [a.]  abs(f2[relevant_bool_indices[:]]) is non zero
    // [b.]  f2[i if zeros1or2[i]] is zero

    // testing [a.]
    #if ASSERT_USED
    {
        // Why this test?
        assert(scalar_shape[0] == n);
        vectorized_scalar f2_(scalar_shape);
        // But it is already done. It can ebe found in f2
        object->eval_implicit(best_result_x, &f2_);  // ??  //why?
        bool everything_alright = true;
        assert(zeros1or2.size() == n);
        for (int i = 0, e = n; i < e; ++i) {
            if (zeros1or2[i]) {
                // why was this done?
                // bool ok = f2_[i] <= ROOT_TOLERANCE; // f(best_result_x[i]);

                bool ok = std::abs(f2_[i]) <= ROOT_TOLERANCE; // f(best_result_x[i]);
                everything_alright = everything_alright && ok;

                assert (f2[i] == f2_[i]);

            }
        }
    }
    #endif

    /*
    // template<>
    // void lookup_<func>(f1, zeros1or2)
    #if ASSERT_USED
    {
        bool everything_alright = true;
        for (int i = 0, e = zeros1or2.size(); i < e; ++i) {
            // what??
            // int j = zeros1or2[i];
            // bool ok = std::abs(f1[j]) <= ROOT_TOLERANCE;

            int j = zeros1or2[i];
            bool ok = std::abs(f1[j]) <= ROOT_TOLERANCE;
            everything_alright = everything_alright && ok;
        }
    }
    #endif
    */

    /*
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
    */

    // testing [b.]
    #if ASSERT_USED
    {
      bool everything_alright = true;
      for (int i = 0, e = relevants_bool_indices_size; i < e; ++i) {
          int j = relevants_bool_indices[i];
          bool ok = std::abs(f2[j]) > ROOT_TOLERANCE;
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
        auto j = relevants_bool_indices[i];
        x1_relevant[i][0] = centroids[j][0];
        x1_relevant[i][1] = centroids[j][1];
        x1_relevant[i][2] = centroids[j][2];

        x2_relevant[i][0] = best_result_x[j][0];
        x2_relevant[i][1] = best_result_x[j][1];
        x2_relevant[i][2] = best_result_x[j][2];

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

        //  f2 = -0.15016 * f1 = -0.00199347        tol=0.001
        for (int i=0; i < m; i++) {
            REAL mult = f2_relevants[i] * f1_relevants[i];
            if (!  (mult <= - ROOT_TOLERANCE*ROOT_TOLERANCE)) {
                clog << mult <<" = " << f2_relevants[i] << " * " << f1_relevants[i] << " tol=" << ROOT_TOLERANCE << " i:[" << i << "]"<< endl;
            }
            // if (0)
            assert(mult <= + ROOT_TOLERANCE*ROOT_TOLERANCE);
            // assert(mult <= + ROOT_TOLERANCE);  // In Python code it is like this.
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
    // note: f2_relevants[] is not valid anymore, because it is now swapped.
    clog << "swapped: " << ctr << endl;

    #if ASSERT_USED
        assert(test_if_points_are_inside(x2_relevant, *object, ROOT_TOLERANCE, true));
        assert(test_if_points_are_outside(x1_relevant, *object, ROOT_TOLERANCE, true));
    #endif

    vectorized_vect  x_bisect(x1_relevant_shape);
    // calling the vectorized bisection
    bisection(object, x_bisect, x1_relevant, x2_relevant, ROOT_TOLERANCE);
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

    // now put them (all of those RELEVEANT ones) into the results array.
    // What about those which were not successful? Answer: They are not "relevant".
    // not Relevant: either the conjugate was not found, or t was already zero (zeros1or2).


    // Shold not be needed if it is not resized:
    assert(relevants_bool_indices.size() == x_bisect.shape()[0]);

    // clog << "1-centroids_output.begin() " << static_cast<void*>(&centroids_output) << " != " << " centroids.begin():" << static_cast<const void*>(&centroids) << std::endl;
    bool same_address1 = centroids_output.begin() != centroids.begin();


    // Copying centroids into centroids_output is necessary, so that when neigher zeros1or2[.] nor relevants_bool_indices[.],
    // ,  ... wait, what?  best_result_x[] is ignored. But it's fine. That was before the bisecion.
    // relevant[] : pased through bisection, hance, correct.
    // zeros1or2[]: not passed through bisection, but correct.
    // rest: either non-successful conjugate or non-successful bisection? But there will be no unsuccessful bisection. It always converges.
    // todo: fix the python code to sepaate the vairables.

    centroids_output = centroids;  // copy
    // clog << "centroids_output.begin() " << static_cast<void*>(centroids_output.begin()) << " != " << " centroids.begin():" << static_cast<void*>(centroids.begin()) << std::endl;
    // clog << "2-centroids_output.begin() " << static_cast<void*>(&centroids_output) << " != " << " centroids.begin():" << static_cast<const void*>(&centroids) << std::endl;
    // assert(centroids_output.begin() != centroids.begin() );
    bool same_address2 = centroids_output.begin() != centroids.begin();
    // Check reference-ness (view-ness)
    assert(same_address1 && same_address2  || !same_address1 && !same_address2);

    // changing the values of the centroids
    assert(m == relevants_bool_indices.size());
    assert(m == relevants_bool_indices_size);
    for (int i=0; i < relevants_bool_indices_size; i++) {
        auto j = relevants_bool_indices[i];
        centroids_output[j][0] = x_bisect[i][0];
        centroids_output[j][1] = x_bisect[i][1];
        centroids_output[j][2] = x_bisect[i][2];
    }

    assert(n == zeros1or2.size());
    for (int i=0; i < n; i++) {
        if (zeros1or2[i]) {
            centroids_output[i][0] = best_result_x[i][0];
            centroids_output[i][1] = best_result_x[i][1];
            centroids_output[i][2] = best_result_x[i][2];
        }
    }
    /*
    // For debugging purpose:
    for (int i=0; i < centroids_output.shape()[0]; i++) {
        centroids_output[i][0] += rand01();
        centroids_output[i][1] += rand01();
        centroids_output[i][2] += rand01();
    }
    */

    // What about the rest of centroids_output?
    // They are copied above.

    auto NN3 = centroids.shape()[0];
    auto NN4 = centroids_output.shape()[0];
    assert(NN1 == NN3);
    assert(NN2 == NN4);
    assert(NN1 == NN2);

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


vectorized_vect  vects2vects(const std::vector<REAL>& result_verts) {
    boost::array<vectorized_vect::index, 2> verts_shape = { (vectorized_vect::index)(result_verts.size()/3) , 3 };
    vectorized_vect  verts(verts_shape);

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
    return verts;
}

void set_vectorverts_from_vectorised_verts(std::vector<REAL>& result_verts, const vectorized_vect & verts) {
    auto n = verts.shape()[0];
    for (int i=0; i < n; i++) {
        result_verts[i*3+0] = verts[i][0];
        result_verts[i*3+1] = verts[i][1];
        result_verts[i*3+2] = verts[i][2];
    }
}

boost::multi_array<int, 2> copy_faces_from_vectorfaces(const std::vector<int> & mesh_faces) {
    // boost::multi_array<int, 2> faces = ;

    vectorized_vect::index  num_faces = static_cast<vectorized_vect::index>(mesh_faces.size()/3);

    boost::array<vectorized_vect::index, 2> faces_shape = { num_faces , 3 };
    boost::multi_array<int, 2> faces(faces_shape);

    int output_faces = 0;
    auto i_f = mesh_faces.begin();
    auto e_f = mesh_faces.end();
    for ( ; i_f != e_f; i_f++, output_faces++) {
        faces[output_faces][0] = (*i_f);
        i_f++;
        faces[output_faces][1] = (*i_f);
        i_f++;
        faces[output_faces][2] = (*i_f);
    }

    return faces;
}


// rename: project_points_on_surface. (output_points)
void centroids_projection(mp5_implicit::implicit_function* object, std::vector<REAL>& result_verts, const std::vector<int>& mesh_faces, bool enable_qem) {

    vectorized_vect  verts = vects2vects(result_verts);

    vectorized_vect::index  num_faces = static_cast<vectorized_vect::index>(mesh_faces.size()/3);

    boost::multi_array<int, 2> faces = copy_faces_from_vectorfaces(mesh_faces);

    REAL average_edge;
    average_edge = compute_average_edge_length(faces,  verts);

    vectorized_vect_shape  centroids_shape = { num_faces , 3 };
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


    /*
    vectorized_bool_shape  treated_shape = {num_faces*3};
    vectorized_bool  treated(treated_shape);
    // IS this necessary? It was missing.
    for (auto it = treated.begin(), e=treated.end(); it != e; ++it) {
        *it = b_false;
    }
    */
    /*
    verts_t mesh_normals(centroids_shape);
    std::clog << "Error: mesh_normals is not initialised" << endl;
    //todo: initialise mesh_normals

    assert(assert_are_normalised(facet_normals));
    mp5_implicit::set_centers_on_surface(object, centroids, average_edge, mesh_normals, treated);
    */
    verts_t facet_normals = produce_facet_normals(faces, verts, true);
    assert(assert_are_normalised(facet_normals));

    mp5_implicit::set_centers_on_surface(object, centroids, average_edge, facet_normals, centroids);


    /*
    if (STORE_POINTSETS) {
    verts_t ps2 = centroids;
    if (VERBOSE_QEM)
        clog << "3.centroids: " << ps2.shape()[0] << "x" << ps2.shape()[1] << " : " << ps2[0][0] << std::endl;
    point_set_set.emplace(std::make_pair(std::string("post_p_centroids"), ps2));
    }
    */
    // LOG_POINTSET("post_p_centroids", centroids);
    STORE_POINTSET("post_p_centroids", centroids);



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
    vertex_neighbours_list = make_neighbour_faces_of_vertex(faces, verts.shape()[0]);

    vectorized_vect  centroid_gradients(centroids_shape);

    compute_centroid_gradient(centroids, centroid_gradients, object);

    /* ***********************
    QEM
    ************************ */
    if (STORE_POINTSETS)
    {
    verts_t ps1 = verts;
    if (VERBOSE_QEM)
        clog << "1.verts: " << ps1.shape()[0] << "x" << ps1.shape()[1] << " : " << ps1[0][0] << std::endl;
    point_set_set.emplace(std::make_pair(std::string("pre_qem_verts"), ps1));
    }

    if (enable_qem) {
        std::clog << "Going for QEM:" << std::endl;
        vertex_apply_qem(&verts, faces, centroids, vertex_neighbours_list, centroid_gradients);

        if (STORE_POINTSETS)
        {
        verts_t ps1 = verts;
        if (VERBOSE_QEM)
            clog << "2.verts: " << ps1.shape()[0] << "x" << ps1.shape()[1] << " : " << ps1[0][0] << std::endl;
        point_set_set.emplace(std::make_pair(std::string("post_qem_verts"), ps1));
        }

        // if APPLY_QEM_RESULT
        set_vectorverts_from_vectorised_verts(result_verts, verts);
    } else {
        std::clog << "qem skipped because you asked for it." << std::endl;
    }

}

}  // namespace
