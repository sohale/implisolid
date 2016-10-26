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

#include <cstddef>   // for std::nullptr only

#include "boost/multi_array.hpp"
#include "boost/array.hpp"

// #include "../js_iteration_2/basic_data_structures.hpp"

#include "../js_iteration_2/vectorised_algorithms/make_random_pm1.hpp"
#include "../js_iteration_2/vectorised_algorithms/add_inplace.hpp"
#include "../js_iteration_2/vectorised_algorithms/cross_product.hpp"
#include "../js_iteration_2/vectorised_algorithms/normalise_inplace.hpp"
#include "../js_iteration_2/vectorised_algorithms/assert_are_normalised.hpp"
using mp5_implicit::vectorised_algorithms::assert_are_normalised;
using mp5_implicit::vectorised_algorithms::first_not_normalised;


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

#include "../js_iteration_2/v2v_f2f.hpp"

#include "../js_iteration_2/subdivision/do_subdivision.hpp"


using namespace std;

namespace mp5_implicit {


REAL compute_average_edge_length(const vectorized_faces& faces, const vectorized_vect& verts) {
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


/*void compute_centroids(const vectorized_faces& faces, const vectorized_vect& verts, vectorized_vect& centroids) {
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
vectorized_vect make_random_pm1(vindex_t n, int dims, REAL amplitude) {
    //vectorized_vect result {n, dims};
    vectorized_vect result {boost::extents[n][dims]};
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
    const vectorized_vect & last_dxc,
    //outputs
    vectorized_vect & dxc_output, std::string & name_output, std::vector<REAL> & alpha_list_output
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
        vectorized_vect perturb = vectorised_algorithms::make_random_pm1(count, 3, R);
        vectorised_algorithms::add_inplace(perturb, directions_basedon_gradient);
        vectorized_vect z = vectorized_vect(vector_shape);
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

        // bug! fixed (but incorrectly)
        //?????
        //dxc_output = facet_normals_directions; // copy. checks size.
        // bug2 fixed:
        dxc_output = z; // copy. checks size.

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
        vectorized_vect mesh_normals_modifiable = facet_normals_directions;
        replace_zero_normals_with_gaussian_random(mesh_normals_modifiable);
        vectorized_vect z2 = vectorized_vect(vector_shape);
        assert(z2.shape()[0] == mesh_normals_modifiable.shape()[0]);
        assert(z2.shape()[1] == mesh_normals_modifiable.shape()[1]);
        //z = ???????????????????
        //set_vector_from(z, dxc_output); // last dxc_output
        //vectorized_vect z = dxc_output;  // copy
        const vectorized_vect& z = last_dxc; // last dxc_output
        // print_vector("mesh_normals_modifiable", mesh_normals_modifiable, 100);
        // print_vector("z", z, 10);
        vectorised_algorithms::cross_product(mesh_normals_modifiable, z, z2);
        // print_vector("z2", z2, 10);


        vectorised_algorithms::normalise_inplace(z2, mp5_implicit::CONFIG_C::center_projection::min_gradient_len);
        // print_vector("z2 after normalisation", z2, 10);

        // FAILS
        #if ASSERT_USED
            if (!assert_are_normalised(z2)) {
                vindex_t  i = first_not_normalised(z2);
                clog << " assersion failed : "
                    << "z2: " << z2[i] << "  = "
                    << "z: " << z[i] << "  (x) "
                    << "mesh_normals_modifiable: " << mesh_normals_modifiable[i] << "  "
                    << std::endl;
            }
        #endif
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
      const vectorized_vect & centroids,
      const REAL average_edge,
      //nones_map
      const vectorized_vect & facet_normals_directions,
      //vectorized_bool& treated,
      vectorized_vect & centroids_output) {

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

    const vectorized_vect& X = centroids;

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
    vectorized_vect  g_a(g_a_shape);
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
        for (int i = 0; i < best_result_x.shape()[0]; ++i) {
            best_result_x[i][0] = MAGIC_VALUE_FOR_DEBUG;
            best_result_x[i][1] = MAGIC_VALUE_FOR_DEBUG;
            best_result_x[i][2] = MAGIC_VALUE_FOR_DEBUG;
        }
    #else
        for (int i = 0; i < best_result_x.shape()[0]; ++i) {
            best_result_x[i][0] = -10.0;
            best_result_x[i][1] = -10.0;
            best_result_x[i][2] = -10.0;
        }
    #endif

    // todo: reverse this. active_indices =  still_nonsuccess_indices
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
    int still_nonsuccess_indices_effective_size = active_indices.shape()[0];
    still_nonsuccess_indices = active_indices;  // copy
    assert(std::begin(still_nonsuccess_indices) != std::begin(active_indices));  // assure it's a proper copy

    /*
    // (vectorized_bool::value_type)(b_false)
    //set_array_to_value<vectorized_bool::value_type, b_false>(already_success_bool);
    //set_array_to_value<vectorized_bool::value_type, b_false>(success);
    */

    vectorized_bool  already_success_bool(scalar_shape);
    set_all_array_elements_to_a_boolean_value(already_success_bool, b_false);

    vectorized_bool  success(scalar_shape);
    set_all_array_elements_to_a_boolean_value(success, b_false);

    //already_success_bool must be all false here.

    assert(USE_MESH_NORMALS);
    assert(assert_are_normalised(facet_normals_directions));


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
        // Uses mesh normals, cross products, random vectors, unit axes.
        cases = std::vector<int>{0, 1, 2, 3, 4, 5, 6};
    } else {
        // Gradients only
        cases = std::vector<int>{0};
    }

    vectorized_vect dxc = vectorized_vect(vector_shape);

    // todo: move alpha_list here
    for (auto iter_type : cases) {

        std::vector<REAL> alpha_list1;
        // always copy
        //vectorized_vect dxc = vectorized_vect(vector_shape);
        std::string iter_type_name;

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
            dxc,
            // outputs
            dxc, iter_type_name, alpha_list1
        );
        #if DEBUG_VERBOSE
        print_out_array_of_vectors(dxc, dxc.shape()[0], "dx");
        print_alpha_list(alpha_list1, "Alphas");
        #endif
        assert(alpha_list1.size() > 0);

        // ********************
        //  ALPHA LOOP
        // ********************
        // input: dxc = search direction at centroid points

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
            REAL c = (length_factor * alpha) * 4.0 / 4.0;
            assert(X.shape()[0] == dxc.shape()[0]);


            super_quiet = true;
            if (!super_quiet)
            std::clog << "alpha: " << alpha << " c:" << c << " length_factor:" << length_factor << std::endl;

            // todo: ONLY ACTIVE INDICES

            // for (int j=0; j < n; j++) {
            for (int j=0; j < X.shape()[0]; j++) {
                // But we need this only for active points
                //No we already have all of dxc, so no need to prune.
                // todo: only on chosen ones. ( fixme: )  set_h1_half --> h1_half_of_selected_points
                x1_half[j][0] = X[j][0] + c * dxc[j][0];
                x1_half[j][1] = X[j][1] + c * dxc[j][1];
                x1_half[j][2] = X[j][2] + c * dxc[j][2];
            }

            #if DEBUG_VERBOSE
            print_out_array_of_vectors(x1_half, x1_half.shape()[0], "x1_half");
            #endif

            // simply: active_indices = still_nonsuccess_indices
            // clog << still_nonsuccess_indices_effective_size << "==" << still_nonsuccess_indices.shape()[0] << std::endl;
            assert(still_nonsuccess_indices_effective_size == still_nonsuccess_indices.shape()[0]);
            active_indices.resize(boost::extents[still_nonsuccess_indices.shape()[0]]);

            active_indices_size = still_nonsuccess_indices_effective_size;
            active_indices = still_nonsuccess_indices;  // just copy all contents. Tobe simplified later.

            #if ASSERT_USED
            // Check if C++ actually copies!
            for (int i=0; i < active_indices.size(); i++) {
                assert(active_indices [i] == still_nonsuccess_indices[i]);
            }
            #endif

            vectorized_scalar f_a = eval_implicit_on_selected_points_indexed(object, x1_half, active_indices, active_indices.shape()[0]);

            //KEEPS repeating forever for 2 points only
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


            set_all_array_elements_to_a_boolean_value(success, b_false);  // todo: unit test.

            for (int j=0; j < active_indices.shape()[0]; j++) {
                assert(signs_a.shape()[0] == active_indices_size);
                assert(signs_a.shape()[0] == active_indices.shape()[0]);
                // unnecessary intermediate variable success0 removed
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

            array_of_indices  new_success_indices(scalar_shape);
            int new_success_indices___n_s = 0;
            //todo: refactoring
            //array_of_indices_struct  new_success_indices(scalar_shape);

            still_nonsuccess_indices_effective_size = 0;
            for (int k=0; k < active_indices.shape()[0]; k++) {
                /* **********************
                *  Major BUG DETECTED here:
                ************************ */
                auto j = active_indices[k];  // BUG!!!!

                if (success[j] && !already_success_bool[j]) {
                    new_success_indices[new_success_indices___n_s] = j;
                    new_success_indices___n_s ++;
                }
                if (!success[j] && !already_success_bool[j]) {
                    // still_nonsuccess_indices is next round's active_indices
                    still_nonsuccess_indices[still_nonsuccess_indices_effective_size] = j;
                    still_nonsuccess_indices_effective_size ++;
                }
            }

            #if ASSERT_USED
            // DEBUG THINGY
            {
                for (int j=0; j < active_indices.shape()[0]; j++) {
                    if(already_success_bool[j])
                    if (best_result_x[j][0] == MAGIC_VALUE_FOR_DEBUG || (best_result_x[j][0] == 0 && best_result_x[j][1] == 0 && best_result_x[j][2] == 0)) {
                        // happens
                        clog << "[a1]is absolute zero. j=" << j << " " << best_result_x[j][0] << std::endl;
                        abort();
                    }
                }
            }
            #endif

            // not really necessary
            still_nonsuccess_indices.resize(boost::extents[still_nonsuccess_indices_effective_size]);

            #if DEBUG_VERBOSE
            print_out_array_of_indices(still_nonsuccess_indices, still_nonsuccess_indices_effective_size, "still_nonsuccess_indices");

            print_out_array_of_indices(new_success_indices, new_success_indices___n_s, "new_success_indices");

            print_out_array_of_bools(already_success_bool, already_success_bool.shape()[0], "already_success_bool");
            print_out_array_of_bools(success, already_success_bool.shape()[0], "success");

            print_out_array_of_signs(signs_a, signs_a.shape()[0], "signs_a");
            print_out_array_of_signs(signs_c, signs_c.shape()[0], "signs_c");

            print_out_array_of_indices(active_indices, active_indices_size, "active_indices");
            #endif

            // active_indices is NOT up to date at this point.


            // note: new_success_indices[] can be empty

            for (int j=0; j < new_success_indices___n_s; j++) {
                auto k = new_success_indices[j];
                // todo: re-factor
                best_result_x[k][0] = x1_half[k][0];
                best_result_x[k][1] = x1_half[k][1];
                best_result_x[k][2] = x1_half[k][2];
            }

            // already_success_bool = already_success_bool OR success
            // Monotonically increasing (i.e. once true => forever true).
            // already_success_bool Accumulates success.
            assert(already_success_bool.shape()[0] == n);
            assert(success.shape()[0] == n);
            for (int j=0; j < n; j++) {
                if (success[j]) {  // it should be new_success? ****
                    #if ASSERT_USED
                    // DEBUG THINGY
                    {
                        // THIS DOES FAIL.
                        if (best_result_x[j][0] == MAGIC_VALUE_FOR_DEBUG || (best_result_x[j][0] == 0 && best_result_x[j][1] == 0 && best_result_x[j][2] == 0)) {
                            clog << "signs_c[j] size " << signs_c.size() << " j=" << j << std::endl;
                            clog << "success[j] size " << success.size() << " j=" << j << std::endl;
                            clog << "already_success_bool[j] size test" << already_success_bool.size() << " j=" << j << std::endl;
                            clog << "IN-POINT" <<
                                X[j][0] << "," << X[j][1] << "," << X[j][2] << " .  " <<
                                " DX=" << dxc[j][0] << "," << dxc[j][1] << "," << dxc[j][2] << " .  " <<
                                " alpha = " << alpha << "c=" << c <<
                                " already_success_bool=" << already_success_bool[j] <<
                                " success=" << success[j] <<
                                // " signs_a=" << signs_a[j] <<
                                " signs_c[j]=" << signs_c[j] <<

                                "            " << best_result_x[j][0] << std::endl;
                            clog << "[2]is absolute zero. j=" << j << " " << best_result_x[j][0] << std::endl;
                            abort();
                        }
                    }
                    #endif

                    already_success_bool[j] = b_true;

                }
            }

            // invariant:
            //     if already_success_bool[j]   <=>  best_result_x[j] is not empty

            clog << "[" << counter << "](+" << new_success_indices___n_s << ")" << still_nonsuccess_indices_effective_size << " ";

            #if ASSERT_USED
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
            #endif

            #if ASSERT_USED
                // fails
                assert(test_if_conjugate_opposite_signs_indexed(
                    object, X, best_result_x, already_success_bool, ROOT_TOLERANCE ));
                // X === centroids
            #endif


            // Code that is Different to bisection_vectorized5_()

            if (still_nonsuccess_indices_effective_size == 0) {
                break;
            }
        }  // alphas loop

        if (still_nonsuccess_indices_effective_size == 0) {
            break;
        }

    }  // iter_type loop

    // Now conjugeates are found for all points.
    // Except for still_nonsuccess_indices
    // best_result_x, centroids are conjugates.

    /////////////////////////////////////////////////////////////////////////////////////

    // todo: still_nonsuccess_indices__effective_size
    assert(still_nonsuccess_indices_effective_size == still_nonsuccess_indices.size());

    /* The points that failed to find their conjugates, just copy the original point (centroids).*/
    assign_vects_chosen_by_fancy_indexing(best_result_x, centroids, still_nonsuccess_indices, still_nonsuccess_indices_effective_size);  // A = B[C];

    // *************************************

    // vectorized_scalar f1(scalar_shape);
    const auto & f1 = fc_a;
    // xa1 <->  = x0_v3  == centroids =?= X    // todo: use X
    // todo: assert f1 == F(xa1)
    const auto & xa2 = best_result_x;

    vectorized_scalar f2(scalar_shape);
    object->eval_implicit(xa2, &f2);

    //****************************

    vectorized_bool  zeros2_bool(scalar_shape);
    // bool_find_zero_scalars(zeros2_bool, f2, ROOT_TOLERANCE);  //refactor like this
    for (int i=0; i < n; i++) {
        zeros2_bool[i] = std::abs(f2[i]) <= ROOT_TOLERANCE;
    }

    vectorized_bool  zeros1_bool(scalar_shape);
    // todo: turn this into a canonical function
    // bool_find_zero_scalars(zeros1_bool, f1, ROOT_TOLERANCE);
    for (int i=0; i < n; i++) {
        zeros1_bool[i] = std::abs(f1[i]) <= ROOT_TOLERANCE;
    }

    // todo: turn this into a canonical function
    // best_result_x = x0_v3[zeros1_bool]
    set_a_b_if_c(best_result_x, centroids, zeros1_bool);


    vectorized_bool  zeros1or2(scalar_shape);
    // todo: turn this into a canonical function
    for (int i=0; i < n; i++) {
        zeros1or2[i] = zeros2_bool[i] || zeros1_bool[i];
    }

    array_of_indices  relevants_bool_indices(scalar_shape);
    // bug fixed. Initialisation was missing.
    int relevants_bool_indices_size = 0;
    for (int i=0; i < n; i++) {
        if (already_success_bool[i] && !zeros1or2[i]) {   // another bug detected here. "!" was missing
            relevants_bool_indices[relevants_bool_indices_size] = i;
            relevants_bool_indices_size ++;
        }
    }

    relevants_bool_indices.resize(boost::extents[relevants_bool_indices_size]);
    int m = relevants_bool_indices_size; //relevants_bool_indices.shape()[0];
    assert(m == relevants_bool_indices_size);

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



    // Checks if all x1 & x2 points have srinctly opposite signs (i.e. non zero)
    #if ASSERT_USED
    {
        // no need to re-evaluate f1
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
        // fails!: assert(x1_relevant.shape() == x2_relevant.shape());   # C++ syntax
        assert(x1_relevant.shape()[0] == x2_relevant.shape()[0]);
        assert(x1_relevant.shape()[1] == x2_relevant.shape()[1]);
        assert(f1_relevants.size() == f2_relevants.size());

        // assert np.all(f1_relevants*f2_relevants <= +THRESHOLD_zero_interval)

        //  example failure: f2 = -0.15016 * f1 = -0.00199347        tol=0.001
        for (int i=0; i < m; i++) {
            REAL mult = f2_relevants[i] * f1_relevants[i];
            if (!  (mult <= - ROOT_TOLERANCE*ROOT_TOLERANCE)) {
                clog << mult <<" = " << f2_relevants[i] << " * " << f1_relevants[i] << " tol=" << ROOT_TOLERANCE << " i:[" << i << "]"<< endl;
            }
            assert(mult <= + ROOT_TOLERANCE*ROOT_TOLERANCE);
            // todo(sohail): fix the python code:
            // In Python code it is like this:
            // assert(mult <= + ROOT_TOLERANCE);
        }
    #endif

    // Swap: Conjugate points may be in the reverse order.
    // They must be (x1=outside, x2=inside)
    int ctr = 0;
    for (int i=0; i < m; i++) {
        // If x2 is inside, swap it. => x1 has to be outside.
        if (f2_relevants[i] < -ROOT_TOLERANCE) {  // outside
            std::swap(x2_relevant[i][0], x1_relevant[i][0]);
            std::swap(x2_relevant[i][1], x1_relevant[i][1]);
            std::swap(x2_relevant[i][2], x1_relevant[i][2]);
            // assert(x2_relevant[i][0] != x1_relevant[i][0]);
            // should we assert that they are not exactly equal? There used to be a bug that sometimes they were exactly equal
            ctr++;
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

    // now put them (all of those RELEVEANT ones) into the results array.
    // What about those which were not successful? Answer: They are not "relevant".
    // not Relevant: either the conjugate was not found, or t was already zero (zeros1or2).


    // Shold not be needed if it is not resized:
    assert(relevants_bool_indices.size() == x_bisect.shape()[0]);

    bool same_address1 = centroids_output.begin() != centroids.begin();


    // Copying centroids into centroids_output is necessary, so that when neigher zeros1or2[.] nor relevants_bool_indices[.],
    // ,  ... wait, what?  best_result_x[] is ignored. But it's fine. That was before the bisecion.
    // relevant[] : pased through bisection, hance, correct.
    // zeros1or2[]: not passed through bisection, but correct.
    // rest: either non-successful conjugate or non-successful bisection? But there will be no unsuccessful bisection. It always converges.
    // todo: fix the python code to sepaate the vairables.

    centroids_output = centroids;  // copy
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
    add_noise_to_vectors(centroids_output, 0.1);

    void add_noise_to_vectors(auto & centroids_output, const REAL noise_amplitude) {
    for (int i=0; i < centroids_output.shape()[0]; i++) {
        centroids_output[i][0] += rand01() * noise_amplitude;
        centroids_output[i][1] += rand01() * noise_amplitude;
        centroids_output[i][2] += rand01() * noise_amplitude;
    }
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

// #include "v2v_f2f.hpp"

// rename: project_points_on_surface. (output_points)
void centroids_projection(mp5_implicit::implicit_function* object, std::vector<REAL>& result_verts, const std::vector<vertexindex_type>& mesh_faces, bool enable_qem) {

    // clog << "mesh_faces.size():" << mesh_faces.size() << std::endl;

    vectorized_vect  verts = vects2vects(result_verts);
    vectorized_faces  faces = copy_faces_from_vectorfaces(mesh_faces);

    REAL average_edge;
    average_edge = compute_average_edge_length(faces,  verts);

    vectorized_vect::index  num_faces = static_cast<vectorized_vect::index>(mesh_faces.size()/3);
    vectorized_vect_shape  centroids_shape = { num_faces , 3 };
    vectorized_vect centroids(centroids_shape);

    //if (VERBOSE_QEM)
    //    clog << "1.centroids: " << centroids.shape()[0] << "x" << centroids.shape()[1] << " : "  << centroids[0][0] << std::endl;

    compute_centroids(faces, verts, centroids);

    STORE_POINTSET("pre_p_centroids", centroids);

    // bug resolved: treated was not initialised
    // bug revolved: mesh_normals (facet_normals) was not initialised.

    vectorized_vect facet_normals = produce_facet_normals(faces, verts, true);
    assert(assert_are_normalised(facet_normals));

    mp5_implicit::set_centers_on_surface(object, centroids, average_edge, facet_normals, centroids);

    STORE_POINTSET("post_p_centroids", centroids);

    /*
    clog << "post_p_centroids" << " :" << centroids.shape()[0] << std::endl;
    clog << "faces" << " :" << faces.shape()[0] << std::endl;
    clog << "mesh_faces" << " :" << mesh_faces.size() << std::endl;
    */


    std::vector< std::vector<faceindex_type>> face_neighbours_of_vertex_list;
    face_neighbours_of_vertex_list = make_neighbour_faces_of_vertex(faces, verts.shape()[0]);

    vectorized_vect  centroid_gradients(centroids_shape);

    compute_centroid_gradient(centroids, centroid_gradients, object);

    // fine_tune_mesh(verts);

    /* ***********************
    QEM
    ************************ */

    STORE_POINTSET("pre_qem_verts", verts);

    if (enable_qem) {
        std::clog << "Going for QEM:" << std::endl;
        // array_of_indices  ranks {num_faces}; // will not work
        array_of_indices  ranks {boost::extents[num_faces]};
        REAL maximum_qem_displacement = average_edge;

        vertex_apply_qem(&verts, faces, centroids, face_neighbours_of_vertex_list, centroid_gradients, &ranks, maximum_qem_displacement);

        if (STORE_POINTSETS)
        {

        STORE_POINTSET("post_qem_verts", verts);

        // To visualise the ranks, the qem verts pointset is removed except for rank==ONLY_RANK.
        int ONLY_RANK = -1;
        for (int i=0;i < verts.shape()[0]; i++) {
            if (ONLY_RANK >= 0)
            if (ranks[i] != ONLY_RANK) {
                for (int d=0;d<3;d++) {
                    point_set_set["post_qem_verts"][i][d] = -1000;
                    point_set_set["pre_qem_verts"][i][d] = -1000;
                }
            }
        }

        }

        // if APPLY_QEM_RESULT
        set_vectorverts_from_vectorised_verts(result_verts, verts);

    } else {
        std::clog << "QEM skipped because you asked for it." << std::endl;
        // todo(Sohail): Check whether we always update the verts (and faces) even if QEM is desabled.

        // should we??
        // set_vectorverts_from_vectorised_verts(result_verts, verts);
    }
}


void my_subdiv_(
        std::vector<REAL>& given_verts,
        std::vector<vertexindex_type>& given_faces
    ) {
    //todo: vector_to_multi_array(given_verts);
    //todo: vector_to_multi_array(given_faces);
    vectorized_vect  verts = vects2vects(given_verts);
    vectorized_faces  faces = copy_faces_from_vectorfaces(given_faces);

    cout << "given_verts.size" << given_verts.size()/3 << std::endl;
    cout << "given_faces.size" << given_faces.size()/3 << std::endl;


    // do_subdivision()
    std::set<faceindex_type>  which_facets_set;
    // which_facets_set.insert(0);
    for (int i = 0; i < faces.shape()[0]; ++i) {
        which_facets_set.insert(i);
    }

    auto fv2 = ::mp5_implicit::subdivision::subdivide_given_faces (faces, verts, which_facets_set);


    cout << "done" << std::endl;
    auto& f = std::get<0>(fv2);
    cout << "f.size" << f.shape()[0] << std::endl;

    faces.resize(boost::extents[f.shape()[0]][f.shape()[1]]);
    faces = f;

    // TODO(Sohail): Too much copying (moving data).
    auto& v = std::get<1>(fv2);
    cout << "v.size" << v.shape()[0] << std::endl;
    verts.resize(boost::extents[v.shape()[0]][v.shape()[1]]);
    verts = v;

    cout << "now randomizing" << std::endl;

    randomize_verts(verts, 0.3);

    cout << "convertings back" << std::endl;

    /* todo(Sohail): unify the types. vector<> and multi_array. remove the former, or only-once in this function.
     * Why is it called mesh_faces?
    */
    set_vectorverts_from_vectorised_verts(given_verts, verts);
    // mesh_faces = copy_faces_from_vectorfaces(faces);
    set_vectorfaces_from_vectorised_faces(given_faces, faces);

    cout << "V:" << given_verts.size()/3 << " , F:" << given_faces.size()/3 << std::endl;
    cout << "subdiv done" << std::endl;

}

}  // namespace


