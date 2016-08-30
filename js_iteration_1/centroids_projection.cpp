/*
Copyright 2016 MyMiniFactory Ltd.
*/

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

using namespace std;
using namespace mp5_implicit;

#define ASSERTS 1

REAL compute_average_edge_length(const faces_t& faces, const verts_t& verts) {
    int nfaces = faces.shape()[0];
    REAL edge_length;
    for (int j=0; j < nfaces; j++) {
        edge_length += norm_2(verts[faces[j][0]][0] - verts[faces[j][1]][0], verts[faces[j][0]][1] - verts[faces[j][1]][1], verts[faces[j][0]][2] - verts[faces[j][1]][2]);
        edge_length += norm_2(verts[faces[j][0]][0] - verts[faces[j][2]][0], verts[faces[j][0]][1] - verts[faces[j][2]][1], verts[faces[j][0]][2] - verts[faces[j][2]][2]);
        edge_length += norm_2(verts[faces[j][2]][0] - verts[faces[j][1]][0], verts[faces[j][2]][1] - verts[faces[j][1]][1], verts[faces[j][2]][2] - verts[faces[j][1]][2]);
    }
    return edge_length/(3.*nfaces);
}


void compute_centroids(const faces_t& faces, const verts_t& verts, verts_t& centroids) {
    int nt = faces.shape()[0];
    for (int j=0; j < nt; j++) {
        int f0 = faces[j][0];
        int f1 = faces[j][1];
        int f2 = faces[j][2];
        for (int di=0; di < 3; di++) {
                centroids[j][di] = (verts[f0][di] + verts[f1][di] + verts[f2][di])/3.;

        }
    }
}

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

template<typename T>
void here(T arg) {
    std::cout << arg << std::endl << std::flush;
}

inline bool test_points_sign(verts_t& x_vectorized, const mp5_implicit::implicit_function& object, REAL ROOT_TOLERANCE, REAL sign) {
    assert(sign == +1 || sign == -1); // || sign == 0.0);

    cout << "test_points_sign" << endl;

    int n = x_vectorized.shape()[0];
    // int n = v_arr.size(); //[0];
    boost::array<int, 1> v1_shape = {n};
    vectorized_scalar v_arr(v1_shape);  // x_vectorized.shape());

    object.eval_implicit(x_vectorized, &v_arr);

    // int n = v_arr.size(); //[0];
    bool everything_alright = true;
    for (int i=0; i < n; i++) {
        auto s1 = my_sign(v_arr[i], ROOT_TOLERANCE);
        bool ok = (s1 * sign >  0 + 0.0001);
        if (!ok) {
            cout <<  v_arr[i] << " " << s1 << " " << sign << " :" << i << endl;
        }
        everything_alright = everything_alright && ok;
    }
    cout << "bool" << everything_alright << endl;
    return everything_alright;
}

inline bool test_if_points_are_inside(verts_t& x2_vectorized, const mp5_implicit::implicit_function& object, REAL ROOT_TOLERANCE) {
    return test_points_sign(x2_vectorized, object, ROOT_TOLERANCE, +1);
}
inline bool test_if_points_are_outside(verts_t& x1_vectorized, const mp5_implicit::implicit_function& object, REAL ROOT_TOLERANCE) {
    return test_points_sign(x1_vectorized, object, ROOT_TOLERANCE, -1);
}


inline void check_bisection_input_signs(
    vectorized_scalar& v1_arr,
    vectorized_scalar& v2_arr,
    REAL ROOT_TOLERANCE,
    int n,
    int active_count,
    int iteration
  ) {
    #if ASSERT_USED
        here("a1. iteration"+std::to_string(iteration));

        bool assert1 = true;
        for (int i=0; i < active_count; i++) {
            auto s1 = my_sign(v1_arr[i], ROOT_TOLERANCE);
            auto s2 = my_sign(v2_arr[i], ROOT_TOLERANCE);
            // assert(s1* s2 < 0 - EPS);
            bool ok = (s1* s2 < 0 - EPS);
            // if (!ok) {
            //     std::cout << "["<<i<<"]"<< s1 << " " << s2 << " v1:" << v1_arr[i] << " v2:" << v2_arr[i] << endl;
            // }
            assert1 = assert1 && ok;
        }
        if (!assert1) cout << "";
        assert(assert1);

        here("a2");

        bool assert2 = true;
        for (int i=0; i < active_count; i++) {
            bool ok = v1_arr[i] < 0 - ROOT_TOLERANCE;
            assert2 = assert2 && ok;
        }
        here("a3");

        assert(assert2);
        here("a4");

        assert(active_count <= n);
        here("a5");

        assert(true);
    #endif
}

/*
    vectorized bisection.
    x1_arr: points outside
    x2_arr: points inside

    All x1 must be strictly outside, and x2 must be inside the solid.
*/
void bisection(
    mp5_implicit::implicit_function* object,
    verts_t& res_x_arr,
    verts_t& x1_arr,
    verts_t& x2_arr,
    REAL ROOT_TOLERANCE,
    boost::multi_array<bool, 1>& treated) {

    // initilization step
    int n = x1_arr.shape()[0];
    assert(x2_arr.shape()[0] == n);
    assert(x1_arr.shape()[1] == 3);
    assert(x2_arr.shape()[1] == 3);

    assert(res_x_arr.shape()[0] == n);
    assert(res_x_arr.shape()[1] == 3);

    // implicit function of the two arrays
    boost::array<int, 1> v1_shape = {n};

    #if ASSERT_USED
        vectorized_scalar v1_arr(v1_shape);   // todo: rename to v1
        vectorized_scalar v2_arr(v1_shape);

        object->eval_implicit(x1_arr, &v1_arr);
        object->eval_implicit(x2_arr, &v2_arr);


        for (int i=0; i < n; i++) {
            res_x_arr[i][0] = -10000.0;
            res_x_arr[i][1] = -10000.0;
            res_x_arr[i][2] = -10000.0;
        }

    #endif

    boost::multi_array<int, 1> active_indices(v1_shape);

    for (int i=0; i < n; i++) {
        active_indices[i] = i;
    }
    int active_indices_size = n;

    for (int i=0; i < n; i++) {
        treated[i] = false;
    }

    int active_count = n;
    int solved_count = 0;

    boost::array<int, 2> x_mid_shape = {n, 3};
    boost::multi_array<REAL, 2> x_mid(x_mid_shape);
    /*
    #if ASSERT_USED
        for (int i = 0; i < n; i++) {
            for (int d = 0; d < 3; d++) {
                x_mid[i][d] = 1.0;
            }
        }
    #endif
    */

    vectorized_scalar v_mid(v1_shape); // implicit function for x_mid
    /*
    #if ASSERT_USED
        for (int i = 0; i < n; i++) {
            v_mid[i] = 0.0;
        }
    #endif
    */

    vectorized_scalar abs_(v1_shape); // absolute value of the implicit function


    // array of indices
    boost::multi_array<vindex_t, 1> indices_boundary(v1_shape);
    boost::multi_array<vindex_t, 1> indices_outside(v1_shape);
    boost::multi_array<vindex_t, 1> indices_inside(v1_shape);
    boost::multi_array<vindex_t, 1> indices_eitherside(v1_shape);
    boost::multi_array<vindex_t, 1> which_zeroed(v1_shape);

    int iteration = 1;

    // loop
    while (true) {
        here("start");

        /* Checks if all v1 and v2 have opposite signs and are not zero*/
        #if ASSERT_USED

        check_bisection_input_signs(v1_arr, v2_arr, ROOT_TOLERANCE, n, active_count, iteration);
        #endif

        here("2");

        // mean of (x1, x2)
        for (int i=0; i < active_count; i++) {
            x_mid[i][0] = (x1_arr[i][0] + x2_arr[i][0]) / 2.;
            x_mid[i][1] = (x1_arr[i][1] + x2_arr[i][1]) / 2.;
            x_mid[i][2] = (x1_arr[i][2] + x2_arr[i][2]) / 2.;
        }

        here("3");

        // *************************************
        // Fix me
        object->eval_implicit(x_mid, &v_mid); // no. only first ones**
        //
        // *************************************

        here("4");

        assert(active_indices_size == active_count);


        here("5");

        for (int i=0; i < active_count; i++) {
            abs_[i] = ABS(v_mid[i]);
        }
        // int abs_size = active_count;

        here("6");

        // imcrementing the size of the indices arrays
        int indices_boundary_size = 0;
        for (int i=0; i < active_count; i++) {
            if (abs_[i] <= ROOT_TOLERANCE) {
                indices_boundary[indices_boundary_size] = i;
                indices_boundary_size ++;
            }
        }

        here("7");

        int i_e = 0;
        for (int i=0; i < active_count; i++) {
            if (abs_[i] > ROOT_TOLERANCE) {
                indices_eitherside[i_e] = i;
                i_e ++;
            }
        }
        int indices_eitherside_size = i_e;

        int i_o = 0;
        for (int i=0; i < active_count; i++) {
            if (v_mid[i] < -ROOT_TOLERANCE) {
                indices_outside[i_o] = i;
                i_o ++;
            }
        }
        int indices_outside_size = i_o;

        int i_i = 0;
        for (int i=0; i < active_count; i++) {
            if (v_mid[i] > +ROOT_TOLERANCE) {
                indices_inside[i_i] = i;
                i_i ++;
            }
        }
        int indices_inside_size = i_i;

        here("8");

        assert(indices_boundary_size + indices_inside_size + indices_outside_size == active_count);

        here("9");

        assert(indices_eitherside_size + indices_boundary_size == active_count);

        here("10");

////////////////
        // which_zeroed = active_indices[indices_boundary]
        for (int i = 0; i < indices_boundary_size; ++i) {
            which_zeroed[i] = active_indices[indices_boundary[i]];
            treated[active_indices[indices_boundary[i]]]=true;
        }
        const int which_zeroed_size = indices_boundary_size;

        int found_count = indices_boundary_size;
        solved_count += found_count;

        assert(active_count - found_count + solved_count == n);

        // copy into the result, the x_mid that solved the equation.
        for (int i = 0; i < indices_boundary_size; ++i) {
            int w = which_zeroed[i];
            int b = indices_boundary[i];
            res_x_arr[w][0] = x_mid[b][0];
            res_x_arr[w][1] = x_mid[b][1];
            res_x_arr[w][2] = x_mid[b][2];
        }

        // assert indices_boundary[:] < active_count
        #if ASSERT_USED
        {
            bool assertok = true;
            for (int i = 0; i < indices_boundary_size; ++i) {
                bool ok = indices_boundary[i] < active_count;
                assertok = assertok && ok;
            }
            assert(assertok);
        }
        #endif
        //**********************************************
        //* Code review of "bisection" done up to here *
        //**********************************************

        // changing the values of x2 and x1
        for (int i=0; i < i_i; i++) {
            #if ASSERT_USED
                    v2_arr[indices_inside[i]] = v_mid[indices_inside[i]];
            #endif
            x2_arr[indices_inside[i]][0] = x_mid[indices_inside[i]][0];
            x2_arr[indices_inside[i]][1] = x_mid[indices_inside[i]][1];
            x2_arr[indices_inside[i]][2] = x_mid[indices_inside[i]][2];
        }

        for (int i=0; i < i_o; i++) {
            #if ASSERT_USED
                v1_arr[indices_outside[i]] = v_mid[indices_outside[i]];
            #endif
            x1_arr[indices_outside[i]][0] = x_mid[indices_outside[i]][0];
            x1_arr[indices_outside[i]][1] = x_mid[indices_outside[i]][1];
            x1_arr[indices_outside[i]][2] = x_mid[indices_outside[i]][2];
        }

        // next round

        for (int i=0; i < i_e; i++) {
            active_indices[i] = active_indices[indices_eitherside[i]];
        }

        indices_boundary.resize(boost::extents[i_e]);
        indices_eitherside.resize(boost::extents[i_e]);
        indices_outside.resize(boost::extents[i_e]);
        indices_inside.resize(boost::extents[i_e]);
        which_zeroed.resize(boost::extents[i_e]);
        active_indices.resize(boost::extents[i_e]);

        active_count = active_count - found_count;

        iteration += 1;

        for (int i=0; i < active_count; i++) {
            #if ASSERT_USED
                v1_arr[i] = v1_arr[indices_eitherside[i]];
                v2_arr[i] = v2_arr[indices_eitherside[i]];
            #endif

            x1_arr[i][0] = x1_arr[indices_eitherside[i]][0];
            x1_arr[i][1] = x1_arr[indices_eitherside[i]][1];
            x1_arr[i][2] = x1_arr[indices_eitherside[i]][2];

            x2_arr[i][0] = x2_arr[indices_eitherside[i]][0];
            x2_arr[i][1] = x2_arr[indices_eitherside[i]][1];
            x2_arr[i][2] = x2_arr[indices_eitherside[i]][2];
        }

        if (active_indices.shape()[0] == 0 || iteration==10) {
            cout << "projection treated this much points" << endl;
            cout << solved_count << endl;
            break;
        }

    }

}


std::vector<REAL> make_alpha_list(REAL initial_step_size, REAL unit_of_length, REAL min_step_size, REAL max_dist, int max_iter, bool EXTREME_ALPHA) {
    // Prepare a list of step sizes.
    // vectorized_scalar alpha_list(scalar_shape);
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

// main function
void  set_centers_on_surface(
      mp5_implicit::implicit_function* object,
      verts_t& centroids,
      const REAL average_edge,
      boost::multi_array<bool, 1>& treated) {

    // intilization, objects creation

    constexpr REAL min_gradient_len = 0.000001;  // Gradients smaller than this are considered zero.
    constexpr int max_iter = 20;
    constexpr bool USE_MESH_NORMALS = true;
    constexpr bool EXTREME_ALPHA = false;  // keep false

    const verts_t& x = centroids;

    REAL max_dist = average_edge;

    int n = centroids.shape()[0];
    boost::array<int, 1> scalar_shape = {n};
    // auto scalar_shape = boost::extents[n];

    vectorized_scalar fc_a(scalar_shape);
    object->eval_implicit(centroids, &fc_a);

    boost::array<int, 2> g_a_shape = { n , 3 };
    boost::multi_array<REAL, 2> g_a(g_a_shape);
    // applying the centroid gradient

    /*compute_centroid_gradient(x, g_a, object);
    Compute the gradients, and normalise them.
    */
    // now g_a --> g_direction
    {
    constexpr REAL MIN_NORM = 1.0;  // An absolute value: not good.
    object->eval_gradient(x, &g_a);
    for (int i = 0, e = g_a.shape()[0]; i < e; i++) {
        REAL norm = norm_2(g_a[i][0], g_a[i][1], g_a[i][2]);
        if (norm < min_gradient_len) {
            norm = MIN_NORM;
        }
        REAL factor = 1.0 / norm;
        g_a[i][0]=g_a[i][0] / norm;
        g_a[i][1]=g_a[i][1] / norm;
        g_a[i][2]=g_a[i][2] / norm;
        #if ASSERT_USED
            for (int j = 0; j < 3; j++) {
                assert(g_a[i][j] <= 1.);
                assert(g_a[i][j] >= -1.);
            }
        #endif
    }
    }

    // Elements will be modified
    auto& g_direction_a = g_a;
    // g_a = delete;
    // delete g_a ;
    // delte g_a
    // int g_a;
    // void qq;

    // todo: remove: negative_f_c


    vectorized_scalar signs_c(scalar_shape);
    for (int i=0; i < fc_a.shape()[0]; i++) {
        if (fc_a[i] > ROOT_TOLERANCE) {
            signs_c[i] = +1.0;
        } else if (fc_a[i] < -ROOT_TOLERANCE) {
            signs_c[i] = -1.0;
        } else {
            signs_c[i] = 0.0;
        }
    }

    assert(g_direction_a.shape()[0] == fc_a.shape()[0]);

    for (int i = 0, e = g_direction_a.shape()[0]; i < e; i++) {
        if (signs_c[i] < 0.0) {
            // i.e., if fc_a[i] < -ROOT_TOLERANCE
            g_direction_a[i][0] = - g_direction_a[i][0];
            g_direction_a[i][1] = - g_direction_a[i][1];
            g_direction_a[i][2] = - g_direction_a[i][2];
        }
    }


    // Move toward surface: If inside (the f value is positive) move outwards (oppoosite the -gradient), and if outside, move inwards (+gradient).
    boost::array<int, 2> vector_shape = {n, 3};
    boost::multi_array<REAL, 2> dx0_c_grad(vector_shape);

    for (int i=0; i < fc_a.shape()[0]; i++) {
        // Bug fixed: negative sign missing.
        dx0_c_grad[i][0] = - g_direction_a[i][0]*signs_c[i];
        dx0_c_grad[i][1] = - g_direction_a[i][1]*signs_c[i];
        dx0_c_grad[i][2] = - g_direction_a[i][2]*signs_c[i];
    }

    // stepsize happens to be equal to the max_dist.
    REAL initial_step_size = max_dist * 1.0;
    std::vector<REAL> alpha_list = make_alpha_list(initial_step_size, average_edge, 0.001, max_dist, max_iter, EXTREME_ALPHA);

    // THE algorithm

    // array definition
    boost::multi_array<REAL, 2> best_result_x(vector_shape);
    boost::multi_array<REAL, 2> x1_half(vector_shape);
    boost::multi_array<REAL, 2> xa4(vector_shape);

    vectorized_scalar f_a(scalar_shape);
    vectorized_scalar signs_a(scalar_shape);

    // boolean
    boost::multi_array<bool_t, 1> success0(scalar_shape);
    boost::multi_array<bool_t, 1> success(scalar_shape);

    // indices arrays
    boost::multi_array<int, 1> active_indices(scalar_shape);
    boost::multi_array<int, 1> still_nonsuccess_indices(scalar_shape);


    for (int i=0; i < n; i++) {
        active_indices[i] = i;
    }

    int active_count = n;

    still_nonsuccess_indices = active_indices;

    boost::multi_array<bool_t, 1> already_success(scalar_shape);
    boost::multi_array<bool_t, 1> new_success_indices(scalar_shape);


    for (int i=0; i < n; i++) {
        already_success[i] = b_false;
    }

    int counter = -1;
    int s_n_s = 0;
    // main part of the algor
    for (int i=0; i < alpha_list.size(); i++) {
        counter += 1;

        for (int j=0; j < n; j++) {
            x1_half[j][0] = centroids[j][0] + (max_dist*alpha_list[i])*dx0_c_grad[j][0];
            x1_half[j][1] = centroids[j][1] + (max_dist*alpha_list[i])*dx0_c_grad[j][1];
            x1_half[j][2] = centroids[j][2] + (max_dist*alpha_list[i])*dx0_c_grad[j][2];
        }

        active_indices = still_nonsuccess_indices;

        for (int j=0; j < active_indices.shape()[0]; j++) {
            xa4[j][0] = x1_half[active_indices[j]][0];
            xa4[j][1] = x1_half[active_indices[j]][1];
            xa4[j][2] = x1_half[active_indices[j]][2];
        }
        object->eval_implicit(xa4, &f_a);

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


        int n_s = 0;
        s_n_s = 0;
        for (int j=0; j < active_indices.shape()[0]; j++) {
            if (success[j] == b_true && already_success[j] == b_false) {
                new_success_indices[n_s] = j;
                n_s ++;
            } else {
                still_nonsuccess_indices[s_n_s] = j;
                s_n_s ++;
            }
        }

        for (int j=0; j < n_s; j++) {
            best_result_x[new_success_indices[j]][0] = x1_half[new_success_indices[j]][0];
            best_result_x[new_success_indices[j]][1] = x1_half[new_success_indices[j]][1];
            best_result_x[new_success_indices[j]][2] = x1_half[new_success_indices[j]][2];
        }

        for (int j=0; j < n; j++) {
            if (success[j] == b_true) {
                already_success[j] = b_true;
            }
        }
        /* Code that is Different to bisection_vectorized5_()

        */

        /*
        // MOVED BELOW

        if (s_n_s == 0) {
            break;
        }
        */

        active_indices.resize(boost::extents[s_n_s]);
        still_nonsuccess_indices.resize(boost::extents[s_n_s]);

        if (s_n_s == 0) {
            break;
        }

    }

/////////////////////////////////////////////////////////////////////////////////////

    // todo: still_nonsuccess_indices__effective_size
    assert(s_n_s == still_nonsuccess_indices.size());

    for (int i=0; i < s_n_s; i++) {
        auto k = still_nonsuccess_indices[i];
        best_result_x[k][0] = centroids[k][0];
        best_result_x[k][1] = centroids[k][1];
        best_result_x[k][2] = centroids[k][2];
    }


    // vectorized_scalar f1(scalar_shape);
    const auto & f1 = fc_a;

    // boost::multi_array<REAL, 2> xa1(vector_shape);
    // boost::multi_array<REAL, 2> xa2(vector_shape);
    const auto & xa1 = centroids;  // = x0_v3
    const auto & xa2 = best_result_x;

    // todo: assert f1 == F(xa1)

    vectorized_scalar f2(scalar_shape);
    object->eval_implicit(xa2, &f2);

    boost::multi_array<bool_t, 1> zeros2_bool(scalar_shape);
    boost::multi_array<bool_t, 1> zeros1_bool(scalar_shape);
    boost::multi_array<bool_t, 1> zeros1or2(scalar_shape);
    boost::multi_array<int, 1> relevants_bool_indices(scalar_shape);

    cout << "3" << endl;



    //****************************
    // UP TO HERE

    for (int i=0; i < n; i++) {
        zeros2_bool[i] = ABS(f2[i])<= ROOT_TOLERANCE;
        /*

        if (ABS(f2[i])<= ROOT_TOLERANCE) {
            zeros2_bool[i] = b_true;
        } else {
            zeros2_bool[i] = b_false;
        }
        */
    }

    for (int i=0; i < n; i++) {
        zeros1_bool[i] = ABS(f1[i]) <= ROOT_TOLERANCE;
        /*
        if (ABS(f1[i])<= ROOT_TOLERANCE) {
            zeros1_bool[i] = b_true;
        } else {
            zeros1_bool[i] = b_false;
        }
        */
    }

    // best_result_x = x0_v3[zeros1_bool]
    for (int i=0; i < n; i++) {
        // if (ABS(f1[i])<= ROOT_TOLERANCE) {
        if (zeros1_bool[i]) {
            best_result_x[i][0] = centroids[i][0];
            best_result_x[i][1] = centroids[i][1];
            best_result_x[i][2] = centroids[i][2];
        } else {
        }
    }


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
        if (already_success[i] && zeros1or2[i]) {
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
            bool ok = ABS(f1[j]) <= ROOT_TOLERANCE;
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
          bool ok = ABS(f1[j]) > ROOT_TOLERANCE;
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
    boost::multi_array<REAL, 2> x1_relevant(x1_relevant_shape);
    boost::multi_array<REAL, 2> x2_relevant(x1_relevant_shape);

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
    }
    cout << "x1x2" << endl;

    // Now we have x1_relevant and x2_relevant which have the same size.

    object->eval_implicit(x2_relevant, &f2_relevants);

    // Check the signs are opposite
    #if ASSERT_USED
        vectorized_scalar f1_relevants(f1_relevant_shape);
        // object->eval_implicit(x2_relevant, &f1_relevants_);
        object->eval_implicit(x1_relevant, &f1_relevants);

        cout << x1_relevant.shape()[0] << " " << x1_relevant.shape()[1] << " " << x1_relevant.shape() << "  ,  " <<
            " " << x2_relevant.shape()[0] << " " << x2_relevant.shape()[1] << " " <<
            x2_relevant.shape()[0]<<"/"<<x2_relevant.shape()[1]<<"/"<<x2_relevant.shape()[2]<<"/"<<x2_relevant.shape()[3]<<"/"<<x2_relevant.shape()[4]<<"/"<<x2_relevant.shape()[5]<<"/"<<x2_relevant.shape()[6]<<"/"<<x2_relevant.shape()[7]
            << "  " <<
            x1_relevant.shape()[0]<<"/"<<x1_relevant.shape()[1]<<"/"<<x1_relevant.shape()[2]<<"/"<<x1_relevant.shape()[3]<<"/"<<x1_relevant.shape()[4]<<"/"<<x1_relevant.shape()[5]<<"/"<<x1_relevant.shape()[6]<<"/"<<x1_relevant.shape()[7]
            << endl;
        cout << "fff" << endl;
        // fails: assert(x1_relevant.shape() == x2_relevant.shape());
        // assert(x1_relevant.shape() == x2_relevant.shape());
        assert(x1_relevant.shape()[0] == x2_relevant.shape()[0]);
        assert(x1_relevant.shape()[1] == x2_relevant.shape()[1]);

        cout << "ggg" << endl;
        cout << f1_relevants.size() << " " << f2_relevants.size() << endl;
        assert(f1_relevants.size() == f2_relevants.size());
        cout << "hhh" << endl;

        // assert np.all(f1_relevants*f2_relevants <= +THRESHOLD_zero_interval)

        for (int i=0; i < m; i++) {
            REAL mult = f2_relevants[i] * f1_relevants[i];
            if (!  (mult <= - ROOT_TOLERANCE*ROOT_TOLERANCE)) {
                cout << mult <<" = " << f2_relevants[i] << " * " << f1_relevants[i] << " tol=" << ROOT_TOLERANCE << "[" << i << "]"<< endl;
            }
            // if (0)
            assert(mult <= - ROOT_TOLERANCE*ROOT_TOLERANCE);
        }
    #endif

    // Swap
    /*
    REAL temp0;
    REAL temp1;
    REAL temp2;
    */

    cout << "m=" << m << endl;
    int ctr = 0;
    for (int i=0; i < m; i++) {
        // If x2 is inside, swap it. => x1 has to be outside.
        if (f2_relevants[i] < -ROOT_TOLERANCE) {
            // ****************************
            // problem: sometimes they are exactly equal!!
            // if (ctr<10) cout << x2_relevant[i][0] << ", " << x1_relevant[i][0] << " <-> ";
            std::swap(x2_relevant[i][0], x1_relevant[i][0]);
            std::swap(x2_relevant[i][1], x1_relevant[i][1]);
            std::swap(x2_relevant[i][2], x1_relevant[i][2]);
            // if (ctr<10) cout << x2_relevant[i][0] << ", " << x1_relevant[i][0] << endl;

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
        assert(test_if_points_are_inside(x2_relevant, *object, ROOT_TOLERANCE));
        assert(test_if_points_are_outside(x1_relevant, *object, ROOT_TOLERANCE));
    #endif

    boost::multi_array<REAL, 2> x_bisect(x1_relevant_shape);
    // calling the vectorized bisection
    bisection(object, x_bisect, x1_relevant, x2_relevant, ROOT_TOLERANCE, treated);
    // x1_relevant: outside, x2_relevant: inside


    #if ASSERT_USED
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
            ok = ok && ABS(f1_relevants[i]) < ROOT_TOLERANCE;
        }
        assert(ok);
    #endif




    assert(relevants_bool_indices.size() == x_bisect.shape()[0]);

    // changing the values of the centroids
    assert(m == relevants_bool_indices.size());
    for (int i=0; i < m; i++) {
        centroids[relevants_bool_indices[i]][0] = x_bisect[i][0];
        centroids[relevants_bool_indices[i]][1] = x_bisect[i][1];
        centroids[relevants_bool_indices[i]][2] = x_bisect[i][2];
    }

    assert(n == zeros1or2.size());
    for (int i=0; i < n; i++) {
        if (zeros1or2[i]) {
            centroids[i][0] = best_result_x[i][0];
            centroids[i][1] = best_result_x[i][1];
            centroids[i][2] = best_result_x[i][2];
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


std::vector< std::vector<int>> make_neighbour_faces_of_vertex(const verts_t& verts, const faces_t& faces) {
    int nt = faces.shape()[0];
    int vt = verts.shape()[0];
    std::vector< std::vector<int>> neighbour_faces_of_vertex;
    for (int fi=0; fi < vt; fi++) {
        neighbour_faces_of_vertex.push_back(std::vector<int>());
    }
    for (int fi=0; fi < nt; fi++) {
        for (int vi=0; vi < 3; vi++) {
            int v1 = faces[fi][vi];
            neighbour_faces_of_vertex[v1].push_back(fi);
        }
    }

    return neighbour_faces_of_vertex;
}

// get the matrix A and b used in vertex_apply_qem
void get_A_b(const std::vector<int> nai, const verts_t& centroids, const verts_t& centroid_gradients, verts_t* A, vectorized_scalar* b) {

    int m = nai.size();

    boost::array<int, 2> center_array_shape = {m, 3};
    boost::multi_array<REAL, 2> center_array(center_array_shape);
    boost::multi_array<REAL, 2> normals(center_array_shape);

    for (int i=0; i < m; i++) {
        vindex_t cn = nai[i];
        normals[i][0] = centroid_gradients[cn][0];
        normals[i][1] = centroid_gradients[cn][1];
        normals[i][2] = centroid_gradients[cn][2];

        center_array[i][0] = centroids[cn][0];
        center_array[i][1] = centroids[cn][1];
        center_array[i][2] = centroids[cn][2];
    }

    REAL a00;
    REAL a01;
    REAL a02;
    REAL a11;
    REAL a12;
    REAL a22;

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
        a00 = normals[j][0] * normals[j][0];
        a01 = normals[j][0] * normals[j][1];
        a02 = normals[j][0] * normals[j][2];
        a11 = normals[j][1] * normals[j][1];
        a12 = normals[j][1] * normals[j][2];
        a22 = normals[j][2] * normals[j][2];

        (*A)[0][0] += a00;
        (*A)[0][1] += a01;
        (*A)[0][2] += a02;
        (*A)[1][0] += a01;
        (*A)[1][1] += a11;
        (*A)[1][2] += a12;
        (*A)[2][2] += a22;
        (*A)[2][0] += a02;
        (*A)[2][1] += a12;

        REAL cx = center_array[j][0];
        REAL cy = center_array[j][1];
        REAL cz = center_array[j][2];

        (*b)[0] -= a00 * cx + a01 * cy + a02 * cz;
        (*b)[1] -= a01 * cx + a11 * cy + a12 * cz;
        (*b)[2] -= a02 * cx + a12 * cy + a22 * cz;

    }

}


void compute_centroid_gradient(const verts_t& centroids, verts_t& centroid_normals_normalized, implicit_function* gradou) {

    gradou->eval_gradient(centroids, &centroid_normals_normalized);
    for (int i = 0; i < centroid_normals_normalized.shape()[0]; i++) {
        REAL norm = norm_2(centroid_normals_normalized[i][0], centroid_normals_normalized[i][1], centroid_normals_normalized[i][2]);
        assert(norm != 0.);
        for (int j = 0; j < 3; j++) {
            centroid_normals_normalized[i][j]=centroid_normals_normalized[i][j]/norm;
            assert(centroid_normals_normalized[i][j] <= 1.);
            assert(centroid_normals_normalized[i][j] >= -1.);
        }
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

void vertex_apply_qem(
    verts_t* verts, const faces_t faces,
    const verts_t centroids,
    const std::vector< std::vector<int>> vertex_neighbours_list,
    const verts_t centroid_gradients,
    const boost::multi_array<bool, 1>& treated)
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
        cout << vi << endl;
        continue;
    }

    get_A_b(nlist, centroids, centroid_gradients, &A, &b);
    SVD(A, u, s, v); // the SVD

    // in python the SVD values of s are sorted by the svd function, this is a possible workaround
    // (we may need to keep the A=u*s*v equality, which is done this way)
    // assert(np.allclose(A, np.dot(u, np.dot(np.diag(s), v)))) validation assert,
    // also note that u and v are supposed to be unitary
    // so we can add: assert(u.T = u^-1) and assert(v.T == v^-1)

    REAL maxi = max(s[0][0], max(s[1][1], s[2][2]));

    REAL tau = 680;
    int rank = 0;
    if (s[0][0]/maxi < 1./tau) {
        s[0][0] = 0.;
    } else {
        ++rank;
    }
    if (s[1][1]/maxi < 1./tau) {
        s[1][1] = 0.;
    } else {
        ++rank;
    }

    if (s[2][2]/maxi < 1./tau) {
        s[2][2] = 0.;
    } else {
        ++rank;
    }

    // assert s[0] == np.max(s)  asserts that SVD produces descending order eigenvalues
    y[0] = v[0][0]*(*verts)[vi][0] + v[1][0]*(*verts)[vi][1] + v[2][0]*(*verts)[vi][2];
    y[1] = v[0][1]*(*verts)[vi][0] + v[1][1]*(*verts)[vi][1] + v[2][1]*(*verts)[vi][2];
    y[2] = v[0][2]*(*verts)[vi][0] + v[1][2]*(*verts)[vi][1] + v[2][2]*(*verts)[vi][2];



    utb[0] = - u[0][0]*b[0] - u[1][0]*b[1] - u[2][0]*b[2];
    utb[1] = - u[0][1]*b[0] - u[1][1]*b[1] - u[2][1]*b[2];
    utb[2] = - u[0][2]*b[0] - u[1][2]*b[1] - u[2][2]*b[2];

    for (int i=0; i < rank; i++) {
      if (s[i][i] != 0) {
          y[i] = utb[i]/s[i][i];
      } else {
          rank++;
      }
    }

    new_x[0] = v[0][0]*y[0] + v[0][1]*y[1] + v[0][2]*y[2];
    new_x[1] = v[1][0]*y[0] + v[1][1]*y[1] + v[1][2]*y[2];
    new_x[2] = v[2][0]*y[0] + v[2][1]*y[1] + v[2][2]*y[2];


    (*verts)[vi][0] = new_x[0];
    (*verts)[vi][1] = new_x[1];
    (*verts)[vi][2] = new_x[2];
    }


}


void centroids_projection(mp5_implicit::implicit_function* object, std::vector<REAL>& result_verts, const std::vector<int>& result_faces) {
    boost::array<unsigned int, 2> verts_shape = { (unsigned int)result_verts.size()/3 , 3 };
    boost::multi_array<REAL, 2> verts(verts_shape);
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

    boost::array<unsigned int, 2> centroids_shape = { result_faces.size()/3 , 3 };
    boost::multi_array<REAL, 2> centroids(centroids_shape);

    compute_centroids(faces, verts, centroids);
    boost::array<unsigned int, 1> treated_shape = {result_faces.size()};
    boost::multi_array<bool, 1> treated(treated_shape);
    set_centers_on_surface(object, centroids, average_edge, treated);

    std::vector< std::vector<int>> vertex_neighbours_list;
    vertex_neighbours_list = make_neighbour_faces_of_vertex(verts, faces);

    boost::multi_array<REAL, 2> centroid_gradients(centroids_shape);

    compute_centroid_gradient(centroids, centroid_gradients, object);

    vertex_apply_qem(&verts, faces, centroids, vertex_neighbours_list, centroid_gradients, treated);

    for (int i=0; i < verts.shape()[0]; i++) {
        result_verts[i*3+0] = verts[i][0];
        result_verts[i*3+1] = verts[i][1];
        result_verts[i*3+2] = verts[i][2];
    }

}
