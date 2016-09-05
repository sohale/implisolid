

template<typename T>
void here(T arg) {
    std::clog << arg << std::endl << std::flush;
}

inline bool test_points_sign(verts_t& x_vectorized, const mp5_implicit::implicit_function& object, REAL ROOT_TOLERANCE, REAL sign) {
    assert(sign == +1 || sign == -1); // || sign == 0.0);

    clog << "test_points_sign" << endl;

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
            clog <<  v_arr[i] << " " << s1 << " " << sign << " :" << i << endl;
        }
        everything_alright = everything_alright && ok;
    }
    clog << "bool" << everything_alright << endl;
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
            //     std::clog << "["<<i<<"]"<< s1 << " " << s2 << " v1:" << v1_arr[i] << " v2:" << v2_arr[i] << endl;
            // }
            assert1 = assert1 && ok;
        }
        if (!assert1) clog << "";
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
    vectorized_vect  x_mid(x_mid_shape);
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
            clog << "projection treated this much points" << endl;
            clog << solved_count << endl;
            break;
        }

    }

}








