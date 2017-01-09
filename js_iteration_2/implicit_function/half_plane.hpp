#pragma once

#include "implicit_function.hpp"
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
#include "boost/multi_array.hpp"
#include "Eigen/Dense"
#include <iostream>

//#include "../basic_functions.hpp"

namespace mp5_implicit {

using Eigen::Dynamic;


class half_plane : public implicit_functions::transformable_implicit_function {

protected:
    //REAL h;

    Eigen::Matrix<REAL, 4, 4> inv_transf_matrix;
    Eigen::Matrix<REAL, 3, 3> inv_transf_matrix_3_3;
    Eigen::Matrix<REAL, 3, 1> inv_transf_matrix_neg_xyz;


    Eigen::Matrix<REAL, 3, 1> plane_vector;
    Eigen::Matrix<REAL, 3, 1> plane_point;

    virtual bool integrity_invariant() const {return true;};



public:

    static void getHalfPlaneParameters(
        Eigen::Matrix<REAL, 3, 4>& matrix, Eigen::Matrix<REAL, 3, 1>& plane_vector, Eigen::Matrix<REAL, 3, 1>& plane_point,
        const pt::ptree& shapeparams_dict
    ){

        int i = 0;
        int j = 0;
        for (const pt::ptree::value_type &element : shapeparams_dict.get_child("matrix")){



            REAL x = element.second.get_value<float>();

            matrix(i, j) = x;

            std::cout << matrix(i, j) << std::endl;

            // my_assert(j == 3, "there shuold be three points")

            if (j == 3) {
                j = 0;
                i++;
            } else {
                j++;
            }
            if (i==3) {
                break;
            }
        }

        // my_assert(i == 3, "there should be three row");
        // my_assert(j == 2, "last row there should be three col");

        std::cout << "---matrix in half plane parameter---" << std::endl;

        std::cout << matrix << std::endl;
        // my_assert(!(A(0,0)==0 && A(1,0)==0 && A(2,0)==0), "possibly A is not initialised correclt");
        // my_assert(!(B(0,0)==0 && B(1,0)==0 && B(2,0)==0), "possibly B is not initialised correclt");

        std::cout << "getplane_vector" << std::endl;
        int getv_counter = 0;
        for (const pt::ptree::value_type &element : shapeparams_dict.get_child("plane_vector")) {
                //std::clog << "matrix value : " << x << std::endl;
            plane_vector(getv_counter, 0) = element.second.get_value<float>();

            std::cout << element.second.get_value<float>() << std::endl;
            // my_assert(getv_counter<=2, "i should not exceed number 2");
            getv_counter++;
        }

        std::cout << plane_vector << std::endl;

        std::cout << "plane_point" << std::endl;
        getv_counter = 0;
        for (const pt::ptree::value_type &element : shapeparams_dict.get_child("plane_point")) {
                //std::clog << "matrix value : " << x << std::endl;
            plane_point(getv_counter, 0) = element.second.get_value<float>();

            std::cout << element.second.get_value<float>() << std::endl;
            // my_assert(getv_counter<=2, "i should not exceed number 2");
            getv_counter++;
        }

        std::cout << plane_point << std::endl;
    }


    static void getHalfPlaneParametersMatrixOnly(
        Eigen::Matrix<REAL, 3, 4>& matrix, 
        const pt::ptree& shapeparams_dict
    ){

        int i = 0;
        int j = 0;
        for (const pt::ptree::value_type &element : shapeparams_dict.get_child("matrix")){



            REAL x = element.second.get_value<float>();

            matrix(i, j) = x;

            std::cout << matrix(i, j) << std::endl;

            // my_assert(j == 3, "there shuold be three points")

            if (j == 3) {
                j = 0;
                i++;
            } else {
                j++;
            }
            if (i==3) {
                break;
            }
        }

        // my_assert(i == 3, "there should be three row");
        // my_assert(j == 2, "last row there should be three col");

        std::cout << "---matrix in half plane parameter---" << std::endl;

        std::cout << matrix << std::endl;
        // my_assert(!(A(0,0)==0 && A(1,0)==0 && A(2,0)==0), "possibly A is not initialised correclt");
        // my_assert(!(B(0,0)==0 && B(1,0)==0 && B(2,0)==0), "possibly B is not initialised correclt");
    }

    half_plane(Eigen::Matrix<REAL, 3, 4> matrix,
               Eigen::Matrix<REAL, 3, 1> plane_vector,
               Eigen::Matrix<REAL, 3, 1> plane_point) {

        this->plane_vector = plane_vector;
        this->plane_vector = (this->plane_vector.array()/this->plane_vector.norm()).matrix();

        std::cout << "----this->plane_vector----" << std::endl;
        std::cout << this->plane_vector << std::endl;

        assert(abs(this->plane_vector.norm() - 1) < 0.0001);


        this->plane_point = plane_point;

        std::cout << "----this->plane_point----" << std::endl;
        std::cout << this->plane_point << std::endl;

        // How to make sure 12 elements are provided? (how to assert)
        Eigen::Matrix<REAL, 4, 4> matrix_4_4;

        matrix_4_4 << matrix(0, 0), matrix(0, 1), matrix(0, 2), matrix(0, 3),
                      matrix(1, 0), matrix(1, 1), matrix(1, 2), matrix(1, 3),
                      matrix(2, 0), matrix(2, 1), matrix(2, 2), matrix(2, 3),
                      0, 0, 0, 1;

        this->inv_transf_matrix = matrix_4_4.inverse();

        this->inv_transf_matrix_3_3 << inv_transf_matrix(0, 0), inv_transf_matrix(0, 1), inv_transf_matrix(0, 2),
                                      inv_transf_matrix(1, 0), inv_transf_matrix(1, 1), inv_transf_matrix(1, 2),
                                      inv_transf_matrix(2, 0), inv_transf_matrix(2, 1), inv_transf_matrix(2, 2);

        this->inv_transf_matrix_neg_xyz << inv_transf_matrix(0, 3), inv_transf_matrix(1, 3), inv_transf_matrix(2, 3);

    };

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) const {

        Eigen::Matrix<REAL, Eigen::Dynamic, 3> x_eigen_matrix(x.shape()[0], 3);
        x_eigen_matrix = vectorized_vect_to_Eigen_matrix(x);
        eigen_matrix_vector_product(this->inv_transf_matrix_3_3, this->inv_transf_matrix_neg_xyz, x_eigen_matrix);

        Eigen::Matrix<REAL, Eigen::Dynamic, 1> implicitFunctionOutput(x.shape()[0], 1);
        implicitFunctionOutput = (x_eigen_matrix.rowwise() - this->plane_point.transpose()) * this->plane_vector;
        *(output) = Eigen_matrix_to_vectorized_scalar(implicitFunctionOutput);

    };

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {
        // Eigen::Matrix<REAL, Eigen::Dynamic, 3> x_eigen_matrix(x.shape()[0], 3);
        // x_eigen_matrix = vectorized_vect_to_Eigen_matrix(x);
        // eigen_matrix_vector_product(this->inv_transf_matrix_3_3, this->inv_transf_matrix_neg_xyz, x_eigen_matrix);

        vectorized_scalar  implicitFunctionOutput(boost::extents[x.shape()[0]]);
        this->eval_implicit(x, &implicitFunctionOutput); // inplace change

        for (int i=0;i<x.shape()[0];i++) {
            if (implicitFunctionOutput[i] >= 0) { // distance to the plane, in the half plane
                (*(output))[i][0] = -this->plane_vector(0, 0);
                (*(output))[i][1] = -this->plane_vector(1, 0);
                (*(output))[i][2] = -this->plane_vector(2, 0);
            } else { // distance to the plane, not in the half plane
                (*(output))[i][0] = this->plane_vector(0, 0);
                (*(output))[i][1] = this->plane_vector(1, 0);
                (*(output))[i][2] = this->plane_vector(2, 0);
            }

        REAL g0 = (*output)[i][0]; // gx 
        REAL g1 = (*output)[i][1]; // gy
        REAL g2 = (*output)[i][2]; // gz

        (*output)[i][0] = this->inv_transf_matrix(0, 0)*g0 + this->inv_transf_matrix(1, 0)*g1 + this->inv_transf_matrix(2, 0)*g2;
        (*output)[i][1] = this->inv_transf_matrix(0, 1)*g0 + this->inv_transf_matrix(1, 1)*g1 + this->inv_transf_matrix(2, 1)*g2;
        (*output)[i][2] = this->inv_transf_matrix(0, 2)*g0 + this->inv_transf_matrix(1, 2)*g1 + this->inv_transf_matrix(2, 2)*g2;

        }

    };

    // virtual mp5_implicit::bounding_box  get_boundingbox() const {
    //     return mp5_implicit::bounding_box{-0.51, +0.51, -0.51, +0.51, -0.51, +0.51 };
    // };

};
}  // namespace mp5_implicit
