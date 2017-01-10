#pragma once
#include "implicit_function.hpp"
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
#include "Eigen/Dense"
#include <iostream>



boost::multi_array<REAL, 2> Eigen_matrix_to_vectorized_vector(const Matrix<REAL, Dynamic, 3> x)
{   
    boost::multi_array<REAL, 2> boost_matrix;
    for (int j=0;j<x.rows();j++){
        boost_matrix[j][0] = x(j,0);
        boost_matrix[j][1] = x(j,1);
        boost_matrix[j][2] = x(j,2);

    }
    return boost_matrix;
}

namespace mp5_implicit {
namespace implicit_functions {

class inf_top_bot_bound : public transformable_implicit_function {

protected:
    const REAL height = 1.0;
    implicit_functions::screw imp_func = implicit_functions::screw();
    Eigen::Matrix<REAL, 4, 4> inv_transf_matrix;
    Eigen::Matrix<REAL, 3, 3> inv_transf_matrix_3_3;
    Eigen::Matrix<REAL, 3, 1> inv_transf_matrix_neg_xyz;

    unique_ptr<implicit_function> imf_func;


public:

    inf_top_bot_bound(Eigen::Matrix<REAL, 3, 4> matrix)
    {

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

        std::cout << "height\n";
        std::cout << this->height << "\n";

        std::cout << this->inv_transf_matrix << "\n";
        std::cout << this->inv_transf_matrix_neg_xyz << "\n";
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) const {

        const int x_row_number = x.shape()[0];

        // todo inverse x..
        this->imp_func.eval_implicit(x, output);

        Eigen::Matrix<REAL, Eigen::Dynamic, 1> imp_func_output(x_row_number, 1);
        imp_func_output = vectorized_scalar_to_Eigen_matrix(*(output));


        // setup fpr eigen
        Eigen::Matrix<REAL, Eigen::Dynamic, 3> x_eigen_matrix(x_row_number, 3);
        x_eigen_matrix = vectorized_vect_to_Eigen_matrix(x);
        eigen_matrix_vector_product(this->inv_transf_matrix_3_3, this->inv_transf_matrix_neg_xyz, x_eigen_matrix);

        std::cout << "x_eigen_matrix.row(0)" << "\n";
        std::cout << x_eigen_matrix.row(0) << "\n";
        std::cout << x_eigen_matrix.row(1) << "\n";

        Eigen::Matrix<REAL, Eigen::Dynamic, 1> tbb_output(x_row_number, 1);
        // todo: change o.25 to 0.5


        tbb_output = ((x_eigen_matrix.col(2).array() - 0.25).max((x_eigen_matrix.col(2).array() + 0.25) * -1)).matrix();

        std::cout << tbb_output.row(0) << "\n";
        std::cout << tbb_output.row(1) << "\n";


        tbb_output = (imp_func_output.array().min(tbb_output.array()*-1)).matrix();

        // std::cout << "tbb_output" << "\n";
        // std::cout << tbb_output << "\n";

        *(output) = Eigen_matrix_to_vectorized_scalar(tbb_output);

    }


    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        std::cout << "top_bottom_bound new eval_gradient" << std::endl;


        Eigen::Matrix<REAL, Eigen::Dynamic, 3> x_eigen_matrix(x.shape()[0], 3);
        x_eigen_matrix = vectorized_vect_to_Eigen_matrix(x);
        eigen_matrix_vector_product(this->inv_transf_matrix_3_3, this->inv_transf_matrix_neg_xyz, x_eigen_matrix);


        // for (int i=0;i<x.shape()[0];i++) {
        //     if (x_eigen_matrix(i, 2) >= 0.25) { // distance to the plane, in the half plane
        //         (*(output))[i][0] = 0;
        //         (*(output))[i][1] = 0;
        //         (*(output))[i][2] = -1;
        //     } else if (0.25 > x_eigen_matrix(i, 2) >= 0.0) {
        //         (*(output))[i][0] = 0;
        //         (*(output))[i][1] = 0;
        //         (*(output))[i][2] = 1;
        //     } else if (0.0 > x_eigen_matrix(i, 2) >= -0.25) {
        //         (*(output))[i][0] = 0;
        //         (*(output))[i][1] = 0;
        //         (*(output))[i][2] = -1;
        //     } else if (-0.25 > x_eigen_matrix(i, 2)) { // distance to the plane, not in the half plane
        //         (*(output))[i][0] = 0;
        //         (*(output))[i][1] = 0;
        //         (*(output))[i][2] = 1;
        //     } else {
        //         std::cout << "this should never happen" << std::endl;
        //     }


        //     REAL g0 = (*output)[i][0]; // gx 
        //     REAL g1 = (*output)[i][1]; // gy
        //     REAL g2 = (*output)[i][2]; // gz

        //     (*output)[i][0] = this->inv_transf_matrix(0, 0)*g0 + this->inv_transf_matrix(1, 0)*g1 + this->inv_transf_matrix(2, 0)*g2;
        //     (*output)[i][1] = this->inv_transf_matrix(0, 1)*g0 + this->inv_transf_matrix(1, 1)*g1 + this->inv_transf_matrix(2, 1)*g2;
        //     (*output)[i][2] = this->inv_transf_matrix(0, 2)*g0 + this->inv_transf_matrix(1, 2)*g1 + this->inv_transf_matrix(2, 2)*g2;

        // }


        // new code 

        // const int x_row_number = x.shape()[0];

        // implicit value calculation
        // boost::multi_array<REAL, 1> imp_func_value_output;
        // this->imp_func.eval_implicit(x, &imp_func_value_output);

        // Eigen::Matrix<REAL, Eigen::Dynamic, 1> imp_func_output(x_row_number, 1);
        // imp_func_output = vectorized_scalar_to_Eigen_matrix(imp_func_value_output);

        // Eigen::Matrix<REAL, Eigen::Dynamic, 1> tbb_output(x_row_number, 1);
        // tbb_output = ((x_eigen_matrix.col(2).array() - 0.25).max((x_eigen_matrix.col(2).array() + 0.25) * -1)).matrix();

        this->imp_func.eval_gradient(x, output);
        for (int i=0; i < x_eigen_matrix.rows(); i++) {
            if (x_eigen_matrix(i, 2)>0.25) {
                (*output)[i][0] = 0;
                (*output)[i][1] = 0;
                (*output)[i][2] = -1;
            } else if (x_eigen_matrix(i, 2)<-0.25) {
                (*output)[i][0] = 0;
                (*output)[i][1] = 0;
                (*output)[i][2] = 1;
            }

            REAL g0 = (*output)[i][0]; // gx 
            REAL g1 = (*output)[i][1]; // gy
            REAL g2 = (*output)[i][2]; // gz

            (*output)[i][0] = this->inv_transf_matrix(0, 0)*g0 + this->inv_transf_matrix(1, 0)*g1 + this->inv_transf_matrix(2, 0)*g2;
            (*output)[i][1] = this->inv_transf_matrix(0, 1)*g0 + this->inv_transf_matrix(1, 1)*g1 + this->inv_transf_matrix(2, 1)*g2;
            (*output)[i][2] = this->inv_transf_matrix(0, 2)*g0 + this->inv_transf_matrix(1, 2)*g1 + this->inv_transf_matrix(2, 2)*g2;


        }



        return;

    };


    /*
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {


        std::cout << "tt eval_gradient start" << std::endl;

        const int x_row_number = x.shape()[0];
     
        Eigen::Matrix<REAL, Eigen::Dynamic, 3> x_eigen_matrix(x_row_number, 3);
        x_eigen_matrix = vectorized_vect_to_Eigen_matrix(x);
        eigen_matrix_vector_product(this->inv_transf_matrix_3_3, this->inv_transf_matrix_neg_xyz, x_eigen_matrix);

        // implicit value calculation

        std::cout << "tt eval_gradient gate 0" << std::endl;

        boost::multi_array<REAL, 1> imp_func_value_output;
        this->imp_func.eval_implicit(x, &imp_func_value_output);

        std::cout << "imp_func_value_output" << "\n";
        std::cout << imp_func_value_output[0] << "\n";
        std::cout << imp_func_value_output[1] << "\n";
        std::cout << imp_func_value_output[2] << "\n";

        Eigen::Matrix<REAL, Eigen::Dynamic, 1> imp_func_output(x_row_number, 1);
        imp_func_output = vectorized_scalar_to_Eigen_matrix(imp_func_value_output);

        std::cout << "---imp_func_output start" << "\n";
        std::cout << imp_func_output.row(0) << "\n";
        std::cout << imp_func_output.row(1) << "\n";
        std::cout << imp_func_output.row(2) << "\n";
        std::cout << "---imp_func_output end" << "\n";


        Eigen::Matrix<REAL, Eigen::Dynamic, 1> tbb_output(x_row_number, 1);
        // todo: change o.25 to 0.5
        tbb_output = ((x_eigen_matrix.col(2).array() - 0.25).max((x_eigen_matrix.col(2).array() + 0.25) * -1)).matrix();

        // my_assert(imp_func_output.rows() == tbb_output.rows(), "row number diff");


        boost::multi_array<REAL, 2> output_gradient;
        this->imp_func.eval_gradient(x,  &output_gradient);
        Eigen::Matrix<REAL, Eigen::Dynamic, 3> imp_func_gradient(x_row_number, 3);

        imp_func_gradient = vectorized_vect_to_Eigen_matrix(output_gradient);

        Eigen::Matrix<REAL, Eigen::Dynamic, 3> eigen_gradient(x_row_number, 3);
        eigen_gradient.col(0).setZero();
        eigen_gradient.col(1).setZero();

        std::cout << "eval_gradient middle" << std::endl;


        // std::cout << "x_eigen_matrix.col(2) before" << "\n";
        // std::cout << x_eigen_matrix.col(2).transpose() << "\n";

        // eigen_gradient.col(2) = (x_eigen_matrix.col(2).array()>=0.5 && (0>=x_eigen_matrix.col(2).array()>=-0.5)).select(-1, 1);

        for (int i=0;i<x_row_number;i++) {

            if (x_eigen_matrix(i, 2) >= 0.25 || (0>=x_eigen_matrix(i, 2) && x_eigen_matrix(i, 2) >=-0.25)) {
                // cout << x_eigen_matrix(i, 2) << "true\n";
                eigen_gradient(i, 2) = -1;

            } else {
                // cout << x_eigen_matrix(i, 2) << "false\n";

                eigen_gradient(i, 2) = 1;
            }
        }

        // std::cout << "x_eigen_matrix.col(2) after" << "\n";
        // std::cout << eigen_gradient.col(2) << "\n";

        Eigen::Matrix<REAL, Eigen::Dynamic, 3> comparison_output(x_row_number, 3);

        for (int i=0;i<x_row_number;i++) {
            if (imp_func_output(i, 0) > tbb_output(i, 0)) {
                comparison_output.row(i) = eigen_gradient.row(i);
            } else {
                comparison_output.row(i) = imp_func_gradient.row(i);
            }
        }

        
        // std::cout << "comparison_output \n";
        // std::cout << comparison_output;

        for (int i=0;i<x_row_number;i++) {
            REAL g0 = comparison_output(i, 0);
            REAL g1 = comparison_output(i, 1);
            REAL g2 = comparison_output(i, 2);

            //todo: vectorise using eigen
           (*output)[i][0] = this->inv_transf_matrix(0, 0)*g0 + this->inv_transf_matrix(1, 0)*g1 + this->inv_transf_matrix(2, 0)*g2;
           (*output)[i][1] = this->inv_transf_matrix(0, 1)*g0 + this->inv_transf_matrix(1, 1)*g1 + this->inv_transf_matrix(2, 1)*g2;
           (*output)[i][2] = this->inv_transf_matrix(0, 2)*g0 + this->inv_transf_matrix(1, 2)*g1 + this->inv_transf_matrix(2, 2)*g2;
        }
        
        for (int i=0;i<x_row_number;i++) {
                (*output)[i][0] = 0;
                (*output)[i][1] = 1;
                (*output)[i][2] = 1;
        }

        std::cout << "eval_gradient over" << std::endl;

    }
    
    */
};


} // implicit_functions namespace
} // mp5_implicit namespace