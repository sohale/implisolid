#pragma once

#include "implicit_function.hpp"
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"


#include "boost/multi_array.hpp"
#include "Eigen/Dense"
#include <iostream>
#include <random>     /* srand, rand */
#include <math.h>       /* sin */


using Eigen::Matrix;
using Eigen::MatrixXf;
using Eigen::Dynamic;
using Eigen::Array;

const REAL pi = 3.1415926535897;

// Array<REAL, Dynamic, 1> phi(const Array<REAL, Dynamic, 1> x)
// {   
//     // std::cout << x << std::endl;
//     return (2*(x - Eigen::floor(x)) - 1.0).abs()*2.0 - 1.0;

// }

// Array<REAL, Dynamic, 1> phi(const Array<REAL, Dynamic, 1> x)
// {   
//     // std::cout << x << std::endl;
//     return (2*(x - Eigen::floor(x)) - 1.0).abs()*2.0 - 1.0;

// }

// def phi(x):
//     return np.sin(x*2*np.pi);

// Array<REAL, Dynamic, 1> phi(const Array<REAL, Dynamic, 1>& x)
// {   
//     Array<REAL, Dynamic, 1> sin_x(x.rows(),1);
//     for (int i=0;i<x.rows();i++){
//         sin_x(i, 0) = sin(x(i, 0)*2*pi); // math sin
//     }
//     return sin_x;
// }

Array<REAL, Dynamic, 1> phi(const Array<REAL, Dynamic, 1>& x)
{   

    /*

    shape likes

      |-| |-|
      | | | |
    --| | | |-
        | |
        |-|
    */
    Array<REAL, Dynamic, 1> M_shape(x.rows(),1);
    for (int i=0;i<x.rows();i++){

        REAL abs_x = x(i, 0) - std::floor(x(i, 0));

        if (abs_x<0.2){
            M_shape(i, 0) = 0;
        } else if (0.2<=abs_x<0.5){
            M_shape(i, 0) = 1;
        } else if (0.5<=abs_x<0.8){
            M_shape(i, 0) = -1;
        } else if (0.8<=abs_x<=1.0){
            M_shape(i, 0) = 0;
        } else {
            std::cout << "this should not happen..";
        }
    }
    return M_shape;
}


Matrix<REAL, Dynamic, Dynamic> vectorized_vect_to_Eigen_matrix(const vectorized_vect x)
{   
    Matrix<REAL, Dynamic, Dynamic> eigen_matrix(x.shape()[0], x.shape()[1]);
    for (int i=0;i<x.shape()[0];i++){
        for (int j=0;j<x.shape()[1];j++){
            eigen_matrix(i,j) = x[i][j];
        }
    }
    return eigen_matrix;
}


boost::multi_array<REAL, 1> Eigen_matrix_to_vectorized_scalar(const Matrix<REAL, Dynamic, 1> x)
{   
    boost::multi_array<REAL, 1> boost_matrix(boost::extents[x.rows()]);
    for (int j=0;j<x.rows();j++){
        boost_matrix[j] = x(j,0);
    }
    return boost_matrix;
}



Matrix<REAL, Dynamic, 1> implicitFunction(const Matrix<REAL, 3, 1> A,
                      const Matrix<REAL, 3, 1> w,
                      const Matrix<REAL, 3, 3> UVW_inv,
                      const REAL slen,
                      const REAL r0,
                      const REAL delta,
                      const REAL twist_rate,
                      const REAL phi0,
                      const Matrix<REAL, Dynamic, 3>& x)
{   
    int num_points = x.rows();
    const Matrix<REAL, Dynamic, 3> aa(num_points, 3);
    const Matrix<REAL, 1, 3> A_transpose = A.transpose();
    // cout<< "tiger" << endl;
    // cout<< x << endl;
    // cout<< A_transpose << endl;

    Matrix<REAL, Dynamic, 1> t(num_points, 1);
    t = (x.rowwise() - A_transpose)*w; // (recenter) * w where w == 0,0,1 or 0,0,-1 is basically getting the z/height value

    Matrix<REAL, 1, Dynamic> t_transpose(1, num_points);
    t_transpose = t.transpose();

    Matrix<REAL, 3, Dynamic> p(3, num_points);
    p = (w*t_transpose).colwise() + A;

    Matrix<REAL, 3, Dynamic> ab(3, num_points);

    std::cout << "-------UVW_inv-------" << std::endl;
    std::cout << UVW_inv << std::endl;
    // std::cout << (x.transpose() - p) << std::endl;

    ab = UVW_inv * (x.transpose() - p); // ?? this should map the local coordinate polar


    Matrix<REAL, 3, Dynamic> example(3, num_points);

// [[  1.20000000e+01   0.00000000e+00  -1.00000000e+00   0.00000000e+00
//     1.33100000e+03]
//  [  2.30000000e+01   0.00000000e+00   1.00000000e+00   0.00000000e+00
//     1.22100000e+03]
//  [  0.00000000e+00   1.11022302e-16   0.00000000e+00  -8.32667268e-17
//     0.00000000e+00]]

    // example << 1.20000000e+01, 0, -1, 0, 1.33100000e+03, 
    //            2.30000000e+01, 0.00000000e+00, 1.00000000e+00, 0.00000000e+00, 1.22100000e+03,
    //            0.00000000e+00, -1.11022302e-16, 0.00000000e+00, 8.32667268e-17, 0.00000000e+00;

    // std::cout << UVW_inv*example << std::endl;

    // std::cout << "----------ab------------" << std::endl;
    // std::cout << ab << std::endl;
    // std::cout << "----------ab------------" << std::endl;

    // REAL theta_array[ab.cols()];

    // for (int i=0; i<ab.cols(); i++) {
    //     cout << atan2(ab.row(1)[i], ab.row(2)[i]) << endl;
    //     theta_array[i] = atan2(ab.row(1)[i], ab.row(2)[i]);
    //     cout << theta_array[i] << endl;
    // }

    // std::uniform_real_distribution<REAL> dis(0, 2*pi);
    // std::random_device rd;
    // std::mt19937 gen(rd());

    Matrix<REAL, 1, Dynamic> theta(1, num_points);
    for (int i=0; i<ab.cols(); i++) {
        // std::cout << "----------ab 1 i------------" << std::endl;
        // std::cout << ab(2,i) << std::endl; // always 0/
        // std::cout << ab(1,i) << std::endl; // always 0/
        // std::cout << ab(0, i) << std::endl;
        // std::cout << std::atan2(ab(1,i), ab(0, i)) << std::endl;

        theta(0, i) = std::atan2(ab(1,i), ab(0, i)); // angle of any given point on polar coordinates
    }

    // std::cout << "-----------theta-----------" << std::endl;
    // std::cout << theta << std::endl;

    Matrix<REAL, Dynamic, 1> r(num_points, 1);
    r = (x - p.transpose()).rowwise().norm(); // length to center

    Matrix<REAL, Dynamic, 1> inside_ness(num_points, 1);


    inside_ness = t/slen; // whether the point is inside the range of length of screw
    // cout << inside_ness << endl;
    // inside_ness = 1 - 2*(inside_ness - 0.5).abs()

    Matrix<REAL, Dynamic, 1> ones(num_points, 1);
    ones = MatrixXf::Ones(num_points, 1);

    inside_ness = 1*ones - 2*((inside_ness-0.5*ones).cwiseAbs());

    Matrix<REAL, Dynamic, 1> zeros(num_points, 1);
    zeros = MatrixXf::Zero(num_points, 1);
    inside_ness = (inside_ness.array() > 0).select(ones, zeros);

    REAL pi2 = pi*2; //M_PI from math.h

    Matrix<REAL, Dynamic, 1> screw_ness(num_points, 1);

    // std::cout << "t_transpose.array()" << std::endl;
    // std::cout << t_transpose.array() << std::endl;

    screw_ness = (
        (
            -r.array() + r0 + delta * 
            phi( t_transpose.array()/twist_rate - theta.array()/pi2 )
        )
        * inside_ness.array()).matrix();

    // std::cout << "-------------t and theta array-------------" << std::endl;
    // std::cout << t_transpose.array() << std::endl;
    // std::cout << theta.array()/pi2 << std::endl;

    // screw_ness = 
    //     (
    //         -r.array() + r0 + // scale down the triangle
    //         phi( t_transpose.array()) // theta.array()/pi2 is always 0 or 0.5
    //     ).matrix();

    Matrix<REAL, Dynamic, 1> lidness1(num_points, 1);
    lidness1 = (slen - t.array()).matrix();

    Matrix<REAL, Dynamic, 3> m3(num_points, 3);

    m3 << t, lidness1, screw_ness;

    return m3.rowwise().minCoeff();

};


namespace mp5_implicit {

class screw : public transformable_implicit_function {

protected:
    // unsign for slen, r0??
    REAL slen, r0, delta, twist_rate, phi0; // is the REAL defined here??
    Matrix<REAL, 3, 1> A, w, u, v; // not using Eigen::Vector3d since Vector3d has only type double 
    Matrix<REAL, 3, 3> UVW, UVW_inv;
    REAL x0, y0, z0;
public: 
    screw(Matrix<REAL, 3, 1> A,
          Matrix<REAL, 3, 1> w,
          Matrix<REAL, 3, 1> u,
          REAL slen,
          REAL r0,
          REAL delta,
          REAL twist_rate,
          REAL phi0){

        this->A = A;
        this->w = w;
        this->u = u;
        this->slen = slen;
        this->r0 = r0;
        this->delta = delta;
        this->twist_rate = twist_rate;
        this->phi0 = phi0;

        this->v = u.cross(w); // cross_product in cross_product.hpp
        this->UVW << u, v, w;
        this->UVW_inv << 1, -0, -0, -0, 1, -0, -0, -0, -1;

    }

    screw()
    {

        Matrix<REAL, 3, 1> A, u, v, w; // using Eigen::Vector3d
        REAL slen, r0, delta, twist_rate, phi0;
        Matrix<REAL, 3, 3> UVW, UVW_inv;

        A << 0,0,0.9; // center
        w << 0,0,-1; // vector defined orientation, orthogonal to u, v
        u << 1,0,0; // vector defined orientation, orthogonal to w, v

        slen = 5; // screw length
        r0 = 0.4; // radius of the cylinder
        delta = 0.1; // ??
        twist_rate = 0.4; // number of cycle ~= slen/twist_rate
        phi0 = 0.0; // ?? default value, not used

        this->A = A;
        this->w = w;
        this->u = u;
        this->slen = slen;
        this->r0 = r0;
        this->delta = delta;
        this->twist_rate = twist_rate;
        this->phi0 = phi0;

        this->v = u.cross(w); // vector defined orientation, orthogon to u, w
        this->UVW << u, this->v, w;


        this->UVW_inv = this->UVW.inverse();
        // this->UVW_inv << 1, -0, -0, -0, 1, -0, -0, -0, -1;

        // cout << "------------------------------" << endl;
        // cout << this->v << endl;
        // cout << u.cross(w) << endl;
        // cout << this->v << endl;
        // cout << this->v << endl;
        // cout << this->UVW << endl;
        // cout << UVW.inverse() << endl;
        // cout << this->UVW_inv << endl;

    }


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) const {


        Matrix<REAL, Dynamic, 3> x_eigen_matrix(x.shape()[0], 3);
        x_eigen_matrix = vectorized_vect_to_Eigen_matrix(x);

        // std::cout << x_eigen_matrix << endl;
        // std::cout << "-------------------" << endl;

        Matrix<REAL, Dynamic, 1> implicitFunctionOutput(x.shape()[0], 1);
        implicitFunctionOutput = implicitFunction(this->A, this->w,this->UVW_inv, this->slen, this->r0,
                         this->delta, this->twist_rate, this->phi0, x_eigen_matrix);

        // std::cout << "-----------implicitFunctionOutput--------" << endl;
        // std::cout << implicitFunctionOutput << endl;
        // std::cout << "-----------implicitFunctionOutput--------" << endl;

        *(output) = Eigen_matrix_to_vectorized_scalar(implicitFunctionOutput);
        // for (int j=0;j<(*output).shape()[0];j++){
        //     std::cout << (*output)[j] << " ";
        // };
        // std::cout << "-------------------" << endl;
        // std::cout << endl;
        // std::cout << "-------------------" << endl;

    }

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        for (int j=0;j<(*output).shape()[0];j++){
            (*output)[j][0] = 1;
            (*output)[j][1] = 1;
            (*output)[j][2] = 1;
        };
    }

    // bool integrity_invariant() const {
    //     REAL norm_tol, matrix_inv_tol;
    //     norm_tol = 0.00000001;
    //     matrix_inv_tol = 0.000001;
    // }

    virtual mp5_implicit::bounding_box getboundingbox() const {
        return mp5_implicit::bounding_box{1,2,3,4,5,6};
    }

};

} //namespace