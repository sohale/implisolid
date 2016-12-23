#pragma once

#include "implicit_function.hpp"
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
#include "boost/multi_array.hpp"
#include "Eigen/Dense"
#include <iostream>
// #include <math.h>       /* sin */

namespace pt = boost::property_tree ;

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

Array<REAL, Dynamic, 1> phi(const Array<REAL, Dynamic, 1>& x)
{   
    Array<REAL, Dynamic, 1> sin_x(x.rows(),1);
    for (int i=0;i<x.rows();i++){
        sin_x(i, 0) = std::sin(x(i, 0)*2*pi); // math sin
    }
    return sin_x;
}

// Array<REAL, Dynamic, 1> phi(const Array<REAL, Dynamic, 1>& x)
// {   

//     /*

//     shape likes

//       |-| |-|
//       | | | |
//     --| | | |-
//         | |
//         |-|
//     */
//     Array<REAL, Dynamic, 1> M_shape(x.rows(),1);
//     for (int i=0;i<x.rows();i++){

//         REAL abs_x = x(i, 0) - std::floor(x(i, 0));

//         if (abs_x<0.2){
//             M_shape(i, 0) = 0;
//         } else if (0.2<=abs_x<0.5){
//             M_shape(i, 0) = 1;
//         } else if (0.5<=abs_x<0.8){
//             M_shape(i, 0) = -1;
//         } else if (0.8<=abs_x<=1.0){
//             M_shape(i, 0) = 0;
//         } else {
//             std::cout << "this should not happen..";
//         }
//     }
//     return M_shape;
// }


// not tested in screw.hpp
// bool integrity_invariant(const Matrix<REAL, 3, 1> w, 
//                          const Matrix<REAL, 3, 1> u, 
//                          const Matrix<REAL, 3, 1> v,
//                          const Matrix<REAL, 3, 3> UVW,
//                          const Matrix<REAL, 3, 3> UVW_inv,
//                          REAL slen,
//                          REAL r0){

//     bool sane = true;
//     REAL norm_tol, matrix_inv_tol, numerical_min_length;
//     norm_tol = 0.00000001;
//     matrix_inv_tol = 0.000001;
//     numerical_min_length = 0.1;

//     // u_norm_check = u.norm() - 1.0;
//     // use vectorwise norm
//     sane = sane & (abs(w.norm() - 1.0) < norm_tol);
//     sane = sane & (abs(u.norm() - 1.0) < norm_tol);
//     sane = sane & (abs(v.norm() - 1.0)  < norm_tol);

//     Matrix<REAL, 3, 3> UVW_multiply_UVW_inv = UVW * UVW_inv;
//     Matrix<REAL, 3, 3> UVW_inv_multiply_UVW_multiply = UVW_inv*UVW;

//     // skip sane = sane and self.UVW.shape == (3, 3)
//     Matrix<REAL, 3, 3> eye = Matrix<REAL, 3, 3>::Identity();
//     sane = sane & allclose(UVW_multiply_UVW_inv, eye, matrix_inv_tol);
//     sane = sane & allclose(UVW_inv_multiply_UVW_multiply, eye, matrix_inv_tol);
//     sane = sane & (slen > numerical_min_length);
//     sane = sane & (r0 > numerical_min_length);
//     sane = sane & allclose(v, u.cross(w), norm_tol);

//     return sane;
// };

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

namespace mp5_implicit {

class screw : public transformable_implicit_function {

protected:
    // unsign for slen, r0??
    REAL slen, r0, delta, twist_rate, phi0; // is the REAL defined here??
    Matrix<REAL, 3, 1> A, w, u, v; // not using Eigen::Vector3d since Vector3d has only type double 
    Matrix<REAL, 3, 3> UVW, UVW_inv;
    REAL x0, y0, z0;

    static inline Matrix<REAL, Dynamic, 1> implicitFunction(const Matrix<REAL, 3, 1>& A,
                          const Matrix<REAL, 3, 1>& w,
                          const Matrix<REAL, 3, 3>& UVW_inv,
                          const REAL& slen,
                          const REAL& r0,
                          const REAL& delta,
                          const REAL& twist_rate,
                          const REAL& phi0,
                          const Matrix<REAL, Dynamic, 3>& x)
    {   
        int num_points = x.rows();
        const Matrix<REAL, Dynamic, 3> aa(num_points, 3);
        const Matrix<REAL, 1, 3> A_transpose = A.transpose();

        Matrix<REAL, Dynamic, 1> t(num_points, 1);
        t = (x.rowwise() - A_transpose)*w; // (recenter) * w where w == 0,0,1 or 0,0,-1 is basically getting the z/height value

        Matrix<REAL, 1, Dynamic> t_transpose(1, num_points);
        t_transpose = t.transpose();

        Matrix<REAL, 3, Dynamic> p(3, num_points);
        p = (w*t_transpose).colwise() + A;

        Matrix<REAL, 3, Dynamic> ab(3, num_points);

        ab = UVW_inv * (x.transpose() - p); // ?? this should map the local coordinate polar


        Matrix<REAL, 3, Dynamic> example(3, num_points);

        Matrix<REAL, 1, Dynamic> theta(1, num_points);
        for (int i=0; i<ab.cols(); i++) {

            theta(0, i) = std::atan2(ab(1,i), ab(0, i)); // angle of any given point on polar coordinates
        }

        Matrix<REAL, Dynamic, 1> r(num_points, 1);
        r = (x - p.transpose()).rowwise().norm(); // length to center

        Matrix<REAL, Dynamic, 1> inside_ness(num_points, 1);


        inside_ness = t/slen; // whether the point is inside the range of length of screw

        Matrix<REAL, Dynamic, 1> ones(num_points, 1);
        ones = MatrixXf::Ones(num_points, 1);

        inside_ness = 1*ones - 2*((inside_ness-0.5*ones).cwiseAbs());

        Matrix<REAL, Dynamic, 1> zeros(num_points, 1);
        zeros = MatrixXf::Zero(num_points, 1);
        // inside_ness = (inside_ness.array() > 0).select(ones, zeros);

        inside_ness = ones; // debug

        REAL pi2 = pi*2; //M_PI from math.h

        Matrix<REAL, Dynamic, 1> screw_ness(num_points, 1);


        screw_ness = (
            (
                -r.array() + r0 + delta * 
                phi( t_transpose.array()/twist_rate - theta.array()/pi2 )
            )
            * inside_ness.array()).matrix();

        return screw_ness;

    };

    static inline void gradient(REAL ax, REAL ay, REAL az, REAL delta, REAL phi0, REAL twist_rate, 
                  REAL uvwi00, REAL uvwi01, REAL uvwi02, REAL uvwi10, REAL uvwi11, REAL uvwi12,
                  REAL wx, REAL wy, REAL wz,
                  REAL x, REAL y, REAL z,
                  REAL& dx, REAL& dy, REAL& dz) {

        dx = M_PI*delta*(-((uvwi00*(-std::pow(wx, 2) + 1) - uvwi01*wx*wy - uvwi02*wx*wz)*(-uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)) + (uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(uvwi10*(-std::pow(wx, 2) + 1) - uvwi11*wx*wy - uvwi12*wx*wz)/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)))/M_PI + 2*wx/twist_rate)*cos(M_PI*(2*phi0 - atan2(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wy*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - wx*wz*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z) + (1.0/2.0)*(-2*std::pow(wx, 2) + 2)*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x))/std::sqrt(std::pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + std::pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + std::pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));
        dy = M_PI*delta*(-((uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(-uvwi10*wx*wy + uvwi11*(-std::pow(wy, 2) + 1) - uvwi12*wy*wz)/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)) + (-uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(-uvwi00*wx*wy + uvwi01*(-std::pow(wy, 2) + 1) - uvwi02*wy*wz)/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)))/M_PI + 2*wy/twist_rate)*cos(M_PI*(2*phi0 - atan2(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wy*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - wy*wz*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z) + (1.0/2.0)*(-2*std::pow(wy, 2) + 2)*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y))/std::sqrt(std::pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + std::pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + std::pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));
        dz = M_PI*delta*(-((uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(-uvwi10*wx*wz - uvwi11*wy*wz + uvwi12*(-std::pow(wz, 2) + 1))/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)) + (-uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(-uvwi00*wx*wz - uvwi01*wy*wz + uvwi02*(-std::pow(wz, 2) + 1))/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)))/M_PI + 2*wz/twist_rate)*cos(M_PI*(2*phi0 - atan2(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wz*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - wy*wz*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + (1.0/2.0)*(-2*std::pow(wz, 2) + 2)*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/std::sqrt(std::pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + std::pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + std::pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));

    }

    


public: 

    // screw(Matrix<REAL, 3, 1> A, Matrix<REAL, 3, 1> B, Matrix<REAL, 3, 1> U, 
    //       REAL pitch_len, std::string profile_shape,
    //       REAL inner_diameter, REAL outer_diameter,  std::string end_type) {

    //     this->A = A;
    //     this->slen = (B - A).norm();
    //     this->w = (B-A)/slen;
    //     my_assert(slen != 0, "cannot divided by zero");
    //     this->u = U;
    //     my_assert(end_type != "0", "end_type can only be zero now");
    //     my_assert(profile_shape != "sin" , "profile_shape can only be sin now");
    //     this->r0 = (inner_diameter + outer_diameter)/2;
    //     this->delta = outer_diameter - this->r0;
    //     this->twist_rate = pitch_len;

    //     this->v = this->u.cross(this->w); // vector defined orientation, orthogon to u, w
    //     this->UVW << this->u, this->v, this->w;
    //     this->UVW_inv = this->UVW.inverse();

    // }

    screw()
    {
        std::cout << "------------------------using empty screw constructor----------------------" << std::endl;
        this->A << 0,0,1; // center of bottom
        this->w << 0,0,-1; // vector defined orientation, orthogonal to u, v
        this->u << 1,0,0; // vector defined orientation, orthogonal to u, v
        this->slen = 5; // screw length
        this->r0 = 0.5; // radius of the cylinder
        this->delta = 0.2; // how much does the screw extend out and subtract in since the phi function is between -1 to 1
        this->twist_rate = 0.4; // number of cycle ~= slen/twist_rate
        // this->phi0 = 0.0; // not used

        this->v = (this->u).cross(this->w); // vector defined orientation, orthogon to u, w
        this->UVW << this->u, this->v, this->w;
        this->UVW_inv = (this->UVW).inverse();
    }

    screw(Matrix<REAL, 4, 4> matrix, REAL pitch_len, std::string profile, 
          std::string end_type, REAL delta_ratio, Matrix<REAL, 3, 1> v)
    {   


        std::cout << matrix << endl;
        std::cout << pitch_len << endl;
        std::cout << profile << endl;
        std::cout << end_type << endl;
        std::cout << delta_ratio << endl;
        std::cout << v << endl; 

        std::cout << "------------------------using matrix screw constructor: v defined----------------------" << std::endl;
        // screw current cannot do non-uniform 

        // consider u as transformation matrix * [0,1,0],
        // consider v as transformation matrix * [1,0,0],
        // consider w as transformation matrix * [0,0,1],

        // this is a problem since if the screw is stretched, the 
        // matrix.col(0).norm() (length of u) and matrix.col(1).norm() (length of v)

        this->u << matrix(0, 0), matrix(1, 0), matrix(2, 0); // first three element from first column
        this->u = (this->u.array()/this->u.norm()).matrix();

        this->w << matrix(0, 2), matrix(1, 2), matrix(2, 2); // first three element from first column
        this->slen = this->w.norm();
        this->w = (this->w.array()/this->w.norm()).matrix();

        this->v = v; // if v is defined or if v is not defined
        REAL outer_diameter = this->v.norm(); 
        this->v = (this->v.array()/this->v.norm()).matrix();


        this->A << matrix(0, 3), matrix(1, 3), matrix(2, 3);
        this->A  = (this->A.array() - this->slen/2).matrix();
        
        REAL inner_diameter = outer_diameter/delta_ratio;
        this->r0 = inner_diameter/2;
        this->delta = (outer_diameter/2 - inner_diameter/2);
        this->twist_rate = pitch_len;

        this->UVW << this->u, this->v, this->w;
        this->UVW_inv = this->UVW.inverse();

        this->phi0 = 0.0; // not used
    }


    // }
    // screw(Matrix<REAL, 3, 4> matrix, REAL pitch_len, std::string profile, 
    //       std::string end_type, REAL delta_ratio)
    // : screw()
    // {   

    //     // screw current cannot do non-uniform 

    //     this->u = matrix.column(0)/matrix.column(0).norm();
    //     this->w = matrix.column(2)/matrix.column(2).norm();
    //     this->v = v;
    //     this->A = matrix.column(3);

    //     this->slen = w.norm();
        
    //     REAL outer_diameter = u.norm();
    //     REAL inner_diameter = outer_diameter/delta_ratio;
    //     this->r0 = (inner_diameter + outer_diameter)/2;
    //     this->delta = outer_diameter - this->r0;
    //     this->twist_rate = pitch_len;

    //     this->UVW << this->u, this->v, this->w;
    //     this->UVW_inv = this->UVW.inverse();


    // }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) const {

        Matrix<REAL, Dynamic, 3> x_eigen_matrix(x.shape()[0], 3);
        x_eigen_matrix = vectorized_vect_to_Eigen_matrix(x);

        Matrix<REAL, Dynamic, 1> implicitFunctionOutput(x.shape()[0], 1);
        implicitFunctionOutput = implicitFunction(this->A, this->w,this->UVW_inv, this->slen, this->r0,
                         this->delta, this->twist_rate, this->phi0, x_eigen_matrix);

        *(output) = Eigen_matrix_to_vectorized_scalar(implicitFunctionOutput);

    }

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        for (int j=0;j<(*output).shape()[0];j++){

            gradient(this->A(0,0), this->A(1,0), this->A(2,0), 
                     this->delta, this->phi0,
                     this->twist_rate,
                     this->UVW(0,0), this->UVW(0,1), this->UVW(0,2),
                     this->UVW(1,0), this->UVW(1,1), this->UVW(1,2),
                     this->w(0, 0), this->w(1, 0), this->w(2, 0),
                     x[j][0], x[j][1], x[j][2],
                     (*(output))[j][0], (*(output))[j][1], (*(output))[j][2]
                     );
        };
    }

    // bool integrity_invariant() const {
    // }

    virtual mp5_implicit::bounding_box getboundingbox() const {
        return mp5_implicit::bounding_box{1,2,3,4,5,6};
    }

    static void getScrewParameters(
        Matrix<REAL, 4, 4>& matrix, REAL& pitch_len, std::string& profile, 
        std::string& end_type, REAL& delta_ratio, Matrix<REAL, 3, 1>& v,
        const pt::ptree& shapeparams_dict
    ){

        int i = 0;
        int j = 0;

        std::cout << "--- here ---" << std::endl;


        for (const pt::ptree::value_type &element : shapeparams_dict.get_child("matrix")){

            std::cout << "--- in for loop ---" << std::endl;

            std::cout << element.second.get_value<float>() << std::endl;

            std::cout << "---i---" << std::endl;
            std::cout << i << std::endl;
            std::cout << "---j---" << std::endl;
            std::cout << j << std::endl;

            std::cout << "---m i j ---" << std::endl;

            REAL x = element.second.get_value<float>();

            matrix(i, j) = x;

            std::cout << matrix(i, j) << std::endl;

            // my_assert(j == 3, "there shuold be three points");
            if (j == 3) {
                j = 0;
                i++;
            } else {
                j++;
            }
        }

        // my_assert(i == 3, "there should be three row");
        // my_assert(j == 2, "last row there should be three col");

        std::cout << "---matrix in screw parameter---" << std::endl;

        std::cout << matrix << std::endl;
        // my_assert(!(A(0,0)==0 && A(1,0)==0 && A(2,0)==0), "possibly A is not initialised correclt");
        // my_assert(!(B(0,0)==0 && B(1,0)==0 && B(2,0)==0), "possibly B is not initialised correclt");

        std::cout << "getV" << std::endl;
        int getv_counter = 0;
        for (const pt::ptree::value_type &element : shapeparams_dict.get_child("v")) {
                //std::clog << "matrix value : " << x << std::endl;
            v(getv_counter, 0) = element.second.get_value<float>();

            std::cout << element.second.get_value<float>() << std::endl;
            // my_assert(getv_counter<=2, "i should not exceed number 2");
            getv_counter++;
        }
        std::cout << "getV" << std::endl;

        std::cout << "getpitch_len" << std::endl;
        std::cout << shapeparams_dict.get<REAL>("pitch") << std::endl;
        pitch_len = shapeparams_dict.get<REAL>("pitch");
        std::cout << pitch_len << std::endl;
        std::cout << "getpitch_len" << std::endl;


        std::cout << "get-profile" << std::endl;
        profile = shapeparams_dict.get<std::string>("profile");
        std::cout << profile << std::endl;
        std::cout << "get-profile" << std::endl;

        // okay 

        std::cout << "delta_ratio" << std::endl;
        std::cout << shapeparams_dict.get<REAL>("delta_ratio") << std::endl;
        delta_ratio = shapeparams_dict.get<REAL>("delta_ratio"); // this line has global affects destroy the screw 
        std::cout << delta_ratio << std::endl;
        std::cout << "delta_ratio" << std::endl;

        // od = shapeparams_dict.get<REAL>("diameter-outer"); // this line has global affects destroy the screw 

        // std::cout << "--------------------------------------" << std::endl;
        // std::cout << shapeparams_dict.get<REAL>("diameter-inner") << std::endl;
        // std::cout << shapeparams_dict.get<REAL>("diameter-outer") << std::endl;
        // std::cout << "--------------------------------------" << std::endl;


        std::cout << "get-end-typer" << std::endl;
        end_type = shapeparams_dict.get<std::string>("end_type");

        std::cout << end_type << std::endl;

        std::cout << "get-end-typer" << std::endl;
    }

};

}; //namespace