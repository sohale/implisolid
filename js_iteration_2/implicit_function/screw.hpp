#pragma once

#include "implicit_function.hpp"
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"


#include "boost/multi_array.hpp"
#include "Eigen/Dense"
#include <iostream>


using Eigen::Matrix;
using Eigen::MatrixXf;
using Eigen::Dynamic;
using Eigen::Array;

#include <math.h>       /* sin */


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


Array<REAL, Dynamic, 1> phi(const Array<REAL, Dynamic, 1>& x)
{   
    Array<REAL, Dynamic, 1> sin_x(x.rows(),1);
    for (int i=0;i<x.rows();i++){
        sin_x(i, 0) = sin(x(i, 0)*2*pi); // math sin
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
    // std::cout << "-------x-------" << std::endl;
    // std::cout << x << std::endl;

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


    // std::cout << (x.transpose() - p) << std::endl;

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
    inside_ness = (inside_ness.array() > 0).select(ones, zeros);

    // inside_ness = ones; // debug

    REAL pi2 = pi*2; //M_PI from math.h

    Matrix<REAL, Dynamic, 1> screw_ness(num_points, 1);


    screw_ness = (
        (
            -r.array() + r0 + delta * 
            phi( t_transpose.array()/twist_rate - theta.array()/pi2 )
        )
        * inside_ness.array()).matrix();

    // std::cout << "-------screw_ness-------" << std::endl;
    // std::cout << screw_ness << std::endl;

    // Matrix<REAL, Dynamic, 1> lidness1(num_points, 1);
    // lidness1 = (slen - t.array()).matrix();

    // Matrix<REAL, Dynamic, 3> m3(num_points, 3);

    // m3 << t, lidness1, screw_ness;

    // return m3.rowwise().minCoeff();

    // std::cout << "-------------x and screw_ness----------" <<std::endl;
    // for(int i=0;i<=5;i++){
    //     std::cout << x.row(i) <<std::endl;
    //     std::cout << screw_ness.row(i) <<std::endl;
    // }
    return screw_ness;

};

// Matrix<REAL, Dynamic, 1> inside_ness_only(const Matrix<REAL, 3, 1> A,
//                       const Matrix<REAL, 3, 1> w,
//                       const Matrix<REAL, 3, 3> UVW_inv,
//                       const REAL slen,
//                       const REAL r0,
//                       const REAL delta,
//                       const REAL twist_rate,
//                       const REAL phi0,
//                       const Matrix<REAL, Dynamic, 3>& x)
// {
//     int num_points = x.rows();
//     const Matrix<REAL, Dynamic, 3> aa(num_points, 3);
//     const Matrix<REAL, 1, 3> A_transpose = A.transpose();

//     Matrix<REAL, Dynamic, 1> t(num_points, 1);
//     t = (x.rowwise() - A_transpose)*w; // (recenter) * w where w == 0,0,1 or 0,0,-1 is basically getting the z/height value

//     Matrix<REAL, Dynamic, 1> inside_ness(num_points, 1);

//     inside_ness = t/slen; // whether the point is inside the range of length of screw

//     Matrix<REAL, Dynamic, 1> ones(num_points, 1);
//     ones = MatrixXf::Ones(num_points, 1);

//     inside_ness = 1*ones - 2*((inside_ness-0.5*ones).cwiseAbs());

//     Matrix<REAL, Dynamic, 1> zeros(num_points, 1);
//     zeros = MatrixXf::Zero(num_points, 1);
//     inside_ness = (inside_ness.array() > 0).select(ones, zeros);

//     return inside_ness;
// }
void gradient(float ax, float ay, float az, float delta, float phi0, float twist_rate, 
              float uvwi00, float uvwi01, float uvwi02, float uvwi10, float uvwi11, float uvwi12,
              float wx, float wy, float wz,
              float x, float y, float z,
              vectorized_vect* output, int j, bool write) {

    
   // UVM_inv hard coded !!!

   if (write){
       // std::cout << "-------------- start here --------------" <<std::endl;
       // std::cout << ax <<std::endl;
       // std::cout << ay <<std::endl;
       // std::cout << az <<std::endl;
       // std::cout << delta <<std::endl;
       // std::cout << phi0 <<std::endl;
       // std::cout << twist_rate <<std::endl;
       // std::cout << wx <<std::endl;
       // std::cout << wy <<std::endl;
       // std::cout << wz <<std::endl;
       // std::cout << x <<std::endl;
       // std::cout << y <<std::endl;
       // std::cout << z <<std::endl;

       // std::cout << M_PI*delta*(-(-wx*wy*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)) + (-pow(wx, 2) + 1)*(ay + wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) - y)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)))/M_PI + 2*wx/twist_rate)*cos(M_PI*(2*phi0 - atan2(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, -ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wy*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - wx*wz*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z) + (1.0L/2.0L)*(-2*pow(wx, 2) + 2)*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x))/sqrt(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2)) << std::endl;
       // std::cout << M_PI*delta*(-(-wx*wy*(ay + wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) - y)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)) + (-pow(wy, 2) + 1)*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)))/M_PI + 2*wy/twist_rate)*cos(M_PI*(2*phi0 - atan2(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, -ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wy*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - wy*wz*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z) + (1.0L/2.0L)*(-2*pow(wy, 2) + 2)*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y))/sqrt(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2)) << std::endl;
       // std::cout << M_PI*delta*(-(-wx*wz*(ay + wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) - y)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)) - wy*wz*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)))/M_PI + 2*wz/twist_rate)*cos(M_PI*(2*phi0 - atan2(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, -ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wz*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - wy*wz*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + (1.0L/2.0L)*(-2*pow(wz, 2) + 2)*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/sqrt(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2)) << std::endl;
       // std::cout << M_PI <<std::endl;
   }

   (*output)[j][0] = M_PI*delta*(-(-wx*wy*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)) + (-pow(wx, 2) + 1)*(ay + wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) - y)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)))/M_PI + 2*wx/twist_rate)*cos(M_PI*(2*phi0 - atan2(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, -ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wy*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - wx*wz*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z) + (1.0L/2.0L)*(-2*pow(wx, 2) + 2)*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x))/sqrt(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));
   (*output)[j][1] = M_PI*delta*(-(-wx*wy*(ay + wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) - y)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)) + (-pow(wy, 2) + 1)*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)))/M_PI + 2*wy/twist_rate)*cos(M_PI*(2*phi0 - atan2(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, -ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wy*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - wy*wz*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z) + (1.0L/2.0L)*(-2*pow(wy, 2) + 2)*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y))/sqrt(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));
   (*output)[j][2] = M_PI*delta*(-(-wx*wz*(ay + wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) - y)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)) - wy*wz*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2)))/M_PI + 2*wz/twist_rate)*cos(M_PI*(2*phi0 - atan2(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, -ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x)/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wz*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - wy*wz*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + (1.0L/2.0L)*(-2*pow(wz, 2) + 2)*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/sqrt(pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));

   // std::cout << (*output)[j][0] << std::endl;
   // std::cout << (*output)[j][0] << std::endl;
   // std::cout << (*output)[j][0] << std::endl;

    // (*output)[j][0] = M_PI*delta*(-((uvwi00*(-std::pow(wx, 2) + 1) - uvwi01*wx*wy - uvwi02*wx*wz)*(-uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)) + (uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(uvwi10*(-std::pow(wx, 2) + 1) - uvwi11*wx*wy - uvwi12*wx*wz)/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)))/M_PI + 2*wx/twist_rate)*cos(M_PI*(2*phi0 - atan2(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wy*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - wx*wz*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z) + (1.0/2.0)*(-2*std::pow(wx, 2) + 2)*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x))/std::sqrt(std::pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + std::pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + std::pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));
    // (*output)[j][1] = M_PI*delta*(-((uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(-uvwi10*wx*wy + uvwi11*(-std::pow(wy, 2) + 1) - uvwi12*wy*wz)/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)) + (-uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(-uvwi00*wx*wy + uvwi01*(-std::pow(wy, 2) + 1) - uvwi02*wy*wz)/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)))/M_PI + 2*wy/twist_rate)*cos(M_PI*(2*phi0 - atan2(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wy*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - wy*wz*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z) + (1.0/2.0)*(-2*std::pow(wy, 2) + 2)*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y))/std::sqrt(std::pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + std::pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + std::pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));
    // (*output)[j][2] = M_PI*delta*(-((uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(-uvwi10*wx*wz - uvwi11*wy*wz + uvwi12*(-std::pow(wz, 2) + 1))/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)) + (-uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) - uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))*(-uvwi00*wx*wz - uvwi01*wy*wz + uvwi02*(-std::pow(wz, 2) + 1))/(std::pow(uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2) + std::pow(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), 2)))/M_PI + 2*wz/twist_rate)*cos(M_PI*(2*phi0 - atan2(uvwi10*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi11*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi12*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z), uvwi00*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) + uvwi01*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + uvwi02*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/M_PI + 2*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z))/twist_rate)) - (-wx*wz*(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x) - wy*wz*(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y) + (1.0/2.0)*(-2*std::pow(wz, 2) + 2)*(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z))/std::sqrt(std::pow(-ax - wx*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + x, 2) + std::pow(-ay - wy*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + y, 2) + std::pow(-az - wz*(wx*(-ax + x) + wy*(-ay + y) + wz*(-az + z)) + z, 2));


}



// Matrix<REAL, Dynamic, 1> implicitFunction(const Matrix<REAL, 3, 1> A,
//                       const Matrix<REAL, 3, 1> w,
//                       const Matrix<REAL, 3, 3> UVW_inv,
//                       const REAL slen,
//                       const REAL r0,
//                       const REAL delta,
//                       const REAL twist_rate,
//                       const REAL phi0,
//                       const Matrix<REAL, Dynamic, 3>& x)
// {
//     int num_points = x.rows();
//     Array<REAL, Dynamic, 1> x_square(num_points, 1);
//     x_square = x.col(0).array();
//     Array<REAL, Dynamic, 1> y_square(num_points, 1);
//     y_square = x.col(1).array();
//     Array<REAL, Dynamic, 1> z_square(num_points, 1);
//     z_square = x.col(2).array();

//     // std::cout << "-------------------------tiger debug-------------------" <<  std::endl;
//     // std::cout << (-(1/x_square.array()).square()-(1/y_square.array()).square()-(1/z_square.array()).square()+1).matrix() << endl;

//     // std::cout << (-x_square.array()-y_square.array()-z_square.array()+9).matrix() << endl;


//     return ((1/x_square.array()).square()+(1/y_square.array()).square()+(1/z_square.array()).square()-2).matrix();
// }

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

        A << 0,0,1; // center
        w << 0,0,-1; // vector defined orientation, orthogonal to u, v
        u << 1,0,0; // vector defined orientation, orthogonal to w, v

        slen = 5; // screw length
        r0 = 0.5; // radius of the cylinder
        delta = 0.2; // how much does the screw extend out and subtract in since the phi function is between -1 to 1
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

        *(output) = Eigen_matrix_to_vectorized_scalar(implicitFunctionOutput);
        // for (int j=0;j<(*output).shape()[0];j++){
        //     std::cout << (*output)[j] << " ";
        // };
        // std::cout << "-------------------" << endl;
        // std::cout << endl;
        // std::cout << "-------------------" << endl;

    }

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        // void gradient(REAL ax, REAL ay, REAL az, 
        //               REAL delta, REAL phi0,
        //               REAL twist_rate,
        //               REAL uvwi00, REAL uvwi01, REAL uvwi02, REAL uvwi10, REAL uvwi11, REAL uvwi12,
        //               REAL wx, REAL wy, REAL wz,
        //               REAL x, REAL y, REAL z,
        //               REAL *f, REAL *g, REAL *h) {


        // Matrix<REAL, Dynamic, 3> x_eigen_matrix(x.shape()[0], 3);
        // x_eigen_matrix = vectorized_vect_to_Eigen_matrix(x);

        // Matrix<REAL, Dynamic, 1> inside_ness(x.shape()[0], 1);
        // inside_ness = inside_ness_only(this->A, this->w,this->UVW_inv, this->slen, this->r0,
        //                  this->delta, this->twist_rate, this->phi0, x_eigen_matrix);

        std::cout << "result---------------" << std::endl;

        for (int j=0;j<(*output).shape()[0];j++){
            // std::cout << "-----------------tiger debug -------------------" << std::endl;
            // std::cout << (*output)[j] << std::endl;

            if(j<=5){
                // std::cout << this->A << std::endl;
                // std::cout << this->delta << std::endl;
                // std::cout << this->phi0 << std::endl;
                // std::cout << this->twist_rate << std::endl;
                // std::cout << this->UVW << std::endl;
                // std::cout << this->w(1, 0) << std::endl;
                // std::cout << this->w(2, 0) << std::endl;
                // std::cout << this->w(2, 0) << std::endl;

                // std::cout << "x j" << std::endl;

                // std::cout << x[j][0]  << std::endl;
                // std::cout << x[j][1]  << std::endl;
                // std::cout << x[j][2]  << std::endl;

                // std::cout << "output j" << std::endl;
                // std::cout << (*output)[j] << std::endl;

                // gradient(this->A(0,0), this->A(1,0), this->A(2,0), 
                //          this->delta, this->phi0,
                //          this->twist_rate,
                //          this->UVW(0,0), this->UVW(0,1), this->UVW(0,2),
                //          this->UVW(1,0), this->UVW(1,1), this->UVW(1,2),
                //          this->w(0, 0), this->w(1, 0), this->w(2, 0),
                //          x[j][0], x[j][1], x[j][2],
                //          output, j, true
                //          );

            }

            gradient(this->A(0,0), this->A(1,0), this->A(2,0), 
                     this->delta, this->phi0,
                     this->twist_rate,
                     this->UVW(0,0), this->UVW(0,1), this->UVW(0,2),
                     this->UVW(1,0), this->UVW(1,1), this->UVW(1,2),
                     this->w(0, 0), this->w(1, 0), this->w(2, 0),
                     x[j][0], x[j][1], x[j][2],
                     output, j, false
                     );

            if(j <= 5) {
                std::cout << x[j] << std::endl;
                std::cout << (*output)[j] << std::endl;
            }   

            // cout << "------------------ result -------------------" << endl;
            // std::cout << x[j] << std::endl;
            // std::cout << (*output)[j] << std::endl;
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