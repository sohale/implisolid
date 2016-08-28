
/**
 * File:
 * 		tetrahedron.hpp
 * Description:
 * 		Defines class tetrahedron, an implicit object with its implicit function
 *   	and implicit gradient.
 */

#pragma once
#include "basic_data_structures.hpp"

/* part of the namespace mp5_implicit */
namespace mp5_implicit {

class tetrahedron : public transformable_implicit_function {

protected:
    // it's not a points, but a coefficients of planes
    // a*x + b*y + c*z + d = 0
    // [a, b, c, d]
    array2d p{boost::extents[4][4]};

    void calculatePlaneCoefficients(
        const REAL x1, const REAL y1, const REAL z1,
        const REAL x2, const REAL y2, const REAL z2,
        const REAL x3, const REAL y3, const REAL z3,
        REAL& a, REAL& b, REAL& c, REAL& d) {
        a = y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2;
        b = x1 * z3 - x1 * z2 + x2 * z1 - x2 * z3 - x3 * z1 + x3 * z2;
        c = x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2;
        d = x1 * y3 * z2 - x1 * y2 * z3 + x2 * y1 * z3 - x2 * y3 * z1 - x3 * y1 * z2 + x3 * y2 * z1;
    }

    void tetrahedronConstructor(std::vector<boost::array<REAL,3>> points, REAL matrix12[12]) {
        // we need to apply matrix to all 4 points
        loger << "We are in (real) tetrahedron constructor" << std::endl;

        vectorized_vect p_c{boost::extents[4][3]};

        p_c[0][0] = points[0][0];
        p_c[0][1] = points[0][1];
        p_c[0][2] = points[0][2];

        p_c[1][0] = points[1][0];
        p_c[1][1] = points[1][1];
        p_c[1][2] = points[1][2];

        p_c[2][0] = points[2][0];
        p_c[2][1] = points[2][1];
        p_c[2][2] = points[2][2];

        p_c[3][0] = points[3][0];
        p_c[3][1] = points[3][1];
        p_c[3][2] = points[3][2];

        loger << "Point 1" << std::endl;
        loger << "x " << p_c[0][0] << std::endl;
        loger << "y " << p_c[0][1] << std::endl;
        loger << "z " << p_c[0][2] << std::endl;

        loger << "Point 2" << std::endl;
        loger << "x " << p_c[1][0] << std::endl;
        loger << "y " << p_c[1][1] << std::endl;
        loger << "z " << p_c[1][2] << std::endl;

        loger << "Point 3" << std::endl;
        loger << "x " << p_c[2][0] << std::endl;
        loger << "y " << p_c[2][1] << std::endl;
        loger << "z " << p_c[2][2] << std::endl;

        loger << "Point 4" << std::endl;
        loger << "x " << p_c[3][0] << std::endl;
        loger << "y " << p_c[3][1] << std::endl;
        loger << "z " << p_c[3][2] << std::endl;


        //apply matrix to points
        matrix_vector_product(matrix12, p_c);

        loger << "Point after matrix applying" << std::endl;


        loger << "Point 1" << std::endl;
        loger << "x " << p_c[0][0] << std::endl;
        loger << "y " << p_c[0][1] << std::endl;
        loger << "z " << p_c[0][2] << std::endl;

        loger << "Point 2" << std::endl;
        loger << "x " << p_c[1][0] << std::endl;
        loger << "y " << p_c[1][1] << std::endl;
        loger << "z " << p_c[1][2] << std::endl;

        loger << "Point 3" << std::endl;
        loger << "x " << p_c[2][0] << std::endl;
        loger << "y " << p_c[2][1] << std::endl;
        loger << "z " << p_c[2][2] << std::endl;

        loger << "Point 4" << std::endl;
        loger << "x " << p_c[3][0] << std::endl;
        loger << "y " << p_c[3][1] << std::endl;
        loger << "z " << p_c[3][2] << std::endl;

        // calculate plane coefficients
        // first plane, based on points: 1, 2, 3
        this->calculatePlaneCoefficients(
        p_c[0][0], p_c[0][1], p_c[0][2],
        p_c[1][0], p_c[1][1], p_c[1][2],
        p_c[2][0], p_c[2][1], p_c[2][2],
        this->p[0][0], this->p[0][1], this->p[0][2], this->p[0][3]);

        // second plane, based on points: 2, 3, 4
        this->calculatePlaneCoefficients(
        p_c[1][0], p_c[1][1], p_c[1][2],
        p_c[2][0], p_c[2][1], p_c[2][2],
        p_c[3][0], p_c[3][1], p_c[3][2],
        this->p[1][0], this->p[1][1], this->p[1][2], this->p[1][3]);

        // third plane, based on points: 1, 3, 4
        this->calculatePlaneCoefficients(
        p_c[0][0], p_c[0][1], p_c[0][2],
        p_c[2][0], p_c[2][1], p_c[2][2],
        p_c[3][0], p_c[3][1], p_c[3][2],
        this->p[2][0], this->p[2][1], this->p[2][2], this->p[2][3]);

        // fourth plane, based on points: 1, 2, 4
        this->calculatePlaneCoefficients(
        p_c[0][0], p_c[0][1], p_c[0][2],
        p_c[1][0], p_c[1][1], p_c[1][2],
        p_c[3][0], p_c[3][1], p_c[3][2],
        this->p[3][0], this->p[3][1], this->p[3][2], this->p[3][3]);

        loger << "Plane 1" << std::endl;
        loger << "a "  << this->p[0][0] << std::endl;
        loger << "b "  << this->p[0][1] << std::endl;
        loger << "c "  << this->p[0][2] << std::endl;
        loger << "d "  << this->p[0][3] << std::endl;

        loger << "Plane 2" << std::endl;
        loger << "a "  << this->p[1][0] << std::endl;
        loger << "b "  << this->p[1][1] << std::endl;
        loger << "c "  << this->p[1][2] << std::endl;
        loger << "d "  << this->p[1][3] << std::endl;


        loger << "Plane 3" << std::endl;
        loger << "a "  << this->p[2][0] << std::endl;
        loger << "b "  << this->p[2][1] << std::endl;
        loger << "c "  << this->p[2][2] << std::endl;
        loger << "d "  << this->p[2][3] << std::endl;


        loger << "Plane 4" << std::endl;
        loger << "a "  << this->p[3][0] << std::endl;
        loger << "b "  << this->p[3][1] << std::endl;
        loger << "c "  << this->p[3][2] << std::endl;
        loger << "d "  << this->p[3][3] << std::endl;


        // calculate sign of each plane

        REAL plane_1_sign = sign(
        this->p[0][0] * p_c[3][0] +
        this->p[0][1] * p_c[3][1] +
        this->p[0][2] * p_c[3][2] +
        this->p[0][3],
        ROOT_TOLERANCE
        );

        REAL plane_2_sign = sign(
        this->p[1][0] * p_c[0][0] +
        this->p[1][1] * p_c[0][1] +
        this->p[1][2] * p_c[0][2] +
        this->p[1][3],
        ROOT_TOLERANCE
        );

        REAL plane_3_sign = sign(
        this->p[2][0] * p_c[1][0] +
        this->p[2][1] * p_c[1][1] +
        this->p[2][2] * p_c[1][2] +
        this->p[2][3],
        ROOT_TOLERANCE
        );

        REAL plane_4_sign = sign(
        this->p[3][0] * p_c[2][0] +
        this->p[3][1] * p_c[2][1] +
        this->p[3][2] * p_c[2][2] +
        this->p[3][3],
        ROOT_TOLERANCE
        );

        loger << "Plane 1 sign:" << std::endl;
        loger << "sign "  << plane_1_sign << std::endl;

        loger << "Plane 2 sign:" << std::endl;
        loger << "sign "  << plane_2_sign << std::endl;

        loger << "Plane 3 sign:" << std::endl;
        loger << "sign "  << plane_3_sign << std::endl;

        loger << "Plane 4 sign:" << std::endl;
        loger << "sign "  << plane_4_sign << std::endl;


        // apply sign to each plane

        this->p[0][0] = plane_1_sign * this->p[0][0];
        this->p[0][1] = plane_1_sign * this->p[0][1];
        this->p[0][2] = plane_1_sign * this->p[0][2];
        this->p[0][3] = plane_1_sign * this->p[0][3];

        this->p[1][0] = plane_2_sign * this->p[1][0];
        this->p[1][1] = plane_2_sign * this->p[1][1];
        this->p[1][2] = plane_2_sign * this->p[1][2];
        this->p[1][3] = plane_2_sign * this->p[1][3];

        this->p[2][0] = plane_3_sign * this->p[2][0];
        this->p[2][1] = plane_3_sign * this->p[2][1];
        this->p[2][2] = plane_3_sign * this->p[2][2];
        this->p[2][3] = plane_3_sign * this->p[2][3];

        this->p[3][0] = plane_4_sign * this->p[3][0];
        this->p[3][1] = plane_4_sign * this->p[3][1];
        this->p[3][2] = plane_4_sign * this->p[3][2];
        this->p[3][3] = plane_4_sign * this->p[3][3];


        loger << "Plane after applying sign" << std::endl;

        loger << "Plane 1" << std::endl;
        loger << "a "  << this->p[0][0] << std::endl;
        loger << "b "  << this->p[0][1] << std::endl;
        loger << "c "  << this->p[0][2] << std::endl;
        loger << "d "  << this->p[0][3] << std::endl;

        loger << "Plane 2" << std::endl;
        loger << "a "  << this->p[1][0] << std::endl;
        loger << "b "  << this->p[1][1] << std::endl;
        loger << "c "  << this->p[1][2] << std::endl;
        loger << "d "  << this->p[1][3] << std::endl;


        loger << "Plane 3" << std::endl;
        loger << "a "  << this->p[2][0] << std::endl;
        loger << "b "  << this->p[2][1] << std::endl;
        loger << "c "  << this->p[2][2] << std::endl;
        loger << "d "  << this->p[2][3] << std::endl;


        loger << "Plane 4" << std::endl;
        loger << "a "  << this->p[3][0] << std::endl;
        loger << "b "  << this->p[3][1] << std::endl;
        loger << "c "  << this->p[3][2] << std::endl;
        loger << "d "  << this->p[3][3] << std::endl;


    }

public:

    tetrahedron(std::vector<boost::array<REAL,3>> points, REAL matrix12[12]) {
        this->tetrahedronConstructor(points, matrix12);

    }

    //constructor for testing
    tetrahedron(REAL matrix12[12]) {
        std::vector<boost::array<REAL,3>> points;
        loger << "We are in tetrahedron constructor" << std::endl;


        boost::array<REAL,3> p1;
        boost::array<REAL,3> p2;
        boost::array<REAL,3> p3;
        boost::array<REAL,3> p4;

        p1[0] = 0.; p1[1] = 0.; p1[2] = 10.;
        p2[0] = 0.; p2[1] = 0.; p2[2] = 0.;
        p3[0] = 0.; p3[1] = 10.; p3[2] = 0.;
        p4[0] = 10.; p4[1] = 0.; p4[2] = 0.;

        points.push_back(p1);
        points.push_back(p2);
        points.push_back(p3);
        points.push_back(p4);

        this->tetrahedronConstructor(points, matrix12);
        
    }

    virtual void eval_implicit(const vectorized_vect & x, vectorized_scalar * f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        int output_ctr = 0;

        auto i = x.begin();
        auto e = x.end();

        REAL i1;
        REAL i2;
        REAL i3;

        for(; i<e; i++, output_ctr++) {
          i1 = (*i)[0];
          i2 = (*i)[1];
          i3 = (*i)[2];
          (*f_output)[output_ctr] = min(
            this->p[0][0] * i1 + this->p[0][1] * i2 + this->p[0][2] * i3 + this->p[0][3],
            min(this->p[1][0] * i1 + this->p[1][1] * i2 + this->p[1][2] * i3 + this->p[1][3],
            min(this->p[2][0] * i1 + this->p[2][1] * i2 + this->p[2][2] * i3 + this->p[2][3],
                this->p[3][0] * i1 + this->p[3][1] * i2 + this->p[3][2] * i3 + this->p[3][3]
          )));
        }
    }

    virtual void eval_gradient(const vectorized_vect & x, vectorized_vect * output) const {

        int output_ctr = 0;

        auto i = x.begin();
        auto e = x.end();

        REAL i1;
        REAL i2;
        REAL i3;

        for(; i < e; i++, output_ctr++) {
        i1 = (*i)[0];
        i2 = (*i)[1];
        i3 = (*i)[2];
        int index = 0;

        REAL min = this->p[index][0] * i1 + this->p[index][1] * i2 + this->p[index][2] * i3 + this->p[index][3];
        REAL possible_min = 0;
        for (int i=1; i<4; i++) {
            possible_min = this->p[i][0] * i1 + this->p[i][1] * i2 + this->p[i][2] * i3 + this->p[i][3];
          if (possible_min < min){
            index = i;
            min = possible_min;
          }
        }

        (*output)[output_ctr][0] = this->p[index][0];
        (*output)[output_ctr][1] = this->p[index][1];
        (*output)[output_ctr][2] = this->p[index][2];

      }

    }

    bool integrity_invariant() const {
        /// this function should make sure all of p points are separate and have at least  a minimum distance.
        /*
            Note: currently not implemented
        */
        return true;
    }

    virtual mp5_implicit::bounding_box get_boundingbox() const {
        // testing, not implemented yet
        return mp5_implicit::bounding_box{-4, 4, -4, 4, -4,4};
    }

};
}
