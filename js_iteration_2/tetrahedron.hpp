
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
    // Array of points. Matrix is already used.
    vectorized_vect p{boost::extents[4][3]};

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

    void savePoints(std::vector<boost::array<REAL,3>> points, REAL matrix12[12]) {
        // we need to apply matrix to all 4 points
        loger << "We are in (real) tetrahedron constructor" << std::endl;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                this->p[i][j] = points[i][j];
            }
        }

        loger << "Point 1" << std::endl;
        loger << "x " << this->p[0][0] << std::endl;
        loger << "y " << this->p[0][1] << std::endl;
        loger << "z " << this->p[0][2] << std::endl;

        loger << "Point 2" << std::endl;
        loger << "x " << this->p[1][0] << std::endl;
        loger << "y " << this->p[1][1] << std::endl;
        loger << "z " << this->p[1][2] << std::endl;

        loger << "Point 3" << std::endl;
        loger << "x " << this->p[2][0] << std::endl;
        loger << "y " << this->p[2][1] << std::endl;
        loger << "z " << this->p[2][2] << std::endl;

        loger << "Point 4" << std::endl;
        loger << "x " << this->p[3][0] << std::endl;
        loger << "y " << this->p[3][1] << std::endl;
        loger << "z " << this->p[3][2] << std::endl;


        // apply matrix to points
        matrix_vector_product(matrix12, this->points);

        loger << "Point after matrix applying" << std::endl;

        loger << "Point 1" << std::endl;
        loger << "x " << this->p[0][0] << std::endl;
        loger << "y " << this->p[0][1] << std::endl;
        loger << "z " << this->p[0][2] << std::endl;

        loger << "Point 2" << std::endl;
        loger << "x " << this->p[1][0] << std::endl;
        loger << "y " << this->p[1][1] << std::endl;
        loger << "z " << this->p[1][2] << std::endl;

        loger << "Point 3" << std::endl;
        loger << "x " << this->p[2][0] << std::endl;
        loger << "y " << this->p[2][1] << std::endl;
        loger << "z " << this->p[2][2] << std::endl;

        loger << "Point 4" << std::endl;
        loger << "x " << this->p[3][0] << std::endl;
        loger << "y " << this->p[3][1] << std::endl;
        loger << "z " << this->p[3][2] << std::endl;  
    }

    vectorized_vect getPlanes() {
        vectorized_vect planes{boost::extents[4][4]};

        // first plane, based on points: 1, 2, 3
        this->calculatePlaneCoefficients(
        this->p[0][0], this->p[0][1], this->p[0][2],
        this->p[1][0], this->p[1][1], this->p[1][2],
        this->p[2][0], this->p[2][1], this->p[2][2],
        planes[0][0], planes[0][1], planes[0][2], planes[0][3]);

        // second plane, based on points: 2, 3, 4
        this->calculatePlaneCoefficients(
        this->p[1][0], this->p[1][1], this->p[1][2],
        this->p[2][0], this->p[2][1], this->p[2][2],
        this->p[3][0], this->p[3][1], this->p[3][2],
        planes[1][0], planes[1][1], planes[1][2], planes[1][3]);

        // third plane, based on points: 1, 3, 4
        this->calculatePlaneCoefficients(
        this->p[0][0], this->p[0][1], this->p[0][2],
        this->p[2][0], this->p[2][1], this->p[2][2],
        this->p[3][0], this->p[3][1], this->p[3][2],
        planes[2][0], planes[2][1], planes[2][2], planes[2][3]);

        // fourth plane, based on points: 1, 2, 4
        this->calculatePlaneCoefficients(
        this->p[0][0], this->p[0][1], this->p[0][2],
        this->p[1][0], this->p[1][1], this->p[1][2],
        this->p[3][0], this->p[3][1], this->p[3][2],
        planes[3][0], planes[3][1], planes[3][2], planes[3][3]);

        loger << "Plane 1" << std::endl;
        loger << "a "  << planes[0][0] << std::endl;
        loger << "b "  << planes[0][1] << std::endl;
        loger << "c "  << planes[0][2] << std::endl;
        loger << "d "  << planes[0][3] << std::endl;

        loger << "Plane 2" << std::endl;
        loger << "a "  << planes[1][0] << std::endl;
        loger << "b "  << planes[1][1] << std::endl;
        loger << "c "  << planes[1][2] << std::endl;
        loger << "d "  << planes[1][3] << std::endl;


        loger << "Plane 3" << std::endl;
        loger << "a "  << planes[2][0] << std::endl;
        loger << "b "  << planes[2][1] << std::endl;
        loger << "c "  << planes[2][2] << std::endl;
        loger << "d "  << planes[2][3] << std::endl;


        loger << "Plane 4" << std::endl;
        loger << "a "  << planes[3][0] << std::endl;
        loger << "b "  << planes[3][1] << std::endl;
        loger << "c "  << planes[3][2] << std::endl;
        loger << "d "  << planes[3][3] << std::endl;
        

        // calculate sign of each plane

        REAL plane_1_sign = sign(
        planes[0][0] * this->p[3][0] +
        planes[0][1] * this->p[3][1] +
        planes[0][2] * this->p[3][2] +
        planes[0][3],
        ROOT_TOLERANCE
        );

        REAL plane_2_sign = sign(
        planes[1][0] * this->p[0][0] +
        planes[1][1] * this->p[0][1] +
        planes[1][2] * this->p[0][2] +
        planes[1][3],
        ROOT_TOLERANCE
        );

        REAL plane_3_sign = sign(
        planes[2][0] * this->p[1][0] +
        planes[2][1] * this->p[1][1] +
        planes[2][2] * this->p[1][2] +
        planes[2][3],
        ROOT_TOLERANCE
        );

        REAL plane_4_sign = sign(
        planes[3][0] * this->p[2][0] +
        planes[3][1] * this->p[2][1] +
        planes[3][2] * this->p[2][2] +
        planes[3][3],
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

        planes[0][0] *= plane_1_sign;
        planes[0][1] *= plane_1_sign;
        planes[0][2] *= plane_1_sign;
        planes[0][3] *= plane_1_sign;

        planes[1][0] *= plane_2_sign;
        planes[1][1] *= plane_2_sign;
        planes[1][2] *= plane_2_sign;
        planes[1][3] *= plane_2_sign;

        planes[2][0] *= plane_3_sign;
        planes[2][1] *= plane_3_sign;
        planes[2][2] *= plane_3_sign;
        planes[2][3] *= plane_3_sign;

        planes[3][0] *= plane_4_sign;
        planes[3][1] *= plane_4_sign;
        planes[3][2] *= plane_4_sign;
        planes[3][3] *= plane_4_sign;


        loger << "Plane after applying sign" << std::endl;

        loger << "Plane 1" << std::endl;
        loger << "a "  << planes[0][0] << std::endl;
        loger << "b "  << planes[0][1] << std::endl;
        loger << "c "  << planes[0][2] << std::endl;
        loger << "d "  << planes[0][3] << std::endl;

        loger << "Plane 2" << std::endl;
        loger << "a "  << planes[1][0] << std::endl;
        loger << "b "  << planes[1][1] << std::endl;
        loger << "c "  << planes[1][2] << std::endl;
        loger << "d "  << planes[1][3] << std::endl;


        loger << "Plane 3" << std::endl;
        loger << "a "  << planes[2][0] << std::endl;
        loger << "b "  << planes[2][1] << std::endl;
        loger << "c "  << planes[2][2] << std::endl;
        loger << "d "  << planes[2][3] << std::endl;


        loger << "Plane 4" << std::endl;
        loger << "a "  << planes[3][0] << std::endl;
        loger << "b "  << planes[3][1] << std::endl;
        loger << "c "  << planes[3][2] << std::endl;
        loger << "d "  << planes[3][3] << std::endl;

        return planes;
    }

public:

    tetrahedron(std::vector<boost::array<REAL,3>> points, REAL matrix12[12]) {
        this->savePoints(points, matrix12);
    }

    //constructor for testing
    tetrahedron(REAL matrix12[12]) {
        std::vector<boost::array<REAL,3>> points;
        
        loger << "We are in tetrahedron constructor(default)" << std::endl;

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

        this->savePoints(points, matrix12);
    }

    virtual void eval_implicit(const vectorized_vect & x, vectorized_scalar * f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        vectorized_vect planes = this->getPlanes();

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
                planes[0][0] * i1 + planes[0][1] * i2 + planes[0][2] * i3 + planes[0][3],
            min(planes[1][0] * i1 + planes[1][1] * i2 + planes[1][2] * i3 + planes[1][3],
            min(planes[2][0] * i1 + planes[2][1] * i2 + planes[2][2] * i3 + planes[2][3],
                planes[3][0] * i1 + planes[3][1] * i2 + planes[3][2] * i3 + planes[3][3]
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

        vectorized_vect planes = this->getPlanes();

        for(; i < e; i++, output_ctr++) {
        i1 = (*i)[0];
        i2 = (*i)[1];
        i3 = (*i)[2];
        int index = 0;

        REAL min = planes[index][0] * i1 + planes[index][1] * i2 + planes[index][2] * i3 + planes[index][3];
        REAL possible_min = 0;
        for (int i=1; i<4; i++) {
            possible_min = planes[i][0] * i1 + planes[i][1] * i2 + planes[i][2] * i3 + planes[i][3];
          if (possible_min < min){
            index = i;
            min = possible_min;
          }
        }

        (*output)[output_ctr][0] = planes[index][0];
        (*output)[output_ctr][1] = planes[index][1];
        (*output)[output_ctr][2] = planes[index][2];

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
