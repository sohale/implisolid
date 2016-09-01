
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
        REAL& a, REAL& b, REAL& c, REAL& d) const {
        
        a = y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2;
        b = x1 * z3 - x1 * z2 + x2 * z1 - x2 * z3 - x3 * z1 + x3 * z2;
        c = x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2;
        d = x1 * y3 * z2 - x1 * y2 * z3 + x2 * y1 * z3 - x2 * y3 * z1 - x3 * y1 * z2 + x3 * y2 * z1;
    
    }

    vectorized_vect getPlanes() const {
        vectorized_vect planes{boost::extents[4][4]};

        // first plane, based on points: 2, 3, 4
        this->calculatePlaneCoefficients(
            this->p[1][0], this->p[1][1], this->p[1][2],
            this->p[2][0], this->p[2][1], this->p[2][2],
            this->p[3][0], this->p[3][1], this->p[3][2],
            planes[0][0], planes[0][1], planes[0][2], planes[0][3]
        );

        // second plane, based on points: 1, 3, 4
        this->calculatePlaneCoefficients(
            this->p[0][0], this->p[0][1], this->p[0][2],
            this->p[2][0], this->p[2][1], this->p[2][2],
            this->p[3][0], this->p[3][1], this->p[3][2],
            planes[1][0], planes[1][1], planes[1][2], planes[1][3]
        );

        // third plane, based on points: 1, 2, 4
        this->calculatePlaneCoefficients(
            this->p[0][0], this->p[0][1], this->p[0][2],
            this->p[1][0], this->p[1][1], this->p[1][2],
            this->p[3][0], this->p[3][1], this->p[3][2],
            planes[2][0], planes[2][1], planes[2][2], planes[2][3]
        );

        // fourth plane, based on points: 1, 2, 3
        this->calculatePlaneCoefficients(
            this->p[0][0], this->p[0][1], this->p[0][2],
            this->p[1][0], this->p[1][1], this->p[1][2],
            this->p[2][0], this->p[2][1], this->p[2][2],
            planes[3][0], planes[3][1], planes[3][2], planes[3][3]
        );

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

        return planes;
    }

public:

    tetrahedron(std::vector<boost::array<REAL,3>> points, REAL matrix12[12]) {
        // we need to apply matrix to all 4 points
        
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                this->p[i][j] = points[i][j];
            }
        }

        // apply matrix to points
        matrix_vector_product(matrix12, this->p);
    }

    virtual void eval_implicit(const vectorized_vect & x, vectorized_scalar * f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        vectorized_vect planes = getPlanes();

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

        vectorized_vect planes = getPlanes();

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

        bool integrity = true;

        int i;
        int j;

        REAL d;

        // check if the points are not too close to each other

        for (i = 0; i < 3 & integrity; i++) {
            for (j = i + 1; j < 4 & integrity; j++) {
                d = sqrt(
                    (this->p[i][0] - this->p[j][0]) * (this->p[i][0] - this->p[j][0]) +
                    (this->p[i][1] - this->p[j][1]) * (this->p[i][1] - this->p[j][1]) +
                    (this->p[i][2] - this->p[j][2]) * (this->p[i][2] - this->p[j][2])
                );

                if (d < MIN_PRINTABLE_LENGTH) {
                    loger << " Points are too close: " << d << std::endl;

                    loger << "( " << this->p[i][0] << ", " << this->p[i][1] << ", "  << this->p[i][1] << " )"<< std::endl;
                    loger << "( " << this->p[j][0] << ", " << this->p[j][1] << ", "  << this->p[j][1] << " )"<< std::endl;

                    integrity = false;
                }
            }
        }

        if (!integrity) 
            return integrity;

        // TODO: check if we can build plane through each 3 points

        // check if distance from each point to opposite face is not too small
        // because first plane is based at points 2, 3, 4, second at 1,3,4 etc
        // we use it, so index of plane, is the index of the point which is not at that plane

        // d = (a * x0 + b * y0 + c * z0 + d) / sqrt(a * a + b * b + c * c)

        planes = getPlanes();

        for (i = 0; i < 4 & integrity; i++) {
            d = planes[i][0] * this->p[i][0] + 
                planes[i][1] * this->p[i][1] + 
                planes[i][2] * this->p[i][2] + 
                planes[i][3];

            d /= sqrt(
                planes[i][0] * planes[i][0] +
                planes[i][1] * planes[i][1] +
                planes[i][2] * planes[i][2]
            );

            if (d < MIN_PRINTABLE_LENGTH) {
                loger << " Point " << i + 1 << " are too close to opposite plane :" << d << std::endl;

                integrity = false;
            }
        }

        return integrity;
    }

    virtual mp5_implicit::bounding_box get_boundingbox() const {

        REAL maxX = this->p[0][0];
        REAL minX = this->p[0][0];
        REAL maxY = this->p[0][1];
        REAL minY = this->p[0][1];
        REAL maxZ = this->p[0][2];
        REAL minZ = this->p[0][2];

        for (int i = 1; i < 4; i++) {
            maxX = this->p[i][0] > maxX ? this->p[i][0] : maxX;
            minX = this->p[i][0] < minX ? this->p[i][0] : minX;

            maxY = this->p[i][1] > maxY ? this->p[i][1] : maxY;
            minY = this->p[i][1] < minY ? this->p[i][1] : minY;

            maxZ = this->p[i][2] > maxZ ? this->p[i][2] : maxZ;
            minZ = this->p[i][2] < minZ ? this->p[i][2] : minZ;
        }

        return mp5_implicit::bounding_box{minX, maxX, minY, maxY, minZ, maxZ};
    }

};
}
