#pragma once
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"

namespace mp5_implicit {
namespace implicit_functions {


//  Nov  2 18:36 2016
//  Dec  23 9-10pm 2016

// To reproduce the object in the demo by Henrik Rydgård/@greggman/alteredq//
/* AlteredQualia's Marching Cubes, C++ version, based on AQ's code based on Henrik Rydgård and @greggman.
 * https://github.com/WebGLSamples/WebGLSamples.github.io/blob/master/blob/marching_cubes.js
 *
 * Based on alteredq's version  https://github.com/mrdoob/three.js/blob/master/examples/js/MarchingCubes.js
 *
 * Port of greggman's ThreeD version of marching cubes to Three.js
 * http://webglsamples.googlecode.com/hg/blob/blob.html
*/
class meta_ball_Rydgård : public transformable_implicit_function {
   /**
   The size is 1 x 1 x 1, and the boundaries of a metaballs object is [-0.5,+0.5]^3 but sometimes the metaballs go beyond this range.
   */

    int DUMP_VALUES = 0;  // used for debugging
    constexpr static bool DEBUG_VERBOSE_ = false;

 public:
    meta_ball_Rydgård (REAL matrix12[12], int num_blobs, REAL time, REAL scale)
    {
        std::cout << "Metaballs Constructor: (";
        std::cout << " blobs:" << num_blobs << ", time:" << time << ", scale:" << scale << ")  matrix: [";
        for (int i=0; i<12; i++) {
            if (i != 0)
                std::cout << ",";
            std::cout << matrix12[i];
        }
        std::cout << "]"  << std::endl;

        meta_ball_Rydgård::global_counter = 0;

        //this->balls.resize(0);
        //this->planes.resize(0);
        my_assert(this->balls.size() == 0, "needs to be empty to begin with");
        my_assert(this->planes.size() == 0, "needs to be empty to begin with");
        std::cout << "bp:" << this->balls.size() << this->planes.size() << std::endl;

        int numblobs = num_blobs;  // default: 4
        for (int ball_i = 0; ball_i < numblobs; ball_i++) {
            REAL x0 = 0.5, y0 = 0.5, z0 = 0.5;
            REAL D = 1;
            REAL ballx = sin(ball_i + 1.26 * time * (1.03 + 0.5*cos(0.21 * ball_i))) * 0.27 * D + 0.5   - x0;
            REAL bally = std::abs(cos(ball_i + 1.12 * time * cos(1.22 + 0.1424 * ball_i))) * 0.77 * D   - y0;  // dip into the floor
            REAL ballz = cos(ball_i + 1.32 * time * 0.1*sin((0.92 + 0.53 * ball_i))) * 0.27 * D + 0.5   - z0;
            REAL subtract = 12;
            REAL strength = 1.2 / ((sqrt(numblobs)- 1) / 4 + 1);
            this->addBall(ballx, bally, ballz, strength, subtract, scale);
        }

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            transf_matrix[i] = matrix12[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");


    }

    static int global_counter;

    class ball_t {
        REAL ballx, bally, ballz, strength, subtract, scale;
    public:
        ball_t (REAL ballx, REAL bally, REAL ballz, REAL strength, REAL subtract, REAL scale)
        : ballx(ballx), bally(bally), ballz(ballz), strength(strength), subtract(subtract), scale(scale)
        {
            if (meta_ball_Rydgård::DEBUG_VERBOSE_) {
                std::cout << "ball_t: Constructor:(" <<
                this->ballx << "," << bally << "," << ballz << "), strength=" << strength << ",subtract=" << subtract << ",scale=" << scale
                << std::endl;
            }
        }

        inline void implicit_ball(REAL realx, REAL realy, REAL realz, REAL & output) const {
            REAL fx = realx - this->ballx;
            REAL fx2 = fx * fx;
            REAL fy = realy - this->bally;
            REAL fy2 = fy * fy;
            REAL fz = realz - this->ballz;
            REAL fz2 = fz * fz;
            REAL h = ( (REAL)0.000001 + fx2 + fy2 + fz2 );
            REAL hinv = 1.0 / h;
            REAL val = strength * hinv - this->subtract;
            //if ( val > 0.0 )
                output += val / 100;

            if (meta_ball_Rydgård::DEBUG_VERBOSE_) {
                if( meta_ball_Rydgård::global_counter < 100 ){
                    std::cout << "implicit_ball: " <<
                        realx <<","<< realy <<","<< realz <<","<< output <<","<< val
                        << std::endl;
                }
                meta_ball_Rydgård::global_counter++;
            }
        }
        inline void gradient_ball(REAL realx, REAL realy, REAL realz, REAL & outputx, REAL & outputy, REAL & outputz) const {
            REAL fz = realz - this->ballz;
            REAL fz2 = fz * fz;
            REAL fy = realy - this->bally;
            REAL fy2 = fy * fy;
            REAL fx = realx - this->ballx;
            REAL fx2 = fx * fx;
            REAL h = 0.000001 + fx2 + fy2 + fz2;
            REAL hinv = 1.0 / h;
            REAL val = strength * hinv - this->subtract;
            //if ( val > 0.0 )
            {
                // output += val / 100;
                REAL gradx = this->strength * (
                    // y^-1  --> (-1) * (y^-2) * y'
                    // (0.001+xyz**2)^-1  --> (-1) * (((0.001+xyz**2))^-2) * (0.001+xyz**2)'
                    // --> -((0.001+xyz**2)^-2) * (2*xyz)
                    // --> -2 * x * (h^-2)
                    -2 * fx * hinv * hinv
                ) / 100;
                REAL grady = this->strength * ( -2 * fy * hinv * hinv ) / 100;
                REAL gradz = this->strength * ( -2 * fz * hinv * hinv ) / 100;

                outputx += gradx;
                outputy += grady;
                outputz += gradz;

                /*
                std::cout << "implicit_ball: " <<
                    realx <<","<< realy <<","<< realz <<": "<< outputx <<","<< outputy <<","<< outputz <<","<< val
                    << std::endl;
                */
            }
        }
        mp5_implicit::bounding_box  get_boundingbox() const {
            REAL radius_ = std::sqrt(this->strength / this->subtract);
            //REAL radius = VOXELS * radius_;
            REAL factor = REAL(1.0) / this->scale;
            REAL
                zs_ = this->ballz * factor,
                ys_ = this->bally * factor,
                xs_ = this->ballx * factor;
            /*REAL
                zs = zs_ * VOXELS,
                ys = ys_ * VOXELS,
                xs = xs_ * VOXELS;
            */

            return mp5_implicit::bounding_box {
                (xs_ - radius_), (xs_ + radius_), 
                (ys_ - radius_), (ys_ + radius_), 
                (zs_ - radius_), (zs_ + radius_)
            };

            /*
            int min_zi = floor( zs_ * VOXELS - VOXELS * radius_ ); if ( min_zi < 1 ) min_zi = 1;
            int max_zi = floor( zs_ * VOXELS + VOXELS * radius_ ); if ( max_zi > VOXELS - 1 ) max_zi = VOXELS - 1;
            int min_yi = floor( ys_ * VOXELS - VOXELS * radius_ ); if ( min_yi < 1 ) min_yi = 1;
            int max_yi = floor( ys_ * VOXELS + VOXELS * radius_ ); if ( max_yi > VOXELS - 1 ) max_yi = VOXELS - 1;
            int min_xi = floor( xs_ * VOXELS - VOXELS * radius_ ); if ( min_xi < 1  ) min_xi = 1;
            int max_xi = floor( xs_ * VOXELS + VOXELS * radius_ ); if ( max_xi > VOXELS - 1 ) max_xi = VOXELS - 1;

            return mp5_implicit::bounding_box{-****max_size, max_size, -max_size, max_size, -max_size, max_size};
            */

        }
    };

    class plane_t {
    
        REAL strength;
        REAL subtract;
    public:
        char xyz;
    public:
        plane_t (REAL strength, REAL subtract, char xyz) : strength(strength), subtract(subtract), xyz(xyz)
        {}

        inline void implicitX(REAL realx, REAL realy, REAL realz, REAL & output) const {
            REAL fx = realx - 0.5;
            REAL x2 = fx * fx;
            REAL val = this->strength / (REAL)( 0.0001 + x2 ) - this->subtract;
            if ( val > 0.0 ) {
                output += val;
            }
        }
        inline void gradientX(REAL realx, REAL realy, REAL realz, REAL & outputx, REAL & outputy, REAL & outputz) const {
            REAL fx = realx - 0.5;
            REAL x2 = fx * fx;
            REAL h = 0.0001 + x2;
            REAL hinv = static_cast<REAL>(1) / h;
            // REAL val = this->strength / h - this->subtract;
            REAL val = this->strength * hinv - this->subtract;
            if ( val > 0.0 ) {
                //output += val;
                REAL gradx = this->strength * (
                    // y^-1  --> (-1) * (y^-2) * y'
                    // (0.001+x**2)^-1  --> (-1) * (((0.001+x**2))^-2) * (0.001+x**2)'
                    // --> -((0.001+x**2)^-2) * (2*x)
                    // --> -2 * x * (h^-2)
                    -2 * fx * hinv * hinv
                );
                outputx += gradx;
            }
        }

        inline void implicitY(REAL realx, REAL realy, REAL realz, REAL & output) const {
            REAL fy = realy - 0.5;
            REAL y2 = fy * fy;
            REAL val = strength / (REAL)( 0.0001 + y2 ) - subtract;
            if ( val > 0.0 ) {
                output += val;
            }
        }
        inline void gradientY(REAL realx, REAL realy, REAL realz, REAL & outputx, REAL & outputy, REAL & outputz) const {
            REAL fy = realy - 0.5;
            REAL y2 = fy * fy;
            REAL h = 0.0001 + y2;
            REAL hinv = static_cast<REAL>(1) / h;
            REAL val = strength * hinv - subtract;
            if ( val > 0.0 ) {
                // output += val;
                REAL grady = this->strength * (-2) * fy * hinv * hinv;
                outputy += grady;
            }
        }

        inline void implicitZ(REAL realx, REAL realy, REAL realz, REAL & output) const {
            //realz = z / (REAL)VOXELS;
            REAL fz = realz - 0.5;
            REAL z2 = fz * fz;
            REAL val = strength / (REAL)( 0.0001 + z2 ) - subtract;
            if ( val > 0.0 ) {
                output += val;
            }
        }

        inline void gradientZ(REAL realx, REAL realy, REAL realz, REAL & outputx, REAL & outputy, REAL & outputz) const {
            //realz = z / (REAL)VOXELS;
            REAL fz = realz - 0.5;
            REAL z2 = fz * fz;
            REAL h = 0.0001 + z2;
            REAL hinv = static_cast<REAL>(1) / h;
            REAL val = this->strength * hinv - subtract;
            if ( val > 0.0 ) {
                // output += val;
                REAL gradz = this->strength * (-2) * fz * hinv * hinv;
                outputz += gradz;
            }
        }
    };




    std::vector<ball_t> balls;
    std::vector<plane_t> planes;

    void addBall( REAL ballx, REAL bally, REAL ballz, REAL strength, REAL subtract, REAL scale) {
        ball_t b {ballx, bally, ballz, strength, subtract, scale};
        balls.push_back(b);
    }
    void addPlaneX( REAL strength, REAL subtract ) {
        plane_t p {strength, subtract, 'x'};
        planes.push_back(p);
    }
    void addPlaneY( REAL strength, REAL subtract ) {
        plane_t p (strength, subtract, 'y');
        planes.push_back(p);
    }
    void addPlaneZ( REAL strength, REAL subtract ) {
        plane_t p (strength, subtract, 'z');
        planes.push_back(p);
    }







    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        {
            int output_ctr=0;
            for(auto i = x.begin(), e = x.end(); i < e; i++, output_ctr++){
                (*f_output)[output_ctr] = 0.0; //-0.00001;
            }
        }
        const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix, x);

        for( const ball_t & b : balls ) {
            int output_ctr=0;
            for(auto i = local_x.begin(), e = local_x.end(); i < e; i++, output_ctr++){
                const REAL x = (*i)[0];
                const REAL y = (*i)[1];
                const REAL z = (*i)[2];
                b.implicit_ball(x, y, z, (*f_output)[output_ctr]);
            }
        }

        for( const plane_t & p : planes ) {
            switch (p.xyz) {
            case 'x':{
                int output_ctr = 0;
                for(auto i = local_x.begin(), e = local_x.end(); i < e; i++, output_ctr++){
                    const REAL x = (*i)[0];
                    const REAL y = (*i)[1];
                    const REAL z = (*i)[2];
                    p.implicitX(x, y, z, (*f_output)[output_ctr]);
                }}
                break;
            case 'y':{
                int output_ctr = 0;
                for(auto i = local_x.begin(), e = local_x.end(); i < e; i++, output_ctr++){
                    const REAL x = (*i)[0];
                    const REAL y = (*i)[1];
                    const REAL z = (*i)[2];
                    p.implicitY(x, y, z, (*f_output)[output_ctr]);
                }}
                break;
            case 'z':{
                int output_ctr = 0;
                for(auto i = local_x.begin(), e = local_x.end(); i < e; i++, output_ctr++){
                    const REAL x = (*i)[0];
                    const REAL y = (*i)[1];
                    const REAL z = (*i)[2];
                    p.implicitZ(x, y, z, (*f_output)[output_ctr]);
                }}
                break;
            default:
                std::cout << "Error" << std::endl;
            }
        }

        if (meta_ball_Rydgård::DEBUG_VERBOSE_) {
            if (DUMP_VALUES > 0) {
                std::cout << "eval i " ;
                int output_ctr = 0;
                for(auto i = local_x.begin(), e = local_x.end(); i < e; i++, output_ctr++){
                    std::cout << (*f_output)[output_ctr] << " ";
                    if (output_ctr >= DUMP_VALUES) break;
                }
                std::cout << std::endl;
            }
        }
    }


    // boost::sub_array<float, 2U - 1>
    // boost::multi_array<REAL, 1>
    // static vectorized_vect sample;
    // todo: use this everywhere, move to transformable_implicit_function
    inline void transform_and_store_gradient( /*decltype(sample[0])  **/vectorized_vect::iterator g, const REAL gx, const REAL gy, const REAL gz) const {
        // M_inv.transpose dot grad
        (*g)[0] = this->inv_transf_matrix[0] * gx + this->inv_transf_matrix[4] * gy + this->inv_transf_matrix[ 8] * gz;
        (*g)[1] = this->inv_transf_matrix[1] * gx + this->inv_transf_matrix[5] * gy + this->inv_transf_matrix[ 9] * gz;
        (*g)[2] = this->inv_transf_matrix[2] * gx + this->inv_transf_matrix[6] * gy + this->inv_transf_matrix[10] * gz;
    }

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix, x);

        int output_ctr=0;
        //for(auto i = local_x.begin(), e = local_x.end(); i != e; i++, output_ctr++){
        auto oi = output->begin();
        for(auto i = local_x.begin(), e = local_x.end(); i != e; ++i, ++oi ) {
            const REAL x = (*i)[0];
            const REAL y = (*i)[1];
            const REAL z = (*i)[2];

            /*REAL gx = 0;
            REAL gy = 0;
            REAL gz = 0;
            b.implicit_ball(x, y, z, (*f_output)[output_ctr]);*/
            vectorized_vect g { boost::extents[1][3] };
            g[0][0] = 0;
            g[0][1] = 0;
            g[0][2] = 0;
            auto gp = g.begin();

            for( const ball_t & b : balls ) {
                    b.gradient_ball(x, y, z,  (*gp)[0], (*gp)[1], (*gp)[2] );
            }

            for( const plane_t & p : planes ) {
                switch (p.xyz) {
                case 'x':
                        p.gradientX(x, y, z,  (*gp)[0], (*gp)[1], (*gp)[2] );
                    break;
                case 'y':
                        p.gradientY(x, y, z,  (*gp)[0], (*gp)[1], (*gp)[2] );
                    break;
                case 'z':
                        p.gradientZ(x, y, z,  (*gp)[0], (*gp)[1], (*gp)[2] );
                    break;
                //default:
                    //std::cout << "Error" << std::endl;
                }
            }

            transform_and_store_gradient( /*&(((*output)[output_ctr]))*/ oi,  g[0][0], g[0][1], g[0][2] );
        }

        if (meta_ball_Rydgård::DEBUG_VERBOSE_) {
            if (DUMP_VALUES > 0) {
                std::cout << "eval g " ;
                int output_ctr = 0;
                for(auto i = (*output).begin(), e = (*output).end(); i < e; i++) {
                    std::cout << (*i)[0] <<","<< (*i)[1] <<","<< (*i)[2] << " ";
                    if (output_ctr >= DUMP_VALUES) break;
                    output_ctr++;
                }
                std::cout << std::endl;
            }
        }


    }
    bool integrity_invariant() const {
          return true;
    }

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        mp5_implicit::bounding_box result;
        bool first = true;
        for ( const ball_t & b : balls ) {
            mp5_implicit::bounding_box bb =  b.get_boundingbox();
            if (first)
                result = bb;
            else
                result.grow(bb);
            first = false;
        }
        my_assert(!first, "metaballs needs at least one ball added. zero balls.");
        return result;
    }
};




/*
// void MarchingCubes::addBall( REAL ballx, REAL bally, REAL ballz,  REAL strength, REAL subtract, REAL scale)

inline
void MarchingCubes::addBall(const ball_t & b)
{
    // Solves this equation:
    // 1.0 / (0.000001 + radius^2) * strength - subtract = 0
    auto VOXELS = this->resolution;

    REAL radius_ = std::sqrt(b.strength / b.subtract);
    REAL radius = VOXELS * radius_;
    REAL factor = REAL(1.0) / b.scale;
    REAL
        zs_ = b.ballz * factor,
        ys_ = b.bally * factor,
        xs_ = b.ballx * factor;
    REAL
        zs = zs_ * VOXELS,
        ys = ys_ * VOXELS,
        xs = xs_ * VOXELS;

        return mp5_implicit::bounding_box {
            xs_ * VOXELS - VOXELS * radius_, xs_ * VOXELS + VOXELS * radius_, 
            ys_ * VOXELS - VOXELS * radius_, ys_ * VOXELS + VOXELS * radius_, 
            zs_ * VOXELS - VOXELS * radius_, zs_ * VOXELS + VOXELS * radius_
        };


    int min_zi = floor( zs_ * VOXELS - VOXELS * radius_ ); if ( min_zi < 1 ) min_zi = 1;
    int max_zi = floor( zs_ * VOXELS + VOXELS * radius_ ); if ( max_zi > VOXELS - 1 ) max_zi = VOXELS - 1;
    int min_yi = floor( ys_ * VOXELS - VOXELS * radius_ ); if ( min_yi < 1 ) min_yi = 1;
    int max_yi = floor( ys_ * VOXELS + VOXELS * radius_ ); if ( max_yi > VOXELS - 1 ) max_yi = VOXELS - 1;
    int min_xi = floor( xs_ * VOXELS - VOXELS * radius_ ); if ( min_xi < 1  ) min_xi = 1;
    int max_xi = floor( xs_ * VOXELS + VOXELS * radius_ ); if ( max_xi > VOXELS - 1 ) max_xi = VOXELS - 1;


    // Don't polygonize_cube in the outer layer because normals aren't
    // well-defined there.

    REAL realx, realy, realz;
    //int x, y, z;
    //REAL val;
    //int y_offset, z_offset;

    for ( int z = min_zi; z < max_zi; z++ ) {
        for ( int y = min_yi; y < max_yi; y++ ) {
            for ( int x = min_xi; x < max_xi; x++ ) {

            int index = 0;

            // z_offset = this->size2 * z,
            //y_offset = z_offset + VOXELS * y;
            // int index = y_offset + x;

            REAL realx = x / (REAL)VOXELS;
            REAL realy = y / (REAL)VOXELS;
            REAL realz = z / (REAL)VOXELS;

            b.implicit(realx, realy, realz, this->field[ index ]);
        }
    }
}

void MarchingCubes::addPlaneX(REAL strength, REAL subtract ) {
    int cxy;

    // cache attribute lookups
    int yd = this->yd;
    auto VOXELS = this->resolution;
    int VOXELS = VOXELS;
    int zd = this->zd;
    array1d& field = this->field;
    REAL dist = VOXELS * sqrt(strength / (REAL)subtract);

    if ( dist > VOXELS ) dist = VOXELS;
    for ( int x = 0; x < dist; x++ ) {
        for ( int y = 0; y < VOXELS; y++ ) {
            cxy = x + y * yd;
            for ( int z = 0; z < VOXELS; z++ ) {
                REAL realx = x / (REAL)VOXELS;

                auto index = zd * z + cxy;
                b.implicit(realx, realy, realz, this->field[ index ]);
            }
        }
    }
}


void MarchingCubes::addPlaneY(REAL strength, REAL subtract ) {
    //int x, y, z;
    REAL y2;
    REAL val;
    REAL realy;
    int cy;
    int cxy;

    // cache attribute lookups
    auto VOXELS = this->resolution;
    int VOXELS = VOXELS;
    int yd = this->yd;
    int zd = this->zd;
    array1d& field = this->field;
    REAL dist = VOXELS * sqrt(strength / subtract);

    if ( dist > VOXELS ) dist = VOXELS;

    for (int y = 0; y < dist; y++ ) {
        REAL realy = y / (REAL)VOXELS;

            cy = y * yd;
            for (int x = 0; x < VOXELS; x++ ) {
                cxy = cy + x;
                for (int z = 0; z < VOXELS; z++ )
                    auto index = zd * z + cxy;

        REAL y2 = realy * realy;
        REAL val = strength / (REAL)( 0.0001 + y2 ) - subtract;
        if ( val > 0.0 ) {
                    field[ index ] += val;
            }
        }
    }
}

void MarchingCubes::addPlaneZ( REAL strength, REAL subtract )
{
    //int x, y, z;
    REAL z2, val, realz;
    int cz, cyz;

    // cache attribute lookups
    auto VOXELS = this->resolution;
    int VOXELS = VOXELS;
    int yd = this->yd;
    int zd = this->zd;
    array1d& field = this->field;
    REAL dist = VOXELS * sqrt( strength / subtract );

    if ( dist > VOXELS ) dist = VOXELS;
    for (int z = 0; z < dist; z++ ) {
            cz = zd * z;
            for (int y = 0; y < VOXELS; y++ ) {
                cyz = cz + y * yd;
                for (int x = 0; x < VOXELS; x++ )
                    auto index = cyz + x;
        realz = z / (REAL)VOXELS;
        z2 = realz * realz;
        val = strength / (REAL)( 0.0001 + z2 ) - subtract;
        if ( val > 0.0 ) {
                    field[ index ] += val;
            }
        }
    }
}
*/

}  // namespace implicit_functions
}  // namespace mp5_implicit
