
#pragma once

#include <cmath>

/*,
================================================================
=           Configuration Parameters                           =
================================================================
*/

namespace mp5_implicit {
namespace vectorised_algorithms {
    //constants: (no change)

    /*
    ALL_CLOSE_EPS:
    Used for check if the norm is 1.0
    */
    constexpr REAL ALL_CLOSE_EPS = 0.000001;  // used for comparing the
}
}

/**
 * Configuration parameters
 * --------------------------
 * Contains global configuration options for tolerances, epsilon and other
 * parameters
 *
 * See: implicit_config.py
 */

const REAL ROOT_TOLERANCE = 0.001;  //
const REAL MIN_PRINTABLE_LENGTH = 0.01 / 2.0;

namespace mp5_implicit {

constexpr static const REAL MILLIMETER = 1.0;
constexpr static const REAL MICROMETER = MILLIMETER / 1000.0;
constexpr static const REAL NANOMETER = MICROMETER / 1000.0;

struct CONFIG_C {


    const REAL MIN_PRINTABLE_LENGTH = 0.01 / 2.0;
    const REAL MIN_THICKNESS = 0.2;  // for FDM

    constexpr static REAL MIN_AREA = (30 * NANOMETER)*(30 * NANOMETER); //
        // std::pow(30 * NANOMETER, 2); //0.000000001 == (30* NANOMETER)^2;    // 0.00001;
    //static const REAL MIN_AREA = 0.00001;  // 0.000000001;
    /* Also used for facet normal vector being zero */
    constexpr static REAL MIN_NORMAL_LEN = 0.1 * NANOMETER; //0.0000001; //  # 0.000001 = I millions of millimeter = 1 nanometer  #TH_N


    struct center_projection {
        // Gradients smaller than this are considered zero.
        constexpr static REAL min_gradient_len = 0.000001;
    };

} CONFIG;

}

/*
//namespace mp5_implicit {
struct CONFIG {

    struct polygonizer{
        static REAL ROOT_TOLERANCE = 0.001;  //
        static REAL MIN_PRINTABLE_LENGTH = 0.01;  //

        struct center_projection{
            //set_centers_on_surface
            static REAL min_gradient_len = 0.000001;  // Gradients smaller than this are considered zero.
            static int max_iter = 20;
            static bool USE_MESH_NORMALS = true;
            static bool EXTREME_ALPHA = false;  // Whether use alpha that exceeeds max_dim
        };
    };
};
//}
*/


constexpr bool VERBOSE_QEM = false;

// constexpr bool STORE_POINTSETS = true;
// Now: see pointset_set.hpp

