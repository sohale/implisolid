#include "gtest/gtest.h"
#include "../subdivision.hpp"
#include "../mesh_algorithms.hpp"
#include "faces_test_tools.hpp"

#include "../configs.hpp"

using mp5_implicit::CONFIG_C;

vectorized_faces make_example_1234() {
    vectorized_faces_shape
        shape = vectorized_faces_shape{6, 3};

    vectorized_faces faces {shape};

    std::vector<std::vector<int>> f = {{1,2,4}, {3,2,4}, {1,4,5}, {2,3,1}, {7,1,5}, {0,1,2}};
    for(int i=0; i < f.size(); ++i) {
        for(int j=0; j < 3; ++j) {
            faces[i][j] = f[i][j];
        }
    }

    /*
    vectorized_faces
        faces //vectorized_faces{shape};
        //faces=
            {{1,2,4}, {3,2,4}, {1,4,5}, {2,3,1}, {7,1,5}, {0,1,2}};
    */
    return faces;
}


using mp5_implicit::encode_edge__sort;

TEST(Subdivision_1to2, a) {
    // subdivide_1to2();

    vectorized_faces facets = make_example_1234();
    std::set<edge_pair_type> edges_with_1_side = {
        encode_edge__sort(1, 2, CONFIG_C::edgecode_base),
        encode_edge__sort(2, 4, CONFIG_C::edgecode_base),
        encode_edge__sort(1, 4, CONFIG_C::edgecode_base),  // not there
        // Why is it converting 1,7 ?
        encode_edge__sort(99,100, CONFIG_C::edgecode_base)
    };
    // if edge e=(1,7) is in the map midpoint_map , but it is not in the map set edges_with_1_side, why does it replace it?

    std::map<edge_pair_type, vectorized_vect::index> midpoint_map;

    midpoint_map[encode_edge__sort(1, 2, CONFIG_C::edgecode_base)] = 9;
    midpoint_map[encode_edge__sort(2, 4, CONFIG_C::edgecode_base)] = 10;
    midpoint_map[encode_edge__sort(1, 7, CONFIG_C::edgecode_base)] = 12;
    midpoint_map[encode_edge__sort(99, 100, CONFIG_C::edgecode_base)] = 9;

    // for any i, midpoint_map[edges_with_1_side[i]] must exist
    bool careful_for_twosides=true;

    vectorized_faces result= subdivide_1to2(facets, edges_with_1_side, midpoint_map, careful_for_twosides);

    EXPECT_TRUE( 1 != 2 );

}

