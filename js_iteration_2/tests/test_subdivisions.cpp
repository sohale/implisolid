#include "gtest/gtest.h"
#include "../subdivision/subdiv_1to2.hpp"
#include "../subdivision/subdiv_1to4.hpp"
#include "../mesh_algorithms.hpp"
#include "faces_test_tools.hpp"
#include "../configs.hpp"

#include "../v2v_f2f.hpp"

using mp5_implicit::CONFIG_C;
using mp5_implicit::easy_edge;

// using mp5_implicit::subdivide_multiple_facets_1to4;


vectorized_faces make_example_1234() {

    std::vector<std::vector<vertexindex_type>> f = {{1,2,4}, {3,2,4}, {1,4,5}, {2,3,1}, {7,1,5}, {0,1,3}};
    vectorized_faces faces = f2f(f);

    return faces;
}

auto testcase_square() {

    auto vv = std::vector<std::vector<REAL>>{
        {0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}
    };
    vectorized_vect v = v2v(vv, 1.0 * 10.0 / 4.0 );

    vectorized_faces f = f2f(std::vector<std::vector<vertexindex_type>>{
        {0, 1, 2}, {0, 2, 3}
    });

    return std::pair<vectorized_vect, vectorized_faces>(v,f);
}


using mp5_implicit::encode_edge__sort;

TEST(Subdivision_1to2, a) {
    // subdivide_1to2();

    vectorized_faces facets = make_example_1234();
    std::set<edge_pair_type> edges_with_1_side = {
        encode_edge__sort(1, 2, CONFIG_C::edgecode_base),
        encode_edge__sort(3, 4, CONFIG_C::edgecode_base),
        //encode_edge__sort(1, 4, CONFIG_C::edgecode_base),  // the mappnig is not there
        encode_edge__sort(1, 7, CONFIG_C::edgecode_base),
        // Why is it converting 1,7 ?
        encode_edge__sort(99,100, CONFIG_C::edgecode_base)
    };
    // if edge e=(1,7)   is in the map midpoint_map , but it is not in the map set edges_with_1_side, why does it replace it?

    std::map<edge_pair_type, vectorized_vect::index> midpoint_map;

    midpoint_map[encode_edge__sort(1, 2, CONFIG_C::edgecode_base)] = 9;
    midpoint_map[encode_edge__sort(2, 4, CONFIG_C::edgecode_base)] = 10;
    midpoint_map[encode_edge__sort(1, 7, CONFIG_C::edgecode_base)] = 12;
    midpoint_map[encode_edge__sort(2, 3, CONFIG_C::edgecode_base)] = 999;  // not used
    midpoint_map[encode_edge__sort(99, 100, CONFIG_C::edgecode_base)] = 9;
    midpoint_map[encode_edge__sort(4, 3, CONFIG_C::edgecode_base)] = 13;

    // for any i, midpoint_map[edges_with_1_side[i]] must exist
    bool careful_for_twosides=true;

    vectorized_faces result = subdivide_1to2(facets, edges_with_1_side, midpoint_map, careful_for_twosides);

    EXPECT_TRUE( 1 != 2 );

}

TEST(Subdivision_1to2, square) {

    auto vf = testcase_square();
    auto faces = vf.second;
    // {0, 1, 2}, {0, 2, 3}


    std::set<edge_pair_type> edges_to_subdivide = {
        encode_edge__sort(1, 2, CONFIG_C::edgecode_base),
    };

    std::map<edge_pair_type, vectorized_vect::index> midpoint_map;
    midpoint_map[encode_edge__sort(1, 2, CONFIG_C::edgecode_base)] = 99;

    bool careful_for_twosides=true;

    vectorized_faces result = subdivide_1to2(faces, edges_to_subdivide, midpoint_map, careful_for_twosides);
}


TEST(Subdivision_1to4, square) {

    auto vf = testcase_square();
    auto faces = vf.second;
    auto verts = vf.first;

    /*
    std::set<edge_pair_type> edges_to_subdivide = {
        easy_edge(1, 2),
    };
    */
    std::set<faceindex_type> triangles_to_subdivide {0};

    std::map<edge_pair_type, vectorized_vect::index> midpoint_map;

    std::cout << "a" << std::endl;

    midpoint_map[easy_edge(1, 2)] = 99;

    std::cout << "b" << std::endl;

    bool careful_for_twosides=true;

    auto result = subdivide_multiple_facets_1to4 (
        faces, verts, triangles_to_subdivide, midpoint_map);
}
