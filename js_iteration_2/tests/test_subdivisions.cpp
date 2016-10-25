#include "gtest/gtest.h"
#include "../subdivision/subdiv_1to2.hpp"
#include "../subdivision/subdiv_1to4.hpp"
#include "../mesh_algorithms.hpp"
#include "faces_test_tools.hpp"
#include "../configs.hpp"

#include "../v2v_f2f.hpp"

using mp5_implicit::CONFIG_C;
using mp5_implicit::easy_edge;
using mp5_implicit::subdivision::midpointmap_type;

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

auto testcase_triangle() {
    auto vv = std::vector<std::vector<REAL>>{
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}
    };
    vectorized_vect v = v2v(vv, 1.0 );   // 2.0

    vectorized_faces f = f2f(std::vector<std::vector<vertexindex_type>>{
        {0, 1, 2}
    });

    return std::pair<vectorized_vect, vectorized_faces>(v,f);
}

using mp5_implicit::encode_edge__sort;

using mp5_implicit::subdivision::subdivide_1to2;

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

    midpointmap_type  midpoint_map;

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

    midpointmap_type midpoint_map;
    midpoint_map[encode_edge__sort(1, 2, CONFIG_C::edgecode_base)] = 99;

    bool careful_for_twosides=true;

    vectorized_faces result = subdivide_1to2(faces, edges_to_subdivide, midpoint_map, careful_for_twosides);
}





/**
 * @brief      {Prints a faces array in the columnar format}
 * @param[in]  mark   index of the item that is highlighted. Useful for separating the last index of the original faces. Use -1 to avoid this demarcation.
 */
void print_faces (const vectorized_faces& faces, int mark) {
    cout << "#faces: " << faces.shape()[0] << std::endl;
    for (int i = 0; i < faces.shape()[0]; ++i ) {
        cout << i <<": ";
        for (int j=0; j < 3; ++j) {
            cout << " " << faces[i][j];
        }
        if (i == mark )
            cout << " <--- ";
        cout << std::endl;
    }
}

bool check_verts_faces (const vectorized_faces& faces, const vectorized_vect& verts) {
    bool ok = true;
    for (int i = 0; i < faces.shape()[0]; ++i ) {
        for (int j=0; j < 3; ++j) {
            ok = ok && faces[i][j] < verts.shape()[0];
        }
    }
    return ok;
}


template <typename T>
std::tuple<T,T,T>  sort_triple(T a, T b, T c) {
    // boost::multi_array<T, 1> triangle_face_
    typename boost::multi_array<T, 1> triangle_face // = triangle_face_;
            {boost::extents[3]};
    triangle_face[0] = a;
    triangle_face[1] = b;
    triangle_face[2] = c;
    typedef typename boost::multi_array<T, 1>::index idxt;
    idxt nd = triangle_face.shape()[0];
    assert(nd==3);
    // cout << "SORT.1" << std::endl;
    for (idxt i = 0; i < nd-1; ++i) {
        // cout << "SORT.i:" << i << std::endl;
        for (idxt j = 0; j <= i; ++j) {
            // cout << "SORT.j:" << j <<  " >?? " << j+1 << "  values: " << triangle_face[j] << " ? " << triangle_face[j+1] <<  std::endl;

            if (triangle_face[j] > triangle_face[j+1]) {
                // cout << "swapping " << j << "<->" << j+1 << std::endl;
                std::swap(triangle_face[j], triangle_face[j+1]);
            }
        }
    }
    for (idxt i = 0; i < nd-1; ++i) {
        assert (triangle_face[i] <= triangle_face[i+1]);
    }
    assert (nd == 3);
    typename std::tuple<T,T,T> s;
    std::get<0>(s) = triangle_face[0];
    std::get<1>(s) = triangle_face[1];
    std::get<2>(s) = triangle_face[2];
    return s;
}

/*
template <typename T>
std::tuple<T,T,T>  sort_triple(boost::multi_array<T, 1> triangle_face_) {
    typename boost::multi_array<T, 1> triangle_face = triangle_face_;
    typedef typename boost::multi_array<T, 1>::index idxt;
    idxt nd = triangle_face.shape()[0];
    for (idxt i = 0; i < nd; ++i) {
        for (idxt j = 0; j < i; ++j) {
            if (triangle_face[j] > triangle_face[i]) {
                std::swap(triangle_face[i], triangle_face[j]);
            }
        }
    }
    assert (nd == 3);
    typename std::tuple<T,T,T> s;
    std::get<0>(s) = triangle_face[0];
    std::get<1>(s) = triangle_face[1];
    std::get<2>(s) = triangle_face[2];
}
*/

// not used
vectorized_faces sort_faces(const vectorized_faces& faces) {
    auto nf = faces.shape()[0];
    vectorized_faces  faces_idx {boost::extents[nf][3]};
    for (int i = 0; i < faces.shape()[0]; ++i ) {
        auto t3 = sort_triple(faces[i][0], faces[i][1], faces[i][2]);
        // auto t3 = sort_triple(faces[i]);
        faces_idx[i][0] = std::get<0>(t3);
        faces_idx[i][1] = std::get<1>(t3);
        faces_idx[i][2] = std::get<2>(t3);
    }
    return faces_idx;
}

boost::multi_array<long, 1> faces_indices(const vectorized_faces& faces) {
    const long B = 10000;
    auto nf = faces.shape()[0];
    // cout << "===1.1" << std::endl;
    boost::multi_array<long, 1>  faces_idx {boost::extents[nf]};
    // cout << "===1.2" << std::endl;
    for (int i = 0; i < faces.shape()[0]; ++i ) {
        // cout << "===1.3" << std::endl;
        auto t3 = sort_triple(faces[i][0], faces[i][1], faces[i][2]);
        // cout << "===1.3.a" << std::endl;
        long l0 = std::get<0>(t3);
        long l1 = std::get<1>(t3);
        long l2 = std::get<2>(t3);
        faces_idx[i] = l0 + l1*B + l2 * B*B;
    }
    //  cout << "===1.4" << std::endl;
    assert(faces_idx.shape()[0] == faces.shape()[0]);
    return faces_idx;
}

bool check_faces_equality (const vectorized_faces& faces1, const vectorized_faces& faces2 ) {
    const bool verbose_about_mismatches = false;
    // auto cfaces = sort_faces(faces);
    //cout << "-.1==" << std::endl;
    auto faces_indices1 = faces_indices(faces1);
    //cout << "-.2==" << std::endl;
    auto faces_indices2 = faces_indices(faces2);
    std::sort(faces_indices1.begin(), faces_indices1.end());
    std::sort(faces_indices2.begin(), faces_indices2.end());


    //cout << "-.1-" << std::endl;
    bool ok = true;
    ok = ok && faces1.shape()[0] == faces2.shape()[0];
    //cout << "-.2-" << std::endl;
    ok = ok && faces1.shape()[1] == faces2.shape()[1];
    //cout << "-.3-" << std::endl;
    ok = ok && faces_indices1.shape()[0] == faces_indices2.shape()[0];
    //cout << "*" << std::endl;
    for (int i = 0; i < std::min(faces_indices1.shape()[0],faces_indices2.shape()[0]); ++i ) {
        //cout << "i*" << std::endl;
        ok = ok && (faces_indices1[i] == faces_indices2[i]);
        if (verbose_about_mismatches) {
            if (!(faces_indices1[i] == faces_indices2[i])) {
                cout << faces_indices1[i] << " != " << faces_indices2[i] << std::endl;
            }
        }
    }
    return ok;
}

TEST(Subdivision_1to4, square) {

    // cout << "=============================" << std::endl;
    // cout << ">>>>>>>>>>>>>>>>>>>>>" << std::endl;

    auto vf = testcase_square();
    auto faces = vf.second;
    auto verts = vf.first;
    print_faces(faces, faces.shape()[0]-1);

    /*
    std::set<edge_pair_type> edges_to_subdivide = {
        easy_edge(1, 2),
    };
    */
    std::set<faceindex_type> triangles_to_subdivide {0};

    midpointmap_type  midpoint_map;


    midpoint_map[easy_edge(1, 2)] = 99;

    std::cout << "b" << std::endl;

    //bool careful_for_twosides=true;

    auto result = subdivide_multiple_facets_1to4 (
        faces, verts, triangles_to_subdivide, midpoint_map);

    print_faces(std::get<1>(result), faces.shape()[0]-1);

    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    //cout << "=============================" << std::endl;

}

TEST(Subdivision_Utils, test_check_faces_equality) {

    auto faces1 = testcase_square().second;
    EXPECT_TRUE( check_faces_equality(faces1, faces1) );

    auto faces2 = testcase_triangle().second;
    EXPECT_TRUE( check_faces_equality(faces2, faces2) );

    EXPECT_FALSE( check_faces_equality(faces1, faces2) );

    auto faces3 =  f2f(std::vector<std::vector<vertexindex_type>> {
            {0, 3, 4}, {1,3,99}, {2,4,99}, {3,4,99-1}
        });
    EXPECT_FALSE( check_faces_equality(faces3, faces1) );
    EXPECT_FALSE( check_faces_equality(faces3, faces2) );

}

#define FACES_LITERAL(a, brackets) (f2f(std::vector<std::vector<vertexindex_type>> brackets))

TEST(Subdivision_1to4, triangle) {

    cout << ">>>>>>>>>>>>>>>>>>>>>" << std::endl;
    auto vf = testcase_triangle();
    auto faces_old = vf.second;
    auto verts = vf.first;
    // print_faces(faces_old, -1);
    std::set<faceindex_type> triangles_to_subdivide {0};

    midpointmap_type  midpoint_map;
    midpoint_map[easy_edge(1, 2)] = 99;
    //bool careful_for_twosides = true;
    assert(check_verts_faces(faces_old, verts));
    auto result = subdivide_multiple_facets_1to4 (
        faces_old, verts, triangles_to_subdivide, midpoint_map);
    auto result_faces = std::get<1>(result);
    // print_faces(result_faces, -1);



    EXPECT_FALSE( check_faces_equality(result_faces,
        f2f(std::vector<std::vector<vertexindex_type>> {
            {0, 3, 4}, {1,3,99}, {2,4,99}, {3,4,99-1}
        })
        //FACES_LITERAL({
        //    {0, 3, 4}, {1,3,99}, {2,4,99}, {3,4,99-1}
        //})
        //
        //FACES_LITERAL(a, {{0, 3, 4}, {1,3,99}, {2,4,99}, {3,4,99-1}})
    ) );


    vectorized_faces  desired = f2f(std::vector<std::vector<vertexindex_type>>{
        {0, 3, 4}, {1,3,99}, {2,4,99}, {3,4,99}
    });
    EXPECT_TRUE( check_faces_equality(desired, result_faces) );


    // {{0, 1, 2}} ->  {{4, 3, 5},  {0, 5, 3},  {1, 3, 4},  {2, 4, 5}}
    auto f1 = f2f(std::vector<std::vector<vertexindex_type>> {{0, 1, 2}} );
    auto f2 = f2f(std::vector<std::vector<vertexindex_type>> {{4, 3, 5},  {0, 5, 3},  {1, 3, 4},  {2, 4, 5}} );
    //EXPECT_TRUE( check_faces_equality(f2, std::get<1>(subdivide_multiple_facets_1to4 (f1, verts, triangles_to_subdivide, midpoint_map)) ));
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;


}

/*
    Cannot make it a loop in C++. How do you update a multi_array?
         faces_old = std::move(result_faces); ? Doesn't work
*/
std::pair<vectorized_faces, vectorized_vect> subdivide_recursively(
    const std::pair<vectorized_faces, vectorized_vect>& fv,
    int fi,
    int n,
    midpointmap_type& midpoint_map
) {
    if (n==0) {
        return fv;
    }

    std::set<faceindex_type> triangles_to_subdivide {0};  // one face only
    if (bool all_faces = true) {
        for (int i = 0; i < fv.first.shape()[0]; ++i) {  // all faces
            triangles_to_subdivide.insert(i);
        }
    }

    cout << "Going to subdivide the follinwg faces: ";
    for (auto i = triangles_to_subdivide.begin(); i != triangles_to_subdivide.end(); ++i) {
        cout << *i << ", ";
    }
    cout << std::endl;

    auto result_vf_tuple = subdivide_multiple_facets_1to4 (
        fv.first, fv.second, triangles_to_subdivide, midpoint_map);
    auto result_verts = std::get<0>(result_vf_tuple);
    auto result_faces = std::get<1>(result_vf_tuple);
    auto result_vf_pair = std::make_pair(result_faces, result_verts);

    if (bool verbose = false) {
        print_faces(fv.first, fv.first.shape()[0]-1);
        cout << " --> " << std::endl;
        print_faces(result_vf_pair.first, fv.first.shape()[0]-1);
        cout << "faces: " <<  result_vf_pair.second.shape()[0] << std::endl;
        cout << "vertices: " <<  result_vf_pair.first.shape()[0] << std::endl;
    }

    return subdivide_recursively(result_vf_pair, fi, n-1, midpoint_map);
}

int serpinski_expected_f(int n) {
    if (n==0)
        return 1;
    return serpinski_expected_f(n-1) * 4;
}

/*
    f(n) = (f==0? 1 ) || (f(n-1) * 4)   // faces of serpinski
    v(n) = (v==0? 3 ) || (v(n-1) * 4 - (3 + p(n-1) ) )  // vertices of serpinski
    p(n) = 3*2**n  // perimeter of the serpinski triangle
*/
int serpinski_expected_v(int n) {
    if (n==0)
        return 3;
    return serpinski_expected_v(n-1) * 4 - ( 3 + 3 * std::pow(2, n-1) );
}

TEST(Subdivision_1to4, randomised_iterative) {

    // to randomise:
    // http://stackoverflow.com/questions/6942273/get-random-element-from-container

    auto vf = testcase_triangle();
    auto faces_old = vf.second;
    auto verts = vf.first;
    midpointmap_type   midpoint_map;
    int NREC = 4;
    auto fv2 = subdivide_recursively(std::make_pair(vf.second, vf.first), 0, NREC, midpoint_map);
    /*
    print_faces(fv2.first, vf.second.shape()[0]-1);
    cout << "vertices: " <<  fv2.second.shape()[0] << std::endl;
    cout << "faces: " <<  fv2.first.shape()[0] << std::endl;
    */

    EXPECT_EQ(serpinski_expected_v(NREC), fv2.second.shape()[0]) << "number of vertices is incorrect";
    EXPECT_EQ(serpinski_expected_f(NREC), fv2.first.shape()[0]) << "number of faces is incorrect";

 }


#include "../subdivision/propagate_subdiv.hpp"
#include "../subdivision/edge_triplets.hpp"

TEST(propagate_subdiv_, trivial_example) {

    auto vf = testcase_triangle();
    auto faces = vf.second;

/*
    midpointmap_type   midpoint_map;
    midpoint_map[easy_edge(1, 2)] = 99;
*/

    std::set<edge_pair_type> requested_edges = {
        easy_edge(1, 2),
        easy_edge(3, 4)
    };

    auto r = mp5_implicit::subdivision::propagate_subdiv (
        faces,
        requested_edges
        // midpoint_map  //we should not need this
    );



}


#include "../subdivision/do_subdivision.hpp"
#include "../object_factory.hpp"

using mp5_implicit::subdivision::do_subdivision;
using mp5_implicit::subdivision::subdivide_given_faces;

/*
TEST(full_subdivision, trivial_example) {

    auto vf = testcase_triangle();
    auto faces = vf.second;
    auto verts = vf.first;

    bool use_metaball;
    //    std::string json = R"( {"type":"root","children":[{"type":"Difference","protected":false,"children":[{"type":"cylinder","displayColor":{"x":0.7675997200783986,"y":0.03892568708507049,"z":0.1754374135888661},"matrix":[35,0,0,0,0,35,0,0,0,0,9,0,0,0,0,1],"index":652818},{"type":"Difference","protected":false,"children":[{"type":"cylinder","displayColor":{"x":0.8122645344236872,"y":0.657334404743416,"z":0.7357336310755096},"matrix":[10,0,0,0,0,10,0,0,0,0,10,0,0,0,0,1],"index":1272174},{"type":"cylinder","displayColor":{"x":0.11421729990684737,"y":0.07562705374348999,"z":0.6324600862122098},"matrix":[10,0,0,0.658889604636343,0,10,0,6.215549332615993,0,0,10,1.3327027659215673e-7,0,0,0,1],"index":2463576}],"initialSize":{"x":1,"y":1,"z":1},"displayColor":{"x":0.6627450980392157,"y":0.4549019607843137,"z":0.7215686274509804},"matrix":[2.381193509886417,0,0,0.3600215429489424,0,2.381193509886417,0,0.5604901669421452,0,0,2.381193509886417,6.9059681360437395,0,0,0,1],"index":413872}],"initialSize":{"x":1,"y":1,"z":1},"displayColor":{"x":0.5529411764705883,"y":0.06666666666666667,"z":0.11764705882352941},"matrix":[1,0,0,0.32938436512727,0,1,0,0.15604124684634,0,0,1,0.000000000000014,0,0,0,1],"index":6565922}]})";
    std::string json = R"({"type":"Difference","protected":false,"children":[{"type":"cylinder","displayColor":{"x":0.7675997200783986,"y":0.03892568708507049,"z":0.1754374135888661},"matrix":[35,0,0,0,0,35,0,0,0,0,9,0,0,0,0,1],"index":652818},{"type":"Difference","protected":false,"children":[{"type":"cylinder","displayColor":{"x":0.8122645344236872,"y":0.657334404743416,"z":0.7357336310755096},"matrix":[10,0,0,0,0,10,0,0,0,0,10,0,0,0,0,1],"index":1272174},{"type":"cylinder","displayColor":{"x":0.11421729990684737,"y":0.07562705374348999,"z":0.6324600862122098},"matrix":[10,0,0,0.658889604636343,0,10,0,6.215549332615993,0,0,10,1.3327027659215673e-7,0,0,0,1],"index":2463576}],"initialSize":{"x":1,"y":1,"z":1},"displayColor":{"x":0.6627450980392157,"y":0.4549019607843137,"z":0.7215686274509804},"matrix":[2.381193509886417,0,0,0.3600215429489424,0,2.381193509886417,0,0.5604901669421452,0,0,2.381193509886417,6.9059681360437395,0,0,0,1],"index":413872}],"initialSize":{"x":1,"y":1,"z":1},"displayColor":{"x":0.5529411764705883,"y":0.06666666666666667,"z":0.11764705882352941},"matrix":[1,0,0,0.32938436512727,0,1,0,0.15604124684634,0,0,1,0.000000000000014,0,0,0,1],"index":6565922})";
    std::cout << json << std::endl;
    implicit_function* implicit_func = object_factory(json, use_metaball, true);

    REAL curvature_epsilon = 0.01;

    auto fv2 = do_subdivision (
        faces,
        verts,
        *implicit_func,
        curvature_epsilon,
        1.0  // randomized_probability = 1.0
    );
    auto f = std::get<0>(fv2);
    auto v = std::get<1>(fv2);
    cout << "new sizes: "<< f.shape()[0] << " " << v.shape()[0] << std::endl;

    gc_objects();
}
*/

void test_no_subdivision() {

}

TEST(full_subdivision, specific_triangles) {

    auto vf = testcase_triangle();
    auto faces = vf.second;
    auto verts = vf.first;

    std::set<faceindex_type>  which_facets_set;
    which_facets_set.insert(0);

    auto fv2 = subdivide_given_faces (
        faces,
        verts,
        which_facets_set
    );
    auto f2 = std::get<0>(fv2);
    auto v2 = std::get<1>(fv2);
    // cout << "new sizes: "<< f2.shape()[0] << " " << v2.shape()[0] << std::endl;
    bool b = check_faces_equality(f2, faces);
    EXPECT_TRUE(b);
    // EXPECT_EQ(verts, v2);
}
