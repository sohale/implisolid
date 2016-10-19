#pragma once

#include <cstddef>   // for std::nullptr only

namespace mp5_implicit {
namespace subdivision {


    template<typename iterator, typename T>
bool check_minimum(iterator begin, iterator end, T minimum_value) {
    bool ok = true;
    for (auto i = begin; i != end; ++i) {
        ok = ok && (*i > minimum_value);
        if (!ok) {
            break;
        }
    }
    return ok;
}

    template<typename value_type>
inline value_type bool_to_1(bool x) {
    return x? static_cast<value_type>(1) : static_cast<value_type>(0);
}

/*
typedef boost::multi_array<vectorized_vect::index, 4>  vectorized_new_faces_type;
typedef boost::multi_array<edge_pair_type, 3>  edgecode_triplets_type;
typedef boost::multi_array<bool_t, 3>  triplet_bool_type;
*/
typedef boost::multi_array<vectorized_vect::index, 2>  vectorized_new_faces_type;
typedef boost::multi_array<edge_pair_type, 2>  edgecode_triplets_type;
typedef boost::multi_array<bool_t, 2>  triplet_bool_type;

// typedef boost::array<edgecode_triplets_type::index, 2>  edgecode_triplets_shape_type;

/**
 * @brief      Makes an edge triplets bool based on a given set of edges.
 *
 * @param[in]  faces                         The facets of the mesh.
 * @param[in]  requested_1side_edgecode_set  A std::set of edges specified by edgecodes.
 * @param[in]  midpoint_map_ptr              only for checking the consistency: whether requested_1side_edgecode_set is a subset of this map.
 *
 * @return     { two arrays, one containing all the edge codes, the other contains 3xF booleans for whether that edge is among the requested edges}
 */
std::tuple<triplet_bool_type, edgecode_triplets_type>
make_edge_triplets_bool (
    const vectorized_faces & faces,
    const std::set<edge_pair_type>& requested_1side_edgecode_set,
    const subdivision::midpointmap_type* midpoint_map_ptr = nullptr
) {
    /* assertions */

    if (midpoint_map_ptr != nullptr) {
    const subdivision::midpointmap_type&  midpoint_map = *midpoint_map_ptr;

    // Assert edgceodes are non-zero. No correct edgecode is 0.
    assert( check_minimum(requested_1side_edgecode_set.begin(), requested_1side_edgecode_set.end(), 1)
        && "Assert edgecodes are non-zero, i.e. no correct edgecode is 0." );

    //"assert requested_1side_edgecode_set is subset of co-range of midpoint_map"
    auto belongs_to_midpoint_map = [&midpoint_map](auto v) {return midpoint_map.find(v) != midpoint_map.end();};  // being_found__in_map_ness
    // Make sure there is a mapping for every in requested_1side_edgecode_set

    assert(std::all_of( requested_1side_edgecode_set.begin(), requested_1side_edgecode_set.end(),
            belongs_to_midpoint_map
        ) && "Make sure there is a mapping for every in requested_1side_edgecode_set"
    );
    }



    //edgecode_triplets_type::shape_type s;


    //  why is this long? :  faces.shape()[0]
    auto
        original_faces_count
            // = static_cast<edgecode_triplets_type::index>(faces.shape()[0]);
            = faces.shape()[0];

    edgecode_triplets_type
        all_edgecodes
            { boost::extents[original_faces_count][3] };
        //python name: all_edges_triples


    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {
        vertexindex_type_ v0 = faces[fi][0];
        vertexindex_type_ v1 = faces[fi][1];
        vertexindex_type_ v2 = faces[fi][2];

        edge_pair_type e0 = encode_edge__sort(v0, v1, CONFIG_C::edgecode_base);
        edge_pair_type e1 = encode_edge__sort(v1, v2, CONFIG_C::edgecode_base);
        edge_pair_type e2 = encode_edge__sort(v2, v0, CONFIG_C::edgecode_base);

        all_edgecodes[fi][0] = e0;
        all_edgecodes[fi][1] = e1;
        all_edgecodes[fi][2] = e2;


        // all_edgecodes[fi][0] = encode_edge(faces[fi][0], faces[fi][1], CONFIG_C::edgecode_base);
        // all_edgecodes[fi][1] = encode_edge(faces[fi][1], faces[fi][2], CONFIG_C::edgecode_base);
        // all_edgecodes[fi][2] = encode_edge(faces[fi][2], faces[fi][0], CONFIG_C::edgecode_base);

    }

    // contains the result of set lookup
    triplet_bool_type
        edge_triplets_bool {boost::extents[original_faces_count][3]};

    /*
    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {
        const auto not_found = midpoint_map.end();
        edge_triplets_bool[fi][0] = (midpoint_map.find(all_edgecodes[fi][0]) != not_found);
        edge_triplets_bool[fi][1] = (midpoint_map.find(all_edgecodes[fi][1]) != not_found);
        edge_triplets_bool[fi][2] = (midpoint_map.find(all_edgecodes[fi][2]) != not_found);
    }
    */
    // todo: rename: requested_1side_edgecode_set -> requested_edges_with_1_side
    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {
        const auto not_found = requested_1side_edgecode_set.end();
        edge_triplets_bool[fi][0] = (requested_1side_edgecode_set.find(all_edgecodes[fi][0]) != not_found);
        edge_triplets_bool[fi][1] = (requested_1side_edgecode_set.find(all_edgecodes[fi][1]) != not_found);
        edge_triplets_bool[fi][2] = (requested_1side_edgecode_set.find(all_edgecodes[fi][2]) != not_found);
    }

    #if ASSERT_USED
    // Make sure at most one edge is requested per triangle.
    auto at_most_one_edge_requested = [](auto triple) {
        return bool_to_1<short>(triple[0]) + bool_to_1<short>(triple[1]) + bool_to_1<short>(triple[2])
            <= 1;
    };  // being_found__in_map_ness

    assert(
        std::all_of( edge_triplets_bool.begin(), edge_triplets_bool.end(),
            at_most_one_edge_requested
        )
        && "Make sure at most one edge is requested per triangle to be subdivided 1->2."
    );
    #endif


    //#if ASSERT_USED
    if (bool debug_print_triplets = true) {
        for (
                edgecode_triplets_type::index fi = 0;
                fi < original_faces_count;
                ++fi )
        {
            cout << fi <<": ";
            for (int j=0; j < 3; ++j) {
                cout << " " << (edge_triplets_bool[fi][j]? "Y" : "-");
            }
            cout << "   ";
            for (int j=0; j < 3; ++j) {
                cout << " " << faces[fi][j];
            }
            cout << std::endl;
        }
    }
    //#endif

    // doe we need to say std::move() explicily?
    return std::make_tuple(edge_triplets_bool, all_edgecodes);
}


}  // namespace subdivision
}  // namespace mp5_implicit
