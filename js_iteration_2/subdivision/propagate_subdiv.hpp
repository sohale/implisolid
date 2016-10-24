#pragma once

#include <vector>
#include <map>
#include <set>

#include "../basic_data_structures.hpp"
#include "../configs.hpp"
#include "../mesh_algorithms.hpp"
#include "../configs.hpp"

#include "edge_triplets.hpp"

namespace mp5_implicit {
namespace subdivision {


auto propagate_subdiv (
    const vectorized_faces &faces,
    //const boost::multi_array<edge_pair_type, 1>& subdivided_edges_codes,  // list of pair tuples?
    const std::set<edge_pair_type>& requested_1side_edgecode_set
    // alternative to subdivided_edges
    // subdivision::midpointmap_type& midpoint_map
) {
    //you can use either midpoint_map or subdivided_edges to get (and update) the edges that are being killed (and subdivided). In that case, all_edgecodes has to be updated also.

    // todo: (performance): we can pass a shared "edge_triplets_bool", and only undate the ones asked in set. Then request set is flushed.
    // subdivided_edges_codes

    auto et = subdivision::make_edge_triplets_bool (
        faces,
        requested_1side_edgecode_set
        // &midpoint_map
    );
    subdivision::triplet_bool_type    edge_triplets_bool  = std::move(std::get<0>(et));
    subdivision::edgecode_triplets_type   all_edgecodes   = std::move(std::get<1>(et));

    #define sides_booleans_Fx3 edge_triplets_bool

    int total_edges = 0;

    auto nfaces = edge_triplets_bool.shape()[0];
    boost::multi_array<unsigned char, 1> numsides {boost::extents[nfaces]};
    {
        auto it = edge_triplets_bool.begin(),
             e = edge_triplets_bool.end();
        auto numsides_iter = numsides.begin();
        for (; it != e; ++it) {
            unsigned char   // bool_t to short int
                c0 = (*it)[0] ? 1 : 0,
                c1 = (*it)[1] ? 1 : 0,
                c2 = (*it)[2] ? 1 : 0;
            *numsides_iter  = c0;
            *numsides_iter += c1;
            *numsides_iter += c2;

            total_edges += c0 + c1 + c2;

            /*
            *numsides_iter  = (*it)[0] ? 1 : 0;
            *numsides_iter += (*it)[1] ? 1 : 0;
            *numsides_iter += (*it)[2] ? 1 : 0;
            */
            ++numsides_iter;
        }
    }

    {
        // combine 2 and 3
        auto it = numsides.begin(),
             e = numsides.end();
        for (; it != e; ++it) {
            if (*it == 3) {
                *it = 2;
            }
        }
    }

    boost::multi_array<edge_pair_type, 1>  edges_need_subdivision {boost::extents[total_edges]};
    {
        // Copy the edges into a flat array, those whole bool is true in edge_triplets_bool.
        auto eit = edges_need_subdivision.begin();
        auto it = edge_triplets_bool.begin(),
             e = edge_triplets_bool.end();
        auto all_edgecodes_iter = all_edgecodes.begin();
        for (; it != e; ++it) {
            if ((*it)[0]) *(eit++) = (*all_edgecodes_iter)[0];
            if ((*it)[1]) *(eit++) = (*all_edgecodes_iter)[1];
            if ((*it)[2]) *(eit++) = (*all_edgecodes_iter)[2];
            ++all_edgecodes_iter;
        }
        assert( eit == edges_need_subdivision.end());
        assert( all_edgecodes_iter == all_edgecodes.end());
    }

    // We only propagate triangles with # subdivided sides 1,2, or 3 (i.e. not including 0 and 4)
    // i.e. {i: numsides[i] = 1,2,3}
    // std::vector<edge_index>
    std::vector<faceindex_type> faces_subset; //amortized
    int capacity = nfaces / 5;
    faces_subset.reserve(capacity);  // Use 75% median

    std::map<int, subdivision::faces_subset_type>  propag_dict;

    // no need for 3!  Insrtead if 1,2,3 ---> 1, "2,3", i.e. 1,2.
    for (int c = 1; c <= 2; ++c ) {
        // Reuse the vector for performance. Reset the size from the previous iteration.
        faces_subset.resize(0);
        for (faceindex_type fi = 0; fi < nfaces; ++fi) {
            if (numsides[fi] == c) {
                faces_subset.push_back(fi);
            }
        }

        subdivision::faces_subset_type   subset_rigid_size {boost::extents[faces_subset.size()]};
        assert (subset_rigid_size.shape()[0] == faces_subset.size());
        std::copy(faces_subset.begin(), faces_subset.end(), subset_rigid_size.begin());

        // propag_dict.insert(std::move(faces_subset));
        // propag_dict.emplace(...std::move(faces_subset));
        // propag_dict.insert(std::make_pair(c, std::move(faces_subset)));
        propag_dict.emplace(std::make_pair(c, std::move(subset_rigid_size)));
        // faces_subset is invalid now.
    }
    // return std::move(propag_dict);  // do we need the move?
    return std::make_tuple(propag_dict, edges_need_subdivision);
}



}  // namespace subdivision
}  // namespace mp5_implicit

