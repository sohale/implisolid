#pragma once

// Each method does the subdivision on multiple triangles.

/*
Two types of subdivision are done here:
1 -> 2
1 -> 4
There are others possible: proportional dubdivisions,

enum class SubdivType {SD12, SD14, };
// SubdivType::SD12
// SD1  SD4

// enum {a,b,c}
*/

/*
Applies subdivision on multiple facets
*/
/*
void subdivide_SD2(faces, requested_1side_edgecode_set, midpoint_map, bool careful_for_twosides=true) {
    ;
}
*/

/*
Sep 23 17:14

// Eigen3
SparseMatrix


http://regexr.com/

var is scoped to the nearest function block and let is scoped to the nearest enclosing block (
!!!!!!!!!
http://stackoverflow.com/questions/762011/let-keyword-vs-var-keyword-in-javascript
*/


#include "../basic_data_structures.hpp"
#include "../configs.hpp"
#include "../mesh_algorithms.hpp"
#include "edge_triplets.hpp"

using mp5_implicit::CONFIG_C;

using mp5_implicit::encode_edge__sort;

using mp5_implicit::subdivision::make_edge_triplets_bool;
using mp5_implicit::subdivision::vectorized_new_faces_type;
using mp5_implicit::subdivision::edgecode_triplets_type;
using mp5_implicit::subdivision::triplet_bool_type;
using mp5_implicit::subdivision::midpointmap_type;



namespace mp5_implicit {
namespace subdivision {


// #if ASSERT_USED
//     // Assert with an error message
//     #define assert_message1(x, error_message) { \
//         if(!(#x)) \
//             {cerr << "Assertion error:" << std::endl << (#error_message);abort();} \
//     }
// #else
//     // #define assert_message(x, error_message) {}
// #endif
//
// #define assert_message(x, m) assert(#x)


/*
bool (setbegin, setend, map)
    bool ok = true;
    for (auto i = begin; i != end; ++i) {
        ok = ok && (map.find(*i) != end);
        if (!ok) {
            break;
        }
    }
    return ok;
}
*/

/*  An alternative streategy for subdivide_1to2():

    Can be done using a map. When we put in the "set", we can provide the original face index too.
    So the job of subdivide_1to2() will be trivial. It may need to do a set intersection.

    Also alternative for deleting: For each new face, we replace the first one in the place of the original one, and only the second one is to be appended.
    Hence, we avoid the cost of deleting, and we add less.

    We can use a std::set of pair<v1,v2> or even tuple<v1,v2,face_index>.

    For now, the following implements the original python implementation.
*/

/*!
*/
auto subdivide_1to2_LRTM(const vectorized_faces & faces,
    const std::set<edge_pair_type>& requested_1side_edgecode_set,
    const midpointmap_type& midpoint_map,
    // const edge_pair_type EdgecodeBase,
    bool careful_for_twosides=true)
{
    //alternative names: edges_with_1_side requested_1side_edgecode_set
    //alternative names: midpoint_map

    // Brief: Forms an aray of  (L,T,T,M) and does lookup in the SET and in the MAP.
    // Lookup and compactify, shift, lookup2, and store all in the LRTM array.

    //Summry of the algorithm:
    // * lookup face indices in the set
    //   - do this by looking up the edges.
    //       - for this: calculate the unique edge codes for each edge.
    // * compactify
    //   - keep the face  indices (as new faces)
    //   - keep the side/edge indices
    // * lookup edge indices in the map
    // * shift the triples of each (new) face
    // * make two faces for each "new" face:
    //   - comprised of : the third vertex, the mapped vertex, and two other sides.
    //   - array of 4 elements! v_shared, v_left, v_right, v_mapped (i.e. v_new)
    // * apply addition and deletion.

    // Make an array of ids for edges (edges_triplets). This makes it easy to make requested_faces_indices.
    // For each edge in the mesh, find whether it belongs to the set requested_1side_edgecode_set or not. (the set of requested edges = requested_1side_edgecode_set). A boolean for each edge: Fx3x<bool>
    // - There will be at most one true value for each face.
    // We make a vector of faces indices (relevant_faces_indices) based on that booleans array.
    // Each face will an index (0..2) that specified which side belongs there. Call it relevant_sides.
    // relevant_faces_sides is the compactified version of relevant_sides is necessary.
    // relevant_faces_indices and relevant_faces_sides together paired (another version of relevant_sides is necessary, which is compacted).
    // has M elements.
    // We use an array of the size M, which contains the mapped values.
    // make two new faces arrays that comibne: 1-the mapped index, 2-indices of the original faces of the third vertex, (by shifting them) 3- indices of the original faces: one for each of the remaining vertices (by shifting them).

    // Another direction of movement can be based on the fact that in that algorithm, each edge appears only in one triangle. Becasue the other triangles is already subdivided, probably by subdivide1->4.
    //

    // std::tuple<triplet_bool_type, edgecode_triplets_type>
    auto et = make_edge_triplets_bool (
        faces,
        requested_1side_edgecode_set,
        &midpoint_map
    );
    triplet_bool_type   edge_triplets_bool = std::move(std::get<0>(et));
    edgecode_triplets_type  all_edgecodes   = std::move(std::get<1>(et));

    auto original_faces_count = faces.shape()[0];

    /*
    // Why does this lookup the set in the map?
    // why even set lookup?
    boost::multi_array<unsigned char, 1>  set_lookup_side_plus1 {boost::extents[original_faces_count]};

    // todo (performance): The below method uses lookup in the set as a map/set (log(N) access time). However, since not only the set is sorted, but also all_edgecodes[] are almost sorted, we May be able to find membership of our item more efficiently.
    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {

        const auto not_found = midpoint_map.end();

        // only one side will match
        unsigned char found_side = 0; // zero => not found
        for (int side = 0; side < 3; ++side ) {
            // Problem: It first should check if it's in the SET. This doesnt guarantee it is actually requested to be midpoint-mapped, unless it was in the set.
            auto ll = all_edgecodes[fi][side];
            if (midpoint_map.find(ll) != not_found) {  // slow
                found_side = side + 1;  // non-zero => found
                break;  // skip the other slow parts
            }
        }
        set_lookup_side_plus1[fi] = found_side;
    }
    */

    // Not sure.
    // set_lookup_side_plus1 contains which side was true after the lookup. It also indicates whether each edge was found in the set lookup.
        boost::multi_array<unsigned char, 1>  set_lookup_side_plus1 {boost::extents[original_faces_count]};
        for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {
            unsigned char found_side = 0; // zero => not found
            for (int side = 0; side < 3; ++side ) {
                if (edge_triplets_bool[fi][side]) {
                    found_side = side + 1;  // non-zero => found
                    break;  // skip the other slow parts
                }
            }
            set_lookup_side_plus1[fi] = found_side;
        }

    /*!
    set_lookup_side_plus1[] values:
        0: no subdivision necessary,
        1:subdiv first side 2 ([0]),
        2: subdiv side 2 ([1]),
        3 subdiv side 3 ([2]).
    */


    // next: 1-make the faces in an array
    // Then, later, or at the same time, circular-shift, and get the 4 indices.


    /*
    boost::multi_array<vectorized_faces::index, 1>::size_type preallocated_compactified_maximum_size =
        original_faces_count; // set_lookup_side_plus1.shape()[0] == original_faces_count
    boost::multi_array<vectorized_faces::index, 1>  compactified_faces_indices{boost::extents[preallocated_compactified_maximum_size]};
    boost::multi_array<vectorized_faces::index, 1>::size_type compactified_faces_indices_effective_size = 0;
    boost::multi_array<unsigned char, 1>  compactified_side{boost::extents[preallocated_compactified_maximum_size]};
    */

    // for more memory efficient but less speed effciency, use vector<vectorized_faces::index> (and use capacity)
    // but now it can be less time-efficient
    std::vector<vectorized_faces::index> compactified_faces_indices;  // {0} will leave one element, (0) will be empty.
    auto capacity = requested_1side_edgecode_set.size();  // should be enough given the "single_tiangle edges" condition.
    compactified_faces_indices.reserve(capacity);
    boost::multi_array<vectorized_faces::index, 1>::size_type compactified_faces_indices_effective_size = 0;


    cout << "compactified_faces_indices.size() " << compactified_faces_indices.size() << std::endl;
    assert(compactified_faces_indices.size() == 0);
    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {
        // cout << fi << ":" << std::flush;
        // cout << set_lookup_side_plus1[fi] << " = " << static_cast<int>(set_lookup_side_plus1[fi]) << " " << std::endl << std::flush;;
        if (set_lookup_side_plus1[fi]) {
            compactified_faces_indices.push_back(fi);
            assert(compactified_faces_indices [compactified_faces_indices_effective_size] == fi);
            ++ compactified_faces_indices_effective_size;
        }
    }
    // consolidate compactified_faces_indices_effective_size; i.e. make it const

    boost::multi_array<unsigned char, 1>  compactified_side{boost::extents[compactified_faces_indices_effective_size]};

    // since we know the final size, it is more efficient not to use vector
    edgecode_triplets_type::index cfi = 0;
    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {
        if (set_lookup_side_plus1[fi]) {
            compactified_side [cfi] = set_lookup_side_plus1[fi] - 1; // careful, it is unsigned
            ++ cfi;
        }
    }
    assert(compactified_faces_indices_effective_size == cfi);
    // del cfi
    assert(compactified_faces_indices.size() >= compactified_faces_indices_effective_size);
    assert(compactified_side.shape()[0] >= compactified_faces_indices_effective_size);

    // Test all compactified_faces_indices[] are in requested_1side_edgecode_set[]
    #if ASSERT_USED
        for (auto fi : compactified_faces_indices) {
            auto notfound = requested_1side_edgecode_set.end();
            int side = static_cast<int>(set_lookup_side_plus1[fi]) - 1;
            //bool edgebool = edge_triplets_bool[fi][side];
            auto edgecode = all_edgecodes[fi][side];
            //cout << fi << ":  edgecode:" << edgecode << " side:"<< side  << "   ";
            //cout << "(";
            for( int j = 0; j < 3; ++j) {
                //cout <<  all_edgecodes[fi][j] << " ";
            }
            //cout << ")";
            //cout << std::endl;
            assert( requested_1side_edgecode_set.find(edgecode) != notfound );
        }
        //cout << std::endl;
    #endif

    #if ASSERT_USED
        // never test this. Looks like an incorrect assertion. However I use it since it was on the Python side.
        // This test is necessary only if we are sure one side is subdivided already with subdivide1->4, hance, each edge is subdivided only once (as opposed to twice).
        constexpr bool SINGLE_TRIANGLE_ONLY = false;
        if (SINGLE_TRIANGLE_ONLY) {
        std::set<edge_pair_type> unique_edges_found;
        for (auto fi : compactified_faces_indices) {
            int side = static_cast<int>(set_lookup_side_plus1[fi]) - 1;
            auto edgecode = all_edgecodes[fi][side];
            unique_edges_found.insert(edgecode);
        }
        // Coretto
        for (auto ec : unique_edges_found)
            cout << ec << " ";
        cout << std::endl;

        // Why
        cout << compactified_faces_indices_effective_size << std::endl;
        cout << " == " << unique_edges_found.size() << std::endl;
        assert(compactified_faces_indices_effective_size == unique_edges_found.size());
        }
    #endif

    // variables: compactified_*  or newface_*
    typedef vectorized_new_faces_type::index  compactified_index_type;
    vectorized_new_faces_type  compactified_newfaces_LRTM_specs
        {boost::extents[compactified_faces_indices_effective_size][4]};
    /* in fact each row in compactified_newfaces_LRTM_specs[][] is
        an array of tuple: tuple<left,right,third,mid> :

    <L,R,T,M>

    compactified_newfaces_LRTM_specs: 0,1,2  contains edges: 0-1, 1-2, 2-0
    (side,side+1) (side+1,side+2), (side+2, side+0)
    (0,1) (1,2), (2,0)  , where (0,1) is the subsdivided edge. Hence tghe vertices are: (side is subtracted from indices):

        T            2
       /|\          /|\
      / | \        / | \
     /  |  \      /  |  \
    L---M---R    0---M---1

    Edges, as stored in faces after cyclic shift (using side = compactified_side[cfi]) are:
    (LR, RT, TL)


    In array: compactified_newfaces_LRTM_specs:
        T            2
       /|\          /|\
      / | \        / | \
     /  |  \      /  |  \
    L---M---R    0---3---1
    The third elements, compactified_newfaces_LRTM_specs[][3], will contain the indeix of the mid-point (M).

    Hence, compactified_newfaces_LRTM_specs[cfi][:] contains:
    [L,R,T,M]

    */

    boost::multi_array<edge_pair_type, 1>  compactified_splitting_side_edgecode
        {boost::extents[compactified_faces_indices_effective_size]};


    for (
            compactified_index_type cfi = 0;
            cfi < compactified_faces_indices_effective_size;
            ++cfi )
    {
        // v1,v2 (of edge), third v, mapped_v
        vectorized_faces::index fi = compactified_faces_indices[cfi];
        auto side = compactified_side[cfi];
        // First and second nodes of the matched edge.
        auto _L = faces[fi][(side)        ];
        auto _R = faces[fi][(side + 1) % 3];
        auto _T = faces[fi][(side + 2) % 3];
        compactified_newfaces_LRTM_specs[cfi][0] = _L;
        compactified_newfaces_LRTM_specs[cfi][1] = _R;
        // Third: not in the matched edge:
        compactified_newfaces_LRTM_specs[cfi][2] = _T;

        #if ASSERT_USED
            // Fourth: to be looked up from the map
            compactified_newfaces_LRTM_specs[cfi][3] = -1;
        #endif

        // Keep the side (shift)'s code, because we need to look it up in the next step.
        // A bit slow
        compactified_splitting_side_edgecode[cfi] = all_edgecodes[fi][side]; // The LR side.
    }



    for (
            compactified_index_type cfi = 0;
            cfi < compactified_faces_indices_effective_size;
            ++cfi )
    {
        edge_pair_type edge_code = compactified_splitting_side_edgecode[cfi];

        /* Will not work, because a [] access to a map is not read-only.
        assert(midpoint_map.find(edge_code) != midpoint_map.end());
        // will not work because the map midpoint_map[] is const.
        compactified_newfaces_LRTM_specs[cfi][3] = midpoint_map[edge_code];
        */

        // a bit slow
        auto pair_itr = midpoint_map.find(edge_code);
        assert(pair_itr != midpoint_map.end()
            && "The midpoint lookup failed. No lookup for the requested edge."
        );
        // auto pair = *pair_itr;
        compactified_newfaces_LRTM_specs[cfi][3] = pair_itr->second;
    }

    // The LRTM array now contains all information needed for subdivision. (contains vertex indices)
    // But we also will need compactified_faces_indices[]. (contains face indices)
    // So the output is: (LRTM, F) <--> SubDiv1-2


    #if VERBOSE_SUBDIV
    cout << " :  L R T M   -- F" << std::endl;
    for (
            compactified_index_type cfi = 0;
            cfi < compactified_faces_indices_effective_size;
            ++cfi )
    {
        cout << cfi <<": ";
        for (int j=0; j < 4; ++j) {
            cout << " " << compactified_newfaces_LRTM_specs[cfi][j];
        }
        cout << "  \t" << "[" << compactified_faces_indices[cfi] << "]";
        cout << std::endl;
    }
    assert(compactified_faces_indices.size() >= compactified_faces_indices_effective_size);
    #endif

    // note that compactified_faces_indices does not have the right size.

    // uses std::move
    // note: the size of compactified_faces_indices is not right
    return std::make_tuple(compactified_newfaces_LRTM_specs, compactified_faces_indices);

    /*
    return std::tuple<
        boost::multi_array<vectorized_faces::index, 1>,
        vectorized_new_faces_type
    >(compactified_faces_indices, compactified_newfaces_LRTM_specs);
    */

    /*

    vectorized_faces f{boost::extents[5][3]};
    clog << all_edgecodes[0][0];
    clog << edge_triplets_bool[0][0];
    return f;
    */
}

//subdivide_1to2_multiple_facets
vectorized_faces  subdivide_1to2 (const vectorized_faces & faces,
    const std::set<edge_pair_type>& requested_1side_edgecode_set,
    const midpointmap_type& midpoint_map,
    //const edge_pair_type EdgecodeBase,
    bool careful_for_twosides=true)
{
    auto tupl = subdivide_1to2_LRTM(faces, requested_1side_edgecode_set, midpoint_map, careful_for_twosides);

    auto compactified_faces_indices = std::get<1>(tupl);
    auto compactified_newfaces_LRTM_specs = std::get<0>(tupl);

    auto nfaces_old = faces.shape()[0];
    auto nfaces_new = compactified_newfaces_LRTM_specs.shape()[0];
    std::cout << "nfaces_old + nfaces_new = " << nfaces_old <<  " + " << nfaces_new << std::endl;
    vectorized_faces newfaces =  vectorized_faces{boost::extents[nfaces_old + nfaces_new][3]};
    // newfaces = faces;
    std::copy(faces.begin(), faces.end(), newfaces.begin());

    auto newface_ctr = nfaces_old;
    for (int i = 0; i < compactified_newfaces_LRTM_specs.shape()[0]; ++i) {
        const auto _L = compactified_newfaces_LRTM_specs[i][0];
        const auto _R = compactified_newfaces_LRTM_specs[i][1];
        const auto _T = compactified_newfaces_LRTM_specs[i][2];
        const auto _M = compactified_newfaces_LRTM_specs[i][3];

        const auto fi = compactified_faces_indices[i];

        // now add two triangles, LTM and RTM, based on LRTM.

        // replace the new faces: TLM (but in the order: preserve the counter-clockwise order)
        newfaces[fi][0] = _L;
        newfaces[fi][1] = _M;
        newfaces[fi][2] = _T;

        // insert the new faces: RTM
        newfaces[newface_ctr][0] = _R;
        newfaces[newface_ctr][1] = _T;
        newfaces[newface_ctr][2] = _M;
        ++newface_ctr;
    }

    // cout << "NOK " << newfaces.shape()[0] << " " <<  VERBOSE_SUBDIV << std::endl;
    #if VERBOSE_SUBDIV
    // cout << "OK" << std::endl;
    for (int i = 0; i < newfaces.shape()[0]; ++i )
    {
        cout << i <<": ";
        for (int j=0; j < 3; ++j) {
            cout << " " << newfaces[i][j];
        }
        if (i == nfaces_old-1 )
            cout << " <--- last old";
        cout << std::endl;
    }
    #endif

    return newfaces;
}

}  // namespace subdivision
}  // namespace mp5_implicit
