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
void subdivide_SD2(faces, edges_with_1_side, midpoint_map, bool careful_for_twosides=true) {
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


#include "basic_data_structures.hpp"
#include "mesh_algorithms.hpp"

using mp5_implicit::CONFIG_C;

using mp5_implicit::encode_edge__sort;
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
vectorized_faces subdivide_1to2(const vectorized_faces & faces,
    const std::set<edge_pair_type>& edges_with_1_side,
    const std::map<edge_pair_type, vectorized_vect::index>& midpoint_map,
    const edge_pair_type EdgecodeBase,
    bool careful_for_twosides=true)
{
    //alternative names: edges_with_1_side
    //alternative names: midpoint_map

    // Brief: Forms an aray of  (L,T,T,M) and does lookup in the SET and in the MAP

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
    // For each edge in the mesh, find whether it belongs to the set edges_with_1_side or not. (the set of requested edges = edges_with_1_side). A boolean for each edge: Fx3x<bool>
    // - There will be at most one true value for each face.
    // We make a vector of faces indices (relevant_faces_indices) based on that booleans array.
    // Each face will an index (0..2) that specified which side belongs there. Call it relevant_sides.
    // relevant_faces_sides is the compactified version of relevant_sides is necessary.
    // relevant_faces_indices and relevant_faces_sides together paired (another version of relevant_sides is necessary, which is compacted).
    // has M elements.
    // We use an array of the size M, which contains the mapped values.
    // make two new faces arrays that comibne: 1-the mapped index, 2-indices of the original faces of the third vertex, (by shifting them) 3- indices of the original faces: one for each of the remaining vertices (by shifting them).

    /*
    typedef boost::multi_array<vectorized_vect::index, 4>  vectorized_new_faces_type;
    typedef boost::multi_array<edge_pair_type, 3>  edgecode_triplets_type;
    typedef boost::multi_array<bool_t, 3>  triplet_bool_type;
    */
    typedef boost::multi_array<vectorized_vect::index, 2>  vectorized_new_faces_type;
    typedef boost::multi_array<edge_pair_type, 2>  edgecode_triplets_type;
    typedef boost::multi_array<bool_t, 2>  triplet_bool_type;

    // typedef boost::array<edgecode_triplets_type::index, 2>  edgecode_triplets_shape_type;

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

    triplet_bool_type
        edge_triplets_bool {boost::extents[original_faces_count][3]};


    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {
        const auto not_found = midpoint_map.end();
        edge_triplets_bool[fi][0] = (midpoint_map.find(all_edgecodes[fi][0]) != not_found);
        edge_triplets_bool[fi][1] = (midpoint_map.find(all_edgecodes[fi][1]) != not_found);
        edge_triplets_bool[fi][2] = (midpoint_map.find(all_edgecodes[fi][2]) != not_found);
    }


    boost::multi_array<unsigned char, 1>  set_lookup_side {boost::extents[original_faces_count]};

    // todo (performance): The below method uses lookup in the set as a map (log(N) access time). However, since not only the set is sorted, but also all_edgecodes[] are almost sorted, we May be able to find membership of our item more efficiently.
    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {

        const auto not_found = midpoint_map.end();

        // only one side will match
        unsigned char found_side = 0; // zero => not found
        for (int side = 0; side < 3; ++side ) {

            if (midpoint_map.find(all_edgecodes[fi][side]) != not_found) {  // slow
                found_side = side + 1;  // non-zero => found
                break;  // skip the other slow parts
            }
        }
        set_lookup_side[fi] = found_side;
    }
    /*!
    set_lookup_side[] values:
        0: no subdivision necessary,
        1:subdiv first side 2 ([0]),
        2: subdiv side 2 ([1]),
        3 subdiv side 3 ([2]).
    */

    cout << "3434343434" << std::endl << std::flush;

    // next: 1-make the faces in an array
    // Then, later, or at the same time, circular-shift, and get the 4 indices.
    boost::multi_array<vectorized_faces::index, 1>::size_type preallocated_compactified_maximum_size =
        original_faces_count; // set_lookup_side.shape()[0] == original_faces_count
    boost::multi_array<vectorized_faces::index, 1>  compactified_faces_indices{boost::extents[preallocated_compactified_maximum_size]};
    boost::multi_array<vectorized_faces::index, 1>::size_type compactified_faces_indices_effective_size = 0;
    boost::multi_array<unsigned char, 1>  compactified_side{boost::extents[preallocated_compactified_maximum_size]};

    cout << "3.3.3.3.3=." << std::endl << std::flush;
    cout << "3.3.3.3.3=." << std::endl << std::flush;
    cout << "3.3.3.3.3=." << std::endl << std::flush;
    cout << " cc " << std::endl << std::flush;

    cout << std::flush;
    cout << " c = " << original_faces_count << std::endl << std::flush;
    cout << std::flush;

    for (edgecode_triplets_type::index fi = 0; fi < original_faces_count; ++fi ) {
        cout << "here" << std::flush;
        cout << fi << ":" << std::flush;
        cout << set_lookup_side[fi] << " = " << static_cast<int>(set_lookup_side[fi]) << " " << std::endl << std::flush;;
        if (set_lookup_side[fi]) {
            cout << "if=true" << std::endl << std::flush;;
            auto s = set_lookup_side[fi] - 1;  // careful, it is unsigned
            cout << "minus one = " << s << std::endl << std::flush;;

            compactified_side [compactified_faces_indices_effective_size] = set_lookup_side[fi] - 1;
            compactified_faces_indices [compactified_faces_indices_effective_size] = fi;
            ++compactified_faces_indices_effective_size;
        }
    }
    cout << std::endl << std::flush;

    cout << "3333333" << std::endl << std::flush;


    // variables: compactified_*  or newface_*
    typedef vectorized_new_faces_type::index  compactified_index_type;
    vectorized_new_faces_type  compactified_newfaces_specs
        {boost::extents[compactified_faces_indices_effective_size][4]};
    /* in fact each row in compactified_newfaces_specs[][] is
        an array of tuple: tuple<left,right,third,mid> :

    <L,R,T,M>

    compactified_newfaces_specs: 0,1,2  contains edges: 0-1, 1-2, 2-0
    (side,side+1) (side+1,side+2), (side+2, side+0)
    (0,1) (1,2), (2,0)  , where (0,1) is the subsdivided edge. Hence tghe vertices are: (side is subtracted from indices):

        T            2
       /|\          /|\
      / | \        / | \
     /  |  \      /  |  \
    L---M---R    0---M---1

    Edges, as stored in faces after cyclic shift (using side = compactified_side[cfi]) are:
    (LR, RT, TL)


    In array: compactified_newfaces_specs:
        T            2
       /|\          /|\
      / | \        / | \
     /  |  \      /  |  \
    L---M---R    0---3---1
    The third elements, compactified_newfaces_specs[][3], will contain the indeix of the mid-point (M).

    Hence, compactified_newfaces_specs[cfi][:] contains:
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
        auto _L = faces[fi][side];
        auto _R = faces[fi][(side + 1) % 3];
        auto _T = faces[fi][(side + 2) % 3];
        compactified_newfaces_specs[cfi][0] = _L; //faces[fi][(side)        ];
        compactified_newfaces_specs[cfi][1] = _R;
        // Third: not in the matched edge:
        compactified_newfaces_specs[cfi][2] = _T;

        #if ASSERT_USED
            // Fourth: to be looked up from the map
            compactified_newfaces_specs[cfi][3] = -1;
        #endif

        // Keep the side (shift)'s code, because we need to look it up in the next step.
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
        compactified_newfaces_specs[cfi][3] = midpoint_map[edge_code];
        */

        auto pair_itr = midpoint_map.find(edge_code);
        assert(pair_itr != midpoint_map.end());
        // auto pair = *pair_itr;
        compactified_newfaces_specs[cfi][3] = pair_itr->second;

    }

    vectorized_faces f{boost::extents[5][3]};
    clog << all_edgecodes[0][0];
    clog << edge_triplets_bool[0][0];
    return f;
}
