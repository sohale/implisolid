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
        auto q = all_edgecodes[fi][0];
        edge_pair_type e0 = q;
        all_edgecodes[fi][0] = (midpoint_map.find(e0) != midpoint_map.end());
        // all_edgecodes[fi][1] = (midpoint_map.find(all_edgecodes[fi][1]) != midpoint_map.end());
        // all_edgecodes[fi][2] = (midpoint_map.find(all_edgecodes[fi][2]) != midpoint_map.end());
    }


    vectorized_faces f{boost::extents[5][3]};
    clog << all_edgecodes[0][0];
    return f;
}
