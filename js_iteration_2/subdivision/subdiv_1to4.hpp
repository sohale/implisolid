#include "../v2v_f2f.hpp"

/**
 * @brief      { Subdivides the requested triangles into 4 triangles}
 *
 * @param[in]  faces_old               The old faces
 * @param[in]  verts_old               The old vertexes
 * @param[in]  requested_face_indices  set of triagles requested to be subdivided
 * @param[in]  midpoint_map            The map of the edges that are already subdivided
 *
 * @return     { tuple<> of the new faces and new verts }
 */
auto subdivide_multiple_facets(
    const vectorized_faces & faces_old,
    const vectorized_vect & verts_old,
    const std::set<edge_pair_type>& requested_face_indices,
    const std::map<edge_pair_type, vectorized_vect::index>& midpoint_map
    )
{
    /*!
     * Data fate:
     *     input: old_verts
     *     input: old_faces
     * contiguous: new_verts = copy(old_verts) + extra*3  (preacclocate)
     * output: new faces = from vector<>  extra_faces
     * int extra = |requested|*(4-1)
     * ...
     * then vector<> extra_faces, is built up.
     * assert extra_faces's size willbe "extra".
     * then extra_faces is merged with old_faces and are merged into a contiguous new_faces.
     *
     * transformation: "xyz" produces (inserts_into) predubdivision_edges
     *
     * "the concat" := concat(old_faces, extra_faces)
     * "the trimed" := trim(new_verts)
     *
     * output: new verts = "the trimmed"
     * output: new_faces = contiguous from "the concat"   # from content
     * output: vector<> predubdivision_edges  # source = transformation "xyz"
     */
    auto expected_number_of_new_verts = requested_face_indices.size() * (6-3);
    std::vector<REAL> updatable_vertices(
        verts_old.shape()[0] * 3 + 3 * expected_number_of_new_verts);  // 9 elements for each: xyz for 3 points.
    copy_vertices_into_stdvector_incomplete(verts_old, updatable_vertices);  // leaves it half empty

    std::vector<edge_pair_type> predubdivision_edges;

    typedef tuple<vertexindex_type,vertexindex_type,vertexindex_type> face_triple;
    std::vector<face_triple> new_faces;
    // capacity, expected_number_of_new_faces, provisional_new_facets_count
    auto expected_number_of_new_faces = requested_face_indices.size() * (4-1);
    new_faces.reserve(expected_number_of_new_faces);

    // trace_subdivided_facets
    // if (trace_subdivided_facets) trace_subdivided_facets->resize();

    vectorized_vect new_vertices;

    return std::make_tuple(new_vertices, new_faces, predubdivision_edges);
}
