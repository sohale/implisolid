#include <vector>
#include <set>
#include <map>
#include "Eigen/Core"
#include "../v2v_f2f.hpp"


using mp5_implicit::easy_edge;

#define USE_PSDE true
/**
 * @brief      { Subdivides the requested triangles into 4 triangles}
 *
 * @param[in]  old_faces               The old faces
 * @param[in]  old_verts               The old vertexes
 * @param[in]  requested_face_indices  set of triagles requested to be subdivided
 * @param[in]  midpoint_map            The map of the edges that are already subdivided. Updated by this algorithm.
 *
 * @return     { tuple<> of the new faces and new verts }
 * renamed: subdivide_multiple_facets -> subdivide_multiple_facets_1to4
 */
auto subdivide_multiple_facets_1to4 (
    const vectorized_faces & old_faces,
    const vectorized_vect & old_verts,
    const std::set<faceindex_type>& requested_face_indices,
    std::map<edge_pair_type, vectorized_vect::index>& midpoint_map
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
     * transformation: "xyz" produces (inserts_into) presubdivision_edges
     *
     * "the concat" := concat(old_faces, extra_faces)
     * "the trimed" := trim(new_verts)
     *
     * output: new_verts = "the trimmed"
     * output: new_faces = contiguous from "the concat"   # from content
     * output: vector<> presubdivision_edges  # source = transformation "xyz"
     */
    auto expected_number_of_new_verts = requested_face_indices.size() * (6-3);
    std::vector<REAL> updatable_vertices(
        old_verts.shape()[0] * 3 + 3 * expected_number_of_new_verts
    );  // 9 elements for each: xyz for 3 points.
    copy_vertices_into_stdvector_incomplete(old_verts, updatable_vertices);  // leaves it half empty


    typedef tuple<vertexindex_type,vertexindex_type,vertexindex_type> face_triple;

    std::vector<face_triple> new_faces;
    // capacity, expected_number_of_new_faces, provisional_new_facets_count
    auto expected_number_of_new_faces = requested_face_indices.size() * (4-1);
    new_faces.reserve(expected_number_of_new_faces);
    //new_faces_effective_size = old_faces.shape()[0];

    {
        int ctr = 0;
        for (auto f3 : old_faces ) {
            std::get<0>(new_faces[ctr]) = f3[0];
            std::get<1>(new_faces[ctr]) = f3[1];
            std::get<2>(new_faces[ctr]) = f3[2];
            ctr++;
        }
        assert( ctr == old_faces.shape()[0]);
    }


    // trace_subdivided_facets
    // if (trace_subdivided_facets) trace_subdivided_facets->resize();

    auto req_count = requested_face_indices.size();

    #if USE_PSDE
    /**
     *   third output.
     *   Keeps a copy of e0,e1,e2
     */
    // std::vector<edge_pair_type> presubdivision_edges;
    boost::multi_array<edge_pair_type, 2> presubdivision_edges{boost::extents[req_count][3]};
    // presubdivision_edges
    //int presubdivision_edges_counter = 0;
    #endif

    // what does it count? faces or vertices?
    int rcounter = 0;  // i.e. presubdivision_edges_counter
    int fcounter = 0;

    // wasting space with hope of speed!
    vectorized_bool_2d  insert_which_verts {boost::extents[req_count][3]};  // opposite of avoid_which[], simliar to use_which[]
    vectorized_faces  actual_vertexindices_used {boost::extents[req_count][3]};
    vectorized_faces  original_vertexindices_used {boost::extents[req_count][3]};

    boost::multi_array<edge_pair_type, 1> edge_triplet_buffer{boost::extents[3]};  // i.e. e0e1e2

    // vectorized_faces fia
    // will be updated
    vertexindex_type  new_verts_effective_size = old_verts.shape()[0];

    #if USE_ASSERT
        int redundancy_counter = 0;
    #endif

    cout << "requested_face_indices ";
    for (auto fi : requested_face_indices) {
        cout << fi << " ";
    }
    cout << " count: " << requested_face_indices.size() << std::endl;

    for (auto fi : requested_face_indices) {

        /*
        const auto v0 = new_faces[fi][0];
        const auto v1 = new_faces[fi][1];
        const auto v2 = new_faces[fi][2];
        */
        assert(fi < old_faces.shape()[0]);
        const auto v0 = std::get<0>(new_faces[fi]);
        const auto v1 = std::get<1>(new_faces[fi]);
        const auto v2 = std::get<2>(new_faces[fi]);
        //std::cout << "fi = " << fi << std::endl;
        //std::cout << "v012 " << v0 << " " << v1 << " " << v2 << std::endl;

        // fia[++fcounter] = fi;  // no need to store fi
        original_vertexindices_used [fcounter][0] = v0;
        original_vertexindices_used [fcounter][1] = v1;
        original_vertexindices_used [fcounter][2] = v2;
        ++fcounter;

        const edge_pair_type e0 = easy_edge(v0, v1);
        const edge_pair_type e1 = easy_edge(v1, v2);
        const edge_pair_type e2 = easy_edge(v2, v0);

        // The three edges codes
        edge_triplet_buffer[0] = e0;
        edge_triplet_buffer[1] = e1;
        edge_triplet_buffer[2] = e2;

        #if USE_PSDE
            presubdivision_edges[rcounter] = edge_triplet_buffer;
            for (int ai = 0; ai < 3; ++ai) {
                assert (
                    presubdivision_edges[rcounter][ai]
                    == edge_triplet_buffer[ai]);   // should fail
            }

            /*
            presubdivision_edges[rcounter][0] = e0;
            presubdivision_edges[rcounter][1] = e1;
            presubdivision_edges[rcounter][2] = e2;
            */
            // rcounter == presubdivision_edges_counter
            // ++presubdivision_edges_counter;
        #endif

        // new_faces_effective_size <--> last_in_new_faces

        // INCORRECT:
        // or, can be merged
        vertexindex_type newvtx_counter = new_verts_effective_size; // this keeps updated

        for (unsigned char vii = 0; vii < 3; ++vii) {
            cout << "edge_triplet_buffer[vii]=" << edge_triplet_buffer[vii] << " ";
            // auto& midpoint_aslookedup = midpoint_map[edge_triplet_buffer[vii]];  // inserts if it doesnt exist
            // note that midpoint_map is updated inside the loop
            const auto finder = midpoint_map.find(edge_triplet_buffer[vii]);
            const bool edge_in_map = (finder != midpoint_map.end());


            cout << " [vii=" << (int)vii << "] ";
            if (edge_in_map) {
                const vertexindex_type  midpoint_aslookedup = finder->second;
                #if USE_ASSERT
                    ++redundancy_counter;
                #endif

                insert_which_verts [rcounter][vii] = false;
                actual_vertexindices_used [rcounter][vii] = midpoint_aslookedup;

            } else {
                insert_which_verts [rcounter][vii] = true;
                actual_vertexindices_used [rcounter][vii] = newvtx_counter;

                // update the midpoint_map
                midpoint_map [edge_triplet_buffer[vii]] = newvtx_counter;

                cout << " [*] ";

                ++newvtx_counter;  // idx_counter
            }

            auto what = actual_vertexindices_used [rcounter];
        }

        // ++fcounter;

        // was missing. why?
        new_verts_effective_size = newvtx_counter;
        // consolidate it.
    }

    // actual_vertexindices_used ---> used for adding new faces
    // insert_which_verts   ---> used for adding new vertices
    // midpoint_map is updated

    assert(actual_vertexindices_used.shape()[0]
        == insert_which_verts.shape()[0]);

    //vertexindex_type
    //    newvtx_counter = new_verts_effective_size;

    //vectorized_vect new_additional_vertices {newvtx_counter};
    vectorized_vect
        new_vertices {boost::extents[new_verts_effective_size][3]}; // size includes the older vertices

    //boost::multi_array<REAL, 2>
    //    new_3_midpoint_vertices {boost::extents[3][3]};

    //new_3_midpoint_vertices[0][0] = old_verts[0][];
    //new_3_midpoint_vertices[0][0] = old_verts[0][];
    constexpr int NEW_VERT_COUNT3 = 3;
    typedef Eigen::Matrix<REAL, 3, 1> vector3d;
    typedef Eigen::Matrix<REAL, 3, 3> v3x3;
    Eigen::Matrix<REAL, NEW_VERT_COUNT3, 3> new_vert_maker;
    const REAL H_ = 0.5, _F = 1.0, O_ = 0.0;
    new_vert_maker << H_,H_,O_,  O_,H_,H_,  H_,O_,H_;  // simple assignment
    v3x3 triangle_vertices;
    //triangle_vertices << ;
    // new_vert_maker * triangle_vertices

    /*
            O         .
           / \        .
          /___\       .
         M     M      .
        /\     /\     .
       /  \   /  \    .
      /    \ /    \   .
     1------M------2  .

    */

    const int nt = insert_which_verts.shape()[0];
    assert (fcounter == nt);

    cout << "nt:" << nt << std::endl;

    int nv_ctr = 0;  // index just in the new vertices only (probably not needed)
    int nvctr_total = nt;  // global index
    for (int fii = 0; fii < nt ; ++fii) {
        // auto n1 = newvtx_counter;
        //
        //auto fi = fia[fii];
        auto v0 = original_vertexindices_used[fii][0];
        auto v1 = original_vertexindices_used[fii][1];
        auto v2 = original_vertexindices_used[fii][2];

        triangle_vertices <<
            old_verts[v0][0], old_verts[v1][0], old_verts[v2][0],
            old_verts[v0][1], old_verts[v1][1], old_verts[v2][1],
            old_verts[v0][2], old_verts[v1][2], old_verts[v2][2];

        // not all will be used:
        //boost::multi_array<REAL, 2>
        //    new_3_midpoint_vertices {boost::extents[3][3]};
        // new_3_midpoint_vertices  // new_midpoints
        //    new_3_midpoint_vertices;
        v3x3 new_3_midpoint_vertices = new_vert_maker * triangle_vertices;

        for (unsigned char vii = 0; vii < 3; ++vii) {
            if (insert_which_verts[fii][vii]) {
                new_vertices[nvctr_total][0] = new_3_midpoint_vertices(vii, 0);
                new_vertices[nvctr_total][1] = new_3_midpoint_vertices(vii, 1);
                new_vertices[nvctr_total][2] = new_3_midpoint_vertices(vii, 2);

                ++nvctr_total;
                ++nv_ctr;
            }
        }
        // auto n2 = newvtx_counter;

        /*
        for (unsigned char i = n1; i < n2; ++i) {
            new_verts[i] =v456[]
        }
        */
    }
    cout << nv_ctr << " + " <<  nt  << " ==  " << new_vertices.shape()[0] << std::endl;
    //assert(nv_ctr + nt == nvctr_total);  // not needed
    //assert(nv_ctr + nt == new_vertices.shape()[0]);  // no!
    //nv_ctr = number of certices added, which is < nt*3
    assert(nv_ctr <= nt*3);



    return std::make_tuple(new_vertices, new_faces, presubdivision_edges);
}
