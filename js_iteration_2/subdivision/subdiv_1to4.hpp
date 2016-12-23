#pragma once

#include <vector>
#include <set>
#include <map>
#include "Eigen/Core"
#include "../v2v_f2f.hpp"
#include "../mesh_algorithms.hpp"


using mp5_implicit::easy_edge;
using mp5_implicit::subdivision::midpointmap_type;
// bool VERBOSE_SUBDIV = false;


#define USE_PSDE true
/**
 * @brief      { Subdivides the requested triangles into 4 triangles}
 *
 * @param[in]  old_faces               The old faces
 * @param[in]  old_verts               The old vertexes
 * @param[in]  requested_face_indices  set of triagles requested to be subdivided
 * @param[in]  midpoint_map            The map of the edges that are already subdivided. Updated by this algorithm.
 *
 * @return     { tuple<> of the new faces, the new verts, and PSDE }
 * PSDE: presubdivision_edges will contain all the edges
 * of all the requested triangles, as they were before the
 * subdivision. It is a list/array of (pair of vertices). It will be
 * used for ... .
 * Why do we need this although we have midpoint_map?
 * renamed: subdivide_multiple_facets -> subdivide_multiple_facets_1to4
 */
auto subdivide_multiple_facets_1to4 (
    const vectorized_faces & old_faces,
    const vectorized_vect & old_verts,
    const std::set<faceindex_type>& requested_face_indices,
    midpointmap_type& midpoint_map
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

    /*
    auto expected_number_of_new_verts = requested_face_indices.size() * (6-3);
    std::vector<REAL> updatable_vertices(
        old_verts.shape()[0] * 3 + 3 * expected_number_of_new_verts
    );  // 9 elements for each: xyz for 3 points.
    copy_vertices_into_stdvector_incomplete(old_verts, updatable_vertices);  // leaves it half empty
    */



    /*
    typedef tuple<vertexindex_type,vertexindex_type,vertexindex_type> face_triple;

    std::vector<face_triple> new_faces__additional;
    // new_faces__additional.reserve(exact_total_number_of_new_faces);
    new_faces.resize(exact_total_number_of_new_faces);
    expected_number_of_new_faces  ...   requested_face_indices.size(); // * (4-1);
    */
    /*
    // capacity, exact_total_number_of_new_faces, provisional_new_facets_count
    // The exactly size is known
    const auto exact_total_number_of_new_faces = requested_face_indices.size(); // * (4-1);
    vectorized_faces   new_faces  {boost::extents[exact_total_number_of_new_faces][3]};
    */


    //new_faces_effective_size = old_faces.shape()[0];

    /*
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
    */

    // trace_subdivided_facets
    // if (trace_subdivided_facets) trace_subdivided_facets->resize();

    const auto req_count = requested_face_indices.size();

    #if USE_PSDE
    /**
     *   third output.
     *   Keeps a copy of e01,e12,e20
     */

    // Why does it have two dimensions?
    // std::vector<edge_pair_type> presubdivision_edges;
    boost::multi_array<edge_pair_type, 1> presubdivision_edges {boost::extents[req_count * 3]};
    // presubdivision_edges
    //int presubdivision_edges_counter = 0;
    #endif

    // what does it count? faces or vertices?
    int req_fcounter = 0;  // i.e. presubdivision_edges_counter
    int fcounter = 0;

    // wasting space with hope of speed!
    vectorized_bool_2d  insert_which_verts {boost::extents[req_count][3]};  // opposite of avoid_which[], simliar to use_which[]
    vectorized_faces  actual_vertexindices_used {boost::extents[req_count][3]};
    vectorized_faces  original_vertexindices_used {boost::extents[req_count][3]};

    boost::multi_array<edge_pair_type, 1> edge_triplet_buffer{boost::extents[3]};  // i.e. e0e1e2

    // indices of request faces in the original one. i.e. a vectorised version of requested_face_indices for speed.
    boost::multi_array<vectorized_faces::index, 1> fia {boost::extents[ req_count ]};

    // will be updated
    vertexindex_type  new_verts_effective_total_size = old_verts.shape()[0];

    #if USE_ASSERT
        int redundancy_counter = 0;
    #endif

    if (VERBOSE_SUBDIV) {
        cout << "Requested_face_indices {";
        for (auto fi : requested_face_indices) {
            cout << fi << ", ";
        }
        cout << "}. req_count=" << req_count << std::endl;  // req_count
    }


    // req_count
    for (auto fi : requested_face_indices) {

        /*
        const auto v0 = new_faces[fi][0];
        const auto v1 = new_faces[fi][1];
        const auto v2 = new_faces[fi][2];
        */
        assert(fi < old_faces.shape()[0]);
        /*
        const auto v0 = std::get<0>(new_faces[fi]);
        const auto v1 = std::get<1>(new_faces[fi]);
        const auto v2 = std::get<2>(new_faces[fi]);
        */
        const auto v0 = old_faces[fi][0];
        const auto v1 = old_faces[fi][1];
        const auto v2 = old_faces[fi][2];
        //std::cout << "fi = " << fi << std::endl;
        //std::cout << "v012 " << v0 << " " << v1 << " " << v2 << std::endl;

        // fia[++fcounter] = fi;  // no need to store fi
        fia [fcounter] = fi;

        original_vertexindices_used [fcounter][0] = v0;
        original_vertexindices_used [fcounter][1] = v1;
        original_vertexindices_used [fcounter][2] = v2;
        ++fcounter;

        const edge_pair_type e01 = easy_edge(v0, v1);
        const edge_pair_type e12 = easy_edge(v1, v2);
        const edge_pair_type e20 = easy_edge(v2, v0);

        // The three edges codes
        edge_triplet_buffer[0] = e01;
        edge_triplet_buffer[1] = e12;
        edge_triplet_buffer[2] = e20;

        #if USE_PSDE
            presubdivision_edges[req_fcounter*3    ] = edge_triplet_buffer[0];
            presubdivision_edges[req_fcounter*3 + 1] = edge_triplet_buffer[1];
            presubdivision_edges[req_fcounter*3 + 2] = edge_triplet_buffer[2];
            for (int ai = 0; ai < 3; ++ai) {
                assert (
                    presubdivision_edges[req_fcounter*3 + ai]
                    == edge_triplet_buffer[ai]);   // should fail
            }

            /*
            presubdivision_edges[req_fcounter][0] = e01;
            presubdivision_edges[req_fcounter][1] = e12;
            presubdivision_edges[req_fcounter][2] = e20;
            */
            // req_fcounter == presubdivision_edges_counter
            // ++presubdivision_edges_counter;
        #endif

        // new_faces_effective_size <--> last_in_new_faces

        // INCORRECT:
        // or, can be merged
        vertexindex_type newvtx_counter = new_verts_effective_total_size; // this keeps updated

        // 01, 12, 20
        for (unsigned char vii = 0; vii < 3; ++vii) {
            if (VERBOSE_SUBDIV)
            cout << "edge_triplet_buffer[vii]=" << edge_triplet_buffer[vii] << " ";
            // auto& midpoint_aslookedup = midpoint_map[edge_triplet_buffer[vii]];  // inserts if it doesnt exist
            // note that midpoint_map is updated inside the loop
            const auto finder = midpoint_map.find(edge_triplet_buffer[vii]);
            const bool edge_in_map = (finder != midpoint_map.end());


            //cout << " [vii=" << (int)vii << "] ";
            if (edge_in_map) {
                const vertexindex_type  midpoint_aslookedup = finder->second;
                #if USE_ASSERT
                    ++redundancy_counter;
                #endif

                insert_which_verts [req_fcounter][vii] = false;
                actual_vertexindices_used [req_fcounter][vii] = midpoint_aslookedup;

            } else {

                insert_which_verts [req_fcounter][vii] = true;
                actual_vertexindices_used [req_fcounter][vii] = newvtx_counter;
                //actual_vertexindices_used[requested_triangle_index][side]

                // update the midpoint_map
                midpoint_map [edge_triplet_buffer[vii]] = newvtx_counter;

                //cout << " [*] ";
                #if VERBOSE_SUBDIV

                vector<std::string> namemap {"m01", "m12", "m20"};
                cout << "Inserting edge edge_triplet_buffer[" << static_cast<int>(vii) << "]  " <<
                   "into map[" << edge_triplet_buffer[vii] << "] = vertex " <<newvtx_counter << " : " << namemap[vii] << std::endl;
                #endif
                ++newvtx_counter;  // idx_counter
            }

            //auto what = actual_vertexindices_used [req_fcounter];
        }

        ++ req_fcounter;

        assert (fcounter == req_fcounter);

        // was missing. why?
        new_verts_effective_total_size = newvtx_counter;
        // consolidate it.
    }
    assert(fcounter == req_count);
    assert(req_fcounter == req_count);


    // actual_vertexindices_used ---> used for adding new faces
    // insert_which_verts   ---> used for adding new vertices
    // midpoint_map is updated

    assert(actual_vertexindices_used.shape()[0]
        == insert_which_verts.shape()[0]);

    //vertexindex_type
    //    newvtx_counter = new_verts_effective_total_size;

    if (VERBOSE_SUBDIV)
    cout << "-------------------- new_verts_effective_total_size  ============ " << new_verts_effective_total_size << std::endl;
    //vectorized_vect new_additional_vertices {newvtx_counter};
    vectorized_vect
        new_vertices {boost::extents[new_verts_effective_total_size][3]}; // size includes the older vertices

    //#copy old_vertices new_vertices

    for (vectorized_vect::index vi = 0, e = old_verts.shape()[0]; vi < e; ++vi) {
        new_vertices[vi] = old_verts[vi];
    }

    //boost::multi_array<REAL, 2>
    //    new_3_midpoint_vertices {boost::extents[3][3]};

    //new_3_midpoint_vertices[0][0] = old_verts[0][];
    //new_3_midpoint_vertices[0][0] = old_verts[0][];
    constexpr int NEW_VERT_COUNT3 = 3;
    typedef Eigen::Matrix<REAL, 3, 1> vector3d;
    typedef Eigen::Matrix<REAL, 3, 3> v3x3;
    Eigen::Matrix<REAL, NEW_VERT_COUNT3, 3> new_vert_maker;
    const REAL H_ = 0.5, F_ = 1.0, O_ = 0.0;
    // new_vert_maker << H_,H_,O_,  O_,H_,H_,  H_,O_,H_;  // simple assignment
    new_vert_maker << H_,O_,H_,  H_,H_,O_,  O_,H_,H_;  // simple assignment
    //  new_vert_maker = new_vert_maker.transpose();
    // new_vert_maker << F_,O_,O_,  O_,F_,O_,  O_,O_,F_;  // simple assignment
    v3x3 triangle_vertices;
    //triangle_vertices << ;
    // new_vert_maker * triangle_vertices

    {
    const int nt = insert_which_verts.shape()[0];
    assert (req_count == nt);
    }

    int nv_ctr = 0;  // index just in the new vertices only (probably not needed)
    int nvctr_total = old_verts.shape()[0];  // global (total) index
    for (int req_fii = 0; req_fii < req_count ; ++req_fii) {
        // auto n1 = newvtx_counter;
        //
        //auto fi = fia[req_fii];
        auto v0 = original_vertexindices_used[req_fii][0];
        auto v1 = original_vertexindices_used[req_fii][1];
        auto v2 = original_vertexindices_used[req_fii][2];

        triangle_vertices <<
            old_verts[v0][0], old_verts[v1][0], old_verts[v2][0],
            old_verts[v0][1], old_verts[v1][1], old_verts[v2][1],
            old_verts[v0][2], old_verts[v1][2], old_verts[v2][2];
        // triangle_vertices = triangle_vertices.transpose();

        // not all will be used:
        //boost::multi_array<REAL, 2>
        //    new_3_midpoint_vertices {boost::extents[3][3]};
        // new_3_midpoint_vertices  // new_midpoints
        //    new_3_midpoint_vertices;
        // v3x3 new_3_midpoint_vertices = new_vert_maker * triangle_vertices;
        v3x3 new_3_midpoint_vertices = triangle_vertices * new_vert_maker;

        for (unsigned char vii = 0; vii < 3; ++vii) {
            if (insert_which_verts[req_fii][vii]) {
                // bug fixed.
                new_vertices[nvctr_total][0] = new_3_midpoint_vertices(0, vii);
                new_vertices[nvctr_total][1] = new_3_midpoint_vertices(1, vii);
                new_vertices[nvctr_total][2] = new_3_midpoint_vertices(2, vii);

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
    if (VERBOSE_SUBDIV) {
        cout << nv_ctr << " <= " <<  req_count*3  << "?, len(newverts)=" << new_vertices.shape()[0] << std::endl;
    }
    //assert(nv_ctr + nt == nvctr_total);  // not needed
    //assert(nv_ctr + nt == new_vertices.shape()[0]);  // no!
    //nv_ctr = number of certices added, which is < nt*3
    // assert(nv_ctr <= nt*3);
    assert(nv_ctr <= req_count*3);

    if (VERBOSE_SUBDIV) {

        cout << "number of new verices:  " << nv_ctr << std::endl;

        cout <<
            requested_face_indices.size()
            << "*" << "(4-1) + " << old_faces.shape()[0]
            << std::endl;
        cout << "req_count " << req_count << std::endl;
        cout << "old_faces.shape()[0] " << old_faces.shape()[0] << std::endl;
    }


    // assert(old_faces.shape()[0] == nt);
    // assert(old_faces.shape()[0] == req_count);
    assert(requested_face_indices.size() == req_count);

    // The exactly size is known
    const auto exact_total_number_of_new_faces =
        requested_face_indices.size()*(4-1) + old_faces.shape()[0];
    vectorized_faces   new_faces  {boost::extents[exact_total_number_of_new_faces][3]};

    assert(requested_face_indices.size() == req_count);

    {
        int ctr = 0;
        for (auto f3 : old_faces ) {
            new_faces[ctr] = f3;
            ctr++;
        }
        assert( ctr == old_faces.shape()[0]);
    }
    if (VERBOSE_SUBDIV) {
        cout << "new_faces: initial part copied" << std::endl;

        cout << "new_faces: size= "  << new_faces.shape()[0] << std::endl;
        cout << "filled by old_faces: "  << old_faces.shape()[0] << std::endl;
        cout << "expected new faces: "  << (new_faces.shape()[0] - old_faces.shape()[0]) << std::endl;
        cout << "loop size: "  << req_count << std::endl;
    }


    #if ASSERT_USED
        const int old_verts_count = old_verts.shape()[0];
    #endif

    int nf_ctr = 0;
    int nfctr_total_index = old_faces.shape()[0];  // global (total) index
    for (int req_fii = 0; req_fii < req_count ; ++req_fii) {  // request index

        auto fi_originalface = fia [req_fii];  // originalface

        auto v0 = original_vertexindices_used [fi_originalface][0];
        auto v1 = original_vertexindices_used [fi_originalface][1];
        auto v2 = original_vertexindices_used [fi_originalface][2];

        // midpoints
        auto m01 = actual_vertexindices_used[req_fii][0];
        auto m12 = actual_vertexindices_used[req_fii][1];
        auto m20 = actual_vertexindices_used[req_fii][2];

        #if VERBOSE_SUBDIV
            cout << "v0, v1, v2: " << v0 << "," << v1 << "," << v2 << std::endl;
            cout << "m01, m12, m20: " << m01 << "," << m12 << "," << m20 << std::endl;
        #endif

        assert(v0 < old_verts_count);
        assert(v1 < old_verts_count);
        assert(v2 < old_verts_count);

        assert(m01 >= old_verts_count);
        assert(m12 >= old_verts_count);
        assert(m20 >= old_verts_count);

        // old_faces:
        new_faces[fi_originalface][0] = m12;
        new_faces[fi_originalface][1] = m20;
        new_faces[fi_originalface][2] = m01;
        if (VERBOSE_SUBDIV) {
            cout << "new_faces: updated the old face" << fi_originalface << "," <<  std::endl;
        }

        new_faces[nfctr_total_index][0] = v0;
        new_faces[nfctr_total_index][1] = m01;
        new_faces[nfctr_total_index][2] = m20;
        ++nfctr_total_index;
        ++nf_ctr;
        if (VERBOSE_SUBDIV) {
            cout << "new_faces: added" << (nfctr_total_index-1) << " (a)," <<  std::endl;
        }

        new_faces[nfctr_total_index][0] = v1;
        new_faces[nfctr_total_index][1] = m12;
        new_faces[nfctr_total_index][2] = m01;
        ++nfctr_total_index;
        ++nf_ctr;
        if (VERBOSE_SUBDIV) {
            cout << "new_faces: added" << (nfctr_total_index-1) << " (b)," <<  std::endl;
        }

        new_faces[nfctr_total_index][0] = v2;
        new_faces[nfctr_total_index][1] = m20;
        new_faces[nfctr_total_index][2] = m12;
        ++nfctr_total_index;
        ++nf_ctr;
        if (VERBOSE_SUBDIV) {
            cout << "new_faces: added" << (nfctr_total_index-1) << " (c)." <<  std::endl;
        }

        /*

                    v0         .
                    .          .
                   /.\         .
                  /___\        .
                m01   m20      .
                /\     /\      .
               /  \   /  \     .
              /    \ /    \    .
            v1-----m12-----v2  .
        */


    }

    if (VERBOSE_SUBDIV) {
        cout << "nfctr_total_index=" << nfctr_total_index << "   exact_total_number_of_new_faces=" << exact_total_number_of_new_faces<<  std::endl;
    }
    assert( nfctr_total_index == exact_total_number_of_new_faces);

    if (VERBOSE_SUBDIV) {
        cout << "-------------------- ============ " << new_vertices.shape()[0] << std::endl;
    }

    // new_faces.resize(boost::extents[old_faces.shape()[0]][3]);
    return std::make_tuple(new_vertices, new_faces, presubdivision_edges);
}
