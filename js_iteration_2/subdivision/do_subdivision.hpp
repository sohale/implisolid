#pragma once

#include <vector>
#include <map>
#include <set>
#include "../anecdote.hpp"
#include "../implicit_function/implicit_function.hpp"
#include "../basic_functions.hpp"   // for chisle_random_subset only
#include "subdiv_1to4.hpp"
#include "propagate_subdiv.hpp"
#include "compute_curvetures.hpp"
#include "subdiv_1to2.hpp"

using mp5_implicit::implicit_function;
using mp5_implicit::bad_numbers_in_multi_array;
using ::mp5_implicit::subdivision::propagate_subdiv;
using ::mp5_implicit::subdivision::subdivide_1to2;


namespace mp5_implicit {
namespace subdivision {

// based on mp5-private/solidmodeler/implicit/ohtake_belyaev_5.py


auto subdivide_given_faces (
    const vectorized_faces & old_faces,
    const vectorized_vect & old_verts,
    const std::set<faceindex_type>&  which_facets_set
    /*const implicit_function& iobj,
    REAL curvature_epsilon,
    REAL randomized_probability = 1.0*/
) {

    subdivision::midpointmap_type  midpoint_map;  // empty
    assert(midpoint_map.size() == 0);

    // return std::make_tuple
    // (verts2, faces2, presubdivision_edges);

    auto tuple_ = subdivide_multiple_facets_1to4(old_faces, old_verts, which_facets_set, midpoint_map);
    // verts2, faces2, presubdivision_edges = tuple;

    // get<0> is not used!
    vectorized_vect verts2 = std::move(std::get<0>(tuple_));

    vectorized_faces faces2 = std::move(std::get<1>(tuple_));
    //boost::multi_array<edge_pair_type, 1>
    std::set<edge_pair_type>
        presubdivision_edges;  // = std::get<2>(tuple_);
    // todo: presubdivision_edges is not assigned to.
    auto & psde = std::get<2>(tuple_);
    presubdivision_edges.insert(psde.begin(), psde.end());
    std::move(std::get<2>(tuple_));

    // todo: data-fate of this and the map
    std::vector<edge_pair_type> list_edges_with_1_side(0);  // empty
    assert(list_edges_with_1_side.size() == 0);

    while (true) {

        // propag_dict, edges_which_in1 = propagate_subdiv(faces2, presubdivision_edges)

        /*
        std::set<edge_pair_type>
         //boost::multi_array<edge_pair_type, 1>   // ?
           requested_1side_edgecode_set = presubdivision_edges;

           // note: the definition of requested_1side_edgecode_set has changed since the Python version.
        tuple_3= propagate_subdiv ( faces2, requested_1side_edgecode_set);  // why not requested_1side_edgecode_set
        //*/
        // tuple : propag_dict, edges_need_subdivision
        //tuple(std::map<int, faces_subset_type>, boost::multi_array<edge_pair_type, 1>)

        auto tuple_3 = propagate_subdiv ( faces2, presubdivision_edges);

        // faceindices and edgeindices of leftover edges. The faceindices are categorized based on the number of sides.
        // Maybe not all need to be separate.
        // Future improvement: To avoid separating and joining them back we can have have an entry "2or3" and another entry "1", instead of key entried for "1", "2" & "3".

        auto propag_dict = std::move(std::get<0>(tuple_3));
        auto edges_which_in1 = std::move(std::get<1>(tuple_3)); // check what it is

        //list_edges_with_1_side += [std::move(edges_which_in1)]
        list_edges_with_1_side.insert(list_edges_with_1_side.end(), edges_which_in1.begin(), edges_which_in1.end());

        if (propag_dict[2].shape()[0] == 0 && propag_dict[3].shape()[0] == 0) {
            // nothing to propagate
            break;
        }

        /*
        subdivision::faces_subset_type
            facets_with_2_or_3_sides = std::move(propag_dict[2]);
        */
        std::set<faceindex_type> facets_with_2_or_3_sides;
        facets_with_2_or_3_sides.insert(propag_dict[2].begin(), propag_dict[2].end());

        auto& p3 = propag_dict[3];
        assert(p3.begin() == p3.end());
        // facets_with_2_or_3_sides.insert(p3.begin(), p3.end());

        auto tuple_2 = subdivide_multiple_facets_1to4 (
            faces2,
            verts2,
            facets_with_2_or_3_sides,
            midpoint_map);

        // replace
        // replace(verts2_, std::move(std::get<0>(tuple_2)));

        // should we do a move on a ref (LVALUE reference) ?
        auto& verts2_ = std::get<0>(tuple_2);
        auto& faces2_ = std::get<1>(tuple_2);
        auto& old_edges2 = std::get<2>(tuple_2); // std::move()
        //replace them back
        // faces2 = faces2_;
        // verts2 = verts2_;
        verts2.resize(boost::extents[verts2_.shape()[0]][verts2_.shape()[1]]);
        verts2 = verts2_;
        faces2.resize(boost::extents[faces2_.shape()[0]][faces2_.shape()[1]]);
        faces2 = faces2_;

        presubdivision_edges.insert(old_edges2.begin(), old_edges2.end());
        break;
    }
    // Now we are finished with the triangles with 2-3 sides.

    //why not save THAT as vector? (see above)

    std::set<edge_pair_type>  requested_1side_edgecode_set;
    // edges_with_1_side_set
    // requested_1side_edgecode_set

    requested_1side_edgecode_set.insert(list_edges_with_1_side.begin(), list_edges_with_1_side.end());

    auto facets2 = subdivide_1to2(faces2, requested_1side_edgecode_set, midpoint_map, true);

    /*vectorized_faces subdivide_1to2(const vectorized_faces & faces,
    const std::set<edge_pair_type>& requested_1side_edgecode_set,
    const midpointmap_type& midpoint_map,
    const edge_pair_type EdgecodeBase,
    bool careful_for_twosides=true)
    */

    // return facets2;
    return std::make_tuple(facets2, verts2);
}


/**
 * @brief      { subdivide the triangles that have large curvature. }
 *
 * @param[in]  old_faces               The given faces
 * @param[in]  old_verts               The given vertexes
 * @param[in]  iobj                    The iobj
 * @param[in]  curvature_epsilon       The curvature threshold
 * @param[in]  randomized_probability  1.0 => subdivide all, 0.5=> subdivide half, randomly select.
 *
 * @return     { tuple(facets, verts) }
 */
auto do_subdivision (
    const vectorized_faces & old_faces,
    const vectorized_vect & old_verts,
    const implicit_function& iobj,
    REAL curvature_epsilon,
    REAL randomized_probability = 1.0
) {
    REAL  EXPECTED_SUBDIVISION_PERCENTAGE = 10.0;

    assert(! bad_numbers_in_multi_array(old_verts));

    //check NaN in old_verts. Can fail.
    vectrorized_real  curvatures = ::mp5_implicit::subdivision::
        compute_facets_subdivision_curvatures(old_faces, old_verts, iobj);
    //check NaN in curvatures.
    // Convert NaN into 0.0


    /*
    vectrorized_faceindex  which_facets = find_greater_than(curvatures);

    std::transform(
        curvatures.begin(), curvatures.end(),
        which_facets.begin(), [](REAL curvature) { return curvature > curvature_epsilon;});
    */
    int expected_subdivisions = old_faces.shape()[0] * (EXPECTED_SUBDIVISION_PERCENTAGE) / 100.0;
    std::vector<faceindex_type> which_facets(0);  // zero elements, non-zero capacity
    which_facets.reserve(expected_subdivisions);

    /*
    std::transform(
        curvatures.begin(), curvatures.end(),
        std::back_inserter(name_sizes),
        [](REAL& curvature) { if (curvature > curvature_epsilon) return &curvature-begin });
    */
    /*
    for (faceindex_type fi = 0; fi < old_faces.shape()[0]; ++fi) {
        if (curvature_epsilon[fi] > curvature_epsilon) {
            which_facets.push_back(fi);
        }
    }
    */
    {
        faceindex_type  fi = 0;
        auto it = curvatures.begin(), e = curvatures.end();
        for (; it != e; ++it, ++fi) {
            if (*it > curvature_epsilon) {
                which_facets.push_back(fi);
            }
            cout << "  fi=" << fi << std::endl;
        }
        cout << fi << " == " << old_faces.shape()[0] << std::endl;
        cout << "[] " << curvatures.shape()[0] << " == " << old_faces.shape()[0] << std::endl;
        assert (fi == old_faces.shape()[0]);
        assert (curvatures.shape()[0] == old_faces.shape()[0]);
    }

    if (randomized_probability < 1.0) {
        assert(randomized_probability >= 0.0);
        assert(randomized_probability <= 1.0);
        int n0 = which_facets.size();
        int m0 = static_cast<int>(std::ceil(randomized_probability * (REAL)n0 ));

        chisle_random_subset(which_facets, m0);

        /*
        // choose_randomly()
        // which_facets = choose_random_subset(which_facets, randomized_probability);
        which_facets = choose_random_subset(which_facets, m0);
        */
    }

    std::set<faceindex_type>  which_facets_set;
    // which_facets_set.reserve(...);
    for (auto i = which_facets.begin(); i != which_facets.end(); ++i) {
        which_facets_set.insert(*i);
    }
    if (which_facets_set.size() == 0) {
        anecdote("Number of subdivided triangles is zero. Nothing to do. \n");
        anecdote_c(which_facets.begin(), which_facets.end());
    }


    anecdote("Going to subdivide the following face(s): ");
    anecdote_c(which_facets_set.begin(), which_facets_set.end());
    anecdote("\n");

    return subdivide_given_faces (
        old_faces,
        old_verts,
        which_facets_set
    );
}

/*
 * Theorem: a face is never subdivided more than once. A face in the original mesh, if it is subdivided, the resulting faces are never subdivided.
 */

}  // namespace subdivision
}  // namespace mp5_implicit

