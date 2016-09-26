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
void subdivide_SD2(facets, edges_with_1_side, midpoint_map, bool careful_for_twosides=true) {
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


vectorized_faces subdivide_1to2(const vectorized_faces & facets,
    const std::vector<edge_pair_type>& edges_with_1_side,
    std::map<edge_pair_type, vectorized_vect::index> midpoint_map,
    bool careful_for_twosides=true)
{
    vectorized_faces f{};
    return f;
}
