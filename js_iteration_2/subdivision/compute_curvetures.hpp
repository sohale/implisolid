#include "../basic_data_structures.hpp"

namespace mp5_implicit {
namespace subdivision {


vectrorized_real  compute_facets_subdivision_curvatures (
    const vectorized_faces & faces,
    const vectorized_vect & verts,
    const implicit_function& iobj
) {
    vectrorized_real  curvatures{boost::extents[faces.shape()[0]]};
    for (auto i = curvatures.begin(), e = curvatures.end(); i != e; ++i) {
        *i = 0.00000000000001;
    }
    return curvatures;
}



}  // namespace mp5_implicit
}  // namespace subdivision
