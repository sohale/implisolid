namespace mp5_implicit {
namespace vectorised_algorithms {

// THIS WILL NOT WORK!!!
void add_inplace(vectorized_vect A, vectorized_vect B) {
    assert(A.shape()[0] == B.shape()[0]);
    assert(A.shape()[1] == B.shape()[1]);
    for (auto ia = std::begin(A), ea = std::end(A), ib = std::begin(B); ia < ea; ++ia /*, ++ib*/) {
        (*ia)[0] += (*ib)[0];
        (*ia)[1] += (*ib)[1];
        (*ia)[2] += (*ib)[2];
    }
}

}
}
