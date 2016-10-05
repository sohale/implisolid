namespace mp5_implicit {
namespace vectorised_algorithms {

void cross_product(const vectorized_vect& A, const vectorized_vect& B, vectorized_vect &C) {
    assert(A.shape()[0] == B.shape()[0]);
    assert(A.shape()[1] == B.shape()[1]);
    assert(A.shape()[1] == 3);
    assert(B.shape()[1] == 3);
    auto ia = std::cbegin(A);
    auto ea = std::cend(A);
    auto ib = std::cbegin(B);
    auto ic = std::begin(C);
    for (; ia < ea; ++ia, ++ib, ++ic) {
        // todo(sohail): compare performance between: constexpr auto& versus const and direct use
        const auto x1 = (*ia)[0];
        const auto y1 = (*ia)[1];
        const auto z1 = (*ia)[2];

        const auto x2 = (*ib)[0];
        const auto y2 = (*ib)[1];
        const auto z2 = (*ib)[2];

        auto& x = (*ic)[0];
        auto& y = (*ic)[1];
        auto& z = (*ic)[2];

        x = y1 * z2 - z1 * y2;
        y = z1 * x2 - x1 * z2;
        z = x1 * y2 - y1 * x2;
    }
}

}
}
