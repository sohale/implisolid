#pragma once


#include "../my_assert.hpp"

namespace mp5_implicit {
namespace implicit_functions {

inline bool assert_implicit_function_io(const vectorized_vect& x, const vectorized_scalar& output){
    //std::clog << x.shape()[1] << " " << x.shape()[0] << " " << output->shape()[0] << std::endl;

    //my_assert(x.shape()[1] == 3, "Size should be N x 3. Not " ); //+ x.shape()[1]);
    my_assert(x.shape()[1] == 3, ""); //, "Size should be N x 3. Not " << x.shape()[1]);
    my_assert(x.shape()[0] == output.shape()[0], "") ; //, "Sizes don't match. Prepare an output using the same size. Sizes: " << x.shape()[0] << " !== " << output.shape()[0] );
    return true;
}

}  // namespace implicit_functions

} // namespace mp5_implicit
