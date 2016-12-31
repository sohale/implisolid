#pragma once

#include <emscripten.h>

#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
namespace mp5_implicit {
namespace implicit_functions {

/**
    Suggested names: live_js_implicit_function, asmjs_implicit_funciton, etc
*/
class javascript_implicit_function : public transformable_implicit_function {

protected:
    REAL param1;
    int id;  // callback id

public:

    javascript_implicit_function(int id, REAL matrix[12], REAL param1) {
        // matrix[12]: How to make sure 12 elements are provided? (how to assert)
        this->id = id;
        this->param1 = param1;

        // The "delete" is added to ~destructor
        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++)
            transf_matrix[i] = matrix[i];

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");
    }

    ~javascript_implicit_function() {
        delete this->transf_matrix;
        this->transf_matrix = 0;
        delete this->inv_transf_matrix;
        this->inv_transf_matrix = 0;
    }


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {

        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);

        /*
        int output_ctr=0;

        auto e = x_copy.end();
        for(auto i = x_copy.begin(); i<e; i++, output_ctr++){
            (*f_output)[output_ctr] = 1 - norm_squared(((*i)[0]-this->x0)/this->a, ((*i)[1]-this->y0)/this->b, ((*i)[2]-this->z0)/this->c);
        }
        */
        int result = EM_ASM_INT({
            console.log('Calling implicit function callback: id=', $0, "x_ptr=", $1, " count=", $2, "output_pt =", $3, " param1=", $4);
            js_implcit_callback($0,$1,$2,$3,$4); /* id, ptr, count, output_f_ptr*/
            // $0, Module.HEAPF.subarray($1 >> 2, ($1+$2*3)>>2), Module.HEAPF.subarray($3 >> 2, ($3+$2*3)>>2)
            return 1;
        }, this->id, (void*)(&((*(x_copy.begin()))[0])), x_copy.shape()[0], (void*)(f_output->data()), this->param1);
    }

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        vectorized_vect x_copy = x;
        matrix_vector_product(this->inv_transf_matrix, x_copy);

        int result = EM_ASM_INT({
            console.log('Calling gradient function callback: id=', $0, "x_ptr=", $1, " count=", $2, "output_pt =", $3);
            js_gradient_callback($0,$1,$2,$3); /* id, ptr, count, output_f_ptr*/
            // $0, Module.HEAPF.subarray($1 >> 2, ($1+$2*3)>>2), Module.HEAPF.subarray($3 >> 2, ($3+$2*3)>>2)
            return 1;
        }, this->id, &((*(x_copy.begin()))[0]), x_copy.shape()[0], (void*)(output->data()));
        /*
        int output_ctr=0;
        auto e = x_copy.end();
        for(auto i = x_copy.begin(); i < e; i++, output_ctr++) {
            REAL g0 = -2. * ((*i)[0]-this->x0)/a2;
            REAL g1 = -2. * ((*i)[1]-this->y0)/b2;
            REAL g2 = -2. * ((*i)[2]-this->z0)/c2;

            (*output)[output_ctr][0] = g0;
            (*output)[output_ctr][1] = g1;
            (*output)[output_ctr][2] = g2;
        }
        */

        int output_ctr=0;
        auto e = x_copy.end();
        for(auto i = x_copy.begin(); i < e; i++, output_ctr++) {
            // not efficient:
            REAL g0 = (*output)[output_ctr][0];
            REAL g1 = (*output)[output_ctr][1];
            REAL g2 = (*output)[output_ctr][2];

            (*output)[output_ctr][0] = this->inv_transf_matrix[0]*g0 + this->inv_transf_matrix[4]*g1 + this->inv_transf_matrix[8]*g2;
            (*output)[output_ctr][1] = this->inv_transf_matrix[1]*g0 + this->inv_transf_matrix[5]*g1 + this->inv_transf_matrix[9]*g2;
            (*output)[output_ctr][2] = this->inv_transf_matrix[2]*g0 + this->inv_transf_matrix[6]*g1 + this->inv_transf_matrix[10]*g2;
        }
    }

    bool integrity_invariant() const {
        //todo: std::clog << "not implemented. use EM_ASM_()" << std::endl;
        return true;
    }

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        std::clog << "not implemented. use EM_ASM_()" << std::endl;
        return mp5_implicit::bounding_box{-1, 1, -1, 1, -1, 1};
    }
};

}  // namespace implicit_functions
}  // namespace mp5_implicit
