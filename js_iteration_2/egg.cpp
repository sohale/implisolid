#pragma once

namespace mp5_implicit {

class egg : public implicit_function {

protected:
    REAL r;

public:
    egg(REAL r){
        this->r = r;
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL r2 = squared(this->r);

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*f_output)[output_ctr] = r2 - norm_squared((*i)[0], (*i)[1], (*i)[2] )*exp(-pow((*i)[2],2));

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {


        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*output)[output_ctr][0] = -2. * (*i)[0];
            (*output)[output_ctr][1] = -2. * (*i)[1];
            (*output)[output_ctr][2] = -(*i)[2]*2.0*exp(-pow((*i)[2],2)) + (pow((*i)[0],2) + pow((*i)[1],2) + pow((*i)[2],2))*2*(*i)[2]*exp(-pow((*i)[2],2));
        }
    }
    bool integrity_invariant() const {
        return true;
    }
};

}
