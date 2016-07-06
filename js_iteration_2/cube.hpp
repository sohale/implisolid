#pragma once

namespace mp5_implicit {

class cube : public implicit_function {

protected:
    REAL r;

public:
    cube(REAL r){
        this->r = r;
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output, REAL grid_real_size)const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            if (((*i)[0]<=r/2.) && ((*i)[1]<=r/2.) && ((*i)[2]<=r/2.) && ((*i)[0]>=-r/2.) && ((*i)[1]>=-r/2.) && ((*i)[2]>=-r/2.)){
              (*f_output)[output_ctr] = 1.;
            }
            else{

              (*f_output)[output_ctr] = - 1.;
            }

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output, REAL grid_real_size) const {

        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i!=e; i++, output_ctr++){
            if ((*i)[0]<=r/2. - 0.05 && (*i)[0]>= -r/2. + 0.05) {
                  (*output)[output_ctr][0] = +1.* (*i)[0];
            }
            else{
                (*output)[output_ctr][0] = -1. * (*i)[0];
            }
            if ((*i)[1]<=r/2. - 0.05 && (*i)[1]>= -r/2. + 0.05){
                  (*output)[output_ctr][1] = +1.*(*i)[1];
            }
            else{
                (*output)[output_ctr][1] = -1.*(*i)[1];
            }
            if ((*i)[2]<=r/2. - 0.05 && (*i)[2]>= -r/2. + 0.05){
                  (*output)[output_ctr][2] = +1.*(*i)[2];
            }
            else{
                (*output)[output_ctr][2] = -1.*(*i)[2];
            }

        }
    }
    bool integrity_invariant() const {
        return true;
    }
};

}
