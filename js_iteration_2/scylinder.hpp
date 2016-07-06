#pragma once

namespace mp5_implicit {

class scylinder : public implicit_function {

protected:
    REAL r;

public:
    scylinder(REAL r){
        this->r = r;
    }


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output, REAL grid_real_size) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL r2 = squared(this->r);

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
          (*f_output)[output_ctr] = -(*i)[0]*(*i)[0] - (*i)[1]*(*i)[1] + r2;
          if((*i)[2]>=0.8 || (*i)[2]<=-0.8 ){
            (*f_output)[output_ctr] = -1.;
          }

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output, REAL grid_real_size) const {

        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        const REAL r2 = squared(this->r);
        for(; i!=e; i++, output_ctr++){

          if((*i)[2]>=0.78 || (*i)[2]<=-0.78){

            (*output)[output_ctr][0] = -2. * (*i)[0];
            (*output)[output_ctr][1] = -2. * (*i)[1];
            (*output)[output_ctr][2] = 0.;

          }
          else {
            (*output)[output_ctr][0] = 0.;
            (*output)[output_ctr][1] = 0.;
            (*output)[output_ctr][2] = -1.;
          }

        }
    }
    bool integrity_invariant() const {
        return true;
    }
};

}
