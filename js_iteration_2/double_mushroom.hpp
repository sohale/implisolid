#pragma once

namespace mp5_implicit {

class double_mushroom : public implicit_function {

protected:
    REAL r;

public:
    double_mushroom(REAL r){
        this->r = r;
        //works with r>3********
    }

;
    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL r2 = squared(this->r);

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*f_output)[output_ctr] = -(pow((*i)[0],2)*r2+pow((*i)[1],2)*r2-pow((*i)[2],2)*r2-1);
            if((*i)[2]>0.8) (*f_output)[output_ctr]=-1;
            if((*i)[2]<-0.8) (*f_output)[output_ctr]=-1;

        }
    }
      virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        const REAL r2 = squared(this->r);
        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*output)[output_ctr][0] = -r2 * (*i)[0];
            (*output)[output_ctr][1] = -r2 * (*i)[1];
            (*output)[output_ctr][2] = r2 * (*i)[2];
            if((*i)[2]<0.8 && (*i)[2]>0.78 && pow((*i)[0],2)+pow((*i)[1],2)<0.9) {(*output)[output_ctr][0] = 0;
                        (*output)[output_ctr][1] = 0;
                        (*output)[output_ctr][2] = -1;}
            if((*i)[2]>0.80) {(*output)[output_ctr][0] = 0;
                        (*output)[output_ctr][1] = 0;
                        (*output)[output_ctr][2] = -1;}
            if((*i)[2]>-0.8 && (*i)[2]<-0.78 && pow((*i)[0],2)+pow((*i)[1],2)<0.9) {(*output)[output_ctr][0] = 0;
                        (*output)[output_ctr][1] = 0;
                        (*output)[output_ctr][2] = 1;}
            if((*i)[2]<-0.80) {(*output)[output_ctr][0] = 0;
                        (*output)[output_ctr][1] = 0;
                        (*output)[output_ctr][2] = 1;}
        }
    }
    bool integrity_invariant() const {
        return true;
    }
};

}