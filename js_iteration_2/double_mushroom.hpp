#pragma once

namespace mp5_implicit {

class double_mushroom : public implicit_function {

protected:
    REAL r;
    REAL a; REAL b; REAL c;
    REAL x; REAL y; REAL z;

public:
    double_mushroom(REAL r, REAL a, REAL b, REAL c){
        this->r = r;
        this->a = a;
        this->b = b;
        this->c = c;
        this->x = 0;
        this->y = 0;
        this->z = 0;
        //works with r>3********
    }

    double_mushroom(REAL r, REAL a, REAL b, REAL c, REAL x, REAL y, REAL z){
        this->r = r;
        this->a = a;
        this->b = b;
        this->c = c;
        this->x = x;
        this->y = y;
        this->z = z;
        //works with r>3********
    }

;
    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*f_output)[output_ctr] = -(pow((*i)[0]-this->x,2)/a2+pow((*i)[1]-this->y,2)/b2-pow((*i)[2]-this->z,2)/c2-1);
            if((*i)[2]-this->z > this->r) (*f_output)[output_ctr]=-1;
            if((*i)[2]-this->z < -this->r) (*f_output)[output_ctr]=-1;

        }
    }
      virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);

        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*output)[output_ctr][0] = -2*(*i)[0]/a2;
            (*output)[output_ctr][1] = -2*(*i)[1]/b2;
            (*output)[output_ctr][2] =  2*(*i)[2]/c2;
            // this may be used :  && pow((*i)[0],2)+pow((*i)[1],2)<0.9
            if((*i)[2]-this->z < this->r && (*i)[2]-this->z > this->r-0.05) {
                        (*output)[output_ctr][0] = 0;
                        (*output)[output_ctr][1] = 0;
                        (*output)[output_ctr][2] = -1;}
            if((*i)[2]-this->z > this->r) {
                        (*output)[output_ctr][0] = 0;
                        (*output)[output_ctr][1] = 0;
                        (*output)[output_ctr][2] = -1;}
            if((*i)[2]-this->z > -this->r && (*i)[2]-this->z < -this->r+0.05) {
                        (*output)[output_ctr][0] = 0;
                        (*output)[output_ctr][1] = 0;
                        (*output)[output_ctr][2] = 1;}
            if((*i)[2]-this->z < -this->r) {
                        (*output)[output_ctr][0] = 0;
                        (*output)[output_ctr][1] = 0;
                        (*output)[output_ctr][2] = 1;}
        }
    }
    bool integrity_invariant() const {
        return true;
    }
};

}
