class egg : public implicit_function {

protected:
    REAL r;

public:
    egg(REAL r){
        this->r = r;
    }

    //boost::array<int, 2> big_shape = {{ 10000, 3 }};
    //boost::multi_array<REAL, 2> huge_test =  boost::multi_array<REAL, 2>(big_shape);

    void eval_implicit(const vectorized_vect& x, vectorized_scalar& f_output){
        my_assert(assert_implicit_function_io(x, f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL r2 = squared(this->r);
        //auto i = x.begin();
        int output_ctr=0;
        //const vectorized_vect::iterator
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            //f_output[output_ctr] = (*x)[0];
            (f_output)[output_ctr] = r2 - norm_squared((*i)[0], (*i)[1], (*i)[2] )*exp(-pow((*i)[2],2));

        }
    }
    void eval_gradient(const vectorized_vect& x, vectorized_vect& output){
        //(*output) = x;

        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (output)[output_ctr][0] = -2. * (*i)[0];
            (output)[output_ctr][1] = -2. * (*i)[1];
            (output)[output_ctr][2] = -(*i)[2]*2.0*exp(-pow((*i)[2],2)) + (pow((*i)[0],2) + pow((*i)[1],2) + pow((*i)[2],2))*2*(*i)[2]*exp(-pow((*i)[2],2));
        }
    }
    bool integrity_invariant(){
        return true;
    }
};
