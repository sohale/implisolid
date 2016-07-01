class cube : public implicit_function {

protected:
    REAL r;

public:
    cube(REAL r){
        this->r = r;
    }

    //boost::array<int, 2> big_shape = {{ 10000, 3 }};
    //boost::multi_array<REAL, 2> huge_test =  boost::multi_array<REAL, 2>(big_shape);

    void eval_implicit(const vectorized_vect& x, vectorized_scalar& f_output){
        my_assert(assert_implicit_function_io(x, f_output), "");
        my_assert(this->integrity_invariant(), "");

        //auto i = x.begin();
        int output_ctr=0;
        //const vectorized_vect::iterator
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            if (((*i)[0]<=r/2.) && ((*i)[1]<=r/2.) && ((*i)[2]<=r/2.) && ((*i)[0]>=-r/2.) && ((*i)[1]>=-r/2.) && ((*i)[2]>=-r/2.)){
              (f_output)[output_ctr] = 1.;
            }
            else{

              (f_output)[output_ctr] = - 1.;
            }

        }
    }
    void eval_gradient(const vectorized_vect& x, vectorized_vect& output){
        //(*output) = x;

        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i!=e; i++, output_ctr++){
            if ((*i)[0]<=r/2. - 0.2 && (*i)[0]>= -r/2. + 0.2) {
                  (output)[output_ctr][0] = -1.;
            }
            else{
                (output)[output_ctr][0] = 0.;
            }
            if ((*i)[1]<=r/2. - 0.2 && (*i)[1]>= -r/2. + 0.2){
                  (output)[output_ctr][1] = -1.;
            }
            else{
                (output)[output_ctr][1] = 0.;
            }
            if ((*i)[2]<=r/2. - 0.2 && (*i)[2]>= -r/2. + 0.2){
                  (output)[output_ctr][2] = -1.;
            }
            else{
                (output)[output_ctr][2] = 0.;
            }

        }
    }
    bool integrity_invariant(){
        return true;
    }
};
