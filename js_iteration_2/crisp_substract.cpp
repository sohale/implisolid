class crisp_subtract : public implicit_function {
protected:
    implicit_function a, b;
public:
    crisp_subtract(implicit_function a, implicit_function b) {
        this->a = a;
        this->b = b;
    }

    virtual void eval.implicit(const vectorized_vect& x, vectorized_scalar* f_output){
        vectorized_scalar *f1,*f2;
        this->a.eval_implicit(x, &f1)
        this->b.eval_implicit(x, &f2)
        boost::multi_array<int, 1> bool_mask;
        auto i = x.begin();
        auto e = x.end();
        int output_ctr = 0;
        // create a vectorized truth table of f1 > f2
        for (; i<e; i++, output_ctr++){
            bool_mask[output_ctr] = f1[output_ctr] > f2[output_ctr] ? 1: 0; // bool_mask[i] set to 1 if f1>f2 else 0
            (*f_output)[output_ctr] = f1 * bool_mask[output_ctr] + f2 * (1 - bool_mask[output_ctr])
        }

    }


    virtual void eval_gradient(){

    }
};
