class implicit_function {

public:
    void eval_implicit(const vectorized_vect& x, vectorized_scalar* output){};
    void eval_gradient(const vectorized_vect& x, vectorized_vect* output){};

protected:
     bool integrity_invariant(){return true;};

public:
    ~implicit_function() {};  // ?
};

//trait:
//SignedDistanceImplicitPointwise, PrimitiveBase
