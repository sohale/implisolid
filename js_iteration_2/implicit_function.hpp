class implicit_function {

public:
    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) = 0;
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) = 0;

protected:
    virtual bool integrity_invariant() = 0;

public:
    virtual ~implicit_function() { };  // ?
};

//trait:
//SignedDistanceImplicitPointwise, PrimitiveBase
