#include "../primitives.cpp"
#include "../timer.hpp"

void immutable_test(const vectorized_vect& v){
    // IS multi_array immutable?
    //v[0][0] = 3.0; // Cannot change!  Unlike Java, it protects.
    //Using the keyword: const in the "body":
    //const_reference operator[](index idx) const {}
}

vectorized_vect make_test_vector(REAL x, REAL y, REAL z) {
    int nsize = 1;
    boost::array<int, 2> values_shape = {{ nsize, 3 }};
    vectorized_vect  values (values_shape);
    values[0][0] = x;
    values[0][1] = y;
    values[0][2] = z;
    return values;

    //vectorized_vect v = vectorized_vect(make_shape_1d({{}}));
    //return vectorized_vect;
}

//make_vectorized_real
auto make_shape_2d(int n, int m) {
    boost::array<int, 2> values_shape = {{ n, m }};
    return values_shape;
}

//unit_sphere s(2.0);

//Sphere s();

// Only use for test
/*
vectorized_x  v = make_test_vector(1,2,3);
vectorized_scalar  f = s.func(v);
vectorized_vect  g = s.grad(v);
*/
//vectorized_scalar  f = s.implicit_func(x);





//vectorized_scalar  f = s.implicit_func(x); //Does it change the contents?
/*
boost::array<int, 2>  x3s = make_shape_2d(num_pts, 3);
vectorized_vect  xx = vectorized_vect(x3s);
*/


bool assertion_side_effect() {
    std::clog << "Side effect: Assertions are activated." << std::endl;
    return true;
}

void test_sphere_one_point(){

    const int num_pts = 1;

    auto sf = make_shape_1d(num_pts);
    vectorized_scalar  f = vectorized_scalar(sf);

    vectorized_vect x = make_test_vector(1.5,2.5,3.5);
    unit_sphere s(2.0);
    //s.implicit_func(x, f);
    s.eval_implicit(x, &f);

    std::clog << f[0] << std::endl;

    implicit_function* s2 = new unit_sphere(2.0);
    s2->eval_implicit(x, &f);

    std::clog << f[0] << std::endl;

    std::clog << "going to delete" << std::endl;

    delete s2;

    // s2 -> ~implicit_function()
}


vectorized_vect  make_empty_x(const int nsize){
    auto sf = make_shape_1d(nsize);
    //vectorized_scalar  f = vectorized_scalar(sf);

    boost::array<int, 2> values_shape = {{ nsize, 3 }};
    vectorized_vect  vectors (values_shape);
    return vectors;
}

//Don't do this. For educational purpose only. This will crash the system.
vectorized_vect&  make_empty_x_2__dontuse(const int nsize){
    //primitives_test.cpp:100:12: warning: reference to stack memory associated with local variable 'values' returned [-Wreturn-stack-address]
    auto sf = make_shape_1d(nsize);
    //vectorized_scalar  f = vectorized_scalar(sf);

    boost::array<int, 2> values_shape = {{ nsize, 3 }};
    vectorized_vect  values (values_shape);
    return values;
}

void  make_empty_x_inplace(const int nsize, vectorized_vect& output){
    auto sf = make_shape_1d(nsize);
    //vectorized_scalar  f = vectorized_scalar(sf);

    boost::array<int, 2> values_shape = {{ nsize, 3 }};
    vectorized_vect  values (values_shape);
    output = values;
}

void test_three_types_of_return_alloc() {

    const int nsize = 10000;
    const int REPEATS = 100;

    std::clog << "make_empty_x";
    timer t1;
    for(int i=0;i<REPEATS;i++) {
        vectorized_vect  x = make_empty_x(nsize);
    }
    t1.stop();

    std::clog << "make_empty_x_2";
    timer t3;
    for(int i=0;i<REPEATS;i++) {
        vectorized_vect   x= make_empty_x_2(nsize);
    }
    t3.stop();

    std::clog << "make_empty_x_inplace";
    timer t2;
    for(int i=0;i<REPEATS;i++) {
        vectorized_vect   x= make_empty_x(nsize);
        //make_empty_x_inplace(nsize);
        make_empty_x_inplace(nsize, x);
    }
    t2.stop();
    /*
make_empty_x
 execution duration: 35.318 msec
make_empty_x_2
 execution duration: 44.6931 msec
make_empty_x_inplace
 execution duration: 3423.99 msec


Using compiler optimisation:

make_empty_x
 execution duration: 3.20714 msec
make_empty_x_2
 execution duration: 134.21 msec
make_empty_x_inplace
 execution duration: 34.6406 msec
!!

After caching:

make_empty_x
 execution duration: 3.01159 msec
make_empty_x_2
 execution duration: 30.6404 msec
make_empty_x_inplace
 execution duration: 28.1285 msec

*/
}




void test_memoryleak_sphere(){

    const int nsize = 10000;

    auto sf = make_shape_1d(nsize);
    vectorized_scalar  f = vectorized_scalar(sf);

    vectorized_vect x = make_empty_x(nsize);

    unit_sphere s(2.0);

    s.eval_implicit(x, &f);

    std::clog << f[0] << std::endl;

    timer t;
    for(int i=0;i<1000*10;i++){
        implicit_function* s2 = new unit_sphere(2.0);
        s2->eval_implicit(x, &f);
        delete s2;
    }
    t.stop();
    //100000 x size(10000), took 5.764 sec. (1 Billion evaluations!)

    // s2 -> ~implicit_function()
}




void grad_test(){
    ;
}

int main() {
    if(0) {
        // Don't call following. It can crash the system.
        test_three_types_of_return_alloc();
    }

    my_assert(assertion_side_effect(), "");
    //immutable_test(x);

    test_sphere_one_point();

    test_memoryleak_sphere();

    std::clog << "Good bye." << std::endl;
}
