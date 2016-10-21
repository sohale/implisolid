#include <dlib/optimization.h>

using namespace dlib;

/*
    Egg
*/

const double a = 1., b = 2., c = 3.; // radius_x, radius_y, radius_z
const double center_x = 0., center_y = 0., center_z = 0.; // center of egg

typedef matrix<double,0,1> column_vector;

inline double norm_squared(double x, double y, double z) {
    return x*x + y*y + z*z;
}

inline double squared(double x) {
    return x*x;
}

double egg_eval_implicit (const column_vector& m)
{
    const double x = m(0); 
    const double y = m(1);
    const double z = m(2);

    double egg_i = 1 - norm_squared((x - center_x)/a, (y - center_y)/b, (z - center_z)/c);
       
    return egg_i;
}

const column_vector egg_eval_gradient  (const column_vector& m)
{
    const double x = m(0);
    const double y = m(1);
    const double z = m(2);

    const double a2 = squared(a);
    const double b2 = squared(b);
    const double c2 = squared(c);

    column_vector egg_g(3);

    egg_g(0) = -2. * (x - center_x)/a2;
    egg_g(1) = -2. * (y - center_y)/b2;
    egg_g(2) = -2. * (z - center_z)/c2;

    return egg_g;
}
