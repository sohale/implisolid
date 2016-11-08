#include <dlib/optimization.h>

using namespace dlib;

/*
    Egg
*/

typedef matrix<double,0,1> column_vector;

inline double norm_squared(double x, double y, double z) {
    return x*x + y*y + z*z;
}

inline double squared(double x) {
    return x*x;
}


class Egg_Abstract {
public:
    double a, b, c; // radius_x, radius_y, radius_z
    double center_x, center_y, center_z; // center of egg

    Egg_Abstract(const column_vector& p) {
        a = p(0);
        b = p(1);
        c = p(2);
        center_x = p(3);
        center_y = p(4);
        center_z = p(5);
    }
};


class Egg_Implicit : public Egg_Abstract{
public:
    Egg_Implicit(const column_vector& p) 
    : Egg_Abstract(p) {
        
    }

    double operator() ( const column_vector& m) const
    {
        const double x = m(0); 
        const double y = m(1);
        const double z = m(2);

        double egg_i = 1 - norm_squared((x - center_x)/a, (y - center_y)/b, (z - center_z)/c);
           
        return egg_i;       
    }
};

class Egg_Gradient : public Egg_Abstract{
public:
    Egg_Gradient(const column_vector& p) 
    : Egg_Abstract(p) {
        
    }

    column_vector operator() ( const column_vector& m) const
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
};