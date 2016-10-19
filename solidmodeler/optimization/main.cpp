#include <dlib/optimization.h>
#include <iostream>
#include <chrono>

using namespace std;
using namespace dlib;

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<float> fsec;

/*
    Optimization tests for egg
*/

const double a = 1., b = 2., c = 3.; // radius_x, radius_y, radius_z
const double center_x = 0., center_y = 0., center_z = 0.; // center of egg

const double huge_number = 10000000.;

typedef matrix<double,0,1> column_vector;

inline double norm_squared(double x, double y, double z) {
    return x*x + y*y + z*z;
}

inline double squared(double x) {
    return x*x;
}

inline double my_abs(double x) {
    if (x < 0) {
        return -x;
    }
    return x;
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

double min_z_egg_eval_implicit (const column_vector& m)
{
    const double x = m(0); 
    const double y = m(1);
    const double z = m(2);

    double egg_i = egg_eval_implicit(m);
    egg_i = my_abs(egg_i);

    double min_z_egg_i = z + huge_number * egg_i;
       
    return min_z_egg_i;
}

const column_vector min_z_egg_eval_gradient  (const column_vector& m)
{
    const double x = m(0);
    const double y = m(1);
    const double z = m(2);

    column_vector min_z_egg_g(3);

    double egg_i = egg_eval_implicit(m); 
    double sign = 1.;
    
    if (egg_i < 0) {
        sign = -1.;
    }

    const column_vector egg_g = egg_eval_gradient(m);

    min_z_egg_g(0) = sign * huge_number * egg_g(0);
    min_z_egg_g(1) = sign * huge_number * egg_g(1);
    min_z_egg_g(2) = 1 + sign * huge_number * egg_g(2);

    return min_z_egg_g;
}


int main()
{
    try
    {

        column_vector starting_point(3);
        
        starting_point = 0, 0, 0;

        auto start = Time::now();
        
        find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                               objective_delta_stop_strategy(1e-7),
                                               min_z_egg_eval_implicit, starting_point, -1);

        auto end = Time::now();
        fsec duration = end - start;

        std::cout << "Duration : " << duration.count() << endl;

        cout << starting_point << endl;

        starting_point = 0, 0, 0;

        start = Time::now();

        find_min(bfgs_search_strategy(),
                 objective_delta_stop_strategy(1e-7),
                 min_z_egg_eval_implicit, min_z_egg_eval_gradient, starting_point, -1);


        end = Time::now();
        duration = end - start;

        std::cout << "Duration : " << duration.count() << endl;

        cout << starting_point << endl;
    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;
    }
}

