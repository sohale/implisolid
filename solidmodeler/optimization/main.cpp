#include <dlib/optimization.h>
#include <iostream>

using namespace std;
using namespace dlib;

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

    double res = 1 - norm_squared((x - center_x)/a, (y - center_y)/b, (z - center_z)/c);
       
    return res;
}

const column_vector egg_eval_gradient  (const column_vector& m)
{
    const double x = m(0);
    const double y = m(1);
    const double z = m(2);

    const double a2 = squared(a);
    const double b2 = squared(b);
    const double c2 = squared(c);

    column_vector res(3);

    res(0) = -2. * (x - center_x)/a2;
    res(1) = -2. * (y - center_y)/b2;
    res(2) = -2. * (z - center_z)/c2;

    return res;
}

double min_z_egg_eval_implicit (const column_vector& m)
{
    const double x = m(0); 
    const double y = m(1);
    const double z = m(2);

    double res = egg_eval_implicit(m);
    res = my_abs(res);

    double res_z = z + huge_number * res;
       
    return res_z;
}

const column_vector min_z_egg_eval_gradient  (const column_vector& m)
{
    const double x = m(0);
    const double y = m(1);
    const double z = m(2);

    column_vector res_grad_z(3);

    double res = egg_eval_implicit(m); 
    double sign = 1.;
    
    if (res < 0) {
        sign = -1.;
    }

    const column_vector res_grad = egg_eval_gradient(m);

    res_grad_z(0) = sign * huge_number * res_grad(0);
    res_grad_z(1) = sign * huge_number * res_grad(1);
    res_grad_z(2) = 1 + sign * huge_number * res_grad(2);

    return res_grad_z;
}


int main()
{
    try
    {

        column_vector starting_point(3);
        /*
        starting_point = 1, 0, 0;

        cout << egg_eval_implicit(starting_point) << endl;
        cout << egg_eval_gradient(starting_point) << endl;
*/
        starting_point = 0, 0, 0;
/*
        cout << egg_eval_implicit(starting_point) << endl;
        cout << egg_eval_gradient(starting_point) << endl;


        starting_point = 1, 1, 1;

        cout << egg_eval_implicit(starting_point) << endl;
        cout << egg_eval_gradient(starting_point) << endl;
*/
        find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                               objective_delta_stop_strategy(1e-7),
                                               min_z_egg_eval_implicit, starting_point, -1);

        cout << starting_point << endl;

        starting_point = 0, 0, 0;
        find_min(bfgs_search_strategy(),
                 objective_delta_stop_strategy(1e-7),
                 min_z_egg_eval_implicit, min_z_egg_eval_gradient, starting_point, -1);

        cout << starting_point << endl;
    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;
    }
}

