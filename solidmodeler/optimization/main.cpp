#include <dlib/optimization.h>
#include <iostream>
#include <chrono>
#include "egg.cpp"

using namespace std;
using namespace dlib;

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<float> fsec;

/*
    Optimization tests for egg
*/

const double huge_number = 10000000.;

inline double my_abs(double x) {
    if (x < 0) {
        return -x;
    }
    return x;
}

class Egg_Implicit_Min_Z : public Egg_Implicit{
public:
    Egg_Implicit_Min_Z(const column_vector& p) 
    : Egg_Implicit(p)
    {
        ;
    }

    double operator() ( const column_vector& m) const
    {
        const double x = m(0); 
        const double y = m(1);
        const double z = m(2);

        double egg_i = Egg_Implicit::operator()(m);
        egg_i = my_abs(egg_i);

        double min_z_egg_i = z + huge_number * egg_i;
           
        return min_z_egg_i;      
    }
};

class Egg_Gradient_Min_Z : public Egg_Gradient {
public:
    Egg_Gradient_Min_Z(const column_vector& p) 
    : Egg_Gradient(p) {
        
    }

    const column_vector operator() ( const column_vector& m) const
    {
        column_vector p(6);
        p = a, b, c, center_x, center_y, center_z;

        const double x = m(0);
        const double y = m(1);
        const double z = m(2);

        column_vector min_z_egg_g(3);
        Egg_Implicit i(p);
        double egg_i = i(m); 
        double sign = 1.;
        
        if (egg_i < 0) {
            sign = -1.;
        }

        const column_vector egg_g = Egg_Gradient::operator()(m);

        min_z_egg_g(0) = sign * huge_number * egg_g(0);
        min_z_egg_g(1) = sign * huge_number * egg_g(1);
        min_z_egg_g(2) = 1 + sign * huge_number * egg_g(2);

        return min_z_egg_g;    
    }
};

int main()
{
    try
    {
        column_vector p(6);
        p = 1, 2, 3, 0, 0, 0;

        Egg_Implicit_Min_Z egg_i(p);
        Egg_Gradient_Min_Z egg_g(p);
        
        column_vector starting_point(3);
        
        starting_point = 0, 0, 0;

        auto start = Time::now();
        find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                               objective_delta_stop_strategy(1e-7),
                                               egg_i, starting_point, -1);

        auto end = Time::now();
        fsec duration = end - start;
        std::cout << "Duration : " << duration.count() << endl;
        cout << starting_point << endl;

        starting_point = 0, 0, 0;

        start = Time::now();
        find_min(bfgs_search_strategy(),
                 objective_delta_stop_strategy(1e-7),
                 egg_i, egg_g, starting_point, -1);
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

