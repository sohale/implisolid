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
    Lagrangian (add 1 variable for 1 constrain)
    Make optimization for 4 variables
*/

class Egg_Implicit_Min_Z_Lagrangian : public Egg_Implicit{
public:
    Egg_Implicit_Min_Z_Lagrangian(const column_vector& p) 
    : Egg_Implicit(p)
    {
        ;
    }

    double operator() ( const column_vector& m) const
    {
        const double x = m(0); 
        const double y = m(1);
        const double z = m(2);
        const double d = m(3);

        double egg_i = Egg_Implicit::operator()(m);

        double min_z_egg_i = z + d * squared(egg_i);
           
        return min_z_egg_i;      
    }
};

class Egg_Gradient_Min_Z_Lagrangian : public Egg_Gradient {
public:
    Egg_Gradient_Min_Z_Lagrangian(const column_vector& p) 
    : Egg_Gradient(p) {
        
    }

    const column_vector operator() ( const column_vector& m) const
    {
        column_vector p(6);
        p = a, b, c, center_x, center_y, center_z;

        const double x = m(0);
        const double y = m(1);
        const double z = m(2);
        const double d = m(3);

        column_vector min_z_egg_g(4);
        Egg_Implicit i(p);
        double egg_i = i(m); 

        const column_vector egg_g = Egg_Gradient::operator()(m);

        min_z_egg_g(0) = 2 * d * egg_g(0);
        min_z_egg_g(1) = 2 * d * egg_g(1);
        min_z_egg_g(2) = 1 + 2 * d * egg_g(2);
        min_z_egg_g(2) = d * squared(egg_i);

        return min_z_egg_g;    
    }
};

int main()
{
    try
    {
        column_vector p(6);
        p = 1, 2, 3, 0, 0, 0;

        Egg_Implicit_Min_Z_Lagrangian egg_i(p);
        Egg_Gradient_Min_Z_Lagrangian egg_g(p);
        
        column_vector starting_point(4);
        
        starting_point = 0.1, 0.1, 0, 10;

        auto start = Time::now();
        find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                               objective_delta_stop_strategy(1e-7),
                                               egg_i, starting_point, -1);

        auto end = Time::now();
        fsec duration = end - start;
        std::cout << "Duration : " << duration.count() << endl;
        cout << starting_point << endl;

        starting_point = 0.1, 0.1, 0, 10;

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

