#include <iostream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <map>
#include <string>

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"


#include "../../js_iteration_2/basic_data_structures.hpp"
#include "../../js_iteration_2/basic_functions.hpp"
#include "../../js_iteration_2/configs.hpp"
#include "../../js_iteration_2/implicit_function/primitives.hpp"

#include <nlopt.h>

using namespace std;
using namespace mp5_implicit;


double func_min_z(unsigned n, const double *m, double *grad, void *my_func_data)
{
    if (grad) {
        grad[0] = 0;
        grad[1] = 0;
        grad[2] = 1;
    }         
    return m[2];       
}

 double my_constraint_eq(unsigned n, const double *m, double *grad, void *data)
{
    transformable_implicit_function *f = (transformable_implicit_function *) data;

    vectorized_vect points(boost::extents[1][3]);
    vectorized_vect gradient(boost::extents[1][3]);
	vectorized_scalar result(boost::extents[1]);

	points[0][0] = m[0];
	points[0][1] = m[1];
	points[0][2] = m[2];

	f->eval_implicit(points, &result);
    
    if (grad) {
    	f->eval_gradient(points, &gradient);
    
        grad[0] = gradient[0][0];
        grad[1] = gradient[0][1];
        grad[2] = gradient[0][2];

    }

    return result[0];
 }

 double getRandomFromRange()
 {
 	//from -5 to 5
 	double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
 	r = 10*r - 5;
 	return r;
 }

int main()
{
	egg my_egg(1, 2, 3);

	double lb[3] = { -10, -10, -10 }; 
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LD_SLSQP, 3); 
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, func_min_z, NULL);

	nlopt_add_equality_constraint(opt, my_constraint_eq, &my_egg, 1e-4);

	nlopt_set_xtol_rel(opt, 1e-4);
	nlopt_set_maxeval(opt, 100);
	nlopt_set_xtol_abs1(opt, 1e-4);

	double m[3] = { 1, -1, 0 };  
	double minf; 
	
	double correctZ = -3;
	int counter = 0;
	for(int i = 0; i < 10; i++) {

		m[0] = getRandomFromRange();
		m[1] = getRandomFromRange();
		m[2] = getRandomFromRange();

		std::cout << m[0] <<  " : " << m[1]  <<  " : " << m[2] << std::endl;
		
		if (nlopt_optimize(opt, m, &minf) < 0) {
		    printf("nlopt failed!\n");
		    counter ++;
		}
		else {
			if (abs(correctZ - m[2]) > 1e-4) {
				counter ++;
				std::cout << "Violating " << m[2] << std::endl;
			} else {
				std::cout << "All good " << m[2] << std::endl;
			}		
		}
	}

	std::cout << "All violation " << counter << std::endl;
	
	nlopt_destroy(opt);

	return 0;
}
