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
    implicit_function *f = (implicit_function *) data;

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

 double getRandom(REAL deviation)
 {
 	double r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
 	r = 2*deviation*r - deviation;
 	return r;
 }

inline REAL my_min(REAL a, REAL b) {

	if (a > b) {
		return b;
	} else {
		return a;
	}
}

inline REAL my_abs(REAL a) {
	a = a > 0 ? a : -a;
	return a;
}



REAL find_min_z(
	implicit_function& f,  
	int random_starting_point_count, 
	REAL random_starting_point_standard_deviation) {

	REAL min_z = 1000000;

	double lb[3] = { -10, -10, -10 }; 
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LD_SLSQP, 3); 
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, func_min_z, NULL);

	nlopt_add_equality_constraint(opt, my_constraint_eq, &f, 1e-4);

	nlopt_set_xtol_rel(opt, 1e-4);
	nlopt_set_maxeval(opt, 100);
	nlopt_set_xtol_abs1(opt, 1e-4);

	double m[3] = { 0, 0, 0 };  
	double minf; 

	vectorized_vect points(boost::extents[1][3]);
	vectorized_scalar result(boost::extents[1]);

	for(int i = 0; i < random_starting_point_count; i++) {

		m[0] = getRandom(random_starting_point_standard_deviation);
		m[1] = getRandom(random_starting_point_standard_deviation);
		m[2] = getRandom(random_starting_point_standard_deviation);

		std::cout << m[0] <<  " : " << m[1]  <<  " : " << m[2] << std::endl;
		
		if (nlopt_optimize(opt, m, &minf) < 0) {
		    printf("nlopt failed!\n");
		}
		else {
			// checking for f(x,y,z) = 0 for result point

			points[0][0] = m[0];			
			points[0][1] = m[1];
			points[0][2] = m[2];

			f.eval_implicit(points, &result);

			if (my_abs(result[0]) < 1e-2) {			
				cout << "Found min z:" << minf << endl;
				min_z = my_min(min_z, minf);		
			} else {
				cout << "Not close enought, min z:" << minf << " but constaints = "<< result[0] << endl;
			}
		}
	}

	nlopt_destroy(opt);

	cout << "Total Found min z:" << min_z << endl;
	return min_z;
}


void egg_test() 
{
	cout << "=================================" << endl;
	cout << "In egg" << endl;
	
	egg obj(1, 2, 4);
	find_min_z(obj, 5, 5);
	cout << "=================================" << endl;
}

void tetrahedron_test() 
{
	cout << "=================================" << endl;
	cout << "In tetrahedron" << endl;
	
	vector<boost::array<REAL,3>> points;
	
	boost::array<REAL,3> p1;
	boost::array<REAL,3> p2;
	boost::array<REAL,3> p3;
	boost::array<REAL,3> p4;

	p1[0] = 0;
	p1[1] = 0;
	p1[2] = 1;
	
	p2[0] = 4;
	p2[1] = 0;
	p2[2] = 3;
	
	p3[0] = 0;
	p3[1] = 4;
	p3[2] = 2;
	
	p4[0] = 0;
	p4[1] = 0;
	p4[2] = 4;

	points.push_back(p1);
	points.push_back(p2);
	points.push_back(p3);
	points.push_back(p4);

	REAL matrix12[12];
	
	matrix12[0] = 1;
	matrix12[1] = 0;
	matrix12[2] = 0;
	matrix12[3] = 0;
	matrix12[4] = 0;
	matrix12[5] = 1;
	matrix12[6] = 0;
	matrix12[7] = 0;
	matrix12[8] = 0;
	matrix12[9] = 0;
	matrix12[10] = 1;
	matrix12[11] = 0;

	tetrahedron obj(points, matrix12);

	find_min_z(obj, 5, 5);
	cout << "=================================" << endl;
}

void heart_test() 
{
	cout << "=================================" << endl;
	cout << "In heart" << endl;

	REAL matrix12[12];
	
	matrix12[0] = 1;
	matrix12[1] = 0;
	matrix12[2] = 0;
	matrix12[3] = 0;
	matrix12[4] = 0;
	matrix12[5] = 1;
	matrix12[6] = 0;
	matrix12[7] = 0;
	matrix12[8] = 0;
	matrix12[9] = 0;
	matrix12[10] = 1;
	matrix12[11] = 0;

	heart obj(matrix12);

	find_min_z(obj, 5, 5);
	cout << "=================================" << endl;
}

void torus_test() 
{
	cout << "=================================" << endl;
	cout << "In torus" << endl;

	REAL matrix12[12];
	
	matrix12[0] = 1;
	matrix12[1] = 0;
	matrix12[2] = 0;
	matrix12[3] = 0;
	matrix12[4] = 0;
	matrix12[5] = 1;
	matrix12[6] = 0;
	matrix12[7] = 0;
	matrix12[8] = 0;
	matrix12[9] = 0;
	matrix12[10] = 1;
	matrix12[11] = 0;

	torus obj(matrix12);

	find_min_z(obj, 5, 5);
	cout << "=================================" << endl;
}

void subtraction_test() 
{
	cout << "=================================" << endl;
	cout << "In subtraction test, egg- cube " << endl;

	REAL matrix12[12];
	
	matrix12[0] = 1;
	matrix12[1] = 0;
	matrix12[2] = 0;
	matrix12[3] = 0;
	matrix12[4] = 0;
	matrix12[5] = 1;
	matrix12[6] = 0;
	matrix12[7] = 0;
	matrix12[8] = 0;
	matrix12[9] = 0;
	matrix12[10] = 1;
	matrix12[11] = 0;

	egg my_egg(1, 1, 1);
	cube my_cube(matrix12);
	
	std::vector<implicit_function*> children;
	children.push_back(&my_egg);
	children.push_back(&my_cube);

	transformed_subtract obj(children, matrix12);

	find_min_z(obj, 5, 5);
	cout << "=================================" << endl;
}

void scone_test() 
{
	cout << "=================================" << endl;
	cout << "In scone" << endl;

	REAL matrix12[12];
	
	matrix12[0] = 1;
	matrix12[1] = 0;
	matrix12[2] = 0;
	matrix12[3] = 0;
	matrix12[4] = 0;
	matrix12[5] = 1;
	matrix12[6] = 0;
	matrix12[7] = 0;
	matrix12[8] = 0;
	matrix12[9] = 0;
	matrix12[10] = 1;
	matrix12[11] = 0;

	scone obj(matrix12);

	find_min_z(obj, 5, 5);
	cout << "=================================" << endl;
}

void scylinder_test() 
{
	cout << "=================================" << endl;
	cout << "In scylinder" << endl;

	REAL matrix12[12];
	
	matrix12[0] = 1;
	matrix12[1] = 0;
	matrix12[2] = 0;
	matrix12[3] = 0;
	matrix12[4] = 0;
	matrix12[5] = 1;
	matrix12[6] = 0;
	matrix12[7] = 0;
	matrix12[8] = 0;
	matrix12[9] = 0;
	matrix12[10] = 1;
	matrix12[11] = 0;

	scylinder obj(matrix12);

	find_min_z(obj, 5, 5);
	cout << "=================================" << endl;
}


void cube_test() 
{
	cout << "=================================" << endl;
	cout << "In cube" << endl;

	REAL matrix12[12];
	
	matrix12[0] = 1;
	matrix12[1] = 0;
	matrix12[2] = 0;
	matrix12[3] = 0;
	matrix12[4] = 0;
	matrix12[5] = 1;
	matrix12[6] = 0;
	matrix12[7] = 0;
	matrix12[8] = 0;
	matrix12[9] = 0;
	matrix12[10] = 1;
	matrix12[11] = 0;

	cube obj(matrix12);

	find_min_z(obj, 5, 5);
	cout << "=================================" << endl;
}



int main()
{
	egg_test();
	tetrahedron_test();
	heart_test();
	torus_test();
	subtraction_test();
	scone_test();//failed
	scylinder_test();//failed
	cube_test();

	return 0;
}
