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

#include "shape_optimization_utils.hpp"

using namespace std;
using namespace mp5_implicit;

const int OPTIMIZATION_ITERATIONS = 5;
const double OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION = 0.5;
const double DELTA = 0.01;


void egg_test()
{
	cout << "=================================" << endl;
	cout << "In egg" << endl;
	
	egg obj(1, 2, 4);
	double correct_min_z = -4;

	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

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

	double correct_min_z = 1;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}
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

	double correct_min_z = -1.02;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

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

	double correct_min_z = -0.2;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

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
	
	double correct_min_z = -1;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

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

	double correct_min_z = -0.5;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

	cout << "=================================" << endl;
}

void scylinder_test() 
{
	cout << "=================================" << endl;
	cout << "In scylinder" << endl;

	REAL matrix12[12];
	
	matrix12[0] = 10;
	matrix12[1] = 0;
	matrix12[2] = 0;
	matrix12[3] = 0;
	matrix12[4] = 0;
	matrix12[5] = 10;
	matrix12[6] = 0;
	matrix12[7] = 0;
	matrix12[8] = 0;
	matrix12[9] = 0;
	matrix12[10] = 10;
	matrix12[11] = 0;

	scylinder obj(matrix12);

	double correct_min_z = -5;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

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

	double correct_min_z = -0.5;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

	cout << "=================================" << endl;
}

void union_test() 
{
	cout << "=================================" << endl;
	cout << "In union test, egg + cube " << endl;

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
	
	//move egg
	vectorized_vect direction(boost::extents[1][3]);;
	direction[0][0] = 0;
	direction[0][1] = 0;
	direction[0][2] = -0.5;
	my_egg.move(direction);

	cube my_cube(matrix12);
	
	std::vector<implicit_function*> children;
	children.push_back(&my_egg);
	children.push_back(&my_cube);

	transformed_union obj(children, matrix12);

	double correct_min_z = -1.5;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

	cout << "=================================" << endl;
}

void intersection_test() 
{
	cout << "=================================" << endl;
	cout << "In intersection test, egg * cube " << endl;

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
	
	//move egg down -1
	vectorized_vect direction(boost::extents[1][3]);;
	direction[0][0] = 0;
	direction[0][1] = 0;
	direction[0][2] = -1;
	my_egg.move(direction);

	//move cube up 0.1
	cube my_cube(matrix12);
	direction[0][0] = 0;
	direction[0][1] = 0;
	direction[0][2] = 0.1;
	my_cube.move(direction);

	
	std::vector<implicit_function*> children;
	children.push_back(&my_egg);
	children.push_back(&my_cube);

	transformed_intersection obj(children, matrix12);

	double correct_min_z = -0.4;
	double our_min_z = find_min_z(
		obj, 
		OPTIMIZATION_ITERATIONS, 
		OPTIMIZATION_STARTING_POINT_STANDART_DEVIATION
	);
	double diff = my_abs(correct_min_z - our_min_z);

	if (diff < DELTA) {
		cout << "Test passed!" << endl;
	} else {
		cout << "Test failed!" << endl;
	}

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
	union_test();
	intersection_test();

	return 0;
}