
/*
	File: perf_crisp_subtract.cpp
	Desc: Test the performance of CrispSubtraction
*/

#include "../primitives.cpp"
#include "../../js_iteration_1/timer.hpp"	// provides a timer class
#include "../crisp_subtract.hpp"	// provides crisp_subtract class
#include "../object_factory.hpp" 	// provides object_factory function
/**
 * Function: test_crisp_subtract_performance
 * ------------------------------------------
 * Usage: test_crisp_subtract_performance(nr_points);
 *
 * Implementation notes:
 *
 * The purpose of this function is to time the performance of the evaluation
 * of a CrispSubtract implicit function on a number of points(nr_points).
 *
 * Timings:
 *
 * #Points		| 	Time(ms)
 * ----------------------------
 * 1M			| 300
 * 1M 10 chunks | 137
 */

void
test_crisp_subtract_performance(const int nr_points)
{
	int nr_chunks=1;
	std::clog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::clog << "Number of points: "  << nr_points << "t" << "Number of chunks: " << ( nr_chunks ? nr_chunks: 1 ) <<std::endl;
	auto sf = make_shape_1d(nr_points); //  Create the shape of the output
	vectorized_scalar f = vectorized_scalar(sf);

	mp5_implicit::unit_sphere s1(2.0);
	mp5_implicit::unit_sphere s2(1.0);

	mp5_implicit::CrispSubtract crs = mp5_implicit::CrispSubtract(s1, s2);

	vectorized_vect x = make_empty_x(nr_points);

	//  Fill x with values
		for (int i = 0; i < nr_points; i++)
		{
	 		x[i][0]= 0.0;
			x[i][1]= 0.0;
			x[i][2]= 0.0;
		}

	timer t;
	crs.eval_implicit(x, &f);
	t.stop();
}

/*
*
 * Function: test_crisp_subtract_performance_in_chunks
 * -------------------------------------------------------------
 * Usage: test_crisp_subtract_performance_in_chunks(nr_points)
 *
 * Implementation notes:
 *
 * The purpose of this function is to measure the
 * performance of CrispSubtract(potentially any other CSG operation) by
 * evaluating the produced implicit function on a number of points(nr_points).
 * In contrast with the function test_crisp_subtract_performance, we allocate
 * memory in chunks(e.g when we initialize the input x)
 *
 * Times:
 *
 * #Points		|		Time(ms)
 * ------------------------------
 * #TODO: provide timings by running node perf_crisp_subtract.compiled.js
 * #TODO: figure out the exception cause of the json object passed as a string:
 * 		a) is it corrupted ?
 * 			- See an example usage
 * 			- Ask S
 * 		b) is there a function mismatch?
 *
 */


void
test_crisp_subtract_performance_in_chunks(const int nr_points, const int nr_chunks)
{
	std::clog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
 	int i = 0; // Loop counter
	int chunk_size = nr_points / nr_chunks;

	// Create the CrispSubtract object.

	mp5_implicit::unit_sphere s1(2.0);
	mp5_implicit::unit_sphere s2(1.0);
	mp5_implicit::CrispSubtract crs = mp5_implicit::CrispSubtract(s1, s2);

	// Allocate input and ouput vectors.

	auto sf = make_shape_1d(chunk_size);
	vectorized_scalar f = vectorized_scalar(sf);
	vectorized_vect x = make_empty_x(chunk_size); // this will be ovewritten
												  // nr_chunks times::::

	// Initialize x . This is done only once.
	for (i = 0; i < chunk_size; i++ ){
		x[i][0] = 0.0;
		x[i][1] = 0.0;
		x[i][2] = 0.0;
	};

	timer t1;
	for (int j=0; j < nr_chunks; j++){
		crs.eval_implicit(x, &f);
	}
	t1.stop();
	std::clog << "Number of points: "  << nr_points << "\t" << "Number of chunks: \
	 " << nr_chunks << "\t" << std::endl;
}

/**
 * Function: test_mp5_object
 * Description:
 *
 * This functions measures the performance of the calculation of nr_points points
 * on an implicit function. To retrieve the implicit function we use the
 * object_factory(str, bool, bool). The first argument of object_factory, is of type string, and it represents
 * a json object that specifies the type of the implicit function and properties.
 * The second argument is of no use to this function so can be set to false.
 * Third argument is set to default = true;
 */

/**
 * Performance:
 *
 * #Points		|		Time(ms)
 * 1M 			|			159
 * 100K			|		 	 24
 * 10K			|	      	  7
 */
void
test_mp5_object(int nr_points) {

	std::clog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::string str_json = "{\"printerSettings\":{\"name\":\"test\",\"layerThickness\":0.2,\"emptyLayer\":0,\"infillSpace\":4,\"topThickness\":0.6,\"paramSamples\":75,\"speedRate\":1000,\"circleSpeedRate\":1000,\"temperature\":220,\"inAirSpeed\":7000,\"flowRate\":0.035,\"criticalLength\":35,\"retractionSpeed\":2400,\"retractionLength\":5,\"shellNumber\":3,\"material\":\"PLA 2.85mm\",\"autoZScar\":true,\"zScarGap\":0.5,\"critLayerTime\":6,\"filamentDiameter\":2.85},\"mp5-version\":\"0.3\",\"root\":{\"type\":\"root\",\"children\":[{\"type\":\"Difference\",\"protected\":false,\"children\":[{\"type\":\"cylinder\",\"displayColor\":{\"x\":0.7675997200783986,\"y\":0.03892568708507049,\"z\":0.1754374135888661},\"matrix\":[42.583184547062736,0,0,0,0,49.55270399250958,0,0,0,0,10,0,0,0,0,1],\"index\":652818},{\"type\":\"Difference\",\"protected\":false,\"children\":[{\"type\":\"cylinder\",\"displayColor\":{\"x\":0.8122645344236872,\"y\":0.657334404743416,\"z\":0.7357336310755096},\"matrix\":[10,0,0,0,0,10,0,0,0,0,10,0,0,0,0,1],\"index\":1272174},{\"type\":\"cylinder\",\"displayColor\":{\"x\":0.11421729990684737,\"y\":0.07562705374348999,\"z\":0.6324600862122098},\"matrix\":[10,0,0,0.658889604636343,0,10,0,6.215549332615993,0,0,10,1.3327027659215673e-7,0,0,0,1],\"index\":2463576}],\"initialSize\":{\"x\":1,\"y\":1,\"z\":1},\"displayColor\":{\"x\":0.6627450980392157,\"y\":0.4549019607843137,\"z\":0.7215686274509804},\"matrix\":[2.381193509886417,0,0,0.3600215429489424,0,2.381193509886417,0,0.5604901669421452,0,0,2.381193509886417,6.9059681360437395,0,0,0,1],\"index\":413872}],\"initialSize\":{\"x\":1,\"y\":1,\"z\":1},\"displayColor\":{\"x\":0.5529411764705883,\"y\":0.06666666666666667,\"z\":0.11764705882352941},\"matrix\":[1,0,0,0.32938436512727,0,1,0,0.15604124684634,0,0,1,0.000000000000014,0,0,0,1],\"index\":6565922}]}}";
	bool use_metaball = false, ignore_root_matrix = true;

	std::stringstream shape_json_stream;
    shape_json_stream << str_json ;

	pt::ptree shapeparams_dict;
	pt::ptree root_params;

	pt::read_json(shape_json_stream, shapeparams_dict);
	root_params = shapeparams_dict.get_child("root");


	implicit_function *object; 	//  = object_factory(root_params.get_child("children"), use_metaball, ignore_root_matrix);
	implicit_function * a = NULL;
	implicit_function * b = NULL;

	int count = 0;
	for (pt::ptree::value_type &element : root_params.get_child("children")){
		a = object_factory(element.second, use_metaball, ignore_root_matrix);
		count++;
	}

	auto sf = make_shape_1d(nr_points);
	vectorized_vect x = make_empty_x(nr_points); // This is the input; must be an orthogonal 3D grid, like numpy's mgrid
	vectorized_scalar f = vectorized_scalar(sf); // This allocates memory for the output data structure,
												 // which holds the values of the implicit_function corresponding to x
	// Initialize x
	for (int i = 0; i < nr_points; i++ ){
		x[i][0] = 0.0;
		x[i][1] = 0.0;
		x[i][2] = 0.0;
	};

	timer t1;
			a->eval_implicit(x, &f);
		t1.stop();
}

int main()
{
	test_crisp_subtract_performance(1000000);
	test_crisp_subtract_performance_in_chunks(1000000, 1);
	test_crisp_subtract_performance_in_chunks(1000000, 10);
	test_crisp_subtract_performance_in_chunks(1000000, 20);
	test_crisp_subtract_performance_in_chunks(1000000, 50);
	test_crisp_subtract_performance_in_chunks(1000000, 100);
	test_crisp_subtract_performance_in_chunks(1000000, 200);
	test_crisp_subtract_performance_in_chunks(1000000, 500);
	test_crisp_subtract_performance_in_chunks(1000000, 1000);
	test_mp5_object(10000);
	return 0;
}
