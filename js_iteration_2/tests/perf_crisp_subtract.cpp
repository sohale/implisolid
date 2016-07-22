
/*
	File: perf_crisp_subtract.cpp
	Desc: Test the performance of CrispSubtraction
*/

#include "../primitives.cpp"
#include "../timer.hpp"				// provides a timer class
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
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "Number of points: "  << nr_points << "\t" << "Number of chunks: " << ( nr_chunks ? nr_chunks: 1 ) <<std::endl;
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
 */


void
test_crisp_subtract_performance_in_chunks(const int nr_points, const int nr_chunks)
{
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
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
	std::cout << "Number of points: "  << nr_points << "\t" << "Number of chunks: \
	 " << nr_chunks << "\t" << std::endl;
}

/**
 * Function: test_mp5_object
 * Description:
 * This functions measures the performance of the calculation of nr_points points
 * on an implicit function. To retrieve the implicit function we use the
 * object_factory(str, bool). The first argument of object_factory, is of type string, and it represents
 * a json object that specifies the type of the implicit function and various properties.
 * The second argument is of no use to this function.
 *
 */

#include "json.hpp"
using json = nlohmann::json;

void
test_mp5_object(int nr_points) {

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::string str_json = "{\"printerSettings\":{\"name\":\"test\",\"layerThickness\":0.2,\"emptyLayer\":0,\"infillSpace\":4,\"topThickness\":0.6,\"paramSamples\":75,\"speedRate\":1000,\"circleSpeedRate\":1000,\"temperature\":220,\"inAirSpeed\":7000,\"flowRate\":0.035,\"criticalLength\":35,\"retractionSpeed\":2400,\"retractionLength\":5,\"shellNumber\":3,\"material\":\"PLA 2.85mm\",\"autoZScar\":true,\"zScarGap\":0.5,\"critLayerTime\":6,\"filamentDiameter\":2.85},\"mp5-version\":\"0.3\",\"root\":{\"type\":\"root\",\"children\":[{\"type\":\"Union\",\"protected\":false,\"children\":[{\"type\":\"cylinder\",\"displayColor\":{\"x\":0.9661750055896976,\"y\":0.8085857086202395,\"z\":0.41578037212168595},\"matrix\":[10,0,0,-3.7368777450135866,0,10,0,-1.9559832356144682,0,0,10,1.7323194345664206e-7,0,0,0,1],\"index\":7575510},{\"type\":\"cube\",\"displayColor\":{\"x\":0.23399071378141634,\"y\":0.31584816496653323,\"z\":0.35457351563365425},\"matrix\":[10,0,0,1.867397834493545,0,10,0,-1.7325527119763677,0,0,10,-9.734106853898084e-10,0,0,0,1],\"index\":5587759},{\"type\":\"cylinder\",\"displayColor\":{\"x\":0.43814645627496795,\"y\":0.39556472441055845,\"z\":0.3415798414286939},\"matrix\":[10,0,0,1.8694799105200275,0,10,0,3.688535947590836,0,0,10,-1.7225853365943067e-7,0,0,0,1],\"index\":6657333}],\"initialSize\":{\"x\":1,\"y\":1,\"z\":1},\"displayColor\":{\"x\":0.6470588235294118,\"y\":0.2784313725490196,\"z\":0.5882352941176471},\"matrix\":[1,0,0,82.63768850593796,0,1,0,126.37324151118989,0,0,1,5.000000079354265,0,0,0,1],\"index\":4872526}]}}";
	bool use_metaball = false;
	
	/*****************************************************
	* There is an exception happening in the next line  *
	******************************************************/
	
	implicit_function* object = object_factory(str_json,use_metaball,true);	
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

	object->eval_implicit(x, &f);

	t1.stop();

	evaluate object-->function
	put timer around it
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
	test_mp5_object(1000);
	return 0;
}
