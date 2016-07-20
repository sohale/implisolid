
/*
	File: perf_crisp_subtract.cpp
	Desc: Test the performance of CrispSubtraction
*/



#include "../primitives.cpp"
#include "../timer.hpp"
#include "../crisp_subtract.hpp"

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


/**
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
 *
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
												  // nr_chunks times

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
	std::cout << "Number of points: "  << nr_points << "\t" << "Number of chunks: \
	 " << nr_chunks << "\t" << "Execution duration: " << t1.stop() << std::endl;
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
	return 0;
}
