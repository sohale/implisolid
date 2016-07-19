#include "../primitives.cpp"
#include "../timer.hpp"
#include "../crisp_subtract.hpp"

/*
	File: perf_crisp_subtract.cpp 
	Desc: Test the performance of CrispSubtraction 
*/

/*
	Function: test_crisp_subtract_performance 
	Args: int nr_points
*/
void 
test_crisp_subtract_performance(const int nr_points)
{
	auto sf = make_shape_1d(nr_points); //  Create the shape of the output
	vectorized_scalar f = vectorized_scalar(sf); 

	vectorized_vect x = make_empty_x(nr_points); 

	mp5_implicit::unit_sphere s1(2.0); 	
	mp5_implicit::unit_sphere s2(1.0);

	mp5_implicit::CrispSubtract crs = mp5_implicit::CrispSubtract(s1, s2);  

	//  Fill x with values 
	int i = 0;
	for (int j = 1; j < 50; j++){
		for (;i < nr_points / 50; i++)
		{
	 		x[i][0]= rand();
			x[i][1]= rand();
			x[i][2]= rand();
		}
	}
	timer t; 

	crs.eval_implicit(x, &f); 

	t.stop();
}

int main()
{
	test_crisp_subtract_performance(1000000); 
	return 0; 
}