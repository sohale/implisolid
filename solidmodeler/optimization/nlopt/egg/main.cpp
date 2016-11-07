#include <math.h>
#include <nlopt.h>
#include <stdio.h>
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>

inline double norm_squared(double x, double y, double z) {
    return x*x + y*y + z*z;
}

inline double squared(double x) {
    return x*x;
}

typedef struct {
    double a, b, c, center_x, center_y, center_z;

    void print() {
    	std::cout << a << " "<< b << " "<<c << " "<< center_x << " "<< center_y <<  " "<< center_z << std::endl;
    }
} my_egg_data;

double myfunc(unsigned n, const double *m, double *grad, void *my_func_data)
{
    if (grad) {
        grad[0] = 0;
        grad[1] = 0;
        grad[2] = 1;
    }
         
    return m[2];       
}

 double myconstraint_eq(unsigned n, const double *m, double *grad, void *data)
{
    my_egg_data *d = (my_egg_data *) data;

    double x = m[0]; 
    double y = m[1];
    double z = m[2];

	//std::cout << m[0] <<  " " << m[1]  <<  " " << m[2] << std::endl;
    
    if (grad) {
    	double a2 = squared(d->a);
        double b2 = squared(d->b);
        double c2 = squared(d->c);

        grad[0] = -2. * (x - d->center_x)/a2;
        grad[1] = -2. * (y - d->center_y)/b2;
        grad[2] = -2. * (z - d->center_z)/c2;

        //std::cout << grad[0] <<  " " << grad[1]  <<  " " << grad[2] << std::endl;
    }

    double egg_i = 1 - norm_squared((x - d->center_x)/d->a, (y - d->center_y)/d->b, (z - d->center_z)/d->c);
           
    //std::cout << "constr eq: " << egg_i << std::endl;
	return egg_i;
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

	double lb[3] = { -10, -10, -10 }; /* lower bounds */
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LD_SLSQP, 3); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc, NULL);

	my_egg_data data = {1, 2, 3, 0, 0, 0};

	nlopt_add_equality_constraint(opt, myconstraint_eq, &data, 1e-4);

	nlopt_set_xtol_rel(opt, 1e-4);
	nlopt_set_maxeval(opt, 100);
	nlopt_set_xtol_abs1(opt, 1e-4);

	double m[3] = { 1, -1, 0 };  /* some initial guess */
	double minf; /* the minimum objective value, upon return */
	
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
		    //printf("found minimum at f(%g,%g,%g) = %0.10g\n", m[0], m[1], m[2], minf);
		}
	}

	std::cout << "All violation " << counter << std::endl;
	
	nlopt_destroy(opt);

	return 0;
}
