/**
 * Purpose: Perform timings to compare Boost, Armadillo and native arrays.
 *
 * Author: @michael-karotsieris
 *
 * Notes: It is important to compile using compiler optimization flags
 * 	      in order to get major speed boosts.
 *
 * 	-O2, -O3: level 2 and 3 optimization --> increases execution speed
 *
 * 	-DBOOST_DISABLE_ASSERTS: disables internal asserts for the boost library
 *  -
 */


#include <iostream>
#include <boost/multi_array.hpp>
#include <armadillo>


using namespace std;
int main(int argc, char* argv[])
{
    // declarations
    const int X_SIZE = 1000;
    const int Y_SIZE = 1000;

    const int ITERATIONS = 100;

    clock_t start;
    double duration;


    // create boost array
    {
        typedef boost::multi_array<double, 2> ImageArrayType;
        ImageArrayType boostMatrix(boost::extents[X_SIZE][Y_SIZE]);

        // init starter
        start = clock();

        // computation loop
        for (int i=0; i < ITERATIONS; ++i)
        {
            for (int x=0; x < X_SIZE; ++x)
            {
                for (int y = 0; y < Y_SIZE; ++y)
                {
                    boostMatrix[x][y] = 1.001;
                }
            }
        }
        duration = ( clock() - start) / ((double) CLOCKS_PER_SEC);
        // cout << ("[Boost] Elapsed Time: " << duration <<  "seconds\n";
        cout << "[Boost]: \t" << duration << endl;
    }
    // measure armadillo loop
    {
        arma::mat matArray( X_SIZE, Y_SIZE);

        start = clock();

        for (int i=0; i < ITERATIONS; ++i)
        {
            for (int x=0; x < X_SIZE; ++x)
            {
                for (int y = 0; y < Y_SIZE; ++y)
                {
                    matArray(x,y) = 1.001;
                }
            }
        }

        duration  = (clock() - start) / ((double) CLOCKS_PER_SEC);
        cout << "[Armadillo]: \t" << duration << endl;
    }

    // measure native
    {
        double *nativeMatrix = new double [X_SIZE * Y_SIZE];
        start = clock();

        for (int i=0; i < ITERATIONS; ++i)
        {
            for (int y = 0; y < X_SIZE*Y_SIZE; ++y)
                {
                    nativeMatrix[y] = 1.001;
                }
        }

        duration  = (clock() - start) / ((double) CLOCKS_PER_SEC);
        cout << "[Native]: \t" << duration << endl;
    }
    return 0;
}
