

###################################################################
# This file does a simple profiling on a numpy function so we can #
# compare the running time with a corresponding C++ function.     #
###################################################################

# Notes
# 1. Run  with python -m kernprof -v -l
# 2. If it complains about not recognizing profile , try running with
#   /usr/bin/python -m kernprof -v -l (uses your default python interpreter)
# or install line_profiler and kernprof again 
import numpy as np

@profile
def function_circle(dim=10000):

    # create random nX3 vector
    x = np.random.random((dim, 3))

    # create nX1 vector to store the result
    res = np.zeros((dim, 1))

    # compute result
    for i in range(x.shape[0]):
        res[i] = 2 - x[i, 0]**2 - x[i, 1]**2 - x[i, 2] ** 2

    return res

if __name__ == "__main__":
    function_circle()
