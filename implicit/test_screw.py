import numpy as np
VERBOSE = False
from basic_types import check_vector4_vectorized
FAST_VERSION = True
import math
# from ipdb import set_trace
from utils import optimised_used
from utils import flush


from example_objects import *

def main():
    import numpy as np

    h = 0.00001;

    xa = np.array([[12.0,23.0,2.0],
                   [0,0,-0.4],
                   [-1,1,1],
                   [0, 0, -0.1],
                   [1331,1221,1222]
                  ])

    # xa_plus_h = np.array([[12.0,23.0 + h,2.0],
    #                [0,0 + h,-0.4],
    #                [-1,1 + h,1],
    #                [0, 0 + h, -0.1]])


    screw = make_example_vectorized('screw2');

    print(screw.implicitFunction(xa))
    # print(screw.implicitFunction(xa_plus_h))

    # print((screw.implicitFunction(xa_plus_h) - screw.implicitFunction(xa))/h)


    # print(screw.implicitGradient(xa))

if __name__ == '__main__':
    main()