import sys
sys.path.append("..")
import matplotlib.pyplot as plt

import numpy as np
global STEPSIZE
from example_objects import make_example_vectorized
ifunc = make_example_vectorized("cube_example")
(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2*1.5/1.5  *2. /2.)


global x

def test1():
    f = ifunc.implicitFunction(x)

def test2():
    for i in range(x.shape[0]):
        f = ifunc.implicitFunction(x[i, np.newaxis, :])

def experiment1():

    #asymp = (10.**4, 2.37+03/10.**4 / 100)
    #print asymp
    #asymp_ratio1 = (2.37e+03) /(10. * 4)/100000000.  # 0.592 microseconds
    asymp_ratio1 = 0.59 / (10. ** 6)  # 0.592 microseconds
    #asymp_ratio2 = 0.50 / (10. ** 6)
    #asymp_ratio2 = (7.33606168e+01)/100000.
    asymp_ratio2 = (0.688)/1000./1000.   # 0.68 microseconds + 350 microseconds
    asymp_bias = 0.000350
    print asymp_ratio1, asymp_ratio2, asymp_ratio2/asymp_ratio1, "bias=", asymp_bias



    global x
    na = [1,2,5,10, 100,200,400,600,] #800,  1000, 2000, 10000, 100000]
    tl = []
    for ei in range(len(na)):
        n = na[ei]
        x = (np.random.rand(n, 4)*2-1)*(RANGE_MAX-RANGE_MIN)+RANGE_MIN
        x[:, 3] = 1
        print "init:", ei,

        repeats = max(int(math.ceil(1000/n)), 2)
        import timeit
        t1 = timeit.timeit('test1()', "from __main__ import test1", number=repeats)
        t2 = timeit.timeit('test2()', "from __main__ import test2", number=repeats)
        tl.append((t1/repeats, t2/repeats))
        print "."
    ta = np.array(tl)
    print ta.shape
    print ta * 1000  # in msec
    print "ratio =", ta[:, 1]/ta[:, 0]

    #import matplotlib.pyplot as plt
    #ratio = ta[:, 1]/ta[:, 0]
    #plt.plot(na, ratio, "r*-")
    #ratio = /ta[:, 0]
    naa = np.array(na)

    #g0 = plt.loglog(naa, naa/asymp[0]*asymp[1]/100, "b:", label='asympt')
    g1 = plt.loglog(naa, ta[:, 1], "r*-", label='point-wise')
    g2 = plt.loglog(naa, ta[:, 0], "bs-", label='numpy vectorised')
    g0 = plt.loglog(naa, naa * asymp_ratio1, "b:", label='asympt')
    g01 = plt.loglog(naa, naa * asymp_ratio2+ asymp_bias, "b:", label='asympt')
    plt.legend(loc='upper left', numpoints = 1)
    plt.show()

    #g1 = plt.plot(naa, ta[:, 1], "r*-", label='point-wise')
    g2 = plt.plot(naa, ta[:, 0], "bs-", label='numpy vectorised')
    plt.plot(naa, naa * asymp_ratio1, "b:", label='asympt')
    g0 = plt.plot(naa, naa * asymp_ratio2 + asymp_bias, "b:", label='asympt')
    plt.legend(loc='upper left', numpoints = 1)
    plt.show()

if __name__ == "__main__":

    import math
    experiment1()
