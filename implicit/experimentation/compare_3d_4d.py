import sys

import sys
import math


import numpy as np

@profile
def test1():
    sys.path.append("..")
    import profile_support
    from vtk_mc import vtk_mc as v1

    from basic_types import check_vector4_vectorized, normalize_vector4_vectorized

    from example_objects import make_example_vectorized as e1
    iobj = e1("cube_with_cylinders")  #
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2/4.)

    from stl_tests import make_mc_values_grid as m1
    gridvals = m1(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = v1(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()
    print sys.path
#    del stl_tests


@profile
def test2():
    sys.path.insert(0, '../clean_code')
    import profile_support
    from vtk_mc import vtk_mc as v2
    from basic_functions import check_vector3_vectorized, normalize_vector3_vectorized

    from example_objects import make_example_vectorized as e2

    iobj = e2("cube_with_cylinders")  #

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2/4.)

    from stl_tests import make_mc_values_grid as m2
    gridvals = m2(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE, old=False)
    verts, facets = v2(gridvals, (RANGE_MIN, RANGE_MAX, STEPSIZE))
    print("MC calculated.");sys.stdout.flush()
    del sys.path[0]
    print sys.path

#test1()
test2()
