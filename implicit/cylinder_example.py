
from example_objects import cyl1, cyl2, cage_rods, cyl3, cyl4

# import vectorized


import time

class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start



import sys;


def cyl_test_example1():

    print("started")
    sys.stdout.flush()


    # c = cyl1()
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 1.92 * 0.2 * 10/ 2.0)

    #c, (RANGE_MIN, RANGE_MAX, STEPSIZE) = cyl2()
    ##(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32, +32, 1.92 / 4.0*1.5/ 1.5)
    # #takes 15 sec!

    c, (RANGE_MIN, RANGE_MAX, STEPSIZE) = cyl3()
    ##(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32, +32, 1.92 / 4.0)   #15 sec!  2.5 millions voxels
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32/2, +32/2, 1.92 / 4.0)  # 2.5 sec!

    # spiral cage
    #c, (RANGE_MIN, RANGE_MAX, STEPSIZE) = cyl4()
    #(RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32/2, +32/2, 1.92 / 4.0)

    c, dims = cyl4()
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = dims

    print("made the function")
    sys.stdout.flush()

    from stl_tests import make_mc_mesh_scikit
    with Timer() as t:
        # @profile
        verts, faces = make_mc_mesh_scikit(c, RANGE_MIN, RANGE_MAX, STEPSIZE)
    print("made MC : done within ", t.interval)  # takes 15 sec!!
    sys.stdout.flush()

    ACTUALLY_SAVE = False
    #from stl_tests import optimise_mesh
    #verts = optimise_mesh(verts, faces, c)
    #Does not work yet
    from stl_tests import m2stl_mesh
    m = m2stl_mesh(verts, faces)
    if ACTUALLY_SAVE:
        m.save('stl/cage.stl')  # wow

    from mesh_utils import mesh_invariant
    mesh_invariant(faces)

    from ohtake_surface_projection import display_simple_using_mayavi_

    print("mayavi start")
    sys.stdout.flush()

    display_simple_using_mayavi_([(verts, faces), ], [], minmax=(RANGE_MIN, RANGE_MAX))

    print("mayavi done")
    sys.stdout.flush()


if __name__ == '__main__':
    cyl_test_example1()
