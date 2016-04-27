import numpy as np

from basic_functions import make_vector4, normalize_vector

import simple_blend

import vectorized

#definition of the vectorized objects

def dice(dice_scale):
    dice_size = 2.0*dice_scale
    c = vectorized.UnitCube1(size=dice_size)

    def hole_crisp(c, i, j, k):
        m = np.eye(4)
        dot_size = 0.25 * dice_scale
        m[0, 0] = dot_size
        m[1, 1] = dot_size
        m[2, 2] = dot_size
        distance = (0.5-0.05)*dice_size
        m[0:3, 3] = np.array([distance*i, distance*j, distance*k])
        s1 = vectorized.Ellipsoid(m)
        c = vectorized.CrispSubtract(c, s1)
        return c

    def hole_r(c, i, j, k):
        m = np.eye(4)
        dot_size = 0.25 * dice_size
        m[0, 0] = dot_size
        m[1, 1] = dot_size
        m[2, 2] = dot_size
        distance = (0.5-0.05)*dice_size
        m[0:3, 3] = np.array([distance*i, distance*j, distance*k])
        s1 = vectorized.Ellipsoid(m)
        c = vectorized.RSubtract(c, s1, 1.0)
        return c

    hole = hole_crisp
    #hole = hole_r

    """ 1 """
    c = hole(c, 1, 0, 0)  # 1
    #c = hole(c, 0, 1, 0)  # 2
    c = hole(c, 0, 0, 1)  # 3
    #c = hole(c, -1, 0, 0)  # 4
    c = hole(c, 0, -1, 0)  # 5
    #c = hole(c, 0, 0, -1)  # 6

    """ 2 """
    d = 0.3
    c = hole(c, -d, 1, -d)
    c = hole(c, +d, 1, +d)

    """ 3 """
    d = 0.6
    c = hole(c, -d, -d, 1)
    c = hole(c, +d, +d, 1)

    """ 4 """
    d = 0.4
    c = hole(c, -1, -d, +d)
    c = hole(c, -1, -d, -d)
    c = hole(c, -1, +d, -d)
    c = hole(c, -1, +d, +d)

    """ 5 """
    d = 0.6
    c = hole(c, -d, -1, +d)
    c = hole(c, -d, -1, -d)
    c = hole(c, +d, -1, -d)
    c = hole(c, +d, -1, +d)

    """ 6 """
    d = 0.7
    for i in range(3):
        c = hole(c, d*(i-1), +d, -1)
        c = hole(c, d*(i-1), -d, -1)

    return c



def rdice_vec(scale, rotated=True):

    d = dice(scale)
    return d
    iobj = vectorized.Transformed(d)
    iobj  \
        .move(-0.2 * scale, -0.2 * scale, 0) \
        .resize(0.9)
    iobj.rotate(20, along=make_vector4(1, 1, 1), units="deg")
    return iobj


def rcube_vec(scale, rotated=True):
    d = vectorized.UnitCube1(size=2.0 * scale)

    iobj = vectorized.Transformed(d)
    iobj  \
        .move(-0.2 * scale, -0.2 * scale, 0) \
        .resize(0.9)
    if rotated:
        iobj.rotate(10 * 2, along=make_vector4(1, 1, 1), units="deg")
    return iobj

def udice_vec(scale=1.):
   """ Un-rotated dice """
   return rdice_vec(scale, rotated=False)

def blend_example2(scale=1.):
    # not tested for scale != 1.
    m1 = np.eye(4) * 1.3 * scale
    m1[1, 1] = 0.4 * scale
    m1[1, 2] = 0.4 * scale
    m1[0:3, 3] = [0, 1 * scale, 0]
    m1[3, 3] = 1

    m2 = np.eye(4) * 2
    m2[0:3, 3] = [2 * scale, 0, 0]
    m2[2, 2] = 0.4 * scale
    m2[3, 3] = 1

    a, b = vectorized.Ellipsoid(m1), vectorized.Ellipsoid(m2)
    iobj = vectorized.SimpleBlend(a, b)

    return iobj


def ell_example1(scale):
    # scale not tested
    m1 = np.eye(4) * 1.3 * scale
    m1[1, 1] = 0.4 * scale
    m1[1, 2] = 0.4 * scale
    m1[0:3, 3] = [0, 1 * scale, 0]
    m1[3, 3] = 1

    iobj = vectorized.Ellipsoid(m1)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)

    return iobj


def sphere_example(scale=1.):
    iobj = vectorized.UnitSphere()

    return iobj


def cube_example(scale=1.):
    iobj = vectorized.UnitCube1()

    iobj = vectorized.Transformed(iobj) \
        .move(-0.1 * scale, -0.1 * scale, -0.1 * scale) .resize(3 * scale) \
        .rotate(-20, along=make_vector4(1, 1, 1), units="deg") .move(0.2 * scale, 0, 0)

    return iobj


def first_csg(scale):
    # scale not tested
    m1 = np.eye(4) * 2 * scale
    m1[1, 1] = 0.4 * scale
    m1[0:3, 3] = [0, 0, 0]
    m1[3, 3] = 1

    m2 = np.eye(4) * 1. * scale
    rcenter = np.array([0.5, 0, 0]) * scale
    m2[0:3, 3] = rcenter[0:3]
    m2[3, 3] = 1

    m3 = np.eye(4) * 1.2 * scale
    m3[0:3, 3] = [1.5 * scale, 0, 0]
    m3[3, 3] = 1

    iobj = vectorized.CrispSubtract(vectorized.CrispUnion(vectorized.Ellipsoid(m1), vectorized.Ellipsoid(m2)), vectorized.Ellipsoid(m3))
    return iobj

def bowl_15_holes(scale_ignored):

    big_radius = 3
    big_radius2 = 2.7  # +0.3-0.001
    small_radius = 0.6

    """ Make a sphere """
    m_big = np.eye(4) * big_radius
    m_big[0:3, 3] = [0, 0, 0]
    m_big[3, 3] = 1
    iobj = vectorized.Ellipsoid(m_big)

    """ Cut a sphere out pf it """
    m_big2 = np.eye(4) * big_radius2
    m_big2[0:3, 3] = [0, 0, 0]
    m_big2[3, 3] = 1

    iobj = vectorized.CrispSubtract(iobj, vectorized.Ellipsoid(m_big2))

    """ Cut the top part"""
    m_big3 = np.eye(4) * 10
    m_big3[0:3, 3] = [0, 0, 10.2]
    m_big3[3, 3] = 1
    iobj = vectorized.CrispSubtract(iobj, vectorized.Ellipsoid(m_big3))

    """ Cut 15 small spheres from it """
    for i in range(0, 15):
        m_small = np.eye(4) * small_radius

        if False:
            unsat = True
            while unsat:
                c = make_random_vector(big_radius, 1)[0:3]
                unsat = c[2] > 0
        # th0 = (i/15.0)*(3.1415926536*2)
        # th = 5*th0
        z0 = (i / 15.0) * big_radius
        # z0 = np.sqrt(1-(float(i)/15.0)**2) * big_radius
        # print(z0)
        # th = i * np.pi*2 * 5/ 4.45
        # th =  (np.pi*2) * float(i) * 4.0 / 5.0
        # th =  (np.pi*2) * float(i) * (1.0/5.0 * 1.0/2.0 * 1.0/5.0) * 5
        NN = float(8)
        # th =  (np.pi*2) * float(i) * (1.0/5.0 + 1.0/5.0 * 1.0/2.0 * 1.0/5.0)
        th = (np.pi * 2) * float(i) * (1.0 / NN + 1.0 / NN * 1.0 / NN / 2.0)
        c = make_vector4(np.cos(th) * big_radius, np.sin(th) * big_radius, -z0)
        c = normalize_vector(c) * big_radius
        c[3] = 1

        # print( np.sqrt(np.dot(c,c)) )
        m_small[0:3, 3] = c[0:3]
        m_small[3, 3] = 1
        small_obj = vectorized.Ellipsoid(m_small)
        iobj = vectorized.CrispSubtract(iobj, small_obj)
        # iobj = CrispUnion( iobj, small_obj )

    xa = vectorized.repeat_vect4(10, make_vector4(1, 1, 1))
    ga = iobj.implicitGradient(xa)
    va = iobj.implicitFunction(xa)
    # print(va)
    # print(ga)
    return iobj


def blend_example1(scale_ignored):
    m1 = np.eye(4) * 1
    # m1[1,1] = 0.4
    m1[0:3, 3] = [0, 1, 0]
    m1[3, 3] = 1

    m2 = np.eye(4) * 2
    m2[0:3, 3] = [2.5, 0, 0]
    m2[3, 3] = 1

    iobj = vectorized.SimpleBlend(vectorized.Ellipsoid(m1), vectorized.Ellipsoid(m2))

    return iobj


def blend_example2_discs(scale):
    m1 = np.eye(4) * 1.3 * scale
    m1[1, 1] = 0.4 * scale
    m1[1, 2] = 0.4 * scale
    m1[0:3, 3] = [0, 1 * scale, 0]
    m1[3, 3] = 1

    m2 = np.eye(4) * 2 * scale
    m2[0:3, 3] = [2 * scale, 0, 0]
    m2[2, 2] = 0.4 * scale
    m2[3, 3] = 1

    iobj = vectorized.SimpleBlend(vectorized.Ellipsoid(m1), vectorized.Ellipsoid(m2))
    assert vectorized.SimpleBlend == simple_blend.SimpleBlend

    return iobj

def french_fries_vectorized(scale):
    def rod():
        c = vectorized.UnitCube1()
        m2 = np.eye(4) * scale
        m2[0, 0] = 0.1 * scale
        m2[1, 1] = 0.1 * scale
        m2[0:3, 3] = [+0.1 / 2.0 * 2 * scale, +0.1 / 2.0 * 2 * scale, + 1.0 / 2.0 * 2 * scale]
        m2[3, 3] = 1
        iobj = vectorized.Transformed(c, m2)
        # .move(-0.1/2.0, -0.1/2.0, -1.0/2.0)
        iobj  \
            .move(-0.2 * scale, -0.2 * scale, 0) \
            .resize(2) \
            .rotate(10, along=make_vector4(1, 1, 0), units="deg") \
            .move(0.5 * scale, 0, 0)
        return iobj

    u = None
    for i in range(18):
        c = rod().rotate(-30 * i, along=make_vector4(0, 0, 1), units="deg")
        if u is None:
            u = c
        else:
            u = vectorized.CrispUnion(u, c)
    return u


def rods(scale):
    #the scale needs to be equal to 2 in order to have a nice object
    def rod():
        c = vectorized.UnitCube1()
        m2 = np.eye(4)
        m2[0, 0] = 0.1 * scale
        m2[1, 1] = 0.1 * scale
        m2[2, 2] = scale
        iobj = vectorized.Transformed(c, m2)
        iobj  \
            .move(-0, -0.3*scale, -1.0*scale) \
            .resize(2) \
            .rotate(40, along=make_vector4(1, 1, 0), units="deg") \
            .move(0.5*scale, 0, 0)
        return iobj

    u = None
    for i in range(17):
        c = rod().rotate(-30 * i, along=make_vector4(0, 0, 1), units="deg")
        if u is None:
            u = c
        else:
            u = vectorized.CrispUnion(u, c)
    return u


import screw

def screw1(scale_ignored):
    a = np.array([0, 0, 0])
    w = np.array([0, 0, 1])
    u = np.array([1, 0, 0])
    r0 = 0.5
    delta = 0.2
    lambda_ = 0.4 * 2
    phi0 = 0
    screw_len = 2
    return screw.Screw(a, w, u, screw_len, r0, delta, lambda_, phi0)


def screw2(SCALE=1.):
    # not tested for SCALE != 1.
    a = np.array([0, 0, 0])     # only works with [0,0,0]
    w = np.array([0, 0, -1])
    u = np.array([1, 0, 0])

    r0 = 0.5 * SCALE
    delta = 0.2 * SCALE
    twist_rate = 0.4 * 2 / 2.0 * SCALE
    screw_len = 2 * SCALE
    return screw.Screw(a, w, u, screw_len, r0, delta, twist_rate)   # phi0) , phi0=0


##################################################
# Cylinder
##################################################

from vectorized import SimpleCylinder


def cyl1(scale_ignored):
    SCALE = 8.  # mm

    radius = 0.5 * SCALE
    c_len = 2 * SCALE

    A = make_vector4(0, 0, -c_len / 2.0)  # bug: aa is wrong
    w = make_vector4(0, 1, 0)
    w = w / np.linalg.norm(w[0:3]); w[3] = 1
    u = make_vector4(1, 0, 0)

    c = SimpleCylinder(A, w, u, radius, radius, c_len)
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-16, +32, 1.92 * 0.2 * 10 / 2.0)

    return c, (RANGE_MIN, RANGE_MAX, STEPSIZE)


def rotation_matrix(angle, along, units="rad", around=None):
        import transformations as tf
        if units == "rad":
            pass
        elif units == "deg":
            angle = angle / 360.0 * np.pi * 2.0
        else:
            raise UsageError()

        if around is None:
            # center = 0
            # t = np.eye(4)
            pass
        else:
            check_vector4(around)
            center = around
            # print "center=",center
            t = np.eye(4)
            t[:, 3] = center
            assert t[3, 3] == 1.
            tinv = np.eye(4)
            tinv[:, 3] = -center
            tinv[3, 3] = 1
            # print "t", t
            # print "tinv", tinv

        check_vector4(along)
        rm = tf.rotation_matrix(angle, along[0:3])

        if around is None:
            return rm
        else:
            return np.dot(np.dot(rm, tinv), t)
        # return rm
        # matrix = np.dot(rm , self.matrix)
        # ...
        # return self


def cyl2(scale_ignored):
    SCALE = 8.  # mm

    un = None
    M = 5 + 1
    for i in range(M):
        radius = 0.5 * SCALE
        c_len = 2 * SCALE * 1. / (i + 1)
        a_x = (i - (M - 1.) / 2.) * 1 * SCALE + -0.5 * c_len

        A = make_vector4(a_x, 0, 0)
        w = make_vector4(0, 1, 0)
        # m = np.eye(4)

        center = A  # + make_vector4(0, 0, c_len/2.0)
        center[3] = 1
        R = rotation_matrix(i * 90. / (M - 1.), make_vector4(1, 0, 0), units="deg", around=center)

        # Rinf = R
        # for i in range(15):
        #    Rinf = np.dot(Rinf,Rinf)
        # print "Rinf", Rinf
        # print "w1", w, np.linalg.norm(w[:3])
        # print R
        w = np.dot(R, np.transpose(w))
        # print "w2", w, np.linalg.norm(w[:3])
        # u = make_vector4(0, 0, 1)   #will not work well for 90deg
        u = make_vector4(1, 0, 0)

        def make_uv(w, u):
            w = w / np.linalg.norm(w[0:3]); w[3] = 1
            u = u / np.linalg.norm(u[0:3]); u[3] = 1
            # print w, "u=",u
            assert np.linalg.norm(np.cross(w[:3], u[:3])) > 0.000000001  # cannot be parallel
            v3 = np.cross(w[:3], u[:3])
            v3 = v3 / np.linalg.norm(v3[:3])
            v = make_vector4(v3[0], v3[1], v3[2])
            u3 = np.cross(v[:3], w[:3])
            u3 = u3 / np.linalg.norm(u3[:3])
            u = make_vector4(u3[0], u3[1], u3[2])
            # return w, u
            return u

        u = make_uv(w, u)

        print np.linalg.norm(np.cross(w[:3], u[:3]))
        print np.linalg.norm(np.cross(w[:3], u[:3])) > 0.000000001

        def set4th1(v):
            assert v.shape == (4,)
            v[3] = 1
            return v
        # c = SimpleCylinder(A, w, u, radius, radius, c_len)
        newlen = c_len * 5
        c1 = SimpleCylinder(set4th1(A - 1. * newlen / 2. * w), w, u, radius / 5., radius / 5., newlen)
        delta, twist_rate = radius * 0.2, 2
        from vectorized import Screw
        c = Screw(A[:3], w[:3], u[:3], c_len, radius, delta, twist_rate)

        c = vectorized.CrispUnion(c1, c)

        if un is None:
            un = c
        else:
            un = vectorized.CrispUnion(un, c)
        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32, +32, 1.92 / 4.0 * 1.5 / 1.5)
    return un, (RANGE_MIN, RANGE_MAX, STEPSIZE)


def cage_rods(rod_r, rod_len, cage_r, N):
    import math
    un = None
    for i in range(N):
        th = i / float(N) * np.pi * 2
        x, y = cage_r * math.sin(th), cage_r * math.cos(th)
        A = make_vector4(x, y, -rod_len / 2.)
        w = make_vector4(0, 0, 1)
        u = make_vector4(1, 0, 0)
        c = SimpleCylinder(A, w, u, rod_r, rod_r, rod_len)

        if un is None:
            un = c
        else:
            un = vectorized.CrispUnion(un, c)
    return un


def cyl3(scale):
    # cage
    rod_r = 1.
    cage_r = 10.
    rod_len = 10.
    N = 20
    return cage_rods(rod_r, rod_len, cage_r, N), (-32 / 2, +32 / 2, 1.92 / 4.0)


def cyl4(scale):
    """ Makes a nice cage with spiral bars. Don't change. """
    cage = cage_rods(rod_r=1, rod_len=20, cage_r=10, N=20)
    from twist_z import TwistZ
    t = TwistZ(cage, 0.02)  # cycles per mm
    # 0.06 is too much  0.02 is reasonable
    base_cyl = SimpleCylinder(
        make_vector4(0, 0, -10),  # A
        make_vector4(0, 0, -1),  # w
        make_vector4(1, 0, 0),
        11., 11.,
        1.)
    ifunc = vectorized.CrispUnion(base_cyl, t)

    # (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32, +32, 1.92 / 4.0)   #15 sec!  2.5 millions voxels
    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-32 / 2, +32 / 2, 1.92 / 4.0)  # 2.5 sec!

    return ifunc#, (RANGE_MIN, RANGE_MAX, STEPSIZE)

def cube_with_cylinders(scale):

    SCALE = 2.  # mm
    sz1 = 2.5

    radius = 0.5 * SCALE
    c_len = 2 * SCALE

    A = make_vector4(-c_len/2.0, 0, 0)
    #A = make_vector4(0, 0, c_len / 2.0)  # bug: aa is wrong
    w = make_vector4(1,0 , 0)
    w = w / np.linalg.norm(w[0:3]); w[3] = 1
    u = make_vector4(0, 1, 0)

    cyl = SimpleCylinder(A, w, u, radius, radius, c_len)


    A2 = make_vector4(0, -c_len/2.0, 0)
        #A = make_vector4(0, 0, c_len / 2.0)  # bug: aa is wrong
    w2 = make_vector4(1,0 , 0)
    w2 = w / np.linalg.norm(w[0:3]); w[3] = 1
    u2 = make_vector4(0, 1, 0)

    cyl_2 = SimpleCylinder(A2, u2, w2, radius, radius, c_len)
    cube = vectorized.UnitCube1(size=sz1)
    union = vectorized.CrispSubtract(cube, cyl_2)
    final_object = vectorized.CrispUnion(union,cyl)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-3, +5, 0.2)
    return final_object
    #, (RANGE_MIN, RANGE_MAX, STEPSIZE)

# **************************************************************************************************


""" contains safe functions to call for producing implicit objects for tests. """
# 2 = vectorized only, 3 = a pair of vectorized and non-vectorised are used.
examples = {

    "sphere_example": 2,
    "ell_example1": 2,
    "blend_example2": 2,
    "cube_example": 2,
    "blend_example2_discs": 2,
    "blend_example1": 2,
    "bowl_15_holes": 2,
    "first_csg": 2,
    "french_fries_vectorized": 2,
    "rdice_vec": 2,
    "rcube_vec": 2,
    "screw1": 2,
    "screw2": 2,
    "udice_vec": 2,
    "rods": 2,
    "cyl4": 2,       # spiral cage
    "cube_with_cylinders": 2,

}


def make_example_vectorized(name, scale=1.0):
    assert any(name == s for s in examples)
    assert examples[name] in [2], "Incorrect example type"
    res = globals()[name](scale)
    assert not type(res) is tuple
    return res


def test_creation_of_all_Examples():
    for name in examples:
        scale = 1.
        res = globals()[name](scale)
        # print("OK.")


def get_all_examples(types_list):
    """ types_list i.e. [1] or [2] or [2,3] or [1,3] r [1,2,3] """
    usable_examples = []
    i = 0
    for e in examples:
        if examples[e] in types_list:
            usable_examples += [e]
            # print("e=", e)

            if examples[e] in [2]:
                iobj = make_example_vectorized(e)
                x = vectorized.repeat_vect4(1, make_vector4(0.5, 0.5, 0.5))
                g = iobj.implicitGradient(x)
                v = iobj.implicitFunction(x)

        i += 1
    assert i > 0
    return usable_examples
