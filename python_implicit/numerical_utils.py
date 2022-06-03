import numpy as np
from basic_types import check_vector4, repeat_vect4, make_vector4
import vectorized
import nonvec

def numerical_gradient(iobj, pos0, delta_t=0.01/10.0/10.0, order=5, is_vectorized="unspecified"):
    #0.1 is not enough for delta_t
    assert is_vectorized != "unspecified"
    if is_vectorized:
        check_vector4(pos0)  # incorrect
        #assert False
        #check_vector4_vectorized(pos0)
        assert issubclass(type(iobj), vectorized.ImplicitFunctionVectorized)
    else:
        check_vector4(pos0)
        assert issubclass(type(iobj), nonvec.ImplicitFunctionPointwise)

    m = order  # sample points: -m,...,-1,0,1,2,...,+m

    _VERBOSE = False
    from lib import finite_diff_weights

    sample_points = range(-m, m+1)
    n = m*2+1

    x0 = 0
    findiff_weights = finite_diff_weights.weights(k=1, x0=x0, xs=np.array(sample_points) * delta_t)
    del x0

    pos0_4 = repeat_vect4(1, pos0)
    pos = np.tile(pos0_4, (3*n, 1))
    assert not issubclass(pos.dtype.type, np.integer)

    if pos.shape[0] in [1, long(1718772)]:  # 1718772L
        set_trace()

    dx = repeat_vect4(1, make_vector4(1, 0, 0))
    dy = repeat_vect4(1, make_vector4(0, 1, 0))
    dz = repeat_vect4(1, make_vector4(0, 0, 1))
    dxyz = [dx, dy, dz]

    ci = 0
    for d in range(3):
        for i in sample_points:
            dd = dxyz[d]

            if pos.shape[0] == 1 and ci ==1:
                set_trace()

            pos[ci, :] = pos[ci, :] + (dd * delta_t * float(i))
            #w[ci] = findef(i,n)
            ci += 1

    pos[:, 3] = 1

    if is_vectorized:
        v = iobj.implicitFunction(pos)    # v .shape: (3,11)
    else:
        v = np.zeros((pos.shape[0],))
        for i in range(pos.shape[0]):
            v1 = iobj.implicitFunction(pos[i,:])    # v .shape: (3,11)
            v[i] = v1
    v3 = np.reshape(v, (3, n), order='C')  # v3 .shape: (11,)
    #print( np.diff(v, axis=0) / delta_t )
    #print( np.diff(v3, axis=1) / delta_t )



    #Lipchitz_L
    #Lipschitz constant = Lipchitz_B
    Lipchitz_B = 50
    Lipchitz_beta = 1  # order. Keep it 1
    #H\:older continuous:   |f(h)-f(0)| <= B|h|^beta
    b_h_beta = Lipchitz_B*(np.abs(delta_t)**Lipchitz_beta)
    #print("v3=",v3)
    #print("diff=",np.diff(v3, axis=1) )
    #print(b_h_beta)

    d0 = np.abs( np.diff(v3, axis=1) )
    #lipschitz_condition = d <= b_h_beta
    nonsmooth_ness = d0 / (np.abs(delta_t)**Lipchitz_beta)
    #nonsmooth_ness2 = np.mean(nonsmooth_ness, axis=1)

    d = np.abs( np.diff(v3, n=1, axis=1) ) / np.abs(delta_t)
    d = d - np.tile( np.mean(d, axis=1, keepdims=True), (1, d.shape[1]) )
    d = np.abs(d) / np.abs(delta_t)
    d = d - np.tile( np.mean(d, axis=1, keepdims=True), (1, d.shape[1]) )
    #d = np.abs(d)
    #d = np.abs( np.diff(d, n=1, axis=1) ) / np.abs(delta_t)
    #nonsmooth_ness = d / (np.abs(delta_t)**Lipchitz_beta)
    #if(np.max(np.ravel(d))) > 50:
    #    print("warning")
    #    print(nonsmooth_ness)

    #print(d)
    if(np.max(np.ravel(nonsmooth_ness))) > 100*10:
        print("warning: nonsmooth ",end='')
        #print(nonsmooth_ness)  # lots of zeros and one big value

    """ Calculating the numerical derivative using finite difference (convolution with weights) """
    #convolusion
    grad_cnv = np.dot(v3, findiff_weights)
    #grad_cnv = np.reshape(grad_cnv, (1,3))[:,np.newaxis]
    #grad_cnv = np.concatenate( ( grad_cnv[np.newaxis,:], np.reshape(np.array([1]),(1,1)) ), axis=1)

    def v3_to_v14(v):
        """ Converts shape from (3,) into a (1,4) vector4 """
        assert v.ndim == 1
        return np.concatenate((v[np.newaxis, :], np.reshape(np.array([1]), (1, 1))), axis=1)

    grad_cnv = v3_to_v14(grad_cnv)
    #print("weights: ",findiff_weights)


    #Detecting sharp edges (non-smooth points, i.e. corners and edges and ridges)
    if np.max(np.abs(grad_cnv)) > 100:
        pass
        #print("*******  max(grad) > 100")
        #print(np.abs(grad_cnv))
    #else:
    #    print(np.abs(grad_cnv))

    if _VERBOSE:
        #np.set_printoptions( precision=9 )
        np.set_printoptions(formatter={'all': lambda x: ''+("%2.19f" % (x,))})

    """ Calculating the numerical derivative using 'mean of diff' """
    grad_mean = np.mean(-np.diff(v3, axis=1) / delta_t, axis=1)
    if _VERBOSE:
        print("grad_mean: ", grad_mean)
        print("grad_convolusion: ", grad_cnv)

    if False:
        g = iobj.implicitGradient(pos0_4)
    if _VERBOSE:
        print("grad_analytical: ", g)

        #print( grad_cnv.shape )
        #print( g.shape )
        #print( grad_mean.shape )

        print("Errors:")
        print("conv error: ", g - grad_cnv)
        #Amazing precision: [[0.0000000000001995071 -0.0000000038590481921 0.0000000000000008882  0.0000000000000000000]]

        print("mean error: ", g - v3_to_v14(grad_mean))
        #Terrible error: [-0.262   2.12266  0 ]

        #v3 * findiff_weights

        print("to be continued")

    assert not np.any(np.isnan(grad_cnv.ravel()))
    return grad_cnv

def numerical_gradient_slow_func(self_iobj, x):
    """ Used as a slow implicitGradient(self, x) when too lasy to derive the gradient! see classes Screw, RSubtract"""
    check_vector4_vectorized(x)
    count = x.shape[0]
    g = np.zeros((count,4))
    for i in range(x.shape[0]):
        v = x[i, 0:4]
        #inefficient: not vectorised
        g[i,:] = numerical_gradient(self_iobj, v, is_vectorized=True)
    return g


def optimize_vertex(vert, iobj, radius_of_max_change):
    v = make_vector4(vert[0], vert[1], vert[2])  # inefficient
    check_vector4(v)

    #f = iobj.implicitFunction(v)
    #g = iobj.implicitGradient(v)
    #v += g * np.random.rand() * radius_of_max_change * 2
    #v[3] = 1

    iterations = 2 #14
    for i in range(iterations):

        f = iobj.implicitFunction(v)
        g = iobj.implicitGradient(v)

        tau  = 0.1  # * 3
        a = 1
        z_force = -tau * a * f * g

        v += z_force
        v[3] = 1

    return v[0:3]

def optimize_vertices(verts, iobj, radius_of_max_change):
    for i in range(verts.shape[0]):
        verts[i,:] = optimize_vertex(verts[i,:], iobj, radius_of_max_change)
    return verts

def error_i(verts, iobj, type='sqr'):
    n = verts.shape[0]
    print(n)
    e = np.zeros((n,))
    for i in range(n):
        v = verts[i, :]
        v = make_vector4(v[0], v[1], v[2])  # inefficient
        #check_vector4(v)
        e[i] = iobj.implicitFunction(v)
    #print(np.mean(np.abs(e[i])))
    #print(np.max(np.abs(e[i])))
    if type == 'sqr':
        return np.sqrt(np.mean(e[i]**2))
    elif type == 'abs':
        return np.mean(np.abs(e[i]))
    else:
        raise "wrong type"


def error_i_vectorized(verts, iobj, type='sqr'):
    #assert is iobj is vectorized
    n = verts.shape[0]
    print(n)
    e = np.zeros((n,))
    for i in range(n):
        v = verts[i, :]
        v = make_vector4(v[0], v[1], v[2])  # inefficient
        #check_vector4(v)
        e[i] = iobj.implicitFunction( v[np.newaxis,:] )
    #print(np.mean(np.abs(e[i])))
    #print(np.max(np.abs(e[i])))
    if type == 'sqr':
        return np.sqrt(np.mean(e[i]**2))
    elif type == 'abs':
        return np.mean(np.abs(e[i]))
    else:
        raise "wrong type"


ROOT_TOLERANCE = 0.000001  # root tolerance
def mysign_np(v):
    return np.sign(v) * (np.abs(v) > ROOT_TOLERANCE)

from basic_types import check_vector4_vectorized
def numerical_raycast_bisection_vectorized(iobj, ray_x, ray_n, ROOT_TOLERANCE=ROOT_TOLERANCE):
    """ ray_x must be outside and ray_x+ray_n must be inside the object. Then this function finds points x=ray_x+(lambda)*ray_n where f(x)=0 using the bisection method."""
    check_vector4_vectorized(ray_x)
    check_vector4_vectorized(ray_n)
    assert ray_x.shape[0] == ray_n.shape[0]
    x1_arr = ray_x  # start
    v1_arr = iobj.implicitFunction(x1_arr)
    x2_arr = x1_arr + ray_n * 1.0
    x2_arr[:,3] = 1
    v2_arr = iobj.implicitFunction(x2_arr)

    result_x_arr = np.zeros(ray_x.shape)

    EPS = 0.0000001  # sign

    n = x1_arr.shape[0]
    already_root = np.zeros((n,), dtype=np.int)
    #assert v2_arr * va > 0 - EPS  # greater or equal
    active_indices  = np.arange(0,n)  # mid
    iteration = 1
    while True:
        #print(v1_arr.shape, "*")
        #print(np.vstack((v1_arr,v2_arr)))
        #print(v2_arr * v1_arr )
        #print("aaaaaaaaaaaaaaaa")
        #print(v1_arr)
        #print(v2_arr)
        #print(v2_arr * v1_arr)
        #assert np.all(v2_arr * v1_arr < 0 - EPS)  # greater or equal
        assert np.all(mysign_np(v2_arr) * mysign_np(v1_arr) < 0 - EPS)  # greater or equal

        assert np.all(v1_arr < 0-ROOT_TOLERANCE)
        assert active_indices.shape[0] == x1_arr.shape[0]
        assert active_indices.shape[0] == x2_arr.shape[0]
        x_mid_arr = ( x1_arr + x2_arr ) / 2.0
        x_mid_arr[:,3] = 1
        v_mid_arr = iobj.implicitFunction(x_mid_arr)
        assert v_mid_arr.shape == active_indices.shape
        assert active_indices.ndim == 1

        #flipped_i = v_mid_arr
        #contains the indices
        dif = -v_mid_arr  # assuming x1 is always outside and x2 is inside
        assert dif.shape == active_indices.shape
        boolean_eq = np.abs(dif) <= ROOT_TOLERANCE
        boolean_gt =  dif >  ROOT_TOLERANCE
        boolean_lt =  dif <  -ROOT_TOLERANCE # dif <  ROOT_TOLERANCE
        boolean_neq = np.logical_not( boolean_eq )
        #logical_or
        assert np.all( np.logical_or(boolean_gt, boolean_lt) == np.logical_not(boolean_eq) )
        #print("boolean_gt", boolean_gt)  #t
        #print("boolean_lt", boolean_lt)  #f


        which_zeroed     = active_indices[ boolean_eq ] # new start = mid
        which_flippedAt1 = active_indices[ boolean_gt ] # new end = mid
        which_flippedAt2 = active_indices[ boolean_lt ]
        which_flippedAny = active_indices[ boolean_neq ]

        already_root[which_zeroed] = 1  # iteration
        result_x_arr[which_zeroed,:] = x_mid_arr[boolean_eq,:]


        #x1_arr and x2_arr should have the same size eventually. the boolean_eq should be removed from their indices.
        #the total is np.arange(n)
        v2_arr[boolean_lt] = v_mid_arr[boolean_lt]#[which_flippedAny]
        x2_arr[boolean_lt,:] = x_mid_arr[boolean_lt,:]#[which_flippedAny]   # which_flippedAt2

        #x1_arr and x2_arr both shrink here

        v1_arr[boolean_gt] = v_mid_arr[boolean_gt]#[which_flippedAny]
        x1_arr[boolean_gt,:] = x_mid_arr[boolean_gt,:]#[which_flippedAny]   #which_flippedAt1

        v1_arr = v1_arr[boolean_neq]
        v2_arr = v2_arr[boolean_neq]
        x1_arr = x1_arr[boolean_neq,:]
        x2_arr = x2_arr[boolean_neq,:]
        #print("active_indices = ", active_indices)
        #print("which_flippedAny = ", which_flippedAny)
        #print("active_indices[which_flippedAny] = ", active_indices[boolean_neq])
        active_indices = active_indices[boolean_neq] #which_flippedAt1 || which_flippedAt2
        iteration += 1

        assert x1_arr.shape == x2_arr.shape
        assert v1_arr.shape == v2_arr.shape
        assert active_indices.shape == v1_arr.shape

        #print(active_indices)
        #print(v1_arr, "****", v2_arr)
        #print(boolean_lt)
        #print("*******")

        if len(active_indices) == 0:
            break

    assert len(active_indices) == 0
    #result_x_arr
    v_arr = iobj.implicitFunction(result_x_arr)
    assert np.all(np.abs(v_arr) < ROOT_TOLERANCE)
    return result_x_arr



# def svd3_by_khaled():
#     #from https://gitlab.kitware.com/chunmingchen/MeshSimplifySerial/blob/master/mymath.h
#     #from: http://stackoverflow.com/questions/4372224/fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition
#     # Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
#     # source: http://www.melax.com/diag.html?attredirects=0
#     #void Diagonalize(const Real (&A)[3][3], Real (&Q)[3][3], Real (&D)[3][3])
# {
#     // A must be a symmetric matrix.
#     // returns Q and D such that
#     // Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
#     const int maxsteps=24;  // certainly wont need that many.
#     int k0, k1, k2;
#     Real o[3], m[3];
#     Real q [4] = {0.0,0.0,0.0,1.0};
#     Real jr[4];
#     Real sqw, sqx, sqy, sqz;
#     Real tmp1, tmp2, mq;
#     Real AQ[3][3];
#     Real thet, sgn, t, c;
#     for(int i=0;i < maxsteps;++i)
#     {
#         // quat to matrix
#         sqx      = q[0]*q[0];
#         sqy      = q[1]*q[1];
#         sqz      = q[2]*q[2];
#         sqw      = q[3]*q[3];
#         Q[0][0]  = ( sqx - sqy - sqz + sqw);
#         Q[1][1]  = (-sqx + sqy - sqz + sqw);
#         Q[2][2]  = (-sqx - sqy + sqz + sqw);
#         tmp1     = q[0]*q[1];
#         tmp2     = q[2]*q[3];
#         Q[1][0]  = 2.0 * (tmp1 + tmp2);
#         Q[0][1]  = 2.0 * (tmp1 - tmp2);
#         tmp1     = q[0]*q[2];
#         tmp2     = q[1]*q[3];
#         Q[2][0]  = 2.0 * (tmp1 - tmp2);
#         Q[0][2]  = 2.0 * (tmp1 + tmp2);
#         tmp1     = q[1]*q[2];
#         tmp2     = q[0]*q[3];
#         Q[2][1]  = 2.0 * (tmp1 + tmp2);
#         Q[1][2]  = 2.0 * (tmp1 - tmp2);

#         // AQ = A * Q
#         AQ[0][0] = Q[0][0]*A[0][0]+Q[1][0]*A[0][1]+Q[2][0]*A[0][2];
#         AQ[0][1] = Q[0][1]*A[0][0]+Q[1][1]*A[0][1]+Q[2][1]*A[0][2];
#         AQ[0][2] = Q[0][2]*A[0][0]+Q[1][2]*A[0][1]+Q[2][2]*A[0][2];
#         AQ[1][0] = Q[0][0]*A[0][1]+Q[1][0]*A[1][1]+Q[2][0]*A[1][2];
#         AQ[1][1] = Q[0][1]*A[0][1]+Q[1][1]*A[1][1]+Q[2][1]*A[1][2];
#         AQ[1][2] = Q[0][2]*A[0][1]+Q[1][2]*A[1][1]+Q[2][2]*A[1][2];
#         AQ[2][0] = Q[0][0]*A[0][2]+Q[1][0]*A[1][2]+Q[2][0]*A[2][2];
#         AQ[2][1] = Q[0][1]*A[0][2]+Q[1][1]*A[1][2]+Q[2][1]*A[2][2];
#         AQ[2][2] = Q[0][2]*A[0][2]+Q[1][2]*A[1][2]+Q[2][2]*A[2][2];
#         // D = Qt * AQ
#         D[0][0] = AQ[0][0]*Q[0][0]+AQ[1][0]*Q[1][0]+AQ[2][0]*Q[2][0];
#         D[0][1] = AQ[0][0]*Q[0][1]+AQ[1][0]*Q[1][1]+AQ[2][0]*Q[2][1];
#         D[0][2] = AQ[0][0]*Q[0][2]+AQ[1][0]*Q[1][2]+AQ[2][0]*Q[2][2];
#         D[1][0] = AQ[0][1]*Q[0][0]+AQ[1][1]*Q[1][0]+AQ[2][1]*Q[2][0];
#         D[1][1] = AQ[0][1]*Q[0][1]+AQ[1][1]*Q[1][1]+AQ[2][1]*Q[2][1];
#         D[1][2] = AQ[0][1]*Q[0][2]+AQ[1][1]*Q[1][2]+AQ[2][1]*Q[2][2];
#         D[2][0] = AQ[0][2]*Q[0][0]+AQ[1][2]*Q[1][0]+AQ[2][2]*Q[2][0];
#         D[2][1] = AQ[0][2]*Q[0][1]+AQ[1][2]*Q[1][1]+AQ[2][2]*Q[2][1];
#         D[2][2] = AQ[0][2]*Q[0][2]+AQ[1][2]*Q[1][2]+AQ[2][2]*Q[2][2];
#         o[0]    = D[1][2];
#         o[1]    = D[0][2];
#         o[2]    = D[0][1];
#         m[0]    = fabs(o[0]);
#         m[1]    = fabs(o[1]);
#         m[2]    = fabs(o[2]);

#         k0      = (m[0] > m[1] && m[0] > m[2])?0: (m[1] > m[2])? 1 : 2; // index of largest element of offdiag
#         k1      = (k0+1)%3;
#         k2      = (k0+2)%3;
#         if (o[k0]==0.0)
#         {
#             break;  // diagonal already
#         }
#         thet    = (D[k2][k2]-D[k1][k1])/(2.0*o[k0]);
#         sgn     = (thet > 0.0)?1.0:-1.0;
#         thet   *= sgn; // make it positive
#         t       = sgn /(thet +((thet < 1.E6)?sqrt(thet*thet+1.0):thet)) ; // sign(T)/(|T|+sqrt(T^2+1))
#         c       = 1.0/sqrt(t*t+1.0); //  c= 1/(t^2+1) , t=s/c
#         if(c==1.0)
#         {
#             break;  // no room for improvement - reached machine precision.
#         }
#         jr[0 ]  = jr[1] = jr[2] = jr[3] = 0.0;
#         jr[k0]  = sgn*sqrt((1.0-c)/2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)
#         jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
#         jr[3 ]  = sqrt(1.0f - jr[k0] * jr[k0]);
#         if(jr[3]==1.0)
#         {
#             break; // reached limits of floating point precision
#         }
#         q[0]    = (q[3]*jr[0] + q[0]*jr[3] + q[1]*jr[2] - q[2]*jr[1]);
#         q[1]    = (q[3]*jr[1] - q[0]*jr[2] + q[1]*jr[3] + q[2]*jr[0]);
#         q[2]    = (q[3]*jr[2] + q[0]*jr[1] - q[1]*jr[0] + q[2]*jr[3]);
#         q[3]    = (q[3]*jr[3] - q[0]*jr[0] - q[1]*jr[1] - q[2]*jr[2]);
#         mq      = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
#         q[0]   /= mq;
#         q[1]   /= mq;
#         q[2]   /= mq;
#         q[3]   /= mq;
#     }
# }


#todo: move into mesh_utils


def average_edge_size(verts, faces):
    # slow version
    p = 0.
    for fi in range(faces.shape[0]):
        p += calculate_perimeter(verts, faces, fi)
    return p / (float(faces.shape[0]) * 3.)


def calculate_perimeter(verts, faces, fi):
    e = np.zeros((3,))
    for j in range(3):
        j2 = (j+1) % 3
        assert j2 >= 0
        assert j2 < 3
        v1 = faces[fi, j]
        v2 = faces[fi, j2]
        e[j] = np.linalg.norm(verts[v1, :]-verts[v2, :])
    return np.sum(e)


def my_process_1(verts, faces, iobj):
    #return verts
    import mesh1
    m = mesh1.Mesh_1(faces, verts)
    m.build_centroids()
    m.build_neighbours()
    m.evaluate_centroid_gradients(iobj)

    # m.update_centroids_and_gradients(iobj)
    fi = 0
    mean_edge = calculate_perimeter(verts, faces, fi)
    c = m.centroids[fi, :]
    print(mean_edge)
    from mesh_utils import project_single_point2_ohtake
    p = project_single_point2_ohtake(iobj, start_x=c.reshape(1,4), lambda_val=0.5*mean_edge, max_dist=10)
    print(p)
    print("OK.")
    exit()


    #if not qem_breakdown:
    m.quadratic_optimise_vertices(1)
    m.verts = m.new_verts
    return m.verts

def test_search1():
    import math
    from stl_tests import make_mc_mesh_scikit

    dicesize = 8.
    exname = "udice_vec"
    import example_objects
    iobj = example_objects.make_example_vectorized(exname, dicesize)

    (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-11, +10., 0.8)  # non-spiky!!
    verts, faces = make_mc_mesh_scikit(iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)

    verts = my_process_1(verts, faces, iobj)

    from stl_tests import display_simple_using_mayavi_vf1
    display_simple_using_mayavi_vf1(verts, faces)
    print(np.min(verts.ravel()), np.max(verts.ravel()))

    plot_stlmesh(m)


def cubic_root_x1():

    def c_real_roots1(a, b, c, d):
        import math
        """ Return x1 only.
        todo: three roots. (all real roots)
        https://en.wikipedia.org/wiki/Cubic_function#Root-finding_formula
        """
        if math.fabs(a) < 0.0000000000000001:
            return quadratic_roots(b, c, d)

        q = (3*a*c - b**2) / (9 * a**2)
        r = (9*a*b*c - 27*a**2*d - 2*b**3) / (54*a**3)

        print("q = ",q)
        print("r = ",r)

        delta = q**3 + r**2


        #u1 = 1
        #u2 = (-1 - 1i*math.sqrt(3))/2
        #u3 = (-1 - 1i*math.sqrt(3))/2
        #
        #if delta > 0:
        #    roots = (1,1,1)
        #elif delta ==0:
        #    roots = (2,1)
        #    multiplicity = ()
        #else: #delta<0
        #    roots = (1,2)  #second one is complex

        print("delta = ", delta)

        # here delta is less than zero so we use the second set of equations from the article:

        rho = (-q**3)**0.5

        #incorrect
        #s_real = rho**(1./3.)
        #t_real = rho**(1./3.)

        #corrected: http://stackoverflow.com/questions/1829330/solving-a-cubic-equation
        theta = math.acos(r/rho)
        s_real = rho**(1./3.) * math.cos(+theta/3)
        t_real = rho**(1./3.) * math.cos(-theta/3)


        print("s [real] = ",s_real)
        print("t [real] = ",t_real)

        x1 = s_real + t_real - b / (3. * a)
        return x1

    a = 1.0
    b = 0.0
    c = 0.2 - 1.0
    d = -0.7 * 0.2

    x1 = c_real_roots1(a, b, c, d)

    print("x1 = ", x1)

    print("should be zero: ",a*x1**3+b*x1**2+c*x1+d)

from utils import flush, optimised_used

TEST_ON = not optimised_used()

import sys
def numerical_gradient_vectorized_v1(iobj, x):
    check_vector4_vectorized(x)
    if TEST_ON:
        print("numerical_gradient1",end=''); flush()
        count = x.shape[0]
        g = np.zeros((count, 4))
        for i in range(x.shape[0]):
            v = x[i, 0:4]
            # inefficient: not vectorised
            g[i, :] = numerical_gradient(iobj, v, is_vectorized=True)
        assert not np.any(np.isnan(g), axis=None)

    print("numerical_gradient2",end=''); flush()
    g2 = numerical_gradient_vectorized_v2(iobj, x.copy())
    if TEST_ON:
        assert np.allclose(g, g2)
    print("done"); flush()
    return g2


from lib import finite_diff_weights
def numerical_gradient_vectorized_v2(iobj, pos0, delta_t=0.01/100., order=5):
    """ A proper vectorized implementation. See numerical_gradient() """
    #Note: 0.1 is not enough for delta_t

    check_vector4_vectorized(pos0)
    assert issubclass(type(iobj), vectorized.ImplicitFunctionVectorized)
    assert pos0.ndim == 2
    if pos0.shape[0] == 0:
        return np.zeros((0, 4))

    m = order  # sample points: -m,...,-1,0,1,2,...,+m

    _VERBOSE = False

    sample_points = range(-m, m+1)
    n = m*2+1

    x0 = 0
    findiff_weights = finite_diff_weights.weights(k=1, x0=x0, xs=np.array(sample_points) * delta_t)
    del x0

    assert n < 20
    pos0_4 = pos0[:, np.newaxis, :]
    pos = np.tile(pos0_4, (1, 3*n, 1))
    assert not issubclass(pos.dtype.type, np.integer)

    #if pos.shape[0] in [1, 1718772L]:
    #    set_trace()

    dx = make_vector4(1, 0, 0)[np.newaxis, np.newaxis, :]
    dy = make_vector4(0, 1, 0)[np.newaxis, np.newaxis, :]
    dz = make_vector4(0, 0, 1)[np.newaxis, np.newaxis, :]
    dxyz = [dx, dy, dz]

    ci = 0
    for d in range(3):
        dd = dxyz[d]
        for i in sample_points:
            #if pos.shape[0] == 1 and ci ==1:
            #    set_trace()

            pos[:, ci, :] = pos[:, ci, :] + (dd * (delta_t * float(i)))
            #w[ci] = findef(i,n)
            assert ci < 3*n
            ci += 1

    pos[:, :, 3] = 1

    vsize = pos0.shape[0]

    v = iobj.implicitFunction(pos.reshape((vsize * 3*n), 4))    # v .shape: (3,11)
    v3 = np.reshape(v, (vsize, 3, n), order='C')  # v3 .shape: (11,)


    if True:
        #Lipchitz_L
        Lipchitz_B = 50  # Lipschitz constant
        Lipchitz_beta = 1  # order. Keep it 1
        #H\:older continuous:   |f(h)-f(0)| <= B|h|^beta
        b_h_beta = Lipchitz_B*(np.abs(delta_t)**Lipchitz_beta)
        #lipschitz_condition = d <= b_h_beta

        d0 = np.abs( np.diff(v3, axis=1+1) )
        nonsmooth_ness = d0 / (np.abs(delta_t)**Lipchitz_beta)
        #nonsmooth_ness2 = np.mean(nonsmooth_ness, axis=1)

        del d0, b_h_beta, Lipchitz_beta, Lipchitz_B

        #print nonsmooth_ness.shape
        #set_trace()
        if(np.max(np.ravel(nonsmooth_ness))) > 100*10:
            print("warning: nonsmooth ",end='')
        del nonsmooth_ness

    if False:
        # v3: (vsize x 3 x n)
        #not tested
        d = np.abs( np.diff(v3, n=1, axis=1+1) ) / np.abs(delta_t)
        d = d - np.tile( np.mean(d, axis=1+1, keepdims=True), (1, 1, d.shape[1+1]) )
        d = np.abs(d) / np.abs(delta_t)
        d = d - np.tile( np.mean(d, axis=1+1, keepdims=True), (1, 1, d.shape[1+1]) )
        del d


    """ Calculating the numerical derivative using finite difference (convolution with weights) """
    #convolusion
    grad_cnv = np.dot(v3, findiff_weights)  # "sum product over the last axis of a and the second-to-last of b"

    def v3v_to_v14(v):
        """ Converts shape from (N,3,) into a (N,4) vector4 """
        assert v.ndim == 1+1
        return np.concatenate((v[:, :], np.ones((v.shape[0], 1), dtype=float)), axis=1)

    grad_cnv = v3v_to_v14(grad_cnv)
    assert not np.any(np.isnan(grad_cnv), axis=None)
    return grad_cnv



if __name__ == '__main__':
    test_search1()  # tests my_process_1()
    #cubic_root_x1()
    pass
