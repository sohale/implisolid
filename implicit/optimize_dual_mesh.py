import numpy as np
import mesh_utils
import profile_support
from mayavi import mlab
# from ipdb import set_trace
"""

Optimizing dual mesh
Steps:
Consider the dual mesh formed by triangle centroids.We are
looking for an optimized dual mesh whose vertices are obtained
by projecting the centroids C onto the implicit surface.The projection
P of triangle centroid C onto the surface f(x,y,z) = 0 is estimated as follows

    1. Set P = C
    2. If |f(P)| < tolerance: terminate
        else:
        Q = P
    3. Set R = Q + lambda*d , where d = -f(Q) * gradientf(Q) (vector) and
    lambda is a small constant
    4. if f(Q) * f(R) < 0 : find_root : condition f(P) < tolerance inside [Q,R]
    and terminate
    else :  Set Q = R  and go to step 3 .

"""


class MeshOptimizer(object):

    def __init__(self, centroids=None, triangles=None):
        self.lamda = 0.15
        self.TOLERANCE = 0.001
        self.centroids = centroids
        self.optimized_dual_mesh = []
        self.faces = []
        self.vertices = []
        self.optimized_vertices = []

    def build_neighbours(self):
        self.vertex_neighbours_list = mesh_utils.make_neighbour_faces_of_vertex(self.faces)

    def register_function(self, fnc):
        self.function = fnc

    def register_gradient(self, grd):
        self.gradient = grd

    @profile
    def load_example(self, example_name="ell_example1", res=1):
        """
            res: controls step size of mc grid
            example_name: name of example to build, available names can be
            found in the file example_objects
        """
        print "[Start] Loading example:  %s" % example_name + ' ...'
        from example_objects import make_example_vectorized
        from stl_tests import make_mc_mesh_scikit
        exname = example_name
        self.iobj = make_example_vectorized(exname)
        self.register_function(self.iobj.implicitFunction)
        self.register_gradient(self.iobj.implicitGradient)

        (RANGE_MIN, RANGE_MAX, STEPSIZE) = (-2., 2., 0.1 / res)  # non-spiky!!

        self.vertices, self.faces = make_mc_mesh_scikit(self.iobj, RANGE_MIN, RANGE_MAX, STEPSIZE)
        self.centroids = np.mean(self.vertices[self.faces[:], :], axis=1)
        self.centroids = np.c_[self.centroids, np.ones((len(self.centroids), 1))]
        self.triangles = self.vertices[self.faces]
        self.build_neighbours()
        self.name = exname
        self.lamda = calc_avg_triangle_length(self.triangles) / 2
        print "[End] Loading example:  %s" % example_name + ' .'


    @profile
    def plot_mesh(self, primal_mesh=True, axisOn=True, representation="surface"):
        print "[Start]: Plotting mesh for:  %s" % self.name + ' ...'

        axis_dim = 2
        # Axes
        x = np.linspace(-axis_dim, axis_dim, 2).reshape(2, 1)
        y = np.zeros((2, 1))
        z = np.zeros((2, 1))
        mlab.plot3d(x, y, z, line_width=1, name="x-axis", opacity=0.5)
        mlab.plot3d(y, x, z, line_width=1, name="y-axis", opacity=0.5)
        mlab.plot3d(z, y, x, line_width=1, name="z-axis", opacity=0.5)
        # Axes labels
        # mlab.text3d(0,0,0,"O",scale=.5)
        mlab.text3d(axis_dim, 0, 0, 'X', scale=.1)
        mlab.text3d(0, axis_dim, 0, 'Y', scale=.1)
        mlab.text3d(0, 0, axis_dim, 'Z', scale=.1)
        mlab.text3d(-axis_dim, 0, 0, '-X', scale=.1)
        mlab.text3d(0, -axis_dim, 0, '-Y', scale=.1)
        mlab.text3d(0, 0, -axis_dim, '-Z', scale=.1)

        # Origin
        mlab.points3d(0, 0, 0, color=(0, 0, 0), scale_factor=0.1)
        x = np.zeros(shape=(len(self.vertices), 1))
        y = np.zeros(shape=(len(self.vertices), 1))
        z = np.zeros(shape=(len(self.vertices), 1))
        if primal_mesh:
            for i in range(len(self.vertices)):
                x[i] = self.vertices[i][0]
                y[i] = self.vertices[i][1]
                z[i] = self.vertices[i][2]
            mlab.triangular_mesh(x, y, z, self.faces, representation=representation)
        else:
            for i in range(len(self.vertices)):
                y[i] = self.optimized_vertices[i][1]
                x[i] = self.optimized_vertices[i][0]
                z[i] = self.optimized_vertices[i][2]
            mlab.triangular_mesh(x, y, z, self.faces, representation=representation)
        print "[End]: Plotting mesh for: %s" % self.name + ' .'
        # mlab.show()

    @profile
    def find_bisection_root(self, a, b, val_a=None, val_b=None):
        """
        Bisection method in the line that connects the points a, b .
        Each time we update one of lthe edge points to the median in the line
        """
        a = np.array(a, dtype=np.float32)
        b = np.array(b, dtype=np.float32)
        median = (a + b) / 2.0
        val_a = val_a or self.function(a)[0]
        val_b = val_b or self.function(b)[0]
        val_md = self.function(median)[0]

        if abs(val_a) < self.TOLERANCE:
            return a
        if abs(val_b) < self.TOLERANCE:
            return b
        if abs(val_md) < self.TOLERANCE:
            return median

        if val_md * val_a < 0:
            return self.find_bisection_root(a, median, val_a=val_a)
        elif val_md * val_b < 0:
            return self.find_bisection_root(median, b, val_b=val_b)

    # @profile

    def update_centroids(self):
        self.centroids = np.mean(self.vertices[self.faces[:], :], axis=1)

    def buildAFromNormal(self, point):
        assert not np.any(np.isnan(point)), "there should not be any NaN values"
        point = point.reshape(1, 4)
        n = (self.gradient(point))[0][:3].reshape(3, 1)
        n /= np.linalg.norm(n)
        return (n * np.transpose(n))

    def buildNormalAtPoint(self, point):
        assert not np.any(np.isnan(point)), "there should not be any NaN values"
        point = point.reshape(1, 4)
        n = (self.gradient(point))[0][:3].reshape(3, 1)
        n /= np.linalg.norm(n)
        return n

    @profile
    def minimize_tangent_planes(self, ind):
        print "[Start] minimize_tangent_planes ind: %d" % ind + ' ...'
        neighbors = self.vertex_neighbours_list[ind]
        assert not np.any(np.isnan(neighbors)), "neighbors should not contain NaN values"
        alpha_current = np.zeros(shape=(3, 3))
        p_matrix = np.array([self.optimized_dual_mesh[neighbor] for neighbor in neighbors])
        assert not np.any(np.isnan(p_matrix)), "The matrix of points should not contain NaN values"
        A_total = np.zeros(shape=(3, 3))
        b_total = np.zeros(shape=(3, 1))
        for i, neighbor in enumerate(neighbors):
            p_current = p_matrix[i]
            alpha_current = self.buildAFromNormal(p_current)

            assert alpha_current.shape == (3, 3), "the alpha matrix is 3X3"
            A_total += alpha_current
            b_current = -2 * (np.dot(p_current[0:3], alpha_current)).reshape(3, 1)
            b_total += b_current
            res = np.linalg.lstsq(2 * A_total, -b_total)
        print res
        print "[End] minimize_tangent_planes ind: %d" % ind + ' .'
        return res[0].reshape(3, 1)

    def plot_tangent_planes(self, ind):

        neighbors = self.vertex_neighbours_list[ind]

        print "These are the neighbors of vertex %d" % ind
        print neighbors
        p_matrix = np.array([self.optimized_dual_mesh[neighbor] for neighbor in neighbors])
        assert not np.any(np.isnan(p_matrix)), "The matrix of points should not contain NaN values"
        planes = []
        for i, neighbor in enumerate(neighbors):
            p_current = p_matrix[i]
            n_current = self.buildNormalAtPoint(p_current)
            dim = 1
            planes.append((p_current[:-1], n_current))
        my_colors = [(1, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1), (1, 1, 1)]
        x, y, z = np.mgrid[-dim:dim:10j, -dim:dim:10j, -dim:dim:10j]
        for i, plane in enumerate(planes):
            # set_trace()
            point, normal = plane
            point = point.flatten()
            normal = normal.flatten()

            scalars = normal[0] * x + normal[1] * y + normal[2] * z - np.dot(point, normal)
            mlab.contour3d(x, y, z, scalars, contours=[0], opacity=0.2, color=my_colors[i])

    @profile
    def optimize_centroid(self, centroid, ind):
        print "[Start] optimize_centroid ..."
        #  set_trace()
        print centroid
        centroid = centroid.reshape(1, 4)

        lambda_counter = 0
        projection = centroid
        R = np.ones((1, 4))
        fnc_val = self.function(projection)
        if abs(fnc_val[0]) < self.TOLERANCE:    # if condition is already met return
            return projection
        else:
            Q = projection
        while True:
            lambda_counter += 1
            fnc_val_Q = self.function(Q)[0]
            gradient_Q = self.gradient(Q)[0][:-1]
            R[0][:-1] = Q[0][:-1] + self.lamda * (- fnc_val_Q * gradient_Q) / (abs(fnc_val_Q) * np.linalg.norm(gradient_Q))
            condition = fnc_val_Q * self.function(R)[0]
            # print condition
            if condition < 0:
                try:
                    projection = self.find_bisection_root(Q, R)
                    break
                except:
                    print "find_bisection_root failed"
            elif condition < self.TOLERANCE:
                return Q
            else:
                Q = R.copy()
        # log lambda values for every simulation so we tune it later on .
        # with open('lambda_logs','a') as logfile:
            # print('lambda_val, ' + repr(self.lambda_val) + ', lambda_counter, '+ repr(lambda_counter) + '\n')
        print "[End] optimize_centroid ."

        # set_trace()
        # assert not np.any(np.isnan(projection)), "Projections cant be NaN"
        return projection

    # @profile
    def run(self, calc_proj=True, calc_opt=True, update_centroids=False):
        print "[Start] Run simulation for %s" % self.name
        if calc_proj:
            if update_centroids:
                # set_trace()
                self.update_centroids()

            lngth = len(self.centroids)
            self.optimized_dual_mesh = np.ones(shape=(lngth, 4))
            for i, centroid in enumerate(self.centroids):
                # set_trace()
                self.optimized_dual_mesh[i] = self.optimize_centroid(centroid, i)

            assert not np.any(np.isnan(np.array(self.optimized_dual_mesh)))

        if calc_opt:
            for current in range(len(self.vertices)):
                opt_vertex = self.minimize_tangent_planes(current)
                self.optimized_vertices.append(opt_vertex)
                # print "Stats : "
                # print "Square error of primal mesh: %f " % sum(self)
        print "[End] Run simulation for %s" % self.name
# Util functions


def plot3d_planes(self, ind):
    """

    INPUT :  tuple containing a point and normal vector  in the form ((0,0,1),(0,1,1))
    the function will plot the 3d plane computed by the equation n*(x - P) = 0 , where x = [x y z] is provided .

    """

    # assert type(args) == list
    # print type(planes)
    axis_dim = 2
    # Axes
    x = np.linspace(-axis_dim, axis_dim, 2).reshape(2, 1)
    y = np.zeros((2, 1))
    z = np.zeros((2, 1))
    mlab.plot3d(x, y, z, line_width=1, name="x-axis", opacity=0.5)
    mlab.plot3d(y, x, z, line_width=1, name="y-axis", opacity=0.5)
    mlab.plot3d(z, y, x, line_width=1, name="z-axis", opacity=0.5)
    # Axes labels
    # mlab.text3d(0,0,0,"O",scale=.5)
    mlab.text3d(axis_dim, 0, 0, 'X', scale=.1)
    mlab.text3d(0, axis_dim, 0, 'Y', scale=.1)
    mlab.text3d(0, 0, axis_dim, 'Z', scale=.1)
    mlab.text3d(-axis_dim, 0, 0, '-X', scale=.1)
    mlab.text3d(0, -axis_dim, 0, '-Y', scale=.1)
    mlab.text3d(0, 0, -axis_dim, '-Z', scale=.1)

    # Origin
    mlab.points3d(0, 0, 0, color=(0, 0, 0), scale_factor=0.1)
    dim = 1
    my_colors = [(1, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1), (1, 1, 1)]
    x, y, z = np.mgrid[-dim:dim:10j, -dim:dim:10j, -dim:dim:10j]
    for i, plane in enumerate(planes):
        point, normal = np.array(plane)
        point = point.flatten()
        normal = normal.flatten()

        scalars = normal[0] * x + normal[1] * y + normal[2] * z - np.dot(point, normal)
        mlab.contour3d(x, y, z, scalars, contours=[0], opacity=0.2, color=my_colors[i])


def build_neighbor_list_c2c(centroids, faces, neighbors_list):
    """ This function returns a dictionary with centroid indexes for keys.
        For each key the value is a list containing the centroids that are
        neighbors of the key.
    """

    res = []
    for i, centroid in enumerate(centroids):
        v1, v2, v3 = faces[i]
        n1, n2, n3 = neighbors_list[v1], neighbors_list[v2], neighbors_list[v3]
        partial1 = [item for item in n1 if item in n2]
        partial2 = [item for item in n2 if item in n3]
        partial3 = [item for item in n1 if item in n3]
        n = partial1 + partial2 + partial3
        res.append([item for item in n if n.count(item) == 1])
    return np.array(res)


def calc_avg_triangle_length(triangles):
    """ Calculates average triangle length for the mesh triangles
        triangles :verts[faces] (marching cubes output)) """

    def avg_edge_len(triangle):
        return (np.linalg.norm(triangle[0] - triangle[1]) + \
                np.linalg.norm(triangle[0] - triangle[2]) + \
                np.linalg.norm(triangle[1] - triangle[2])) / 3.
    avg_len = map(avg_edge_len, triangles)
    res = reduce(lambda x, y: x + y, avg_len, 0) / (1. * len(avg_len))
    return res
