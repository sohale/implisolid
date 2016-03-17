import numpy as np

"""
Optimizing dual mesh
Steps:
Consider the dual mesh formed by triangle centroids.We are
looking for an optimized dual mesh whose vertices are obtained
by projecting the centroids C onto the implicit surface.The projection
P of a triangle centroid C onto the surface f(x,y,z) = 0 is estimated as follows:

    1. Set P = C
    2. If |f(P)| < tolerance: terminate
        else:
        Q = P
    3. Set R = Q + lambda*d , where d = -f(Q) * gradientf(Q) (vector) and
    lambda is a small constant
    4. if f(Q) * f(R) < 0 : find_root with condition f(P) < tolerance inside [Q,R]
    and terminate
    else :  Set Q = R  and go to step 3 .

"""

class MeshOptimizer(object):

    def __init__(self, centroids=None, triangles=None):
        self.lambda_val = 0.01
        self.TOLERANCE = 0.0001
        self.centroids = centroids

    def register_function(self, fnc):
        self.function = fnc

    def register_gradient(self, grd):
        self.gradient = grd

    def load_example(self, filename):
        # need to load an object in pieces since we can pass declaration of the implicit function
        from example_objects import make_example_vectorized
        import pickle
        # need to get the implicit and gradient from here
        iobj = make_example_vectoried(filename)
        with open(filename) as fid:
            obj = pickle.load(fid)
        self.centroids = np.c_[ obj['centroids'], np.ones(len(obj['centroids']))]
        self.faces = obj['faces']
        self.vertices = obj['verts']
        self.register_function(iobj.implicitFunction)
        self.register_gradient(iobj.implicitGradient)

    def optimize_centroid(self,centroid):
        # initialize P = C
        centroid = centroid.reshape(1,4)
        lambda_counter = 0
        projection = centroid
        R = np.ones((1,4))

        if abs(self.function(projection)[0]) < self.TOLERANCE: # if condition is already met return
            return projection
        else:
            Q = projection
        while True:

            lambda_counter += 1

            R[0][:-1] = Q[0][:-1] + self.lambda_val * (- 2.0*self.function(Q)[0] * self.gradient(Q)[0][:-1])/(np.linalg.norm(self.gradient(Q)[0][:-1]))
            print self.function(R)[0]
            condition = self.function(Q)[0] * self.function(R)[0]
            if  condition < 0 :
                try:
                    projection = self.find_bisection_root(Q,R)
                    break
                except:
                    print "find_bisection_root failed"
            elif condition < self.TOLERANCE:
                return Q
            else:
                Q = R
        # log lambda values for every simulation so we tune it later on .
        with open('lambda_logs','a') as logfile:
            print('lambda_val, ' + repr(self.lambda_val) + ', lambda_counter, '+ repr(lambda_counter) + '\n')

        return projection



    def find_bisection_root(self, a, b):

        """
        Bisection method in the line that connects the points a, b .
        Each time we update one of lthe edge points to the median in the line
        """
        a  = np.array(a)
        b = np.array(b)
        median = ( a + b ) / 2.0

        val_a = self.function(a)[0]
        val_b = self.function(b)[0]
        val_md = self.function(median)[0]

        if abs(val_a) < self.TOLERANCE:
            return a
        if abs(val_b) < self.TOLERANCE:
            return b
        if abs(val_md) < self.TOLERANCE:
            return median

        if val_md * val_a < 0:
            return self.find_bisection_root(a,median)
        elif val_md * val_b < 0:
            return self.find_bisection_root(median,b)

    def run(self):
        self.dual_mesh = []
        for centroid in self.centroids:
            self.dual_mesh.append(optimize_centroid(centroid))
