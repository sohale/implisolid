from basic_types import check_vector4_vectorized
from implicit_vectorized import ImplicitFunctionVectorized
import numpy as np

class PointCollector(ImplicitFunctionVectorized):
    """ Used for debugging and development purposes"""

    def __init__(self, iobject):
        self.iobject = iobject
        self.eval_collector = []

    def integrity_invariant(self):
        return self.iobject.integrity_invariant()

    def implicitFunction(self, x):
        self.eval_collector.append(x.copy())
        return self.iobject.implicitFunction(x)

    def implicitGradient(self, x):
        #self.grad_collector.append(x.copy())
        return self.iobject.implicitGradient(x)

    def hessian(self, x):
        #self.hessian_collector.append(x.copy())
        return self.iobject.hessian(x)

    def reset(self):
        self.eval_collector = []

    def get_as_array(self):
        a = np.concatenate(tuple(self.eval_collector), axis=0)
        #a = a[::15, :]
        a = a[:5000, :]
        #print a.shape
        if a.shape[0] > 6000:
            raise Exception("Too many points: %d"%a.shape[0])
        return a
"""
        display_simple_using_mayavi_2( [(verts, facets), ],
           pointcloud_list=[point_collector.get_as_array()], pointsizes=[0.05],
           mayavi_wireframe=[False,], opacity=[0.2,],
           gradients_at=None, separate_panels=False, gradients_from_iobj=None,
           minmax=(RANGE_MIN,RANGE_MAX),
           add_noise=[0], noise_added_before_broadcast=True  )
        exit()

"""