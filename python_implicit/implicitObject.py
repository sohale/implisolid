import sympy as sp
import numpy as np
from sympy import symbols, Function, sympify, diff
import numexpr as ne
import types
from vtk_mc import vtk_mc
from python_clean_code.mc_utils import _prepare_grid, make_grid
from mayavi import mlab

"""
    File: ImplicitObject.py
    Description: Defines primitives as implicit functions
    and csg operations (Rvachev functions)

    Implementation Details:
        The classes defined here, depend on  SymPy and Numexpr for symbolic and fast calculations.


    Usage:
            To use this file, first you will have to import it.

            Example:
                    from implicitObject import Heart

                    # Init new obj
                    new_obj = Heart()

                    # Marching Cubes it
                    new_obj.build_marching_cubes()

                    # Show it!

                    new_obj.show(representation="surface") or new_obj.show() which uses wireframe representation as a default.

"""


class ImplicitFunction(Function):
    # Does it need to inherit Function?
    def __init__(self, fnc_expression):
        """

        ImplicitFunction class
        ARGS:

        1. f n c _ e x p r e s s i o n:  is an expression (string or python expression )
        that represents an implicit function

        example: my_implicit_fnc = ImplicitFunction('1-x**2 - y**2 -z**2')

        """

        # init symbols for sympy
        x, y, z = symbols('x y z')

        # str_fnc can be a string or a python expression
        if type(fnc_expression) == str:
            self.formula = sympify(fnc_expression)
        else:
            self.formula = fnc_expression

        # create a valid numpy function for the implicit function
        # import ipdb;ipdb.set_trace()?
        self.fnc_numpy = sp.lambdify((x, y, z), self.formula, modules="numpy")

        # self.fnc_numexpr = sp.lambdify(str_fnc,(x,y,z), modules="numexpr")

        self.gradx, self.grady, self.gradz = [str(sp.diff(self.formula, variable, chop=True)) for variable in [x, y, z]]

    def __call__(self, x, y, z):

        """
        ARGS:

        1. x,y,z : NX1 ndarrays , they can have shape (N,1) or (N,)


        example: points = np.ones((10,3))
        ImplicitFunction(points) will call ImplicitFunction.fnc_numpy on points

        """

        return self.fnc_numpy(x, y, z)

    def gradient(self, x, y, z):

        return np.column_stack((ne.evaluate(self.gradx), ne.evaluate(self.grady), ne.evaluate(self.gradz)))


class Object3D(object):

    """ Base class to inherit from for all objects"""

    def __init__(self, bound_dim=2):
        """
            bound_dim:  an indicator of the size of the object to be used for the marching cubes bounds

        """
        self.triangles = None
        self.vertices = None
        self.faces = None
        self.bound_dim = bound_dim

    def show(self, representation="wireframe"):
        """
        This function shows the triangular_mesh of 3D objects, and the
        helping axes .

        OPTIONS:
            representation="wireframe" || "surface" , default is wireframe
        """
        mlab.figure()
        (RANGE_MIN, RANGE_MAX) = (-5, 5)
        x = np.linspace(RANGE_MIN, RANGE_MAX, 2).reshape(2, 1)
        y = np.zeros((2, 1))
        z = np.zeros((2, 1))

        mlab.plot3d(x, y, z, line_width=3, name="x-axis")
        mlab.plot3d(y, x, z, line_width=3, name="y-axis")
        mlab.plot3d(z, y, x, line_width=3, name="z-axis")

        mlab.text3d(RANGE_MAX, 0, 0, "x", scale=0.3)
        mlab.text3d(0, RANGE_MAX, 0, "y", scale=0.3)
        mlab.text3d(0, 0, RANGE_MAX, "z", scale=0.3)
        mlab.text3d(RANGE_MIN, 0, 0, "-x", scale=0.3)
        mlab.text3d(0, RANGE_MIN, 0, "-y", scale=0.3)
        mlab.text3d(0, 0, RANGE_MIN, "-z", scale=0.3)

        mlab.triangular_mesh([vert[0] for vert in self.vertices],
                             [vert[1] for vert in self.vertices],
                             [vert[2] for vert in self.vertices],
                             self.faces, representation=representation)
        mlab.show()
        #opacity = 0.2 #0.1

    def contour_(self):
        from mayavi import mlab

        """
        Uses mayavi to visualize the formula of object in 3D space
        ARGS
        1.dim, slices: arguments for np.mgrid

        NOTES:

        1. This function assumes that mayavi.mlab is imported:
        if not:
            from mayavi import mlab
        if mayavi is not installed:
            install with pip
        2. This function creates a contour of the mathematical expression of the
        object and not the actual vertices and faces.
        """

        # assert isinstance(self.function.formula, types.FunctionType), "Argument fnc should be a callable"
        dim, slices = 2, 100j
        x, y, z = np.mgrid[-dim:dim:slices, -dim:dim:slices, -dim:dim:slices]
        mlab.figure()
        mlab.contour3d(x, y, z, self.function.fnc_numpy, contours=[0])
        # mlab.show()

    def build_marching_cubes(self):
        minn, maxx, step = -3, 3, 0.05
        rng = np.arange(-3, 3, 0.05)
        x, y, z = np.mgrid[minn:maxx:step, minn:maxx:step, minn:maxx:step]
        xyza = _prepare_grid(rng)
        vgrid_v = self.function.fnc_numpy(xyza[:, 0], xyza[:, 1], xyza[:, 2])
        vgrid = np.reshape(vgrid_v, (len(rng), len(rng), len(rng)), order='C')
        self.vertices, self.faces = vtk_mc(vgrid, (minn, maxx, step))
        from mayavi import mlab
        # mlab.figure()
        # mlab.triangular_mesh([vert[0] for vert in self.vertices],
        #                      [vert[1] for vert in self.vertices],
        #                      [vert[2] for vert in self.vertices], self.faces, representation="wireframe")
        # mlab.show()


class UnitSphere(Object3D):

    """ A unit sphere  """

    def __init__(self, radius=1, center=(0, 0, 0)):
        super(UnitSphere, self).__init__(bound_dim=2.2 * radius)    # set the bound_dim a bit more than the diameter
        # x,y,z = sp.symbols('x y z')
        cx, cy, cz = center
        self.function = ImplicitFunction("{radius}**2 - ((x - {cx}) ** 2) - ((y - {cy}) ** 2) - ((z - {cz}) ** 2)".format(radius=radius, cx=center[0], cy=center[1], cz=center[2]))
        self.function_str = str(self.function)
        # self.gradient = (sp.diff(self.function, x),sp.diff(self.function, y),sp.diff(self.function, y))


class UnitCube(Object3D):

    """ A unit cube """
    #   Implementation Detail: This is a supercube, constructed by using high powers ^20 - ^100

    def __init__(self, sideLen=1):
        super(UnitCube, self).__init__(bound_dim=2.2 * sideLen)     # set the bound_dim a bit more than twice the side length
        x, y, z = symbols('x y z')
        self.function = ImplicitFunction("{sideLen} - (x ** 20 + y ** 20 + z ** 20)** (1/20.)".format(sideLen=sideLen))


class Cone(Object3D):

    """ Cone implicit function """

    def __init__(self, c=1, z0=0):
        super(Cone, self).__init__(bound_dim=4)     # set the bound_dim a bit more than twice the side length
        x, y, z = symbols('x y z')
        self.function = ImplicitFunction("((x**2 + y**2)/({c}**2) - (z - {z0})**4) ".format(c=c, z0=z0))


class Torus(Object3D):

    """ A torus primitive """

    def __init__(self):
        super(Torus, self).__init__(bound_dim=2.2)
        x, y, z = symbols('x y z')
        r, rx, ry, rz = 4, .1, .1, .1
        self.function = ImplicitFunction("1 - ({r} - ((x / {rx}) ** 2 + (y / {ry}) ** 2)**0.5 ) ** 2 - (z / {rz}) ** 2".format(r=r, rx=rx, ry=ry, rz=rz))


class Heart(Object3D):
    """ A heart shaped primitive """
    def __init__(self):
        super(Heart, self).__init__(bound_dim=2.2)
        x, y, z = symbols('x y z')
        self.function = ImplicitFunction("(2*x**2 + y**2 + z**2 -4)**3 - (0.1*x**2 + y**2)*z**3")


class Citrus(Object3D):
    """ """
    def __init__(self):
        super(Citrus, self).__init__(bound_dim=2.2)
        x, y, z = symbols('x y z')
        self.function = ImplicitFunction("x**2 + z**2 - (4*y**3) * (1 - 0.5*y)**3")


class CylinderOpen(Object3D):
    def __init__(self, radius=1):
        super(CylinderOpen, self).__init__(bound_dim=2.2)
        x, y, z = symbols('x y z')
        self.function = ImplicitFunction("x ** 2 + y ** 2 - {radius}**2".format(radius=radius))


class Cylinder(Object3D):
    def __init__(self, cx=1, cy=0.5, cz=0.5, c=1):
        """
        c: controls radius
        cx: controls length
        cy, cz : control y and z axis radius
        """
        super(Cylinder, self).__init__(bound_dim=10)
        x, y, z = symbols('x y z')
        self.function = ImplicitFunction("(x/{cx})**10 + ({cy}*y**2 + {cz}*z**2) - {c}".format(cx=cx, cy=cy, cz=cz, c=c))


class Intersection(Object3D):

    """ Rvachev intersection """
    def __init__(self, obj1, obj2):
        super(Intersection, self).__init__(bound_dim=4)   # bound_dim needs to be extracted from somewhere that makes sense, for now set a default=10
        x, y, z = symbols('x y z')
        self.function = ImplicitFunction("{expr1} + {expr2} - ({expr1}**(2) + {expr2}**2)**(1/2.)".format(expr1=obj1.function.formula, expr2=obj2.function.formula))


class Union(Object3D):
    """ Rvachev Union """
    def __init__(self, obj1, obj2):
        super(Union, self).__init__(bound_dim=4)
        x, y, z = symbols('x y z')
        self.function = ImplicitFunction("{expr1} + {expr2} + ({expr1}**(2) + {expr2}**2)**(1/2.)".format(expr1=obj1.function.formula, expr2=obj2.function.formula))
