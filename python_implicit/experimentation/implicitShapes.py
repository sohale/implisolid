import numpy as np
import numexpr as ne
from circle_square_symbolic_demo import contour_fnc


def torus_fnc(x, y, z, r=1):
    r = 4
    rx = 0.1
    ry = 0.1
    rz = 0.1
    return (1 - (r - ((x / rx) ** 20 + (y / ry) ** 20)**0.05) ** 20 - (z / rz) ** 20)


def torus_circle_ne(x, y, z, r=1):
    r = 2
    rx = 0.4
    ry = 0.4
    rz = 1
    sphere = ne.evaluate('(1 - x ** 2 - y ** 2 - z ** 2)')
    torus = ne.evaluate('(1 - (r - ((x / rx) ** 20 + (y / ry) ** 20)**0.05 ) ** 20 - (z / rz) ** 20 )')
    expr = 'where(sphere<torus, sphere, torus)'  # equivalent to min
    return ne.evaluate(expr)


def torus_two_circles_ne(x, y, z, r=1):
    r = 2
    rx = 0.4
    ry = 0.4
    rz = 1
    sphere1 = ne.evaluate('1 - x ** 2 - y ** 2 - z ** 2')
    torus = ne.evaluate('(1 - (r - ((x / rx) ** 2 + (y / ry) ** 2)**0.05 ) ** 2 - (z / rz) ** 2 )')
    sphere2 = ne.evaluate('0.8 - x ** 2 - y ** 2 - z ** 2')
    expr = ne.evaluate('where(sphere1<torus, sphere1, torus)')  # equivalent to min
    expr2 = ne.evaluate('where(sphere2 < expr, sphere2, expr)')
    return expr2


def torus_ne(x, y, z, r=1):
    r = 2
    rx = 0.4
    ry = 0.4
    rz = 0.4

    torus = ne.evaluate('(1 - (r - ((x / rx) ** 20 + (y / ry) ** 20)**0.05 ) ** 20 - (z / rz) ** 20 )')
    expr = 'where(sphere<torus, sphere, torus)'  # equivalent to min
    return torus
