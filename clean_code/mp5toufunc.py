import vector3
from basic_functions import make_vector4
import  numpy as np

def get_root_node(tree):
    return tree["root"]



def get_fonuky(node):
    type = node["type"]


    switch = {
        "Difference": CrispSubtract,
        "cylinder": SimpleCylinder,
        "root" : root,

    }
    ufunc = switch[type](node)
    return ufunc


def CrispSubtract(node):

    sons = node["children"]
    son1 = get_fonuky(sons[0])
    son2 = get_fonuky(sons[1])
    matrix = make_matrix4(node["matrix"])
    return vector3.Transformed(vector3.CrispSubtract(son1, son2), matrix)

def SimpleCylinder(node):
    matrix = make_matrix4(node["matrix"])
    A = make_vector4(-0.5, 0, 0)
    w = make_vector4(0, 0, 1)
    u = make_vector4(1, 0, 0)
    return vector3.Transformed(vector3.SimpleCylinder(A,w,u,0.5, 0.5, 1), matrix)

def root(node):
    sons = node["children"]
    return get_fonuky(sons[0])

def make_matrix4(array):
    arr = np.array(array).astype(float)
    return np.reshape(arr, (4,4))


#
# class node():
#
#     def __init__(self, type, matrix):
#         self.type = type
#         self.sons =[]
#         self.matrix = matrix
#
#
#     def is_leaf(self):
#         return (self.type == 'Difference' or self.type == "mescouilles")
#
#     def add_son(self, node):
#         self.sons.append(node)


