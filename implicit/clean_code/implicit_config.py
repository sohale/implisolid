VERBOSE = False

TOLERANCE = 0.0000001


INTEGRITY_TOLERANCES_NORM = 0.00000000001


# 3d printing using FDA
# config_FDA={}
# config_FDA["max_iter"] = 20
# config_FDA["numerical_min_length"] = 0.1

class config_FDA(object):

    __slots__ = ('max_iter', 'numerical_min_length')   # high performance

    def __init__(self):
        self.max_iter = 20
        self.numerical_min_length = 0.1


config = config_FDA()
