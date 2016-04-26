import __builtin__

try:
    __builtin__.profile
except AttributeError:
    def profile(func):
        return func
    __builtin__.profile = profile
