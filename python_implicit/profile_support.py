import sys
is_python2 = sys.version_info[0] < 3

if is_python2:
    try:
        __builtin__.profile
    except AttributeError:
        def profile(func):
            return func
        __builtin__.profile = profile

else:
    # python3

    try:
        import builtins
        builtins.profile

    except AttributeError:
        def profile(func):
            return func
        builtins.profile = profile
        #globals()['profile'] = profile # sometimes does not work
        #print(profile)  # does not work
    #print(profile)
