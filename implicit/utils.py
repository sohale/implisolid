
def optimised_used():
    global _optimised_used
    _optimised_used = True
    def side_effect():
        global _optimised_used
        _optimised_used = False
        return True
    assert side_effect()
    #print "optimisation", _optimised_used
    return  _optimised_used


import time

class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

import sys
def flush():
    sys.stdout.flush()
