# Copyright 2016 Sean Patrick Santos
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# See: git@github.com:quantheory/finitediff.git

import numpy as np

__all__ = ['weights']


def weights(k, x0, xs):
    """Calculate weights for the finite difference approximation.
    Arguments:
    k - The k-th derivative will be approximated.
    x0 - The point at which to approximate the derivative.
    xs - The grid points at which the function's value is known.
    This uses the algorithm described in:
    B. Fornberg, "Calculation of weights in finite difference formulas",
    SIAM Review 40 (1998), pp. 685-691.
    """
    # Size of the system.
    n = xs.size
    assert k < n, "need more grid points to calculate this derivative"
    # Measure points relative to x0.
    xs = xs - x0
    # Weight matrix to calculate successive finite difference
    # approximations.
    w = np.zeros((k+1, n))
    # Before starting, we want to pre-compute certain reusable
    # quantities.
    product_ratio = np.ones(n)
    for j in range(n):
        for i in range(j):
            product_ratio[j] *= xs[j] - xs[i]
    product_ratio[1:] = product_ratio[:n-1] / product_ratio[1:]
    # Each iteration of this loop will produce weights for the
    # derivatives that use one point more than the previous iteration.
    for j in range(n - k):
        # 0-th derivative approximation.
        if j == 0:
            # The approximation to the 0th derivative given only one
            # function value is trivially to just use that value.
            w[0, 0] = 1
        else:
            w[0, j] = - xs[j-1] * w[0, j-1] * product_ratio[j]
            for i in range(1, j+1):
                w[0, j-i] = xs[j] * w[0, j-i] / (xs[j] - xs[j-i])
        for m in range(1, k+1):
            # Generate weights for each derivative using the
            # previous one.
            # m is the derivative we are currently working on,
            # and l is the number of points used in this round,
            # minus one.
            l = j+m
            w[m, l] = (m*w[m-1, l-1] - xs[l-1]*w[m, l-1]) * product_ratio[l]
            for i in range(1, l+1):
                w[m, l-i] = (xs[l] * w[m, l-i] - m*w[m-1, l-i]) / (xs[l] - xs[l-i])
    return w[k, :]
