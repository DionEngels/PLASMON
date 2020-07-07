"""Generic interface for least-square minimization."""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.linalg import norm

from scipy.optimize import _minpack, OptimizeResult
from scipy.optimize._numdiff import approx_derivative
from scipy._lib.six import string_types

from scipy.optimize._lsq.common import EPS


TERMINATION_MESSAGES = {
    -1: "Improper input parameters status returned from `leastsq`",
    0: "The maximum number of function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`ftol` termination condition is satisfied.",
    3: "`xtol` termination condition is satisfied.",
    4: "Both `ftol` and `xtol` termination conditions are satisfied."
}


FROM_MINPACK_TO_COMMON = {
    0: -1,  # Improper input parameters from MINPACK.
    1: 2,
    2: 3,
    3: 4,
    4: 1,
    5: 0
    # There are 6, 7, 8 for too small tolerance parameters,
    # but we guard against it by checking ftol, xtol, gtol beforehand.
}


def call_minpack(fun, x0, jac, ftol, xtol, gtol, max_nfev, diag):

    n = x0.size
    epsfcn = EPS

    # Compute MINPACK's `diag`, which is inverse of our `x_scale` and
    # ``x_scale='jac'`` corresponds to ``diag=None``.

    full_output = True
    factor = 100.0

    if max_nfev is None:
        # n squared to account for Jacobian evaluations.
        max_nfev = 100 * n * (n + 1)
    x, info, status = _minpack._lmdif(
        fun, x0, (), full_output, ftol, xtol, gtol,
        max_nfev, epsfcn, factor, diag)

    f = info['fvec']

    J = np.atleast_2d(approx_derivative(fun, x))

    cost = 0.5 * np.dot(f, f)
    g = J.T.dot(f)
    g_norm = norm(g, ord=np.inf)

    nfev = info['nfev']
    njev = info.get('njev', None)

    status = FROM_MINPACK_TO_COMMON[status]
    active_mask = np.zeros_like(x0, dtype=int)

    return OptimizeResult(
        x=x, cost=cost, fun=f, jac=J, grad=g, optimality=g_norm,
        active_mask=active_mask, nfev=nfev, njev=njev, status=status)


def least_squares(
        fun, x0, ftol=1e-8, xtol=1e-8, gtol=1e-8, f_scale=1.0, 
        max_nfev=None, args=(), kwargs={}, cache=None):

    if max_nfev is not None and max_nfev <= 0:
        raise ValueError("`max_nfev` must be None or positive integer.")

    if np.iscomplexobj(x0):
        raise ValueError("`x0` must be real.")

    x0 = np.atleast_1d(x0).astype(float)

    if x0.ndim > 1:
        raise ValueError("`x0` must have at most 1 dimension.")

    x_scale = np.ones(5)

    def fun_wrapped(x):
        nonlocal cache
        
        key = tuple(x) 
        if key in cache:
            return cache[key]
        else:
            result = fun(x, *args, **kwargs)
            cache[key] = result
            return result

    f0 = fun_wrapped(x0)

    if f0.ndim != 1:
        raise ValueError("`fun` must return at most 1-d array_like.")

    if not np.all(np.isfinite(f0)):
        raise ValueError("Residuals are not finite in the initial point.")

    result = call_minpack(fun_wrapped, x0, None, ftol, xtol, gtol,
                          max_nfev, x_scale)


    result.message = TERMINATION_MESSAGES[result.status]
    result.success = result.status > 0

    return result, cache
