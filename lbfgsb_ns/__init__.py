### Note: most of this copy-paste from SciPy's L-BFGS-B wrapper:
### https://github.com/scipy/scipy/blob/master/scipy/optimize/lbfgsb.py

import numpy as np
from lbfgsb_ns import _lbfgsb_ns
from scipy.optimize import OptimizeResult
from scipy.optimize.lbfgsb import LbfgsInvHessProduct
from numpy import array, asarray, float64, int32, zeros
from scipy.optimize._optimize import (MemoizeJac, OptimizeResult,
					   _check_unknown_options, _prepare_scalar_function)
from scipy.optimize._constraints import old_bound_to_new
from scipy.optimize import _lbfgsb

def minimize_lbfgsb_ns(fun, x0, args=(), jac=None, bounds=None,
	disp=None, maxcor=10, ftol=2.2204460492503131e-09,
	gtol=1e-5, eps=1e-8, maxfun=15000, maxiter=15000,
	iprint=-1, callback=None, maxls=20, finite_diff_rel_step=None, taux=1e-3, taud=1e-6, **unknown_options):
	"""
	Minimize a scalar function of one or more variables using the L-BFGS-B-NS
	algorithm.

	Parameters
	----------
	disp : None or int
		If `disp is None` (the default), then the supplied version of `iprint`
		is used. If `disp is not None`, then it overrides the supplied version
		of `iprint` with the behaviour you outlined.
	maxcor : int
		The maximum number of variable metric corrections used to
		define the limited memory matrix. (The limited memory BFGS
		method does not store the full hessian but uses this many terms
		in an approximation to it.)
	ftol : float
		The iteration stops when ``(f^k -
		f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= ftol``.
	gtol : float
		The iteration will stop when ``max{|proj g_i | i = 1, ..., n}
		<= gtol`` where ``pg_i`` is the i-th component of the
		projected gradient.
	eps : float
		Step size used for numerical approximation of the jacobian.
	maxfun : int
		Maximum number of function evaluations.
	maxiter : int
		Maximum number of iterations.
	maxls : int, optional
		Maximum number of line search steps (per iteration). Default is 20.
	taux : float
		Parameter for line search
	taud: float
		Parameter for line search

	Returns
	-------
	res : OptimizeResult
		dict-like object with the optimal values in entry 'x' and other associated info.

	References
	----------
	Henao, W., 2014. An L-BFGS-B-NS optimizer for non-smooth functions. Master's thesis, Courant Institute of Mathematical Science, New York University.
	"""
	_check_unknown_options(unknown_options)
	m = maxcor
	epsilon = eps
	pgtol = gtol
	factr = ftol / np.finfo(float).eps

	x0 = asarray(x0).ravel()
	n, = x0.shape

	if bounds is None:
		bounds = [(None, None)] * n
	if len(bounds) != n:
		raise ValueError('length of x0 != length of bounds')
	# unbounded variables must use None, not +-inf, for optimizer to work properly
	bounds = [(None if l == -np.inf else l, None if u == np.inf else u) for l, u in bounds]
	# LBFGSB is sent 'old-style' bounds, 'new-style' bounds are required by
	# approx_derivative and ScalarFunction
	new_bounds = old_bound_to_new(bounds)

	# check bounds
	if (new_bounds[0] > new_bounds[1]).any():
		raise ValueError("LBFGSB-NS - one of the lower bounds is greater than an upper bound.")

	# initial vector must lie within the bounds. Otherwise ScalarFunction and
	# approx_derivative will cause problems
	x0 = np.clip(x0, new_bounds[0], new_bounds[1])

	if disp is not None:
		if disp == 0:
			iprint = -1
		else:
			iprint = disp

	if jac == True:
		fun = MemoizeJac(fun)
		jac = fun.derivative
	sf = _prepare_scalar_function(fun, x0, jac=jac, args=args, epsilon=eps,
								  bounds=new_bounds,
								  finite_diff_rel_step=finite_diff_rel_step)

	func_and_grad = sf.fun_and_grad

	fortran_int = _lbfgsb.types.intvar.dtype

	nbd = zeros(n, fortran_int)
	low_bnd = zeros(n, float64)
	upper_bnd = zeros(n, float64)
	bounds_map = {(None, None): 0,
				  (1, None): 1,
				  (1, 1): 2,
				  (None, 1): 3}
	for i in range(0, n):
		l, u = bounds[i]
		if l is not None:
			low_bnd[i] = l
			l = 1
		if u is not None:
			upper_bnd[i] = u
			u = 1
		nbd[i] = bounds_map[l, u]

	if not maxls > 0:
		raise ValueError('maxls must be positive.')

	x = array(x0, float64)
	f = array(0.0, float64)
	g = zeros((n,), float64)
	jmax= maxls
	wa = zeros(2*m*n + 5*n + 11*m*m + 8*m+ 2*jmax*n, float64)
	iwa = zeros(3*n, fortran_int)
	task = zeros(1, 'S120')
	csave = zeros(1, 'S120')
	lsave = zeros(4, fortran_int)
	isave = zeros(47, fortran_int)
	dsave = zeros(30, float64)

	task[:] = 'START'

	n_iterations = 0

	nfg = array(0, fortran_int)
	nbisect = array(0, fortran_int)

	taux = array(taux, float64)
	taud = array(taud, float64)
	jmax = array(maxls, fortran_int)
	m = array(m, fortran_int)
	n = array(n, fortran_int)

	while 1:
		_lbfgsb_ns.setulb(
				n = n, 
				m = m, 
				x = x, 
				l = low_bnd, 
				u = upper_bnd, 
				nbd = nbd, 
				f = f, 
				g = g, 
				factr = factr, 
				pgtol = pgtol,  
				wa = wa, 
				iwa = iwa, 
				task = task, 
				iprint = iprint, 
				csave = csave, 
				lsave = lsave,
				isave = isave, 
				dsave = dsave, 
				taux = taux, 
				nfg = nfg,
				jmax = jmax, 
				taud = taud, 
				nbisect = nbisect
			)

		task_str = task.tobytes()
		if task_str.startswith(b'FG'):
			# The minimization routine wants f and g at the current x.
			# Note that interruptions due to maxfun are postponed
			# until the completion of the current minimization iteration.
			# Overwrite f and g:
			f, g = func_and_grad(x)
		elif task_str.startswith(b'NEW_X'):
			# new iteration
			n_iterations += 1
			if callback is not None:
				callback(np.copy(x))

			if n_iterations >= maxiter:
				task[:] = 'STOP: TOTAL NO. of ITERATIONS REACHED LIMIT'
			elif sf.nfev > maxfun:
				task[:] = ('STOP: TOTAL NO. of f AND g EVALUATIONS '
						   'EXCEEDS LIMIT')
		else:
			break

	task_str = task.tobytes().strip(b'\x00').strip()
	if task_str.startswith(b'CONV'):
		warnflag = 0
	elif sf.nfev > maxfun or n_iterations >= maxiter:
		warnflag = 1
	else:
		warnflag = 2

	# These two portions of the workspace are described in the mainlb
	# subroutine in lbfgsb.f. See line 363.
	s = wa[0: m*n].reshape(m, n)
	y = wa[m*n: 2*m*n].reshape(m, n)

	# See lbfgsb.f line 160 for this portion of the workspace.
	# isave(31) = the total number of BFGS updates prior the current iteration;
	n_bfgs_updates = isave[30]

	n_corrs = min(n_bfgs_updates, maxcor)
	hess_inv = LbfgsInvHessProduct(s[:n_corrs], y[:n_corrs])

	return OptimizeResult(fun=f, jac=g, nfev=sf.nfev, njev=sf.ngev,
						  nit=n_iterations, status=warnflag, message=task_str,
						  x=x, success=(warnflag == 0), hess_inv=hess_inv)
