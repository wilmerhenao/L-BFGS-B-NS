# L-BFGS-B-NS

Minimize a differentiable non-smooth function subject to bound constraints. Same as L-BFGS-B, but using a line search strategy that satisfies weak Wolfe conditions instead of strong Wolfe conditions, which aids with non-smoothness.

Code is based on L-BFGS-B (see references at the bottom), and is written in Fortran with a Python wrapper that follows the same API as `scipy.minimize`. See file `DriverRosenbrockp.f90` for a FORTRAN example, and `python_example/lbfgsb_ns.ipynb` for a Python example.

## Python installation

Clone or download the repository and install with:
```
python setup.py install
```

(Requires `numpy`, `scipy`, and `gfortran` or other compiler)


## Notes

The documentation of the FORTRAN code is not up-to-date. Look at the examples for a better idea of the array sizes that it takes as inputs. Note also that:
* If the variables turn to NA it might not stop automatically.
* It might throw final message `'ABNORMAL_TERMINATION_IN_LNSRCH'` when it's actually converged.

## References
* Zhu, C., Byrd, R.H., Lu, P. and Nocedal, J., 1997. Algorithm 778: L-BFGS-B: Fortran subroutines for large-scale bound-constrained optimization. ACM Transactions on Mathematical Software (TOMS), 23(4), pp.550-560.
* Henao, W., 2014. An L-BFGS-B-NS optimizer for non-smooth functions. Master's thesis, Courant Institute of Mathematical Science, New York University.
