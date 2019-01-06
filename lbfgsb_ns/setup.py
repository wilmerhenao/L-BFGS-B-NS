import numpy.distutils
import numpy.distutils.system_info


def configuration(parent_package='', top_path=None):
	from numpy.distutils.misc_util import Configuration
	config = Configuration('lbfgsb_ns', parent_package, top_path)
	lapack = numpy.distutils.system_info.get_info('lapack_opt')
	blas = numpy.distutils.system_info.get_info('blas_opt')

	if len(lapack) == 0:
		raise ValueError("System must have LAPACK library linked to SciPy")

	if len(blas) == 0:
		raise ValueError("System must have BLAS library linked to SciPy")

	dct_keys = set(lapack.keys()).union(set(blas.keys()))
	for k in dct_keys:
		lapack[k] = list(set(lapack[k] + blas[k]))

	sources = ['lbfgsbns.pyf', '../lbfgsbNS.f90', '../timer.f', '../linpack.f']
	config.add_extension('_lbfgsb_ns',
							 sources=sources,
							 **lapack)

	return config
