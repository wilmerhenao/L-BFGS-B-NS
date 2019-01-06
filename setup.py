import numpy.distutils
from numpy.distutils.command.build_ext import build_ext
from numpy.distutils.misc_util import Configuration
import numpy.distutils.system_info
from numpy.distutils.core import setup
# from numpy.distutils.core import setup as setup_np

# config = Configuration('_fortran', "", None)
# config = Configuration('lbfgsb_ns')
# lapack = numpy.distutils.system_info.get_info('lapack_opt')
# blas = numpy.distutils.system_info.get_info('blas_opt')

# if len(lapack) == 0:
# 	raise ValueError("System must have LAPACK library linked to SciPy")

# if len(blas) == 0:
# 	raise ValueError("System must have BLAS library linked to SciPy")

# dct_keys = set(lapack.keys()).union(set(blas.keys()))
# for k in dct_keys:
# 	lapack[k] = list(set(lapack[k] + blas[k]))

# sources = ['lbfgsb_ns/lbfgsbns.pyf', 'lbfgsbNS.f90', 'timer.f', 'linpack.f']
# config.add_extension('_lbfgsb_ns',
#                          sources=sources,
#                          **lapack)

# setup(**configuration(top_path='').todict())
# setup(**config.todict())

# print(config.todict())

# from numpy.distutils.core import setup as np_setup
    # setup(name = "lbfgsb_ns",
    #     author = "Wilmer Henao (wrapper by David Cortes)"
    #     **configuration(top_path='').todict())
    # from setuptools import setup
    # setup(name = "lbfgsb_ns", ext_modules = [configuration])

# config = config.todict()
# print(config)

# config["version"] = "0.1"


# setup(**config)

# from setuptools import setup
# from distutils.extension import Extension

# setup(
# 	name = "lbfgsb_ns",
# 	packages = ["lbfgsb_ns"],
# 	# cmdclass = {'build_ext': build_ext},
# 	ext_modules = [Extension(config["ext_modules"][0])]
# 	# ext_modules = config
# 	)

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py = True,
                       assume_default_configuration = True,
                       delegate_options_to_subpackages = True,
                       quiet = True)

    config.add_subpackage('lbfgsb_ns')

    return config


setup(
	name = "lbfgsb_ns",
	version = "0.1",
	configuration = configuration
	)

