import numpy.distutils
from numpy.distutils.command.build_ext import build_ext
from numpy.distutils.misc_util import Configuration
import numpy.distutils.system_info
from numpy.distutils.core import setup

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

