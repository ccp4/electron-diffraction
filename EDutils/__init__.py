import utils,os
from . import _version
__version__=utils._check_version(os.path.join(os.path.dirname(__file__)),_version.version)
