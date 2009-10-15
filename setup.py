import distutils
from distutils.core import setup, Extension

import sdss3tools
import os

sdss3tools.setup(
        description = "SDSS-3 guider actor.",
        )

