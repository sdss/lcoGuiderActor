import distutils
from distutils.core import setup, Extension

import sdss3tools
import os

sdss3tools.setup(
        ext_modules=[Extension('libguide', 
                               sources=['src/gimg/ipGguide.c', 'src/gimg/gutils.c', 'src/gimg/shUtils.c'],
                               include_dirs=['include'],
                               )],

        description = "SDSS-3 guider actor.",
        )

