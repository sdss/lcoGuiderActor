import distutils
import distutils.command.install as distInstall
from distutils.core import setup, Extension

import sdss3tools
import os
import glob

# jkp: This is a hack, but it should work to get the lib symlink installed.
# there may be a better way, if we can just get everything put into lib automatically.
class my_install(distInstall.install):
    def run(self):
        distInstall.install.run(self)
        build_path = glob.glob(os.path.join(self.install_lib,'build/lib.*'))[0]
        os.symlink(build_path, os.path.join(self.install_lib, 'lib'))

sdss3tools.setup(
        ext_modules=[Extension('libguide', 
                               sources=['src/gimg/ipGguide.c', 'src/gimg/gutils.c', 'src/gimg/shUtils.c'],
                               include_dirs=['include'],
                               )],

        description = "SDSS-3 guider actor.",
        cmdclass = dict(install=my_install),
        )

