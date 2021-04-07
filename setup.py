from setuptools import setup

import versioneer


# Most arguments are set in the `setup.cfg`.
setup(version=versioneer.get_version(), cmdclass=versioneer.get_cmdclass())
