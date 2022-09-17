from numpy.distutils.core import setup, Extension
import re


# Version code from https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
VERSIONFILE = "pyQuickBeam/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


######################
# FORTRAN code setup #
######################
f90_args = ['-fPIC', '-O']

wrappers = [
    Extension(
        'radsim',
        sources=['src/radsim.pyf',
                 'src/math_lib.f90',
                 'src/dsd.f90',
                 'src/array_lib.f90',
                 'src/gases.f90',
                 'src/optical_sphere.f90',
                 'src/optics_lib.f90',
                 'src/zeff.f90'],
        libraries=['math_lib',
                   'array_lib',
                   'm_mrgrnk'],
        extra_f90_compile_args=f90_args),]

setup(
    name='pyQuickBeam',
    version=verstr,
    author='Edward Gryspeerdt',
    author_email='e.gryspeerdt@imperial.ac.uk',
    maintainer='Edward Gryspeerdt',
    maintainer_email='e.gryspeerdt@imperial.ac.uk',
    description='A simple python interface/wrapper for the QuickBeam radar simulator',
    license='Quickbeam license',
    url='https://github.com/EdGrrr/pyQuickBeam',
    install_requires=['numpy'],
    packages=['pyQuickBeam'],
    package_data = {'pyQuickBeam': ['data/*']},
    libraries = [
        ('m_mrgrnk', dict(
            sources=['src/m_mrgrnk.f90'],
            extra_f90_compile_args=f90_args)),
        ('array_lib', dict(
            sources=['src/array_lib.f90',
                     'src/m_mrgrnk.f90'],
            extra_f90_compile_args=f90_args)),
        ('math_lib', dict(
            sources=['src/math_lib.f90',
                     'src/m_mrgrnk.f90'],
            extra_f90_compile_args=f90_args)),
    ],
    ext_modules = wrappers
)
