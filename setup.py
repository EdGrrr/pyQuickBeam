from numpy.distutils.core import setup, Extension

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
    version='0.1',
    author='Edward Gryspeerdt',
    author_email='e.gryspeerdt@imperial.ac.uk',
    maintainer='Edward Gryspeerdt',
    maintainer_email='e.gryspeerdt@imperial.ac.uk',
    description='A simple python interface/wrapper for the QuickBeam radar simulator',
    license='Quickbeam license',
    url='https://github.com/EdGrrr/pyQuickBeam',
    install_requires=['numpy',
                      'matplotlib'],
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
