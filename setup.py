from numpy.distutils.core import setup, Extension

wrappers = [
    Extension('radsim', sources=[
        'radsim.pyf',
        'math_lib.f90',
        'dsd.f90',
        'array_lib.f90',
        'gases.f90',
        'optical_sphere.f90',
        'optics_lib.f90',
        'zeff.f90'], libraries=[
            'math_lib',
            'array_lib',
            'm_mrgrnk']),]
              

f90_args = ['-fPIC', '-O']

setup(
    name='pyQuickBeam',
    version='0.1',
    author='Edward Gryspeerdt',
    author_email='e.gryspeerdt@imperial.ac.uk',
    maintainer='Edward Gryspeerdt',
    maintainer_email='e.gryspeerdt@imperial.ac.uk',
    description='A simple python interface/wrapper for the QuickBeam radar simulator',
    license='BSD (3-clause)',
    url='https://github.com/EdGrrr/pyQuickBeam',
    install_requires=['numpy',
                      'matplotlib'],
    packages=['pyQuickBeam'],
    package_data = {'pyQuickBeam': ['data/*']},
    libraries = [
        ('m_mrgrnk', dict(
            sources=['m_mrgrnk.f90'],
            extra_f90_compile_args=f90_args)),
        ('array_lib', dict(
            sources=['array_lib.f90',
                     'm_mrgrnk.f90'],
            extra_f90_compile_args=f90_args)),
        ('math_lib', dict(
            sources=['math_lib.f90',
                     'm_mrgrnk.f90'],
            extra_f90_compile_args=f90_args)),
    ],
    ext_modules = wrappers
)
