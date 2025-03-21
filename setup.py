from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import sys
import os.path
from io import open


ext_modules = [
    Extension(
        '_phat',
        ['python/_phat.cpp'],
        include_dirs=['include'],
        language='c++',
    ),
]

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'python', 'README.rst'), encoding = 'utf8') as f:
    long_description = f.read()

class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': ['-std=c++11'],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        import pybind11
        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.include_dirs.append(pybind11.get_include())
            ext.include_dirs.append(pybind11.get_include(user=True))
        build_ext.build_extensions(self)

setup(
    name='phat',
    version='1.6.0.post2',
    author='Bryn Keller',
    author_email='bryn.keller@intel.com',
    url='https://bitbucket.org/phat-code/phat',
    description='Python bindings for PHAT, the C++ based Persistent Homology Algorithm Toolbox',
    license = 'LGPL',
    keywords='algebraic-topology PHAT distributed topology persistent-homology',
    long_description=long_description,
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    package_dir={'':'python'},
    py_modules = ['phat'],
    python_requires=">=3.9"
    # packages = find_packages(exclude = ['doc', 'test'])
 )





