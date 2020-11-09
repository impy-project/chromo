"""This setup.py is adopted from iMinuit https://github.com/scikit-hep/iminuit and
https://github.com/pybind/cmake_example/blob/master/setup.py
"""

import os
import re
import sys
import platform
import subprocess
import glob
from os.path import join, dirname, abspath

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class MakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class MakeBuild(build_ext):
    def run(self):
        if platform.system() == "Windows":
            try:
                out = subprocess.check_output(
                    ['mingw32-make.exe', '--version'])
            except OSError:
                raise RuntimeError(
                    "mingw32-make.exe must be installed to build the following extensions: "
                    + ", ".join(e.name for e in self.extensions))
        else:
            try:
                out = subprocess.check_output(['make', '--version'])
            except OSError:
                raise RuntimeError(
                    "Make must be installed to build the following extensions: "
                    + ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        import sysconfig
        import multiprocessing

        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        suffix = sysconfig.get_config_var('EXT_SUFFIX')
        # Some Python 2.7 versions don't define EXT_SUFFIX
        if suffix is None and 'SO' in sysconfig.get_config_vars():
            suffix = sysconfig.get_config_var('SO')

        if platform.system() == "Windows":
            make_command = 'mingw32-make.exe'
            libext = ('.cp' + sysconfig.get_config_var('py_version_nodot') +
                      '-' + sysconfig.get_platform().replace('-', '_') +
                      suffix)

            build_args = [
                'LIB_DIR=' + (extdir + '/impy/lib').replace('/', '\\'),
                'LEXT=' + libext
            ]
        else:
            make_command = 'make'
            libext = suffix
            build_args = ['LIB_DIR=' + extdir + '/impy/lib', 'LEXT=' + libext]

        n_parallel_builds = multiprocessing.cpu_count()

        cfg = 'Debug' if self.debug else 'Release'
        build_args += ['Config=' + cfg]

        make_args = []
        make_args += ['-j' + str(n_parallel_builds)]

        # env = os.environ.copy()
        # env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
        #                                                       self.distribution.get_version())
        print(extdir, make_args, build_args)
        if not os.path.exists(extdir):
            os.makedirs(extdir)
        subprocess.check_call([make_command] + build_args + make_args)


def get_version():
    version = {}
    with open("impy/version.py") as fp:
        exec(fp.read(), version)
    return version['__version__']


__version__ = get_version()

# Require pytest-runner only when running tests
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []


def extract_longdescription():
    this_directory = abspath(dirname(__file__))
    if sys.version_info.major == 3:
        with open(join(this_directory, 'README.md'), encoding='utf-8') as f:
            long_description = f.read()
    else:
        with open(join(this_directory, 'README.md')) as f:
            long_description = f.read()

    skip_marker = "# impy"
    return long_description[long_description.index(skip_marker):].lstrip()


# Data files for interaction models
iamfiles = glob.glob('impy/iamdata/*/*')

setup(
    name='impy',
    version=__version__,
    description='Hadronic Interaction Model interface in PYthon',
    long_description=extract_longdescription(),
    long_description_content_type='text/markdown',
    author='Anatoli Fedynitch',
    maintainer_email='afedynitch@gmail.com',
    license='BSD 3-Clause License',
    url='https://github.com/afedynitch/impy',
    setup_requires=[] + pytest_runner,
    packages=['impy', 'impy.models'],
    data_files=[('', ['LICENSE', 'impy/impy_config.yaml']),
                ('iamdata', iamfiles)],
    include_package_data=True,
    install_requires=['six', 'particletools', 'numpy', 'pyyaml', 'pyhepmc-ng', 'scipy'],
    ext_modules=[MakeExtension('impy_libs')],
    cmdclass=dict(build_ext=MakeBuild),
    zip_safe=False,
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License'
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
    ])