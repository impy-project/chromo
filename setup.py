import os
import re
import sys
import platform
import subprocess
import glob
from os.path import join, dirname, abspath


from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class MakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class MakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['make', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            # cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            # if cmake_version < '3.1.0':
                # raise RuntimeError("CMake >= 3.1.0 is required on Windows")
            raise RuntimeError("impy does not support builds on Windows, yet.")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        import sysconfig
        import multiprocessing

        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        libext = sysconfig.get_config_var("EXT_SUFFIX")
        make_args = ['LIB_DIR=' + extdir,
                     'LEXT=' + libext]
        n_parallel_builds = multiprocessing.cpu_count()

        # cfg = 'Debug' if self.debug else 'Release'
        # build_args = ['--config', cfg]
        build_args = []

        if platform.system() == "Windows":
            # cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            # if sys.maxsize > 2**32:
            #     cmake_args += ['-A', 'x64']
            # build_args += ['--', '/m']
            raise RuntimeError("impy does not support builds on Windows, yet.")
        else:
            # cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['-j' + str(n_parallel_builds)]

        env = os.environ.copy()
        # env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
        #                                                       self.distribution.get_version())
        print(extdir, make_args, build_args)
        if not os.path.exists(extdir):
            os.makedirs(extdir)
        # subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['make'] + make_args + build_args)
        # subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

# This method is adopted from iMinuit https://github.com/scikit-hep/iminuit
# Getting the version number at this point is a bit tricky in Python:
# https://packaging.python.org/en/latest/development.html#single-sourcing-the-version-across-setup-py-and-your-project
# This is one of the recommended methods that works in Python 2 and 3:
def get_version():
    version = {}
    with open("impy/version.py") as fp:
        exec (fp.read(), version)
    return version['__version__']

__version__ = get_version()

# Require pytest-runner only when running tests
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

this_directory = abspath(dirname(__file__))
if sys.version_info.major == 3:
    with open(join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
else:
    with open(join(this_directory, 'README.md')) as f:
        long_description = f.read()

skip_marker = "# impy"
long_description = long_description[long_description.index(skip_marker) :].lstrip()

# Data files for interaction models
iamfiles = glob.glob('iamdata/*')

setup(
    name='impy',
    version=__version__,
    description='Hadronic Interaction Model interface in PYthon',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Anatoli Fedynitch and others',
    maintainer_email='afedynitch@gmail.com',
    license='BSD 3-Clause License',
    url='https://github.com/afedynitch/impy',
    setup_requires=[] + pytest_runner,
    packages=['impy', 'impy.models'],
    package_data={
        'impy': ['impy/impy_config.yaml'] + iamfiles
    },
    install_requires=[
        'six',
        'particletools',
        'numpy',
        'pyyaml',
        'pyhepmc-ng'
    ],
    # ext_modules=[libnrlmsise00, libcorsikaatm],
    ext_modules=[MakeExtension('cmake_example')],
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
    ])