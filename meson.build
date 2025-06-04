project('chromo', 'c', 'cpp', 'fortran',
  version: '2.7.0',
  license: 'MIT',
  meson_version: '>=1.0.0',
  default_options: [
    'buildtype=release',
    'cpp_std=c++14',
    'fortran_std=legacy',
    'warning_level=1',
  ]
)

# Find Python and its development headers
py = import('python').find_installation()
py_dep = py.dependency()

# Find numpy and its headers
numpy_dep = dependency('numpy')
numpy_inc = run_command(py, '-c', 'import numpy; print(numpy.get_include())', check: true).stdout().strip()

# Find f2py include directory
f2py_inc = run_command(py, '-c', 
  'import numpy.f2py; from pathlib import Path; print(Path(numpy.f2py.__file__).parent / "src")',
  check: true).stdout().strip()

# F2py source files
f2py_fortranobject_c = f2py_inc / 'fortranobject.c'

# Add include directories
inc_dirs = include_directories(
  numpy_inc,
  f2py_inc,
  'src/fortran',
)

# Compiler settings
if host_machine.system() != 'windows'
  # Add common compile arguments for Unix-like systems
  add_project_arguments([
    '-Wno-uninitialized',
  ], language: ['c', 'cpp'])
  
  add_project_arguments([
    '-std=legacy',
    '-fno-second-underscore',
    '-cpp',
  ], language: 'fortran')
  
  if host_machine.system() == 'darwin'
    # macOS specific flags for gfortran > 10
    fc = meson.get_compiler('fortran')
    if fc.get_id() == 'gcc' and fc.version().version_compare('>=10.0.0')
      add_project_arguments('-fallow-argument-mismatch', language: 'fortran')
    endif
  else
    add_project_arguments('-Wno-argument-mismatch', language: 'fortran')
  endif
else
  # Windows MinGW flags
  add_project_arguments([
    '-Wno-uninitialized',
    '-std=legacy',
    '-fno-second-underscore',
  ], language: 'fortran')
endif

# Add compile definitions
add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: ['c', 'cpp'])

# Common source files
fortran_dir = 'src/fortran'
logging_source = files('src/fortran/logging.f')
rangen_source = files('src/fortran/rangen.c')
normal_source = files('src/fortran/normal.c')

# Test building just sib21 first to verify f2py integration
sib21_sources = [
  files('src/fortran/sibyll/sibyll_21.f'),
  files('src/fortran/sibyll/sib21aux.f'),
  files('src/fortran/sibyll/sibyll_init.fpp'),
  logging_source,
  rangen_source,
]

# Custom target to generate f2py wrapper for sib21
sib21_f2py_target = custom_target('sib21_f2py',
  output: ['_sib21module.c', '_sib21-f2pywrappers.f', '_sib21.pyf'],
  input: sib21_sources[0],  # Use first source as input for dependency tracking
  command: [
    py, files('generate_f2py.py'),
    '_sib21',
    'sibyll,sibyll_ini,sib_sigma_hp,sib_sigma_hair,sib_sigma_hnuc',
    '@INPUT@',
    '',
    'src/fortran',
    '@OUTDIR@'
  ],
  build_by_default: true
)

# Create the sib21 extension module
py.extension_module('_sib21',
  sources: sib21_sources + sib21_f2py_target + [f2py_fortranobject_c] + [normal_source],
  include_directories: [inc_dirs],
  dependencies: [py_dep, numpy_dep],
  c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', '-DLINUX', '-DSIBYLL_21'],
  fortran_args: ['-DLINUX', '-DSIBYLL_21'],
  install: true,
  subdir: 'chromo/models',
)

message('Simple meson build configuration complete')
