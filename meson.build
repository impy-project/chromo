project('chromo', 'c', 'cpp', 'fortran',
  license: 'BSD-3-Clause',
  version: run_command([
    'python3', meson.project_source_root() + '/scripts/get_version.py'
    ], check: true).stdout().strip(),
  meson_version: '>=1.4.0',
  default_options: [
    'buildtype=release',
    'cpp_std=c++14',
    'fortran_std=legacy',
    'warning_level=1',
    'optimization=3'
  ]
)
add_project_arguments('-Wno-uninitialized', '-fPIC',
                      '-fno-second-underscore', '-Wno-unused', '-Wno-conversion',
                      '-Wno-unused-dummy-argument', '-Wno-character-truncation', 
                      '-Wno-integer-division','-Wno-surprising', '-Wno-tabs',
                      '-Wno-backslash-newline-escape',
                      '-cpp', '-UNDEBUG',
                      language: 'fortran')
add_project_arguments('-fPIC',
                      '-Wno-unused-const-variable',
                      language: 'C')

# Directory for utilities
scripts_dir = meson.project_source_root() + '/scripts/'

# Update submodules if needed
run_command('git', 'submodule', 'update', '--init', '--recursive',
  check: true, capture: true)

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
chromo_inc_dirs = [
  numpy_inc,
  f2py_inc,
  'src/fortran'
]

# Add compile definitions
add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language: ['c', 'cpp'])

# Common source files
fdir = 'src/fortran/'
logging_source = [fdir + 'logging.f']
rangen_source = [fdir + 'rangen.c', fdir + 'rangen.fpp']
normal_source = [fdir + 'normal.c']
chromo_common_functions = 'chromo_openlogfile,chromo_closelogfile,npyrng'
chromo_common_flags = ['-DCHROMO']

# SIBYLL-specific definitions
sibyll_sources_base = [
  fdir + 'sibyll/sibyll_init.fpp',
  logging_source,
  rangen_source,
]
sibyll_common_functions = [
    'sibyll',
    'sibyll_ini',
    'sib_sigma_hp',
    'sib_sigma_hair',
    'sib_sigma_hnuc',
    'int_nuc',
    'decsib',
    'decpar',
    'sibini',
    'sibhep',
    'sib_list',
    'isib_pid2pdg',
    'isib_pdg2pid',
    'pdg_ini'
] + chromo_common_functions

# QGSJET-specific definitions

qgsjetII_common_sources = [
  fdir + 'qgsjet/chromo_qgsII.f',
  logging_source,
  rangen_source
]

qgsjetII_common_functions = [
    'cqgsini',
    'qgsect',
    'qgini',
    'qgconf',
    'qgreg',
    'chepevt',
    'qgcrossc',
    'cqgshh_ha_cs'
] + chromo_common_functions

#UrQMD-specific definitions
urqmd_fdir = fdir / 'urqmd-3.4/sources'
urqmd_sources_f = [
  '1fluid.f', 'bessel.f', 'delpart.f', 'getmass.f', 'hepcmp.f', 'iso.f',
  'numrec.f', 'pythia6409.f', 'siglookup.f', 'upmerge.f', 'addpart.f',
  'blockres.f', 'coload.f', 'detbal.f', 'getspin.f', 'hepnam.f', 'ityp2pdg.f',
  'output.f', 'string.f', 'urqmd.f', 'angdis.f', 'boxprg.f', 'dwidth.f',
  'init.f', 'jdecay2.f', 'paulibl.f', 'saveinfo.f', 'tabinit.f', 'whichres.f',
  'anndec.f', 'cascinit.f', 'dectim.f', 'error.f', 'hepchg.f', 'input.f',
  'make22.f', 'proppot.f', 'scatter.f', 'uhmerge.f', 'urqinit.f'
]
urqmd_sources_f90 = ['CFmax.f90', 'quadri.f90', 'cornelius.f90']
urqmd_chromo_sources = [fdir / 'urqmd-3.4/chromo_urqmd.f']
urqmd_ignore_in_interface_sources = [
  'newpart.f', 'uhmerge.f', 'iso.f'
]
urqmd_all_sources = []
urqmd_interface_sources = []

foreach s : urqmd_sources_f
  urqmd_all_sources += urqmd_fdir / s
  if s not in urqmd_ignore_in_interface_sources
    urqmd_interface_sources += urqmd_fdir / s
  endif
endforeach
foreach s : urqmd_sources_f90
  urqmd_all_sources += urqmd_fdir / s
endforeach

urqmd_all_sources += logging_source + rangen_source
urqmd_interface_sources += logging_source + rangen_source

urqmd_functions = [
  'urqmd', 'init', 'uinit', 'set0', 'params', 'uounit', 'strini', 'loginit',
  'loadwtab', 'norm_init', 'output', 'cascinit', 'nucrad', 'urqini',
  'partname', 'chepevt', 'ptsigtot'
] + chromo_common_functions

# Pythia6-specific definitions
pythia6_sources = [
  fdir / 'pythia6/pythia-6.4.28.f',
  logging_source,
  rangen_source
]
pythia6_functions = [
  'pyinit', 'pyexec', 'pytune', 'pylist', 'pyevnt', 'pyevnw', 'pystat',
  'pyedit', 'pyhepc', 'pychge', 'pycomp', 'pyk'
] + chromo_common_functions

# DPMJET and PHOJET-specific definitions

# Paramterize all targets. The models must have an _ underline prefix
# to avoid conflicts with the Python module names.
models = {
  '_sib21': {
    'sources': sibyll_sources_base + [fdir + 'sibyll/sibyll_21.f', fdir + 'sibyll/sib21aux.f'],
    'functions': ','.join(sibyll_common_functions + ['spgasdev', 'glauber']),
    'flags' : ['-DSIBYLL_21'],
  },
  '_sib23': {
    'sources': sibyll_sources_base + [fdir + 'sibyll/sibyll2.3.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23c00': {
    'sources': sibyll_sources_base + [fdir + 'sibyll/sibyll2.3c00.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23c01': {
    'sources': sibyll_sources_base + [fdir + 'sibyll/sibyll2.3c01.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23c02': {
    'sources': sibyll_sources_base + [fdir + 'sibyll/sibyll2.3c02.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23c03': {
    'sources': sibyll_sources_base + [fdir + 'sibyll/sibyll2.3c03.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23d': {
    'sources': sibyll_sources_base + [fdir + 'sibyll/sibyll2.3d.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_sib23d_star': {
    'sources': sibyll_sources_base + [fdir + 'sibyll/sibyll2.3d-star-p02.f'],
    'functions': ','.join(sibyll_common_functions + ['gasdev', 'sig_had_nuc']),
  },
  '_qgs01': {
    'sources': [
      fdir + 'qgsjet/qgsjet01d.f',
      fdir + 'qgsjet/chromo_qgs1.f',
      logging_source,
      rangen_source
    ],
    'functions': ','.join([
      'cqgsini','sectnu','xxaini','psconf','xxreg','psaini','chepevt','xxfz',
      'cqgshh_ha_cs','crossc'] + chromo_common_functions),
  },
  '_qgsII03': {
    'sources': [
      fdir + 'qgsjet/qgsjet-II-03.f'] + qgsjetII_common_sources,
    'functions': ','.join(qgsjetII_common_functions)
  },
  '_qgsII04': {
    'sources': [
      fdir + 'qgsjet/qgsjet-II-04.f'] + qgsjetII_common_sources,
    'functions': ','.join(qgsjetII_common_functions)
  },
  '_qgsIII': {
    'sources': [
      fdir + 'qgsjet/qgsjet-III.f',
      fdir + 'qgsjet/chromo_qgsIII.f',
      logging_source,
      rangen_source
    ],
    'functions': ','.join(qgsjetII_common_functions + ['qgran'])
  },
  '_urqmd34': {
    'sources': urqmd_all_sources,
    'interface_sources': urqmd_interface_sources, # Adjust if filtering is done elsewhere
    'functions': ','.join(urqmd_functions),
    'include_directories': [fdir / 'urqmd-3.4/sources'],
    'flags': chromo_common_flags,
  },
  '_pythia6': {
    'sources': pythia6_sources,
    'functions': ','.join(pythia6_functions),
    'flags': chromo_common_flags,
  },
}

# Read from pyproject.toml which models to build
enabled_models = run_command([
  py, scripts_dir + 'get_models.py', '--enabled'], 
  check: true).stdout().strip().split('\n')
disabled_models = run_command([
  py, scripts_dir + 'get_models.py', '--disabled'], 
  check: true).stdout().strip().split('\n')

foreach f2py_name, model_info : models
  if f2py_name.substring(1) not in enabled_models
    message('Skipping model: ' + f2py_name + ' (disabled)')
    continue
  endif

  message('Building model: ' + f2py_name)
  output = [f2py_name + 'module.c', f2py_name + '-f2pywrappers.f', f2py_name + '.pyf']
  includes = chromo_inc_dirs + model_info.get('include_directories', [])
  include_args = []
  foreach p : includes
    include_args += '-I' + (meson.project_source_root() / p)
  endforeach

  foreach f : files(model_info.get('interface_sources', model_info['sources']))
    message(f.full_path())
  endforeach
  flags = chromo_common_flags + model_info.get('flags', [])
  message('Flags: ' + ' '.join(flags))

  interface_sources = files(model_info.get('interface_sources', model_info['sources']))
  sources = files(model_info['sources'])

  # Custom target to generate f2py wrapper
  f2py_target = custom_target(f2py_name,
    output: output,
    input: interface_sources,
    command: [
      py, scripts_dir + 'generate_f2py.py',
      f2py_name,
      model_info['functions'],
      meson.global_build_root() / 'meson-logs',
      ' '.join(include_args),
      ' '.join(flags),
      '@OUTDIR@',
      '@INPUT@'
    ],
    build_by_default: true
  )

  foreach f : output
    message('F2PY output: ' + f)
  endforeach
  # Create the sib21 extension module
  py.extension_module(f2py_name,
    sources: sources + f2py_target + [f2py_fortranobject_c] + [normal_source],
    include_directories: include_directories(includes),
    dependencies: [py_dep, numpy_dep],
    c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'] + flags,
    fortran_args: flags,
    install: true,
    subdir: 'chromo/models',
  )
endforeach

message('Multi-model meson build configuration complete with SIBYLL-21, EPOS-LHC, and DPMJET III-19.1')
