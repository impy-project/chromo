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
                      '-cpp', '-UNDEBUG',
                      language: 'fortran')
add_project_arguments('-fPIC',
                      '-Wno-unused-const-variable',
                      language: 'C')

# Optionally use mold linker if available
mold = find_program('mold', required: false)
if mold.found()
  message('Using mold linker')
  add_project_link_arguments('-fuse-ld=mold', language: ['c', 'cpp', 'fortran'])
endif

# Directory for utilities
scripts_dir = meson.project_source_root() + '/scripts/'

# Update submodules if needed
run_command('git', 'submodule', 'update', '--init', '--recursive',
  check: true, capture: true)

# Find Python and its development headers
py = import('python').find_installation(pure: false)
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
]

# QGSJET-specific definitions

qgsjet2_common_sources = [
  fdir + 'qgsjet/chromo_qgsII.f',
]

qgsjet2_common_functions = [
    'cqgsini',
    'qgsect',
    'qgini',
    'qgconf',
    'qgreg',
    'chepevt',
    'qgcrossc',
    'cqgshh_ha_cs'
]

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

urqmd_all_sources += urqmd_chromo_sources
urqmd_interface_sources += urqmd_chromo_sources

urqmd_functions = [
  'urqmd', 'init', 'uinit', 'set0', 'params', 'uounit', 'strini', 'loginit',
  'loadwtab', 'norm_init', 'output', 'cascinit', 'nucrad', 'urqini',
  'partname', 'chepevt', 'ptsigtot'
]

# Pythia6-specific definitions
pythia6_sources = [
  fdir / 'pythia6/pythia-6.4.28.f',
]
pythia6_functions = [
  'pyinit', 'pyexec', 'pytune', 'pylist', 'pyevnt', 'pyevnw', 'pystat',
  'pyedit', 'pyhepc', 'pychge', 'pycomp', 'pyk'
]

# SOPHIA-specific definitions
sophia_sources = [
  fdir / 'sophia/SOPHIA20.f',
  fdir / 'sophia/eventgen.f',
  fdir / 'sophia/sampling.f',
  fdir / 'sophia/inpoutput.f',
  fdir / 'sophia/jetset74dp.f',
  fdir / 'sophia/chromo_sophia.f',
]
sophia_functions = [
  'eventgen',
  'print_event',
  'crossection',
  'initial',
  'icon_pdg_sib',
  'toevt'
]

# DPMJET307-specific definitions
dpmjet307_sources = [
  fdir / 'dpmjet3.0/sources/dpmjet3.0-7.f',
  fdir / 'dpmjet3.0/sources/phojet1.12-36c4.f',
  fdir / 'dpmjet3.0/sources/pythia6115dpm3v1.f',
  fdir / 'dpmjet3.0/sources/user3.0-7.f',
]
dpmjet_phojet_functions = [
  'dt_evtout',
  'dt_eventb',
  'dt_ficonf',
  'dt_flahad',
  'dt_init',
  'dt_initjs',
  'dt_kkinc',
  'dt_phoxs',
  'dt_sigga',
  'dt_siggat',
  'dt_siggp',
  'dt_title',
  'dt_xsglau',
  'dt_xshn',
  'idt_icihad',
  'idt_icihad',
  'idt_ipdgha',
  'ipho_id2pdg',
  'ipho_pdg2id',
  'pho_borncs',
  'pho_event',
  'pho_ghhias',
  'pho_harint',
  'pho_harxto',
  'pho_harxpt',
  'pho_init',
  'pho_mcini',
  'pho_pmass',
  'pho_pname',
  'pho_prevnt',
  'pho_rregpar',
  'pho_setmdl',
  'pho_setpar',
  'pho_setpdf',
  'pho_setpcomb',
  'pho_sregpar',
  'pho_xsect',
  'poevt1',
  'poevt2',
  'pycomp',
  'pytune',
]

# DPMJETIII191-specific definitions
dpmjet191_fdir = fdir / 'dpmjetIII-19.1'
# Use Python script to collect source files and filter excluded patterns
dpmjet191_sources = run_command([
  py, scripts_dir + 'list_dpmjet_sources.py', '--directory', dpmjet191_fdir
], check: true).stdout().strip().split('\n')

# DPMJET193-specific definitions
dpmjet193_fdir = fdir / 'dpmjetIII-19.3'
# Use Python script to collect source files and filter excluded patterns
dpmjet193_sources = run_command([
  py, scripts_dir + 'list_dpmjet_sources.py', '--directory', dpmjet193_fdir
], check: true).stdout().strip().split('\n')

# EPOS-LHC-specific definitions
eposlhc_fdir = fdir / 'epos'
eposlhc_sources = []

# Use Python script to collect EPOS source files and filter excluded patterns
eposlhc_sources = run_command([
  py, scripts_dir + 'list_epos_sources.py', '--directory', eposlhc_fdir
], check: true).stdout().strip().split('\n')

# EPOS interface sources (specific files used for f2py interface)
eposlhc_interface_sources = run_command([
  py, scripts_dir + 'list_epos_sources.py', 
  '--directory', eposlhc_fdir, '--interface-only'
], check: true).stdout().strip().split('\n')

eposlhc_functions = [
  'aaset',
  'ainit',
  'aepos',
  'alist',
  'afinal',
  'hepmcstore',
  'getcharge',
  'idtrafo',
  'initepos',
  'initeposevt',
  'xsection',
  'rmmard',
  'ranfgt',
  'ranfst'
]



# Parametrize all targets. The models must have an _ underline prefix
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
    ],
    'functions': ','.join([
      'cqgsini','sectnu','xxaini','psconf','xxreg','psaini','chepevt','xxfz',
      'cqgshh_ha_cs','crossc']),
  },
  '_qgs2_03': {
    'sources': [
      fdir + 'qgsjet/qgsjet-II-03.f'] + qgsjet2_common_sources,
    'functions': ','.join(qgsjet2_common_functions)
  },
  '_qgs2_04': {
    'sources': [
      fdir + 'qgsjet/qgsjet-II-04.f'] + qgsjet2_common_sources,
    'functions': ','.join(qgsjet2_common_functions)
  },
  '_qgs3': {
    'sources': [
      fdir + 'qgsjet/qgsjet-III.f',
      fdir + 'qgsjet/chromo_qgsIII.f',
    ],
    'functions': ','.join(qgsjet2_common_functions + ['qgran'])
  },
  '_urqmd34': {
    'sources': urqmd_all_sources,
    'interface_sources': urqmd_interface_sources, # Adjust if filtering is done elsewhere
    'functions': ','.join(urqmd_functions),
    'include_directories': [fdir / 'urqmd-3.4/sources'],
  },
  '_pythia6': {
    'sources': pythia6_sources,
    'functions': ','.join(pythia6_functions),
  },
  '_sophia': {
    'sources': sophia_sources,
    'functions': ','.join(sophia_functions),
  },
  '_dpmjet_phojet307': {
    'sources': dpmjet307_sources,
    'functions': ','.join(dpmjet_phojet_functions),
  },
  '_dpmjet_phojet191': {
    'sources': dpmjet191_sources,
    'functions': ','.join(dpmjet_phojet_functions),
    'include_directories': [
      fdir / 'dpmjetIII-19.1/include/phojet',
      fdir / 'dpmjetIII-19.1/include/dpmjet', 
      fdir / 'dpmjetIII-19.1/include/pythia',
      fdir / 'dpmjetIII-19.1/include/flinclude'
    ],
  },
  '_dpmjet_phojet193': {
    'sources': dpmjet193_sources,
    'functions': ','.join(dpmjet_phojet_functions),
    'include_directories': [
      fdir / 'dpmjetIII-19.3/include/phojet',
      fdir / 'dpmjetIII-19.3/include/dpmjet', 
      fdir / 'dpmjetIII-19.3/include/pythia',
      fdir / 'dpmjetIII-19.3/include/flinclude'
    ],
  },
  '_eposlhc': {
    'sources': eposlhc_sources,
    'interface_sources': eposlhc_interface_sources,
    'functions': ','.join(eposlhc_functions),
    'include_directories': [fdir / 'epos/sources'],
  },
}



# Read from pyproject.toml which models to build
enabled_models = run_command([
  py, scripts_dir + 'get_models.py', '--enabled'], 
  check: true).stdout().strip().split('\n')
disabled_models = run_command([
  py, scripts_dir + 'get_models.py', '--disabled'], 
  check: true).stdout().strip().split('\n')


# TODO: Need to somehow let meson know that _dpmjet191 and _phojet191
# share 99% of files except the f2py interface files.

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

  # foreach f : files(model_info.get('interface_sources', model_info['sources']))
  #   message(f.full_path())
  # endforeach
  flags = chromo_common_flags + model_info.get('flags', [])
  # message('Flags: ' + ' '.join(flags))

  interface_sources = files(model_info.get('interface_sources', model_info['sources']) 
    + [logging_source, rangen_source])
  sources = files(model_info['sources'] + [logging_source, rangen_source])

  # Custom target to generate f2py wrapper
  f2py_target = custom_target(f2py_name,
    output: output,
    input: interface_sources,
    command: [
      py, scripts_dir + 'generate_f2py.py',
      f2py_name,
      model_info['functions'] + ',' + ','.join(chromo_common_functions),
      meson.global_build_root() / 'meson-logs',
      ' '.join(include_args),
      ' '.join(flags),
      '@OUTDIR@',
      '@INPUT@'
    ],
    build_by_default: true
  )

  # foreach f : output
  #   message('F2PY output: ' + f)
  # endforeach

  # Create the sib21 extension module
  ext_module = py.extension_module(f2py_name,
    sources: sources + f2py_target + [f2py_fortranobject_c] + [normal_source],
    include_directories: include_directories(includes),
    dependencies: [py_dep, numpy_dep],
    c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'] + flags,
    fortran_args: flags,
    install: true,
    subdir: 'chromo/models',
  )
  

  # For editable installs: copy the extension to the source tree
  # Unfortunately this also copies the extensions to src/chromo/models
  # for non-editable installs, but this is not a real problem.
  custom_target(f2py_name + '_copy_for_editable',
    input: ext_module,
    output: f2py_name + '.so',
    command: ['cp', '@INPUT@', meson.project_source_root() / 'src' / 'chromo' / 'models' / '@OUTPUT@'],
    build_by_default: true,
    install: false
  )
endforeach

# pythia8 About the license and how it affects chromo. Pythia-8 is licensed
# under GPL-v2 or higher, which means GPL-v3 applies. An excerpt from GPL-v3:

# > A compilation of a covered work with other separate and independent > works,
# which are not by their nature extensions of the covered work, > and which are
# not combined with it such as to form a larger program, > in or on a volume of
# a storage or distribution medium, is called an > "aggregate" if the
# compilation and its resulting copyright are not > used to limit the access or
# legal rights of the compilation's users > beyond what the individual works
# permit.  Inclusion of a covered work > in an aggregate does not cause this
# License to apply to the other > parts of the aggregate.
#
# Our legal position is that chromo is an "aggregate". It is not adding
# functionality to Pythia-8 nor is it modifying the original source. The
# Pythia-8 library is redistributed by us under its original license, which is
# allowed. chromo is merely a way to steer Pythia-8 and to provide its output in
# a particular format. The Pythia-8 authors are free to take some of our code
# and make it GPL code, they can do that anyway, but we are still free to
# distribute the code that we write under a less restrictive license.
#
# To strengthen this position, we compile the original Pythia-8 code into a
# shared library to which we only link. This makes clear where the boundary
# between chromo and Pythia-8 is.

if 'pythia8' in enabled_models
  message('Building model: _pythia8')
  cpp_dir = 'src/cpp'
  pybind11_inc = include_directories(cpp_dir / 'pybind11/include')

   if host_machine.system() == 'darwin'
    pythia8_rpath = '@loader_path'
  elif host_machine.system() == 'linux'
    pythia8_rpath = '$ORIGIN'
  else
    pythia8_rpath = ''
  endif

  # Find all Pythia8 source files
  pythia8_source_list = run_command([
    'find', meson.project_source_root() / cpp_dir / 'pythia83/src', 
    '-name', '*.cc', '-type', 'f'
  ], check: true).stdout().strip().split('\n')

  pythia8_sources = files(pythia8_source_list)

  # Create Pythia8 shared library
  libpythia8 = shared_library('pythia8',
    pythia8_sources,
    include_directories: include_directories(cpp_dir / 'pythia83/include'),
    cpp_args: ['-DXMLDIR=""'],  # We don't use this, must set xmldir at runtime
    install: true,
    install_dir: py.get_install_dir() / 'chromo' / 'models'
  )

  # Create the _pythia8 Python extension module
  pythia8_module = py.extension_module('_pythia8',
    cpp_dir / '_pythia8.cpp',
    include_directories: [
      pybind11_inc,
      include_directories(cpp_dir),
      include_directories(cpp_dir / 'pythia83/include')
    ],
    dependencies: [py_dep],
    link_with: libpythia8,
    install: true,
    subdir: 'chromo/models',
    install_rpath: pythia8_rpath,
    build_rpath: pythia8_rpath
  )

  # For editable installs: copy both libpythia8.so and _pythia8.so to source tree
  custom_target('copy_libpythia8_for_editable',
    input: libpythia8,
    output: 'libpythia8_editable.so',
    command: ['cp', '@INPUT@', meson.project_source_root() / 'src' / 'chromo' / 'models' / 'libpythia8.so'],
    build_by_default: true,
    install: false
  )

  custom_target('copy_pythia8_for_editable',
    input: pythia8_module,
    output: '_pythia8_editable.so',
    command: ['cp', '@INPUT@', meson.project_source_root() / 'src' / 'chromo' / 'models' / '_pythia8.so'],
    build_by_default: true,
    install: false
  )
  message('Pythia8 configuration complete')
else
  message('Skipping Pythia8 model (not enabled)')
endif


# Install the pure Python package files
install_subdir('src/chromo', install_dir: py.get_install_dir())

message('CHROMO\'s meson build configuration is complete.')


