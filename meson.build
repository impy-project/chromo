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
# Check if mold found and compiler is gcc
if mold.found() and meson.get_compiler('c').get_id() == 'gcc'
  message('Using mold linker')
  add_project_link_arguments('-fuse-ld=mold', language: ['c', 'cpp', 'fortran'])
endif

cc = meson.get_compiler('cpp')
dl_dep = cc.find_library('dl', required: true)
thread_dep = dependency('threads')

# Directory for utilities
scripts_dir = meson.project_source_root() + '/scripts/'

# Update submodules if needed
run_command('git', 'submodule', 'update', '--init', '--recursive',
  check: false, capture: true)

# Find Python and its development headers
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# Find numpy and its headers
numpy_dep = dependency('numpy')
numpy_inc = numpy_dep.get_variable('cflags').split()[0].replace('-I', '')

f2py_inc    = run_command(py, '-c',
  'import numpy.f2py, pathlib, sys; print((pathlib.Path(numpy.f2py.__file__).parent ))',
  check : true).stdout().strip() / 'src'

# inc_f2py  = include_directories(f2py_inc)

# F2py source files
# f2py_fortranobject_c = f2py_inc / 'fortranobject.c'
f2py_fortranobject_c = configure_file(
                 input  : join_paths(f2py_inc, 'fortranobject.c'),
                 output : 'fortranobject.c',
                 copy   : true) 

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
  # _dpmjet_phojetXXX entries are removed from here and handled separately
  '_eposlhc': {
    'sources': eposlhc_sources,
    'interface_sources': eposlhc_interface_sources,
    'functions': ','.join(eposlhc_functions),
    'include_directories': [fdir / 'epos/sources'],
  },
}

# DPMJET/PHOJET models that will be split into shared library and separate f2py modules
dpmjet_phojet_configs = {
  '307': {
    'sources': dpmjet307_sources,
    'base_model_name': '_dpmjet_phojet307', # For checking if group is enabled
    # No specific include_directories for 307, uses chromo_inc_dirs
  },
  '191': {
    'sources': dpmjet191_sources,
    'base_model_name': '_dpmjet_phojet191',
    'include_directories': [
      fdir / 'dpmjetIII-19.1/include/phojet',
      fdir / 'dpmjetIII-19.1/include/dpmjet',
      fdir / 'dpmjetIII-19.1/include/pythia',
      fdir / 'dpmjetIII-19.1/include/flinclude'
    ],
  },
  '193': {
    'sources': dpmjet193_sources,
    'base_model_name': '_dpmjet_phojet193',
    'include_directories': [
      fdir / 'dpmjetIII-19.3/include/phojet',
      fdir / 'dpmjetIII-19.3/include/dpmjet',
      fdir / 'dpmjetIII-19.3/include/pythia',
      fdir / 'dpmjetIII-19.3/include/flinclude'
    ],
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

# Build all models except DPMJET/PHOJET in a compact loop
foreach f2py_name, model_info : models
  if f2py_name.startswith('_dpmjet_phojet')
    continue
  endif
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
  flags = chromo_common_flags + model_info.get('flags', [])
  interface_sources = files(model_info.get('interface_sources', model_info['sources']) + [logging_source, rangen_source])
  sources = files(model_info['sources'] + [logging_source, rangen_source])
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
  foreach incdir : includes
    message('Adding include directory: ' + incdir)
  endforeach
  ext_module = py.extension_module(f2py_name,
    sources: sources + f2py_target + [f2py_fortranobject_c] + [normal_source],
    include_directories: include_directories(includes),
    dependencies: [py_dep, numpy_dep],
    c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'] + flags,
    fortran_args: flags,
    install: true,
    subdir: 'chromo/models',
  )
endforeach

# Compact DPMJET/PHOJET build loop
foreach version_suffix, model_config : dpmjet_phojet_configs
  dpmjet_module_name_short = 'dpmjet' + version_suffix
  phojet_module_name_short = 'phojet' + version_suffix
  base_model_name_short = model_config['base_model_name'].substring(1)
  should_build_dpmjet_part = dpmjet_module_name_short in enabled_models
  should_build_phojet_part = phojet_module_name_short in enabled_models
  original_group_enabled = base_model_name_short in enabled_models
  if not (should_build_dpmjet_part or should_build_phojet_part or original_group_enabled)
    message('Skipping DPMJET/PHOJET model group for version: ' + version_suffix + ' (all parts disabled)')
    continue
  endif
  if original_group_enabled
    should_build_dpmjet_part = dpmjet_module_name_short in enabled_models or (dpmjet_module_name_short not in disabled_models)
    should_build_phojet_part = phojet_module_name_short in enabled_models or (phojet_module_name_short not in disabled_models)
  endif
  if not (should_build_dpmjet_part or should_build_phojet_part)
    message('Skipping DPMJET/PHOJET model group for version: ' + version_suffix + ' (resolved to no parts enabled)')
    continue
  endif
  message('Building DPMJET/PHOJET model group for version: ' + version_suffix)
  model_sources_str = model_config['sources']
  model_specific_includes = model_config.get('include_directories', [])
  current_includes = chromo_inc_dirs + model_specific_includes
  current_include_args = []
  foreach p : current_includes
    current_include_args += '-I' + (meson.project_source_root() / p)
  endforeach
  current_flags = chromo_common_flags
  model_rpath = ''
  if host_machine.system() == 'darwin'
    model_rpath = '@loader_path'
  elif host_machine.system() == 'linux'
    model_rpath = '$ORIGIN'
  endif
  lib_name = 'dpmjet_phojet_' + version_suffix
  all_sources_for_lib = files(model_sources_str + logging_source + rangen_source + normal_source)
  lib_dpmjet_phojet = static_library(lib_name, # Changed from shared_library
    all_sources_for_lib,
    include_directories: include_directories(current_includes),
    dependencies: [dl_dep, thread_dep], # These might still be needed if the static lib itself uses them
    c_args: ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'] + current_flags,
    fortran_args: current_flags
    # install: true, # Not typically needed for static libraries
    # install_dir: py.get_install_dir() / 'chromo' / 'models' # Not typically needed for static libraries
  )
  common_ext_args = {
    'include_directories': include_directories(current_includes),
    'dependencies': [py_dep, numpy_dep],
    'link_with': [lib_dpmjet_phojet],
    'c_args': ['-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'] + current_flags,
    'install': true,
    'subdir': 'chromo/models',
    'install_rpath': model_rpath, # Restored: rpath for the extension module
    'build_rpath': model_rpath  # Restored: rpath for the extension module
  }

  f2py_interface_sources = files(model_sources_str + logging_source + rangen_source)
  # DPMJET module
  dpmjet_f2py_module_name = '_dpmjet' + version_suffix
  if should_build_dpmjet_part
    message('Building DPMJET module: ' + dpmjet_f2py_module_name)
    dpmjet_functions_to_wrap = dpmjet_phojet_functions + chromo_common_functions.split(',')
    dpmjet_f2py_output_files = [dpmjet_f2py_module_name + 'module.c', dpmjet_f2py_module_name + '-f2pywrappers.f', dpmjet_f2py_module_name + '.pyf']
    dpmjet_f2py_target = custom_target(dpmjet_f2py_module_name + '_f2py',
      output: dpmjet_f2py_output_files,
      input: f2py_interface_sources,
      command: [
        py, scripts_dir + 'generate_f2py.py',
        dpmjet_f2py_module_name,
        ','.join(dpmjet_functions_to_wrap),
        meson.global_build_root() / 'meson-logs',
        ' '.join(current_include_args),
        ' '.join(current_flags),
        '@OUTDIR@',
        '@INPUT@'
      ],
      build_by_default: true
    )
    dpmjet_ext_module = py.extension_module(dpmjet_f2py_module_name,
      [dpmjet_f2py_target, f2py_fortranobject_c],
      kwargs: common_ext_args
    )
    # custom_target(dpmjet_f2py_module_name + '_copy_for_editable',
    #   input: dpmjet_ext_module,
    #   output: dpmjet_f2py_module_name + '.so',
    #   command: ['cp', '@INPUT@', meson.project_source_root() / 'src' / 'chromo' / 'models' / '@OUTPUT@'],
    #   build_by_default: true,
    #   install: false
    # )
  else
    message('Skipping DPMJET module: ' + dpmjet_f2py_module_name + ' (disabled)')
  endif
  # PHOJET module
  phojet_f2py_module_name = '_phojet' + version_suffix
  if should_build_phojet_part
    message('Building PHOJET module: ' + phojet_f2py_module_name)
    phojet_functions_to_wrap = dpmjet_phojet_functions + chromo_common_functions.split(',')
    phojet_f2py_output_files = [phojet_f2py_module_name + 'module.c', phojet_f2py_module_name + '-f2pywrappers.f', phojet_f2py_module_name + '.pyf']
    phojet_f2py_target = custom_target(phojet_f2py_module_name + '_f2py',
      output: phojet_f2py_output_files,
      input: f2py_interface_sources,
      command: [
        py, scripts_dir + 'generate_f2py.py',
        phojet_f2py_module_name,
        ','.join(phojet_functions_to_wrap),
        meson.global_build_root() / 'meson-logs',
        ' '.join(current_include_args),
        ' '.join(current_flags),
        '@OUTDIR@',
        '@INPUT@'
      ],
      build_by_default: true
    )
    phojet_ext_module = py.extension_module(phojet_f2py_module_name,
      [phojet_f2py_target, f2py_fortranobject_c],
      kwargs: common_ext_args
    )
  else
    message('Skipping PHOJET module: ' + phojet_f2py_module_name + ' (disabled)')
  endif
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
    dependencies: [dl_dep, thread_dep],
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

  message('Pythia8 configuration complete')
else
  message('Skipping Pythia8 model (not enabled)')
endif


# Install the pure Python package files, but exclude iamdata
# Install all .py files in src/chromo
install_data('src/chromo/__init__.py', install_dir: py.get_install_dir() / 'chromo')
install_data('src/chromo/__main__.py', install_dir: py.get_install_dir() / 'chromo')
install_data('src/chromo/cli.py', install_dir: py.get_install_dir() / 'chromo')
install_data('src/chromo/common.py', install_dir: py.get_install_dir() / 'chromo')
install_data('src/chromo/constants.py', install_dir: py.get_install_dir() / 'chromo')
install_data('src/chromo/decay_handler.py', install_dir: py.get_install_dir() / 'chromo')
install_data('src/chromo/kinematics.py', install_dir: py.get_install_dir() / 'chromo')
install_data('src/chromo/util.py', install_dir: py.get_install_dir() / 'chromo')
install_data('src/chromo/writer.py', install_dir: py.get_install_dir() / 'chromo')
# Install all subdirs except iamdata
install_subdir('src/chromo/models', install_dir: py.get_install_dir() / 'chromo')

