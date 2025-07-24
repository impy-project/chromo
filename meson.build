project(
  'chromo', ['c', 'cpp', 'fortran'],
  version        : run_command(['python', meson.project_source_root() + '/scripts/get_version.py'],
                               check:true).stdout().strip(),
  license        : 'BSD-3-Clause',
  meson_version  : '>=1.4',
  default_options: [
    'buildtype=release',
    'cpp_std=c++17',           # Pythia 8 likes C++17
    'fortran_std=legacy',
    'warning_level=1',
    'optimization=3',
    'b_staticpic=true',  # static ↔ PIC
  ],
)

# ── Core dependencies ────────────────────────────────────────────────────────────────
py          = import('python').find_installation(pure:false)
py_dep      = py.dependency()
numpy_dep = dependency('numpy')
numpy_inc = numpy_dep.get_variable('cflags').split()[0].replace('-I', '')

# Update git submodules if in source tree
run_command(['git', 'submodule', 'update'], check: false)



# F2PY helpers
f2py_inc = run_command(
  py, '-c',
  'import numpy.f2py, pathlib, sys; print((pathlib.Path(numpy.f2py.__file__).parent / "src"))',
  check:true).stdout().strip()

f2py_obj = configure_file(              # copies fortranobject.c into build dir
  input  : join_paths(f2py_inc, 'fortranobject.c'),
  output : 'fortranobject.c',
  copy   : true)

cc          = meson.get_compiler('c')
dl_dep      = cc.find_library('dl', required : host_machine.system() != 'windows')
thread_dep  = dependency('threads')

# Optional mold linker
if find_program('mold', required:false).found() and cc.get_id() == 'gcc'
 add_project_link_arguments('-fuse-ld=mold', language:['c','cpp','fortran'])
endif

if get_option('buildtype') == 'debug'
  add_global_arguments('-fbacktrace', '-fcheck=all', language : ['fortran']) #'-shared-libasan','-fsanitize=address',
endif

add_project_arguments('-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION', language:['c','cpp'])
add_project_arguments(
  '-Wno-uninitialized','-fno-second-underscore','-Wno-unused',
  '-Wno-conversion','-Wno-unused-dummy-argument','-Wno-character-truncation',
  '-Wno-integer-division','-Wno-surprising','-Wno-tabs','-cpp','-UNDEBUG',
  language:'fortran')
add_project_arguments('-Wno-unused-const-variable', language:'c')

# ── Helper vars ─────────────────────────────────────────────────────────────────────
rpath_tok = host_machine.system() == 'darwin' ? '@loader_path' : '$ORIGIN'

scripts_dir = meson.project_source_root() + '/scripts/'

src_root    = 'src/fortran'
common_inc  = [
  src_root,
  numpy_inc,
  f2py_inc,
]
common_cargs   = []
common_fargs   = ['-DCHROMO']

logging_src = src_root / 'logging.f'
rangen_src  = [src_root/'rangen.c', src_root/'rangen.f']
normal_src  = src_root / 'normal.c'
common_syms = ['chromo_openlogfile','chromo_closelogfile','npyrng']

# ── Symbol tables (unchanged) ───────────────────────────────────────────────────────
sibyll_syms_base  = [
  'sibyll','sibyll_ini','sib_sigma_hp','sib_sigma_hair','sib_sigma_hnuc',
  'int_nuc','decsib','decpar','sibini','sibhep','sib_list',
  'isib_pid2pdg','isib_pdg2pid','pdg_ini', 'sibnuc', 'sigma_nuc_nuc']
qgs_syms_base     = ['cqgsini','qgsect','qgini','qgconf','qgreg','chepevt','qgcrossc','cqgshh_ha_cs']
urqmd_syms        = [
  'urqmd','init','uinit','set0','params','uounit','strini','loginit','loadwtab',
  'norm_init','output','cascinit','nucrad','urqini','partname','chepevt','ptsigtot']
pythia6_syms      = [
  'pyinit','pyexec','pytune','pylist','pyevnt','pyevnw','pystat','pyedit',
  'pyhepc','pychge','pycomp','pyk']
sophia_syms       = ['eventgen','print_event','crossection','initial','icon_pdg_sib','toevt']
dpm_syms          = [
  'dt_evtout','dt_eventb','dt_ficonf','dt_flahad','dt_init','dt_initjs','dt_kkinc',
  'dt_phoxs','dt_sigga','dt_siggat','dt_siggp','dt_title','dt_xsglau','dt_xshn',
  'idt_icihad','ipho_id2pdg','ipho_pdg2id','pho_borncs','pho_event','pho_ghhias',
  'pho_harint','pho_harxto','pho_harxpt','pho_init','pho_mcini','pho_pmass',
  'pho_pname','pho_prevnt','pho_rregpar','pho_setmdl','pho_setpar','pho_setpdf',
  'pho_setpcomb','pho_sregpar','pho_xsect','poevt1','poevt2','pycomp','pytune']
epos_syms         = [
  'aaset','ainit','aepos','alist','afinal','hepmcstore','getcharge','idtrafo',
  'initepos','initeposevt','xsection','rmmard','ranfgt','ranfst']

# ── Helper: build one extension module ---------------------------------------------
# (No custom functions; just plain code within the loop.)

# ── Build catalogue ----------------------------------------------------------------
f = src_root  # shorthand in long lists

models = {
  '_sib21' : {
      'src'   : [f/'sibyll/sibyll_21.f', f/'sibyll/sibyll_init.f', f/'sibyll/sib21aux.f'],
      'funcs' : sibyll_syms_base + ['spgasdev','glauber'] + common_syms,
      'flags' : ['-DSIBYLLSP'],
  },
  '_sib23' : { 'src':[f/'sibyll/sibyll2.3.f',f/'sibyll/sibyll_init.f'],
               'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },
  '_sib23c00' : { 'src':[f/'sibyll/sibyll2.3c00.f',f/'sibyll/sibyll_init.f'],
                  'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },
  '_sib23c01' : { 'src':[f/'sibyll/sibyll2.3c01.f',f/'sibyll/sibyll_init.f'],
                  'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },
  '_sib23c02' : { 'src':[f/'sibyll/sibyll2.3c02.f',f/'sibyll/sibyll_init.f'],
                  'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },
  '_sib23c03' : { 'src':[f/'sibyll/sibyll2.3c03.f',f/'sibyll/sibyll_init.f'],
                  'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },
  '_sib23d'   : { 'src':[f/'sibyll/sibyll2.3d.f',f/'sibyll/sibyll_init.f'],
                  'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },
  '_sib23d_star':{ 'src':[f/'sibyll/sibyll2.3d-star-p02.f',f/'sibyll/sibyll_init.f'],
                   'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },
  '_sib23e'   : { 'src':[f/'sibyll/sibyll2.3e.f',f/'sibyll/sibyll_init.f'],
                  'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },
  '_sib23e_star':{ 'src':[f/'sibyll/sibyll2.3e-star.f',f/'sibyll/sibyll_init.f'],
                   'funcs':sibyll_syms_base+['gasdev','sig_had_nuc']+common_syms },

  '_qgs01' : { 'src':[f/'qgsjet/qgsjet01d.f',f/'qgsjet/chromo_qgs1.f'],
               'funcs':['cqgsini','sectnu','xxaini','psconf','xxreg','psaini',
                        'chepevt','xxfz','cqgshh_ha_cs','crossc']+common_syms },
  '_qgs2_03':{ 'src':[f/'qgsjet/qgsjet-II-03.f',f/'qgsjet/chromo_qgsII.f'],
               'funcs':qgs_syms_base+common_syms },
  '_qgs2_04':{ 'src':[f/'qgsjet/qgsjet-II-04.f',f/'qgsjet/chromo_qgsII.f'],
               'funcs':qgs_syms_base+common_syms },
  '_qgs3'   :{ 'src':[f/'qgsjet/qgsjet-III.f',f/'qgsjet/chromo_qgsIII.f'],
               'funcs':qgs_syms_base+['qgran']+common_syms },

  # URQMD ----------------------------------------------------------------
  '_urqmd34':{
    'src'   : run_command(py, scripts_dir + 'list_urqmd_sources.py',
                          '--directory', f/'urqmd-3.4/sources', check:true)
              .stdout().strip().split('\n'),
    'inc'   : [f/'urqmd-3.4/sources'],
    'iface'   : run_command(py, scripts_dir + 'list_urqmd_sources.py', '--interface-only',
                          '--directory', f/'urqmd-3.4/sources', check:true)
              .stdout().strip().split('\n'),
    'funcs' : urqmd_syms + common_syms,
  },

  '_pythia6':{
    'src'   : [f/'pythia6/pythia-6.4.28.f'],
    'funcs' : pythia6_syms + common_syms,
  },
  '_sophia':{
    'src'   : [f/'sophia/SOPHIA20.f',f/'sophia/eventgen.f',f/'sophia/sampling.f',
               f/'sophia/inpoutput.f',f/'sophia/jetset74dp.f',f/'sophia/chromo_sophia.f'],
    'funcs' : sophia_syms + common_syms,
  },

  # EPOS -----------------------------------------------------------------
  '_eposlhc':{
    'src'   : run_command(py, scripts_dir + 'list_eposlhc_sources.py', '--directory',
                          f/'epos-lhc', check:true).stdout().strip().split('\n'),
    'inc'   : [f/'epos-lhc/sources'],
    'iface' : run_command(py,'scripts/list_eposlhc_sources.py','--directory', f/'epos-lhc',
                      '--interface-only', check:true).stdout().strip().split('\n'),
    'funcs' : epos_syms + common_syms,
    'flags' : ['-O0'],
  },
  '_eposlhcr':{
    'src'   : run_command(py, scripts_dir + 'list_eposlhcr_sources.py', '--directory',
                          f/'epos-lhc-r', check:true).stdout().strip().split('\n'),
    'inc'   : [f/'epos-lhc-r/sources'],
    'iface' : run_command(py,'scripts/list_eposlhcr_sources.py','--directory', f/'epos-lhc-r',
                      '--interface-only', check:true).stdout().strip().split('\n'),
    'funcs' : epos_syms + common_syms,
    'flags' : ['-O0'],
  },
}

# DPMJET/PHOJET shared lib + split front-ends ----------------------------
dpm_cfg = {
 '307': {
   'src'  : [f/'dpmjet3.0/sources/dpmjet3.0-7.f',f/'dpmjet3.0/sources/phojet1.12-36c4.f',
             f/'dpmjet3.0/sources/pythia6115dpm3v1.f',f/'dpmjet3.0/sources/user3.0-7.f'],
   'inc'  : [],
   'split': ['_dpmjet307','_phojet307'],
 },
 '191': {
   'src'  : run_command(py, scripts_dir + 'list_dpmjet_sources.py',
                        '--directory', f/'dpmjetIII-19.1', check:true)
            .stdout().strip().split('\n'),
   'inc'  : [f/'dpmjetIII-19.1/include/phojet',f/'dpmjetIII-19.1/include/dpmjet',
             f/'dpmjetIII-19.1/include/pythia',f/'dpmjetIII-19.1/include/flinclude'],
   'split': ['_dpmjet191','_phojet191'],
 },
 '193': {
   'src'  : run_command(py, scripts_dir + 'list_dpmjet_sources.py',
                        '--directory', f/'dpmjetIII-19.3', check:true)
            .stdout().strip().split('\n'),
   'inc'  : [f/'dpmjetIII-19.3/include/phojet',f/'dpmjetIII-19.3/include/dpmjet',
             f/'dpmjetIII-19.3/include/pythia',f/'dpmjetIII-19.3/include/flinclude'],
   'split': ['_dpmjet193','_phojet193'],
 },
}

# ── selection lists ----------------------------------------------------
enabled_models  = run_command(py,scripts_dir + 'get_models.py','--enabled', check:true).stdout().strip().split('\n')
disabled_models = run_command(py,scripts_dir + 'get_models.py','--disabled',check:true).stdout().strip().split('\n')
message('Enabled models: \n\t' + '\n\t'.join(enabled_models))
message('Disabled models: \n\t' + '\n\t'.join(disabled_models))


# ── loop: regular models ----------------------------------------------
foreach name, cfg : models
  if name.substring(1) not in enabled_models
    message('Skipping disabled model: ' + name)
    continue
  endif

  srcs = files(cfg['src'] + [logging_src] + rangen_src + [normal_src])

  inc_dirs = common_inc
  if 'inc' in cfg
    foreach p : cfg['inc']
      inc_dirs += p
    endforeach
  endif

  # build -I flags string
  inc_flags = []
  foreach d : inc_dirs
    if not d.startswith(meson.project_source_root())
      d = meson.project_source_root() / d
    endif
    inc_flags += '-I' + d
  endforeach
  inc_str = ' '.join(inc_flags)

  flags = common_fargs + cfg.get('flags', [])

  wrap_inputs = files(cfg.get('iface', cfg['src']) + [logging_src] + rangen_src)
  outs  = [name+'module.c', name+'-f2pywrappers.f', name+'.pyf']

  wrap = custom_target(name+'_wrap',
    output : outs,
    input  : wrap_inputs,
    command: [py,scripts_dir + 'generate_f2py.py', name, ','.join(cfg['funcs']),
              meson.global_build_root()/'meson-logs', inc_str, ' '.join(flags),
              '@OUTDIR@','@INPUT@'],
    build_by_default : true)

  # Override optimization for EPOS models
  override_opts = []
  if name in ['_eposlhc', '_eposlhcr']
    override_opts = ['optimization=0']
  endif

  py.extension_module(name,
    sources            : srcs + [wrap, f2py_obj],
    include_directories: include_directories(inc_dirs),
    dependencies       : [py_dep, numpy_dep],
    c_args             : common_cargs + flags,
    fortran_args       : flags,
    override_options   : override_opts,
    subdir             : 'chromo/models',
    install            : true,
    install_rpath      : rpath_tok,
    build_rpath        : rpath_tok)
endforeach

# ── DPMJET/PHOJET ------------------------------------------------------
foreach ver, cfg : dpm_cfg
  model_group = 'dpmjet_phojet'+ver
  if model_group not in enabled_models
    message('Skipping disabled model group: ' + model_group)
    continue
  endif
  base_src = files(cfg['src'] + [logging_src] + rangen_src + [normal_src])

  inc_dirs = common_inc
  if 'inc' in cfg
    foreach p : cfg['inc']
      inc_dirs += p
    endforeach
  endif

  inc_flags = []
  foreach d : inc_dirs
    if not d.startswith(meson.project_source_root())
      d = meson.project_source_root() / d
    endif
    inc_flags += '-I' + d
  endforeach
  inc_str = ' '.join(inc_flags)

  lib = static_library('dpmjet_phojet_'+ver,
          base_src,
          include_directories : include_directories(inc_dirs),
          c_args              : common_cargs,
          fortran_args        : common_fargs)

  foreach front : cfg['split']
    outs = [front+'module.c', front+'-f2pywrappers.f', front+'.pyf']
    wrap = custom_target(front+'_wrap',
      output : outs,
      input  : base_src,
      command: [py,scripts_dir + 'generate_f2py.py', front, ','.join(dpm_syms+common_syms),
                meson.global_build_root()/'meson-logs', inc_str, ' '.join(common_fargs),
                '@OUTDIR@','@INPUT@'],
      build_by_default:true)

    py.extension_module(front,
      sources            : [wrap, f2py_obj],
      include_directories: inc_dirs,
      dependencies       : [py_dep, numpy_dep],
      link_whole         : [lib],
      c_args             : common_cargs,
      fortran_args       : common_fargs,
      subdir             : 'chromo/models',
      install            : true,
      install_rpath      : rpath_tok,
      build_rpath        : rpath_tok)
  endforeach
endforeach

# ── PYTHIA 8 block (unchanged) ----------------------------------------
cpp_dir    = 'src/cpp'
if 'pythia8' in enabled_models
  p8_inc     = include_directories(cpp_dir / 'pythia83/include')
  pybind_inc = include_directories(cpp_dir / 'pybind11/include')

  p8_src_dir = meson.project_source_root() / cpp_dir / 'pythia83/src'
  p8_srcs = run_command(
    'python', 'scripts/list_cc_files.py', p8_src_dir,
    check: true
  ).stdout().strip().split('\n')


  libpy8 = shared_library('pythia8',
            files(p8_srcs),
            include_directories : p8_inc,
            cpp_args            : ['-DXMLDIR=""'],
            dependencies        : [dl_dep,thread_dep],
            install             : true,
            install_dir         : py.get_install_dir()/'chromo'/'models',
            install_rpath       : rpath_tok,
            build_rpath         : rpath_tok)

  py.extension_module('_pythia8',
    cpp_dir/'_pythia8.cpp',
    include_directories : [pybind_inc, include_directories(cpp_dir), p8_inc],
    link_with           : libpy8,
    dependencies        : [py_dep],
    subdir              : 'chromo/models',
    install             : true,
    install_rpath       : rpath_tok,
    build_rpath         : rpath_tok)
endif

# ── Pure-Python files --------------------------------------------------
# Install the pure Python package files, but exclude iamdata
# Install all .py files in src/chromo
py.install_sources([
  'src/chromo/__init__.py',
  'src/chromo/__main__.py',
  'src/chromo/cli.py',
  'src/chromo/common.py',
  'src/chromo/constants.py',
  'src/chromo/decay_handler.py',
  'src/chromo/kinematics.py',
  'src/chromo/util.py',
  'src/chromo/writer.py',
], subdir: 'chromo')
# Install all subdirs except iamdata
install_subdir('src/chromo/models', install_dir: py.get_install_dir() / 'chromo')
#######################################################################################