project('chromo', 'c', 'cpp', 'fortran', version: '0.6.0')

py = import('python').find_installation(pure: false)

# Install Python sources
py_files_str = run_command(py, ['-c', 'import pathlib,sys;print(";".join(str(p) for p in pathlib.Path("src/chromo").rglob("*.py")))'], check: true).stdout().strip()
py_files = py_files_str.split(';')
if py_files_str == ''
  py_files = []
endif
py.install_sources(py_files, subdir: 'chromo')

# Install default models configuration
install_data('default_models.cfg', install_dir: py.get_install_dir(subdir: 'chromo'))

# Build _sib23d extension via f2py
sib_sources = [
  'src/fortran/logging.f',
  'src/fortran/rangen.c',
  'src/fortran/rangen.fpp',
  'src/fortran/normal.c',
  'src/fortran/sibyll/sibyll2.3d.f',
  'src/fortran/sibyll/sibyll_init.fpp',
]

output_name = '_sib23d.so'

sib_ext = custom_target('_sib23d',
  input: sib_sources,
  output: output_name,
  command: [py, files('build_f2py.py'), '@OUTPUT@', '@INPUT@'],
  install: true,
  install_dir: py.get_install_dir(pure: false, subdir: 'chromo/models'),
)


