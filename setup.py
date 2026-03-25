import sys
import os
import platform
import shutil
import subprocess
import tempfile
from setuptools import setup, Extension

long_description= ''
previous_line= ''
with open('README.md') as dfile:
    for line in dfile:
        if '[!' in line: continue
        long_description+= line

libraries= ['m']

def _check_openmp():
    """Return (available, compile_flags, link_flags) for OpenMP support."""
    cc = os.environ.get('CC', 'cc')
    src = '#include <omp.h>\nint main(void){return omp_get_num_threads();}\n'
    tmpdir = tempfile.mkdtemp()
    src_file = os.path.join(tmpdir, 'test_omp.c')
    exe_file = os.path.join(tmpdir, 'test_omp')
    try:
        with open(src_file, 'w') as f:
            f.write(src)
        # Candidate flag sets: try GCC-style first, then Clang/macOS style
        candidates = [
            (['-fopenmp'], ['-lgomp']),                  # GCC / Linux
            (['-Xpreprocessor', '-fopenmp'], ['-lomp']), # Apple clang / macOS
        ]
        for cflags, lflags in candidates:
            result = subprocess.run(
                [cc] + cflags + [src_file, '-o', exe_file] + lflags,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            if result.returncode == 0:
                return True, cflags, lflags
        return False, [], []
    except Exception:
        return False, [], []
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

#Option to forego OpenMP
try:
    openmp_pos = sys.argv.index('--no-openmp')
except ValueError:
    _omp_available, _omp_cflags, _omp_lflags = _check_openmp()
    if _omp_available:
        extra_compile_args = _omp_cflags
        libraries.extend([f[2:] for f in _omp_lflags if f.startswith('-l')])  # strip leading '-l'
    else:
        extra_compile_args = ["-DNO_OMP"]
else:
    del sys.argv[openmp_pos]
    extra_compile_args= ["-DNO_OMP"]

#Option to track coverage
try:
    coverage_pos = sys.argv.index('--coverage')
except ValueError:
    extra_link_args= []
else:
    del sys.argv[coverage_pos]
    extra_compile_args.extend(["-O0","--coverage"])
    extra_link_args= ["--coverage"]

setup(name='wendy',
      version='0.3.dev',
      description='One-dimensional gravitational N-body code',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Jo Bovy',
      author_email='bovy@astro.utoronto.ca',
      license='MIT',
      url='http://github.com/jobovy/wendy',
      package_dir = {'wendy': 'wendy'},
      packages=['wendy'],
      package_data={"": ["README.md","LICENSE"]},
      include_package_data=True,
      install_requires=['numpy','numba'],
      ext_modules=[Extension('wendy_c',
                             sources=['wendy/wendy.c','wendy/bst.c',
                                      'wendy/parallel_sort.c'],
                             libraries=libraries,
                             include_dirs=['wendy/'],
                             extra_compile_args=extra_compile_args,
                             extra_link_args=extra_link_args)])
