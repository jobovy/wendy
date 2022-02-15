import sys
from setuptools import setup
from distutils.core import Extension

long_description= ''
previous_line= ''
with open('README.md') as dfile:
    for line in dfile:
        if '[!' in line: continue
        long_description+= line

libraries= ['m']

#Option to forego OpenMP
try:
    openmp_pos = sys.argv.index('--no-openmp')
except ValueError:
    extra_compile_args = ["-fopenmp"]
    libraries.append('gomp')
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
      package_dir = {'wendy/': ''},
      packages=['wendy'],
      package_data={"": ["README.md","LICENSE"]},
      include_package_data=True,
      install_requires=['numpy!=1.21.*','numba'],
      ext_modules=[Extension('wendy_c',
                             sources=['wendy/wendy.c','wendy/bst.c',
                                      'wendy/parallel_sort.c'],
                             libraries=libraries,
                             include_dirs=['wendy/'],
                             extra_compile_args=extra_compile_args,
                             extra_link_args=extra_link_args)])
