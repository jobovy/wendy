import sys
from setuptools import setup
from distutils.core import Extension
    
#Option to track coverage
try:
    coverage_pos = sys.argv.index('--coverage')
except ValueError:
    extra_compile_args= []
    extra_link_args= []
else:
    del sys.argv[coverage_pos]
    extra_compile_args= ["-O0","--coverage"]
    extra_link_args= ["--coverage"]

setup(name='wendy',
      version='0.1',
      description='One-dimensional gravitational N-body code',
      author='Jo Bovy',
      author_email='bovy@astro.utoronto.ca',
      license='MIT',
      url='http://github.com/jobovy/wendy',
      package_dir = {'wendy/': ''},
      packages=['wendy'],
      package_data={"": ["README.md","LICENSE"]},
      include_package_data=True,
      install_requires=['numpy>=1.7'],
      ext_modules=[Extension('wendy_c',
                             sources=['wendy/wendy.c','wendy/bst.c'],
                             libraries=['m'],
                             include_dirs=['wendy/'],
                             extra_compile_args=extra_compile_args,
                             extra_link_args=extra_link_args)])
