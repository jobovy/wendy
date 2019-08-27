# wendy

A one-dimensional gravitational N-body code. 

[![Build Status](https://travis-ci.org/jobovy/wendy.svg?branch=master)](https://travis-ci.org/jobovy/wendy)
[![Coverage Status](https://codecov.io/gh/jobovy/wendy/branch/master/graph/badge.svg)](https://codecov.io/gh/jobovy/wendy)
[![Binder](http://mybinder.org/badge.svg)](http://beta.mybinder.org/repo/jobovy/wendy)

## Overview

``wendy`` solves the one-dimensional gravitational N-body problem to machine precision with an efficient algorithm [O(log N) / particle-collision]. Alternatively, it can solve the problem with approximate integration, but with exact forces.

## Author

Jo Bovy (University of Toronto): bovy - at - astro - dot - utoronto - dot - ca

## Installation

Install the latest release using
```
pip install wendy
```
or clone/fork/download the repository and install using
```
sudo python setup.py install
```
or locally using
```
python setup.py install --user
```

The behavior of the parallel sorting algorithm used when setting ``sort='parallel'`` in the approximate version of the N-body code (``approx=True``) is controlled by a few compilation-time variables: ``PARALLEL_SERIAL_SORT_SWITCH``, which sets the length of an array below which the serial sort is used, ``PARALLEL_SERIAL_MERGE_SWITCH``, which sets the length of an array below which a serial merge is used (as part of the ``mergesort`` sorting algorithm used), and ``PARALLEL_SORT_NUM_THREADS``, the number of threads used in the parallel sorting algorithm. By default, these are set to ``PARALLEL_SERIAL_MERGE_SWITCH=10000``, ``PARALLEL_SERIAL_MERGE_SWITCH=50000``, and ``PARALLEL_SORT_NUM_THREADS=32``, which appear to work well. Significant speed-ups can be obtained by optimizing these for your system and specific problem. They can be set to different values by running, e.g.,
```
export CFLAGS="$CFLAGS -D PARALLEL_SERIAL_SORT_SWITCH=10 -D PARALLEL_SERIAL_MERGE_SWITCH=10 -D PARALLEL_SORT_NUM_THREADS=2"
```
before compiling the code (if you are trying to change them, make sure to force a re-compilation by removing the ``build/`` directory).

## Usage

Use ``wendy.nbody`` to initialize a generator object for initial *(x,v)* with masses *m*. The generator then returns the state of the system at equally-spaced time intervals:
```
g= wendy.nbody(x,v,m,0.05) # delta t = 0.05
next_x, next_v= next(g) # at t=0.05
next_x, next_v= next(g) # at t=0.10
...
```
The generator initialization with ``wendy.nbody`` has options to (a) solve the problem exactly or not using ``approx=``, (b) include an external harmonic oscillator potential ``omega^2 x^2 / 2`` with ``omega=`` (both for exact and approximate solutions), and (c) include an arbitrary external force ``F(x,t)`` (using ``ext_force=``).

## Examples

You can run these *without* installing ``wendy`` by clicking on [![Binder](http://mybinder.org/badge.svg)](http://beta.mybinder.org/repo/jobovy/wendy) and navigating to the ``examples/`` directory. Note that some of the movies might fail to be rendered on the binder webpage, so you might want to skip those when running the notebooks (or changing the ``subsamp`` input for them).

* Phase mixing and violent relaxation in one dimension: [example notebook](examples/PhaseMixingViolentRelaxation.ipynb) (run locally to see movies, or view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/PhaseMixingViolentRelaxation.ipynb?flush_cache=true))

<img src="https://cloud.githubusercontent.com/assets/1044876/26030657/e29c9efe-3826-11e7-8419-7bf96d565569.gif" width="400"><img src="https://cloud.githubusercontent.com/assets/1044876/26030672/1fafa9bc-3827-11e7-9167-16f10bb40b59.gif" width="400">

* A self-gravitating, sech<sup>2</sup> disk: [example notebook](examples/SelfGravitatingSech2Disk.ipynb) (run locally to see movies, or view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/SelfGravitatingSech2Disk.ipynb?flush_cache=true))

<img src="https://user-images.githubusercontent.com/1044876/26942002-7dce9644-4c4e-11e7-90d0-43beddfd0dd9.gif" width="400"><img src="https://cloud.githubusercontent.com/assets/1044876/26089414/850dd7a2-39cb-11e7-9beb-73dbce2c6c4b.gif" width="400">

* Adiabatic contraction: [example notebook](examples/AdiabaticContraction.ipynb) (run locally to see movies, or view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/AdiabaticContraction.ipynb?flush_cache=true))

<img src="https://user-images.githubusercontent.com/1044876/26941638-2b4eed8e-4c4d-11e7-9804-a63b681b86e6.gif" width="400"><img src="https://user-images.githubusercontent.com/1044876/26941809-d7393b22-4c4d-11e7-97dd-cd3b259aefa7.gif" width="400">

* Adiabatic vs. non-adiabatic energy injection for an exponential disk: [example notebook](examples/AdiabaticVsNonAdiabatic.ipynb) (run locally to see movies, or view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/AdiabaticVsNonAdiabatic.ipynb?flush_cache=true))

<img src="https://user-images.githubusercontent.com/1044876/27014815-3d5fe2de-4ece-11e7-953f-ced11d4993d2.gif" width="400"><img src="https://user-images.githubusercontent.com/1044876/27014825-648498f0-4ece-11e7-8360-cf1f04e13cba.gif" width="400">

* ``wendy`` scaling with particle number: [example notebook](examples/WendyScaling.ipynb) (view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/WendyScaling.ipynb?flush_cache=true))


