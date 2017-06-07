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

Clone/fork/download the repository and install using
```
sudo python setup.py install
```
or locally using
```
python setup.py install --user
```

## Usage

Use ``wendy.nbody`` to initialize a generator object for initial *(x,v)* with masses *m*. The generator then returns the state of the system at equally-spaced time intervals:
```
g= wendy.nbody(x,v,m,0.05) # delta t = 0.05
next_x, next_v= next(g) # at t=0.05
next_x, next_v= next(g) # at t=0.10
...
```
The generator initialization with ``wendy.nbody`` has options to (a) solve the problem exactly or not using ``approx=`` and (b) include an external harmonic oscillator potential ``omega^2 x^2 / 2`` with ``omega=`` (both for exact and approximate solutions).

## Examples

* Phase mixing and violent relaxation in one dimension: [example notebook](examples/PhaseMixingViolentRelaxation.ipynb) (run locally to see movies, or view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/PhaseMixingViolentRelaxation.ipynb?flush_cache=true))

<img src="https://cloud.githubusercontent.com/assets/1044876/26030657/e29c9efe-3826-11e7-8419-7bf96d565569.gif" width="400"><img src="https://cloud.githubusercontent.com/assets/1044876/26030672/1fafa9bc-3827-11e7-9167-16f10bb40b59.gif" width="400">

* A self-gravitating, sech<sup>2</sup> disk: [example notebook](examples/SelfGravitatingSech2Disk.ipynb) (run locally to see movies, or view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/SelfGravitatingSech2Disk.ipynb?flush_cache=true))

<img src="https://cloud.githubusercontent.com/assets/1044876/26089385/4b4f491a-39cb-11e7-9709-278162abe5b7.gif" width="400"><img src="https://cloud.githubusercontent.com/assets/1044876/26089414/850dd7a2-39cb-11e7-9beb-73dbce2c6c4b.gif" width="400">

* ``wendy`` scaling with particle number: [example notebook](examples/WendyScaling.ipynb) (view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/WendyScaling.ipynb?flush_cache=true))


