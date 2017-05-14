# wendy

A one-dimensional gravitational N-body code. 

## Overview

``wendy`` solves the one-dimensional gravitational N-body problem to machine precision with an efficient algorithm [O(log N) / particle-collision].

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
next_x, next_v= g.next() # at t=0.05
next_x, next_v= g.next() # at t=0.10
...
```
A pure Python version of the code is available as ``wendy.nbody_python``.

## Examples

* Phase mixing and violent relaxation in one dimension: [example notebook](examples/PhaseMixingViolentRelaxation.ipynb) (run locally to see movies, or view on [nbviewer](http://nbviewer.jupyter.org/github/jobovy/wendy/blob/master/examples/PhaseMixingViolentRelaxation.ipynb?flush_cache=true))

<img src="https://cloud.githubusercontent.com/assets/1044876/26021056/fc5bf4d4-3754-11e7-9b10-98828c0b9298.gif" width="400"><img src="https://cloud.githubusercontent.com/assets/1044876/26020997/267a5ee6-3754-11e7-94f9-4822a60f5f4d.gif" width="400">

