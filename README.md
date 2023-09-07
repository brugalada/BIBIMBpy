# BIBIMBpy
Backwards Integration Basic Interface Module for Bars

This package is a simple interface to simplify the act of running backward integrations and thus leaving more time for planning and exploring the parameter space.

Backward Integrations are a simple way to generate high resolution renderings of phase-space at an arbitrary time, and for an arbitrary time-evolving potential. The main limitation of this method is that there is no self-gravity envolved so, in a way, are equivalent to test particle simulations. 

The physical principle in which it is based is the following: "_The collisionless Boltzmann equation tells us that the Distribution Function (DF) remains constant along stellar trajectories. Thus, the value f(w,t2) of the DF at some phase-space point w = (x, v) at time t2 is equal to f(w0,0) if w originates from integrating w0 from t = 0 to t = t2._" (extract from Dehnen 2000)

## Documentation

All classes and methods/functions are documented so use the python help() function to find out more.


## Installation

This is a Python3 package (*issues may arise if executed with Python2*).

The required dependencies are:
* [agama](https://github.com/GalacticDynamics-Oxford/Agama)
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [matplotlib](https://matplotlib.org)


To install it:
1. Clone the Github repository or download the source file
2. cd to the directory
3. Run ```python setup.py install```

Make sure to run the latest versions, specially for AGAMA.

## Basic usage

In the folder *TESTS*, you will find examples on how to use the different functions provided.

Also, in the folder *DOCS* there is an schema that depicts the normal flow and intended usages. 

In a few words, the main inputs are:
1. The parameters of the background gravitational potential pre-perturbation
2. The parameters of the final perturbation
3. A description of the mass evolution of the perturbation
4. The pattern speed of the perturbation
5. The grid used to start the simulation (*for now, limited to two variables. The other 4 must be fixed*)

Once all this inputs have been defined, use the functions in the submodule _initialize_ to prepare the variables required by the main module _orbit.runBI_.

The result of _orbit.runBI_ is an array containing the value of the Distribution Function for each particle, at a time prior to the perturbation. It also returns the array of time and the orbits of each particle (in the rotating frame, use _utils.rotating2inertial_ to convert back to the intertial frame).

## Conventions

By default, stars with positive angular momentum will rotate counterclock-wise. In the same manner, a perturbation with a positive pattern speed will also rotate counterclock-wise. This means that, in the backward integration were time decreases, particles will _appear_ to rotate clock-wise, but just because time is moving backwards. 

Also by default, the bar will always start with its long axis aligned with the x-axis (i.e. flat). This means that a particle that starts at positive azimuth is _in front_ of the bar, in the sense that the bar is moving towards it instantaneously. In other words, the Sun would have a negative $\phi_{bar}$ in this convention.


## Attribution

Feel free to use this software and, if you run into any issue, do not hesitate to contact me. 


## Image source
The image used for the logo is not mine. Credit goes to user _su-lin_ @[flickr](https://www.flickr.com/photos/su-lin/4429927455)
