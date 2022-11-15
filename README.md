# BIBIMBpy
 Backwards Integration Basic Interface Module for Bars

 This package is a simple interface to simplify the act of running backward integrations and thus leaving more time for planning and exploring the parameter space.

## Documentation

All classes and methods/functions are documented so use the python help() function to find out more.


## Installation

This is a Python3 package (*issues may arise if executed with Python2*).

The required dependencies are:
* [agama](https://github.com/GalacticDynamics-Oxford/Agama)
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)


To install it:
1. Clone the Github repository or download the source file
2. cd to the directory
3. Run ```python setup.py install```


## Basic usage

In the folder *TUTORIALS*, you will find examples on how to use the different functions provided.

Also, in the folder *DOCS* there is an schema that depicts the normal flow and intended usages. 

In a few words, the main inputs are:
1. The parameters of the background gravitational potential pre-perturbation
2. The parameters of the final perturbation
3. A description of the mass evolution of the perturbation
4. The pattern speed of the perturbation
5. The grid used to start the simulation (*for now, limited to two variables. The other 4 must be fixed*)

Once all this inputs have been defined, use the functions in the submodule _initialize_ to prepare the variables required by the main module _orbit.runBI_.

The result of _orbit.runBI_ is an array containing the value of the Distribution Function for each particle, at a time prior to the perturbation. It also returns the array of time and the orbits of each particle (in the rotating frame, use _utils.rotating2inertial_ to convert back to the intertial frame).


## Attribution

If you make use of this package for your research, please acknowledge the authors in your work.