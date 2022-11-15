import numpy as np
import agama
from .utils import rotating2inertial

def integrate_backwards(particle_ini,pot,number_of_cycles,dyn_time_base,pattern_speed,_trajsize = 2):
    """ 
    particle_ini: initial conditions Nx6
    pot: agama potential
    number_of_cycles: how many dynamical times to integrate for
    dyn_time_base: dynamical time of reference [s*kpc/km]
    pattern_speed: integrate in a rotating reference frame of constant pattern speed [km/s/kpc]
    _trajsize: number of snapshots
    """
    
    trajectories = agama.orbit(ic = particle_ini,potential = pot,timestart = dyn_time_base*number_of_cycles,
                                time=-dyn_time_base*number_of_cycles,lyapunov=False,trajsize=_trajsize,Omega=pattern_speed)
    
    return trajectories

def runBI(particle_ini,pot_timedep,df_gen_func,dyn_time_base,number_of_cycles,pattern_speed,_Npoints = 2):
    """
    Run a backwards integration run in a time dependent potential.

    Inputs:
    - particle_ini: Initial conditions, Nx6 array of positions and velocites at t=tfin.
    - pot_timedep: AGAMA potential instance use for the backwards integration
    - df_gen_func: Arbritrary function that, given a vector Nx6 of positions and velocities in an intertial frame, returns the value of the Distribution Function for each of the particles in pot_base at t<=0, before the perturbation.
    - dyn_time_base: Dynamical timescale. It is meant to be used as the period of the pattern speed, i.e. 2pi/pattern_speed [in units of s*kpc/km]
    - number_of_cycles: Number of dynamical timescales to integrate.
    - pattern_speed: pattern speed of the rotating frame, positive for a prograre rotation [in units of km/s/kpc]
    -_Npoints: Number of snapshots. In case you want to see the whole orbit, you can increase the number of sample taken along the whole orbit.

    Outputs:
    - df_eval: value of the Distribution Function for each particle in particle_ini. 
    - time: time stamps of the snapshots obtained.
    - orbits: position and velocities of all particles in the rotating frame across time. 
    """
    
    #integration
    trajectories = integrate_backwards(particle_ini,pot_timedep,number_of_cycles,dyn_time_base,pattern_speed,_trajsize=_Npoints)
    
    #put particles in inertial frame
    time = trajectories[0,0]
    orbits = trajectories[:,1]
    particle_fin_barframe = np.stack([o[-1] for o in orbits])
    o_inertial = rotating2inertial(time[-1],particle_fin_barframe,pattern_speed)
    
    #evaluate DF
    df_eval = df_gen_func(o_inertial)

    return df_eval,time,orbits