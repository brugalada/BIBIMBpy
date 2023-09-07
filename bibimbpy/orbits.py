import numpy as np
from .utils import rotating2inertial
from .__init__ import agama

def integrate_backwards(particle_ini,pot,t_start,total_time,pattern_speed,_trajsize = 2):
    """ 
    particle_ini: initial conditions Nx6
    pot: agama potential
    t_start: starting time of the integration
    total_time: amount of time to integrate for (positive value for backward integration, negative value for forward)
    pattern_speed: integrate in a rotating reference frame of constant pattern speed [km/s/kpc]
    _trajsize: number of snapshots
    """
    
    trajectories = agama.orbit(ic = particle_ini,potential = pot,timestart = t_start,
                                time=-total_time,lyapunov=False,trajsize=_trajsize,Omega=pattern_speed)
    
    return trajectories

def runBI(particle_ini,pot_timedep,df_gen_func,t_start,total_time,pattern_speed,_Npoints = 2):
    """
    Run a backwards integration run in a time dependent potential.

    Inputs:
    - particle_ini: Initial conditions, Nx6 array of positions and velocites at t=tfin.
    - pot_timedep: AGAMA potential instance use for the backwards integration
    - df_gen_func: Arbritrary function that, given a vector Nx6 of positions and velocities in an intertial frame, returns the value of the Distribution Function for each of the particles in pot_base at t<=0, before the perturbation.
    - t_start: time when the simulation starts
    - total_time: integration time (positive for backwards, negative for foward)
    - pattern_speed: pattern speed of the rotating frame, positive for a prograre rotation [in units of km/s/kpc]
    -_Npoints: Number of snapshots. In case you want to see the whole orbit, you can increase the number of sample taken along the whole orbit.

    Outputs:
    - df_eval: value of the Distribution Function for each particle in particle_ini. 
    - time: time stamps of the snapshots obtained.
    - orbits: position and velocities of all particles in the rotating frame across time. 
    """
    
    #integration
    trajectories = integrate_backwards(particle_ini,pot_timedep,t_start,total_time,pattern_speed,_trajsize=_Npoints)
    
    #put particles in inertial frame
    time = trajectories[0,0]
    orbits = trajectories[:,1]
    particle_fin_barframe = np.stack([o[-1] for o in orbits])
    o_inertial = rotating2inertial(time[-1],particle_fin_barframe,pattern_speed,_t0 = t_start)
    
    #evaluate DF
    df_eval = df_gen_func(o_inertial)

    return df_eval,time,orbits