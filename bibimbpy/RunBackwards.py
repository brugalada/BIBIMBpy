import os
import sys
import numpy as np
import agama
import matplotlib.pyplot as plt
import matplotlib.colors as mcols
import scipy.stats
import pandas as pd
from astropy.coordinates import *
import argparse






def read_inputs():
    parser = argparse.ArgumentParser(description='Backward integration in custom potential using AGAMA.')
    
    parser.add_argument('directory', type=str, help='Folder to store the results')
    parser.add_argument('filename', type=str, help='Address to the potential (.ini) file')
    parser.add_argument('filename_base', type=str, help='Address to the potential (.ini) file')
    parser.add_argument('r0', type=float, help='Galactocentric radius')
    parser.add_argument('phi0', type=float,help='Cylindrical azimuth (w.r.t. the bar, in degrees)')
    parser.add_argument('z0', type=float, help='Height above/below the disc')
    parser.add_argument('vz0', type=float, help='Vertical velocity')
    parser.add_argument('vrmin', type=float, help='Minimum Vr on the grid')
    parser.add_argument('vrmax', type=float, help='Maximum Vr on the grid')
    parser.add_argument('vphimin', type=float, help='Minimum Vphi on the grid: should be negative for prograde orbits')
    parser.add_argument('vphimax', type=float, help='Maximum Vphi on the grid')
    parser.add_argument('Nvr', type=int, help='Number of points in the Vr axis of the grid')
    parser.add_argument('Nvphi', type=int, help='Number of points in the Vphi axis of the grid')
    parser.add_argument('steps_per_dt', type=int, help='Number of snapshots per integration')
    parser.add_argument('number_of_cycles', type=int, help='Number of dynamical times to integrate')
    parser.add_argument('dyn_time_base', type=float, help='Reference dynamical time. Units of second*kpc/km (e.g. 2*pi/pattern_speed).')
    parser.add_argument('pattern_speed', type=float, help='Speed in km/s/kpc of the reference frame. Should be negative for a prograde rotation.')

    args = parser.parse_args()
    
    return args

def read_inputs_file(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            (key, val) = line.split()
            try:
                d[str(key)] = float(val)
            except:
                d[str(key)] = val
    return d


def set_initial_conditions(r0,phi0,z0,vz0,vrmin,vrmax,vphimin,vphimax,Nvr,Nvphi):
    """
    Distances in kpc, angles in degrees and velocities in km/s
    Creates a Nx6 array of initial conditions at a fix location, for a grid of
    Vrs and Vphis
    """
    r0 = r0
    phi0 = np.deg2rad(phi0)
    z0 = z0
    vz0 = vz0

    vr0_ = np.linspace(vrmin,vrmax,Nvr)
    vphi0_ = np.linspace(vphimin,vphimax,Nvphi)

    VR,VPHI = np.meshgrid(vr0_,vphi0_)

    vr0 = VR.flatten()
    vphi0 = VPHI.flatten()

    x0 = r0*np.cos(phi0)*np.ones_like(vr0)
    y0 = r0*np.sin(phi0)*np.ones_like(vr0)
    z0 = z0*np.ones_like(vr0)
    vx0 = vr0*np.cos(phi0)-vphi0*np.sin(phi0)
    vy0 = vr0*np.sin(phi0)+vphi0*np.cos(phi0)
    
    #### CHECK SIGNS AND SHIT! ####
    vz0 = vz0*np.ones_like(vr0)

    return np.column_stack((x0,y0,z0,vx0,vy0,vz0)),vr0,vphi0,vr0_,vphi0_


def integrate_backwards(particle_ini,pot,steps_per_dt,number_of_cycles,dyn_time_base,pattern_speed):
    """ 
    particle_ini: initial conditions Nx6
    pot: agama potential
    steps_per_dt: number of snapshots
    number_of_cycles: how many dynamical times to integrate for
    dyn_time_base: dynamical time of reference
    pattern_speed: integrate in a rotating reference frame of constant pattern speed [km/s/kp]
    """
    trajsize = steps_per_dt*number_of_cycles+1
    
    trajectories = agama.orbit(ic = particle_ini,potential = pot,timestart = dyn_time_base*number_of_cycles,
                                time=-dyn_time_base*number_of_cycles,lyapunov=False,trajsize=trajsize,Omega=pattern_speed)
    
    return trajectories


def inertial2bar_frame(t, o,om):
    return np.column_stack((
    o[:,0] * np.cos(om*t) + o[:,1] * np.sin(om*t),
    o[:,1] * np.cos(om*t) - o[:,0] * np.sin(om*t),
    o[:,2],
    o[:,3] * np.cos(om*t) + o[:,4] * np.sin(om*t),
    o[:,4] * np.cos(om*t) - o[:,3] * np.sin(om*t),
    o[:,5] ))

def bar2inertial_frame(t,o,om):
    return np.column_stack((
    o[:,0] * np.cos(-om*t) + o[:,1] * np.sin(-om*t),
    o[:,1] * np.cos(-om*t) - o[:,0] * np.sin(-om*t),
    o[:,2],
    o[:,3] * np.cos(-om*t) + o[:,4] * np.sin(-om*t),
    o[:,4] * np.cos(-om*t) - o[:,3] * np.sin(-om*t),
    o[:,5] ))

def read_potential(filename):
    return agama.Potential(filename)

def obtain_df(pbase):
    return agama.DistributionFunction(type="QuasiIsothermal", Sigma0 = 1., Rdisk = 2.5, Hdisk = 0.3, 
                                    Rsigmar = 5., sigmar0 = 35./np.exp(-8./5.), potential = pbase)

def plot_result(df_sum,vre,vphie,vrmin,vrmax,vphimin,vphimax,t_,r0,phi0,figname):
    
    VR,VPHI = np.meshgrid(0.5*(vre[1:]+vre[:-1]),0.5*(vphie[1:]+vphie[:-1]))
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.gca()
    cc = ax.pcolormesh(VR,VPHI,df_sum.T,cmap="jet");
    plt.colorbar(cc)
    
    ax.set_xlim(vrmin,vrmax)
    ax.set_ylim(vphimin,vphimax)
    ax.set_xlabel(r"$V_R$ [km$\,$s$^{-1}$]")
    ax.set_ylabel(r"$V_\phi$ [km$\,$s$^{-1}$]")
    ax.set_title(f"Time = {t_} Gyr // R = {r0} kpc // Phi = {phi0} degrees")
    
    fig.savefig(figname)

def store_results(filename,x,y,h):
    np.save(filename+"vr_edges",x)
    np.save(filename+"vphi_edges",y)
    np.save(filename+"dfsum",h)

def main(inputs,_source):

    agama.setUnits(mass=1, length=1, velocity=1)
    
    if _source == "command line":
        r_ic             = inputs.r
        phi_ic           = inputs.phi
        z_ic             = inputs.z
        vz_ic            = inputs.vz
        vrmin            = inputs.vrmin
        vrmax            = inputs.vrmax
        vphimin          = inputs.vphimin
        vphimax          = inputs.vphimax
        Nvr              = inputs.Nvr
        Nvphi            = inputs.Nvphi
        directory         = inputs.directory
        filename         = inputs.filename
        filename_base    = inputs.filename_base
        steps_per_dt     = inputs.steps_per_dt
        number_of_cycles = inputs.number_of_cycles
        dyn_time_base    = inputs.dyn_time_base
        pattern_speed    = inputs.pattern_speed
    elif _source == "file":
        r_ic             = inputs["r"]
        phi_ic           = inputs["phi"]
        z_ic             = inputs["z"]
        vz_ic            = inputs["vz"]
        vrmin            = inputs["vrmin"]
        vrmax            = inputs["vrmax"]
        vphimin          = inputs["vphimin"]
        vphimax          = inputs["vphimax"]
        Nvr              = int(inputs["Nvr"])
        Nvphi            = int(inputs["Nvphi"])
        directory        = inputs["directory"]
        filename         = inputs["filename"]
        filename_base    = inputs["filename_base"]
        steps_per_dt     = int(inputs["steps_per_dt"])
        number_of_cycles = int(inputs["number_of_cycles"])
        dyn_time_base    = inputs["dyn_time_base"]
        pattern_speed    = inputs["pattern_speed"]
    else:
        raise ValueError("Requested source is not an option. Choose 'file' or 'command line'")
    
    #potentials
    pot = read_potential(filename)
    pbase = read_potential(filename_base)
    
    #DF
    df = obtain_df(pbase)

    #action-angle
    af = agama.ActionFinder(pbase)
    
    #intial conditions
    particle_ini,vr0,vphi0,vr0_,vphi0_ = set_initial_conditions(r_ic,phi_ic,z_ic,vz_ic,
                                                vrmin,vrmax,vphimin,vphimax,Nvr,Nvphi)
    
    #integration
    trajectories = integrate_backwards(particle_ini,pot,steps_per_dt,number_of_cycles,dyn_time_base,pattern_speed)
    
    
    #put particles in inertial frame
    if Nvr*Nvphi > 1:
        time_ = trajectories[0,0]
        orbits = trajectories[:,1]
        particle_fin_barframe = np.stack([o[-1] for o in orbits])
        o_inertial = bar2inertial_frame(time_[-1],particle_fin_barframe,pattern_speed)
        
    else:
        
        time_ = trajectories[0]
        orbit = trajectories[1]
        o_inertial = bar2inertial_frame(time_,trajectories[1],pattern_speed)[-1]
        
    #compute actions at t=0
    actions = af(o_inertial)
    
    #evaluate DF
    df_eval = df(actions)
    
    #create velocity histrogram
    df_sum,vre,vphie,_ = scipy.stats.binned_statistic_2d(vr0,vphi0,df_eval,statistic="sum",bins=(vr0_,vphi0_));
    
    savename = directory + f"BackwardsIntegration_R{r_ic}_Phi{phi_ic}_Time{np.round(number_of_cycles*dyn_time_base,2)}_Omega{pattern_speed}"
    #plot it
    figname = savename + ".png"
    plot_result(df_sum,vre,vphie,vrmin,vrmax,vphimin,vphimax,
                np.round(number_of_cycles*dyn_time_base,2),np.round(r_ic,2),np.round(phi_ic,2),figname)
    
    store_results(savename,vre,vphie,df_sum)
    
if __name__ == "__main__":
    inputs = read_inputs()
    main(inputs,"command line")
