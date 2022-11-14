import numpy as np
import os
import agama
from .utils import write_potential, invert_scaling_file

def set_initial_conditions(r,phi,z,vr,vphi,vz):
    """
    Generate an array of initial conditions ready to be feed to the orbit integrator.
    
    NOTES: FOR NOW, ONLY TWO VARIABLES ARE ITERABLE. So, 4 variables have to be numbers and the other two, arrays of arbitrary lenghts N1 and N2. The resulting array will have N1*N2 particles.

    Inputs:
    - r: Galactocentric radius
    - phi: Galactocentric azimuth (in degrees)
    - z: Verticle height
    - vr: Galactocentric radial velocity (positive: outwards, negative: inwards)
    - vphi: Galactocentric azimuthal velocity (negative: prograde, positive: retrograde)
    - vz: Verticle velocity

    Output:
    - particles: N1*N2x6 array of particles
    - var1: array of all used values for the first iterative variable, useful for the final histogram (Galactocentric coordiantes)
    - var2: array of all used values for the second iterative variable, useful for the final histogram (Galactocentric coordiantes)
    """
    #check which variables are iterables
    iter_vars = {}
    final_vars = {}

    #store lengths
    lenghts = []

    #TO-DO: check a better way to do this check
    if isinstance(r,list) or isinstance(r,tuple) or isinstance(r,np.ndarray):
        iter_vars["r"] = r
        lenghts.append(len(r))
    else:
        final_vars["r"] = r

    if isinstance(phi,list) or isinstance(phi,tuple) or isinstance(phi,np.ndarray):
        iter_vars["phi"] = phi
        lenghts.append(len(phi))
    else:
        final_vars["phi"] = phi

    if isinstance(z,list) or isinstance(z,tuple) or isinstance(z,np.ndarray):
        iter_vars["z"] = z
        lenghts.append(len(z))
    else:
        final_vars["z"] = z

    if isinstance(vr,list) or isinstance(vr,tuple) or isinstance(vr,np.ndarray):
        iter_vars["vr"] = vr
        lenghts.append(len(vr))
    else:
        final_vars["vr"] = vr

    if isinstance(vphi,list) or isinstance(vphi,tuple) or isinstance(vphi,np.ndarray):
        iter_vars["vphi"] = vphi
        lenghts.append(len(vphi))
    else:
        final_vars["vphi"] = vphi

    if isinstance(vz,list) or isinstance(vz,tuple) or isinstance(vz,np.ndarray):
        iter_vars["vz"] = vz
        lenghts.append(len(vz))
    else:
        final_vars["vz"] = vz

    #check that all lengths are equal
    ##TO-DO check that all values are the same

    #Cross variables with a meshgrid (accounting for all the cases)
    if len(iter_vars.keys())==2:
        key1,key2 = list(iter_vars.keys())
        aux1,aux2 = list(iter_vars.values())
        var1,var2 = np.meshgrid(aux1, aux2, indexing="ij")
        final_vars[key1] = var1.flatten()
        final_vars[key2] = var2.flatten()
        len_ = len(final_vars[key1])
    else:
        #TO-DO: allow for more variable to iterate through
        raise ValueError("Sorry! More than two variables is not implemented yet!")

    #Generate initial conditions
    x0  = final_vars["r"]*np.cos(final_vars["phi"])*np.ones(len_)
    y0  = final_vars["r"]*np.sin(final_vars["phi"])*np.ones(len_)
    z0  = final_vars["z"]*np.ones(len_)
    vx0 = (final_vars["vr"]*np.cos(final_vars["phi"])-final_vars["vphi"]*np.sin(final_vars["phi"]))*np.ones(len_)
    vy0 = (final_vars["vr"]*np.sin(final_vars["phi"])+final_vars["vphi"]*np.cos(final_vars["phi"]))*np.ones(len_)
    vz0 = final_vars["vz"]*np.ones(len_)

    return np.column_stack((x0,y0,z0,vx0,vy0,vz0))


def generate_TimeDepPot_old(folder_name,file_name,generating_function,times,interpol="false"):
    """
    interpol = true or false
    """

    #Generate individual files for each step of the perturbation
    for i,t in times:
        pot_params_dict = generating_function(t)
        write_potential(pot_params_dict,folder_name+file_name+f"_t{t}.ini")

    #make sure that the files are sorted by time
    bar_files = [f_ for f_ in os.listdir(folder_name) if f_.endswith(".ini")]
    t_snapshot = np.array([float(aux.split("_")[-1][:-4]) for aux in bar_files])
    t_snapshot_sorted = np.sort(t_snapshot)
    bar_files_sorted = [x for _, x in sorted(zip(t_snapshot, bar_files))]

    #generate one file for the final Time Dependent potential
    with open(folder_name+file_name+".ini","w") as f:
        f.write(f"[Potential perturber]\ntype=Evolving\ninterpLinear={interpol}\nTimestamps\n")
        for i,filename in enumerate(bar_files_sorted):
            f.write(f"{str(t_snapshot_sorted[i])} {filename}\n")

    return folder_name+file_name+".ini"


def generate_TimeDepPot(rmin,rmax,**pot_kwargs):
    """
    Generates a potential of a growing perturbation. Starts as a m=0 mode (only mass) and evolves into the final perturbation.
    The perturbation is any arbitrary AGAMA potential and is initialised as one would with AGAMA.
    The other required parameters is a file describing the time evolution of the perturbation and rmin, rmax used for the m=0 expantion.

    Input:
    - rmin: the radius of the innermost nonzero node in the radial grid (for both potential
expansions); zero means automatic determination.
    - rmax: same for the outermost node; zero values mean automatic determination.
    - pot_kwargs: the parameters passed to AGAMA to generate the desired perturbation. Must include the "scale" parameter!
        - scale: address to the scaling file. The expected format of this file is the following:
    #Time Mass_scale Radius_scale
    0 0 1
    0.1 0.5 1
    0.2 1 1
    NOTES: time must be order in increasing order and the separation between values is done with blank spaces.

    Output:
    - perturbation: agama.Potential, time dependant potential of the perturbation
    - m0_static: agama.Potential, only the m=0 component (static, no time dependence)
    """
    #make the timedep part
    pot_pertuber = agama.Potential(**pot_kwargs)
    scaling_file = pot_kwargs["scale"]

    #make the static part
    pot_kwargs_static = pot_kwargs.pop("scale")
    pot_pertuber_static = agama.Potential(**pot_kwargs_static)
    pot_pertuber_m0_static = agama.Potential(type='CylSpline', potential=pot_pertuber_static, mmax=0, rmin=rmin, rmax=rmax)
    pot_pertuber_m0 = agama.Potential(type='CylSpline', potential=pot_pertuber_static, mmax=0, rmin=rmin, rmax=rmax, 
                                        scale=invert_scaling_file(scaling_file))

    return agama.Potential(pot_pertuber,pot_pertuber_m0),pot_pertuber_m0_static

def generate_Pot(base_pot_dict,perturb_pot_dict,_rmin=0,_rmax=20):

    #generate timedep potential
    perturb,m0_mode_stat = generate_TimeDepPot(_rmin,_rmax,**perturb_pot_dict)

    #generate base potential
    pbase_vanilla = agama.Potential(**base_pot_dict)

    