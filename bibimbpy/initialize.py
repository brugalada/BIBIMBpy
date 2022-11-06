import numpy as np

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

    #TO-DO: check a better way to do this check
    if r is isinstance((list,tuple,np.ndarray)):
        iter_vars["r"] = r
    else:
        final_vars["r"] = r

    if phi is isinstance((list,tuple,np.ndarray)):
        iter_vars["phi"] = phi
    else:
        final_vars["phi"] = phi

    if z is isinstance((list,tuple,np.ndarray)):
        iter_vars["z"] = z
    else:
        final_vars["z"] = z

    if vr is isinstance((list,tuple,np.ndarray)):
        iter_vars["vr"] = vr
    else:
        final_vars["vr"] = vr

    if vphi is isinstance((list,tuple,np.ndarray)):
        iter_vars["vphi"] = vphi
    else:
        final_vars["vphi"] = vphi

    if vz is isinstance((list,tuple,np.ndarray)):
        iter_vars["vz"] = vz
    else:
        final_vars["vz"] = vz

    #Cross variables with a meshgrid (accounting for all the cases)
    if len(iter_vars.keys)==2:
        var1,var2 = np.meshgrid(iter_vars.items[0],iter_vars.items[1],indexing="ij")
        final_vars[vars.key[0]] = var1.flatten()
        final_vars[vars.key[1]] = var2.flatten()
    else:
        #TO-DO: allow for more variable to iterate through
        raise ValueError("Sorry! More than two variables is not implemented yet!")

    #Generate initial conditions
    x0 = final_vars["r"]*np.cos(final_vars["phi"])
    y0 = final_vars["r"]*np.sin(final_vars["phi"])
    z0 = final_vars["z"]
    vx0 = final_vars["vr"]*np.cos(final_vars["phi"])-final_vars["vphi"]*np.sin(final_vars["phi"])
    vy0 = final_vars["vr"]*np.sin(final_vars["phi"])+final_vars["vphi"]*np.cos(final_vars["phi"])
    vz0 = final_vars["vz"]

    return np.column_stack((x0,y0,z0,vx0,vy0,vz0)),var1,var2


def generate_TimeDepPot(interpol="false"):
    """
    interpol = true or false
    """

    #Generate individual files for each step of the perturbation
    for i,m in enumerate(bar_mass):
        t = bar_time[i]
        with open("Potentitals/Ferrer_bar/Ferrer_srad5_qy0.15_qz0.035_{}.ini".format(str(m)+"_"+str(t)),"w") as f:
            f.write("[Potential bar]\n")
            f.write(f"type=Ferrers\nmass={m}\nscaleRadius={scaleRadius}\naxisRatioY={axisRatioY}\naxisRatioZ={axisRatioZ}")

    #make sure that the files are sorted by time
    bar_files = [f_ for f_ in os.listdir("Potentitals/Ferrer_bar/") if f_.endswith(".ini") and f_.startswith("Ferrer_srad")]
    t_snapshot = np.array([float(aux.split("_")[-1][:-4]) for aux in bar_files])
    t_snapshot_sorted = np.sort(t_snapshot)
    bar_files_sorted = [x for _, x in sorted(zip(t_snapshot, bar_files))]

    #generate one file for the final Time Dependent potential
    with open(filename_TimeDepPot,"w") as f:
        f.write(f"\n[Potential perturber]\ntype=Evolving\ninterpLinear={interpol}\nTimestamps\n")
        for i,filename in enumerate(bar_files_sorted):
            f.write(f"{str(t_snapshot_sorted[i])} {filename}\n")

    return filename_TimeDepPot