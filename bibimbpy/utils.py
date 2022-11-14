import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcols
import scipy.stats

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

def plot_result(df_sum,xe,ye,xmin,xmax,ymin,ymax,figname):
    #TO-DO: generalise for any two arbitrary coordinates
    X,Y = np.meshgrid(0.5*(xe[1:]+xe[:-1]),0.5*(ye[1:]+ye[:-1]))
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.gca()
    cc = ax.pcolormesh(X,Y,df_sum.T,cmap="jet")
    plt.colorbar(cc)
    
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_xlabel(r"$V_R$ [km$\,$s$^{-1}$]")
    ax.set_ylabel(r"$V_\phi$ [km$\,$s$^{-1}$]")
    #ax.set_title(f"Time = {t_} Gyr // R = {r0} kpc // Phi = {phi0} degrees")
    
    fig.savefig(figname)

def store_results(filename,x,y,h):
    #TO-DO: generalise for any two arbitrary coordinates
    np.save(filename+"vr_edges",x)
    np.save(filename+"vphi_edges",y)
    np.save(filename+"dfsum",h)

def DFhistogram2d(x,y,df_eval,bins):
    """
    Add up the DF values of the particles simulated. X,Y can be anything but should be the two variables that have been explored in the backward integration a t=tfin! That is, the initial conditions, NOT the result of the backwards integration.

    For instance, normally we would use Backwards Integration to generate histograms in Vr-Vphi. In that case, X and Y should be Vr and Vphi (in no particular order).

    Input:
    - x: value of coordinate 1 of phase-space of the particles at t=tfin.
    - y: value of coordinate 2 of phase-space of the particles at t=tfin.
    - df_eval: value of the Distribution Function at t=0 of each particle.
    - bins: iterable containing the bins in x and y to use for the histogram 2D.

    Output:
    - df_sum: Sum of distribution function in each bin.
    """
    
    df_sum,xe,ye,_ = scipy.stats.binned_statistic_2d(x,y,df_eval,statistic="sum",bins=bins)

    return df_sum,xe,ye

def write_potential(pot_params_dict,name):

    pot_str = extract_params(pot_params_dict)

    with (name,"w") as f:
        f.write("[Potential]\n")
        f.write(pot_str)
    
    return None

def extract_params(pot_params_dict):
    pot_str = ""
    for key,param in pot_params_dict.items():
        pot_str += key + " = " + param + "\n"
    return pot_str


def invert_scaling_file(file):
    """
    This function takes a file with the right format (see below) and inverts the second and third column. In other words, the time evolution of the scaling factors is reversed.

    The expected format of the file is the following:
    #Time Mass_scale Radius_scale
    0 0 1
    0.1 0.5 1
    0.2 1 1

    NOTES: time must be order in increasing order and the separation between values is done with blank spaces.
    """
    time, mass, radius = [],[],[]
    file_split = file.split("/")
    directories = ""
    for f in file_split[:-1]: 
        directories += f + "/"
    filename = file_split[-1]
    #print(directories+"reversed_"+filename)
    with open(directories+"reversed_"+filename,"w") as r:
        with open(file,"r") as f:
            for i,line in enumerate(f.readlines()):
                if line.startswith("#"):
                    r.write(line)
                else:
                    line_split = line.split()
                    time.append(line_split[0])
                    mass.append(line_split[1])
                    radius.append(line_split[2])
            for i,t in enumerate(time):
                new_line = t + " " + mass[-i-1] + " " + radius[-i-1] + "\n"
                r.write(new_line)
            r.write(str(float(t)+10) + " " + mass[-i-1] + " " + radius[-i-1] + "\n")

    return directories+"reversed_"+filename

def generate_scaling_file(tf,mode,nodes,filename,_amp=1):
    """
    Generate a file containing the values of the scaling factor as a function of time.

    Input:
    - tf: final time. Can be positive or negative. Regardeless, the file will always be ordered in increasing time.
    - mode: Choose between "linear", "exponential" or "dehnen"
    - nodes: number of points or nodes between 0 and tf
    - filename: full address of the file in which to store the results
    - _amp: final value of the scaling factor (initial value will always be 0) 
    """

    def dehnen(t,tf):
        def aux_func(t,tf):
            return 2*t/tf-1
        aux = aux_func(t,tf)
        return (3/16*aux**5-5/8*aux**3+15/16*aux+1/2)
    
    #create time vector
    t = np.sort(np.linspace(0,tf,nodes))

    #create scaling factors
    if mode=="linear":
        amp = np.linspace(0,1,nodes)
    elif mode=="exponential":
        amp = np.logspace(-10,0,nodes)
    elif mode=="dehnen":
        amp = dehnen(t,tf)

    amp = _amp*amp

    #write to file
    with open(filename,"w") as f:
        f.write("#Time mass_scale radius_scale")
        for i,t_ in enumerate(t):
            new_line = str(t_) + " " + str(amp[i]) + " " + "1" + "\n"
            f.write(new_line)

    return None


def default_df_gen_func(points,pot_base,df):
    #action-angle
    af = agama.ActionFinder(pot_base)

    #compute actions at t=0
    actions = af(points)

    return df(actions)


