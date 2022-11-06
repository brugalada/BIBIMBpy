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
    
    df_sum,_,_,_ = scipy.stats.binned_statistic_2d(x,y,df_eval,statistic="sum",bins=bins)

    return df_sum
