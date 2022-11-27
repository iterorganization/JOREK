"""
Script for reading and plotting the output from diagnostics/poincare.f90
"""
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt

# Read data from HDF5 file
#
# Returned fields:
# r       : R-coordinate
# z       : z-coordinate
# phi     : Toroidal angle [rad]
# psi     : Normalized psi
# iprt    : Marker ID this data point correspond to (1,...,Nmrk)
# iplane  : Poincare plane this data point correspond to (default is phi-tor plane = 1 and Rz-plane = 2 )
# mil     : Mileage, i.e., the distance [m] marker has travelled up to that point
# clength : Connection length at each point or NaN if marker was not lost
#
# To get indices of all points for marker i at plane j, use:
# idx = np.logical_and.reduce([iprt == i, iplane == j])
def read(fn):
    with h5py.File(fn, "r") as f:
        r      = f["r"][:]
        z      = f["z"][:]
        phi    = f["phi"][:]
        psi    = f["psi"][:]
        iprt   = f["iprt"][:]
        iplane = f["pncrid"][:]
        dist   = f["mileage"][:]
        mil    = f["mil"][:]

    # mileage is the connection length from the marker *initial* position whereas
    # mil stores the distance the marker had travelled up to that stored point. Therefore
    # we get the connection length at each Poincare crossing as clength = dist - mil
    prt = np.unique(iprt) - 1 # Marker indices start from 1 but python indexing from 0
    clength = np.zeros(mil.shape)
    for i in prt:
        idx = iprt == i + 1
        clength[idx] = (dist[i] - mil[idx])

        # If mil reaches mileage, then this marker was confined so set connection length to NaN
        if np.amin(clength[idx]) == 0: clength[idx] = np.nan

    return (r, z, phi, psi, iprt, iplane, mil, clength)


def plot(s1, s2, r, z, phi, psi, iprt, iplane, clength=None):
    s1.set_xlabel("R [m]")
    s1.set_ylabel("z [m]")

    s2.set_xlabel(r"$\psi_n$")
    s2.set_ylabel("Toroidal angle [rad]")
    s2.set_ylim(0,2*np.pi)

    
    if clength is not None:
        # Plot with connection length but first plot confined markers separately.
        idx = np.isnan(clength)
        plot(s1, s2, r[idx], z[idx], phi[idx], psi[idx], iprt[idx], iplane[idx], clength=None)

        # Case where all marker are confined
        if( np.sum(idx) == clength.size ):
            print("All markers were confined. Therefore connection length is not shown.")
            return

        idx = ~idx
        r = r[idx]; z = z[idx]; phi = phi[idx]; psi = psi[idx]; iprt = iprt[idx]; iplane = iplane[idx]; clength = clength[idx];
        cmin = np.nanmin(clength)
        cmax = np.nanmax(clength)
        cmap = mpl.cm.get_cmap('Blues').copy()

        idx = iplane == 2
        s1.scatter(r[idx], z[idx], 2, clength[idx], norm=mpl.colors.LogNorm(vmin=cmin,vmax=cmax), cmap=cmap)

        idx = iplane == 1
        h2 = s2.scatter(psi[idx], phi[idx], 2, clength[idx], norm=mpl.colors.LogNorm(vmin=cmin,vmax=cmax), cmap=cmap)

        cax = plt.colorbar(h2, ax=s2, location='right')
        cax.set_label("Connection length [m]")
        
    else:
        # Assign each marker a color when plotting to help in separating field lines
        # from one another. These colors are cycled from a group of six colors.
        colors = mpl.cm.get_cmap('Reds')(np.linspace(0,1,6))

        prt = np.unique(iprt)
        for i in prt:
            idx = np.logical_and.reduce([iprt == i, iplane == 2])
            s1.scatter(r[idx], z[idx], 1, color=colors[np.mod(i,6),:])

        for i in prt:
            idx = np.logical_and.reduce([iprt == i, iplane == 1])
            s2.scatter(psi[idx], phi[idx], 1, color=colors[np.mod(i,6),:])
    
if __name__ == "__main__":
    fn = "poincare.h5"


    cm = 1/2.54
    params = {'legend.fontsize': 12,
              'axes.labelsize':  12,
              'axes.titlesize':  12,
              'xtick.labelsize': 12,
              'ytick.labelsize': 12,
              'font.size' : 12,
              'text.usetex' : True}
    plt.rcParams.update(params)
    
    fig = plt.figure()

    s1 = fig.add_subplot(1,2,1)
    s2 = fig.add_subplot(1,2,2)

    r, z, phi, psi, iprt, iplane, mil, clength = read(fn)

    #clength = None # Uncomment to not display connection length
    plot(s1, s2, r, z, phi, psi, iprt, iplane, clength=clength)
    plt.tight_layout()

    # Uncomment to save figure
    fig.savefig("poincare.png", format="png", dpi=96*2)
    
    plt.show()
