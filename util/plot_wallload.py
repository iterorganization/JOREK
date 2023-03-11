#!/usr/bin/env python
"""
Script for postprocessing and plotting the wall loads evaluated with the particle tracer.

See particles/mod_wall_collision.f90 on how the wall loads are evaluated. This script is used
to plot and postprocess the results. To run this script, you'll need the wall input (HDF5
file containing the wall triangles) and either/both the resulting particle restart file
(where the wall IDs are stored in i_elm field, see mod_wall_collision.f90) and
the wall load file that can be exported with a routine in mod_wall_collision.f90. Place those
files in the folder together with this script, type in the filenames below and just run this
script.

This script needs pyvista to work, but you can comment wherever it is used to plot
some rudimentary results. However, pyvista is required to convert the data to a vtk file
(this is done by default and stored in wallload.vtk).

The 3D plot is by default shown in matplotlib, but you can view the wall loads interactively
by setting showpyvista = True. All in all, you should modify this file to suit your exact
needs.
"""
import numpy as np
import h5py

# Hack to convert e -> x when displaying numbers
def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))

# Function to read wall data
def read_wallinput(wallin, limiter, pyvista):
    with h5py.File(wallin,'r') as h5:
        data = h5['nodes'][:]
        ntriangle = h5['ntriangle'][0]

    verts = [] # Collection of vertices
    faces = [] # Indices for vertices that form each triangle

    # Different coordinates for element triangles' center of mass
    wd    = np.zeros((ntriangle,)) # Distance along the given contour
    wr    = np.zeros((ntriangle,)) # R
    wz    = np.zeros((ntriangle,)) # z
    wphi  = np.zeros((ntriangle,)) # Toroidal angle
    wpol  = np.zeros((ntriangle,)) # Poloidal angle

    warea = np.zeros((ntriangle,)) # Area

    # Slow loop for evaluating all those coordinates
    for i in range(ntriangle):

        # Compute area for each node
        a = data[i*9+3:i*9+6] - data[i*9+0:i*9+3] 
        b = data[i*9+6:i*9+9] - data[i*9+0:i*9+3]
        c = np.cross(a,b)
        warea[i] = np.sqrt( np.sum(c*c) ) / 2

        # Compute (r,z,tor,pol) coordinates for each element center
        x = np.sum(data[i*9+0:i*9+7:3]) / 3
        y = np.sum(data[i*9+1:i*9+8:3]) / 3
        z = np.sum(data[i*9+2:i*9+9:3]) / 3

        wr[i]   = np.sqrt(x**2 + y**2)
        wz[i]   = z
        wphi[i] = np.mod(np.arctan2( y, x ), 2*np.pi)
        wpol[i] = np.arctan2( z - z0, wr[i] - r0 )

        # Compute distance along the contour by finding first the closes segment
        dist = 1e6 # (just some large value)
        idx  = 0
        for j in range(limiter.shape[0]-1):
            dist0 = ( np.abs( (limiter[j+1,0]-limiter[j,0])*(limiter[j,1]-wz[i]) - (limiter[j,0]-wr[i])*(limiter[j+1,1]-limiter[j,1]) )
                      / np.sqrt( (limiter[j+1,0]-limiter[j,0])**2 + (limiter[j+1,1]-limiter[j,1])**2 ) )
            if dist0 < dist:
                dist = dist0
                idx  = j

        # Now sum together all segments before this one
        for j in range(idx):
            wd[i] += np.sqrt( (limiter[j+1,0]-limiter[j,0])**2 + (limiter[j+1,1]-limiter[j,1])**2 )

        # ...and finally at the projected distance along the closest segment
        j = idx
        wd[i] += ( (wr[i] - limiter[j,0]) * (limiter[j+1,0]-limiter[j,0]) + (wz[i] - limiter[j,1]) * (limiter[j+1,1]-limiter[j,1]) ) \
                 / np.sqrt( (limiter[j+1,0]-limiter[j,0])**2 + (limiter[j+1,1]-limiter[j,1])**2 )

        # Collect vertices and faces for pyvista
        verts += [data[i*9+0],data[i*9+1],data[i*9+2]]
        verts += [data[i*9+3],data[i*9+4],data[i*9+5]]
        verts += [data[i*9+6],data[i*9+7],data[i*9+8]]
        faces += [[3, i*3 + 0, i*3 + 1, i*3 + 2]]

    out = [wr, wz, wphi, wpol, warea]
    if limiter is not None:
        out.append(wd)
    
    if pyvista:
        import pyvista as pv
        wallmesh = pv.PolyData(verts, faces)
        wallmesh.cell_data['pload'] = np.zeros((ntriangle,))
        wallmesh.cell_data['eload'] = np.zeros((ntriangle,))
        wallmesh.cell_data['eloadlog10'] = np.zeros((ntriangle,))
        wallmesh.cell_data['iangle'] = np.zeros((ntriangle,))
        out.append(wallmesh)

    return out

def read_loads(wallload, warea, wallmesh=None, fnout=None):
    with h5py.File(wallload,'r') as h5:
        wallid = h5[group]["wallid"][:] - 1 # Fortran indexing starts at 1, python at 0
        eload  = h5[group]["energydepot"][:] / warea[wallid]
        pload  = h5[group]["particledepot"][:] / warea[wallid]
        iangle  = h5[group]["angleofincidence"][:]

    wetted_area = np.sum(warea[wallid])
    peak_pload  = np.amax(pload)
    peak_eload  = np.amax(eload)

    if wallmesh is not None:
        wallmesh.cell_data['pload'][wallid] = pload
        wallmesh.cell_data['eload'][wallid] = eload
        wallmesh.cell_data['eloadlog10'][wallid] = np.log10(eload)
        wallmesh.cell_data['eload'][wallid] = iangle

        if fnout is not None:
            wallmesh.save(fnout)
    
    print("Wetted area: "        + "{:.2e}".format(wetted_area) + " m^2")
    print("Peak particle load: " + "{:.2e}".format(peak_pload)  + " prt/m^2")
    print("Peak energy load: "   + "{:.2e}".format(peak_eload)  + " J/m^2")

    return [wallid, eload, pload, iangle, (wetted_area, peak_pload, peak_eload)]

def read_markers(partout, group, warea):
    with h5py.File(partout,'r') as h5:
        wallid_prt = h5["groups/"+group+"/i_elm"][:]
        r_prt   = h5["groups/"+group+"/x"][:,0]
        z_prt   = h5["groups/"+group+"/x"][:,1]
        phi_prt = 2*np.pi - np.mod(h5["groups/"+group+"/x"][:,2], 2*np.pi) # Change direction of phi from JOREK coords to right-handed
        pol_prt = np.arctan2( z_prt - z0, r_prt - r0 )
        weight  = h5["groups/"+group+"/weight"][:]

        lost = wallid_prt < 0
        wallid_prt = -wallid_prt - 1

        wettedid = np.unique(wallid_prt[lost])
        wetted_area_prt = np.sum(warea[wettedid])

        prtdepot_prt = np.zeros((wettedid.size,))
        for i,e in enumerate(wettedid):
            prtdepot_prt[i] = np.sum(weight[e == wallid_prt])

        pload_prt      = prtdepot_prt / warea[wettedid]
        peak_pload_prt = np.amax(pload_prt)

        print("\nBased on particle restart file:")
        print("Wetted area: "        + "{:.2e}".format(wetted_area_prt) + " m^2")
        print("Peak particle load: " + "{:.2e}".format(peak_pload_prt)  + " prt/m^2")

    return [r_prt, z_prt, phi_prt, pol_prt, lost, pload_prt, wettedid, pload_prt, (wetted_area_prt, peak_pload_prt)]

def plot_wallcontour(ax, wr, wz, r0, z0, wallmesh=None, phi=0):
    if wallmesh is None:
        ax.scatter( wr, wz, 1, 'black' )

    else:
        import pyvista as pv

        planemesh   = pv.Plane(center=(r0*np.cos(phi),r0*np.sin(phi),z0),
                               direction=(-np.sin(phi),np.cos(phi),0), i_size=r0*2,j_size=10).triangulate()
        wallcontour = wallmesh.intersection(planemesh)[0]
        pts = wallcontour.points
        idx = wallcontour.lines

        i = 0
        while i < idx.size:
            nk = idx[i]
            xyz = pts[idx[i+1:i+1+nk],:]
            ax.plot(np.sqrt(xyz[:,0]**2+xyz[:,1]**2), xyz[:,2], color='black')
            i = i + nk + 1

    ax.set_xlabel(r'R [m]')
    ax.set_ylabel(r'z [m]')
    ax.set_aspect('equal', 'box')
    

def plotrz_prt(ax, r_prt, z_prt, lost):
    h1 = ax.scatter(r_prt[lost], z_prt[lost], 1, 'red')
    h2 = ax.scatter(r_prt[~lost], z_prt[~lost], 1, 'blue')
    lh = ax.legend([h1, h2], ["Lost", "Confined"])

    ax.set_xlabel(r'R [m]')
    ax.set_ylabel(r'z [m]')
    ax.set_aspect('equal', 'box')

    return (h1, h2, lh)
    
def plotrz(ax, wr, wz, wallid):
    h1 = ax.scatter(wr[wallid], wz[wallid], 1, 'red')

    ax.set_xlabel(r'R [m]')
    ax.set_ylabel(r'z [m]')
    ax.set_aspect('equal', 'box')

    return (h1)

def plotphipol_prt(ax, phi_prt, pol_prt, lost):
    h1 = ax.scatter(phi_prt[lost]*180/np.pi, pol_prt[lost]*180/np.pi, 1, 'red')

    ax.set_xlabel(r'Tor [deg]')
    ax.set_ylabel(r'Pol [deg]')

    ax.set_xlim([0, 360])
    ax.set_ylim([-180, 180])

    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([-180, -90, 0, 90, 180])

    return(h1)

def plotphipol(ax, wphi, wpol, wallid):
    h1 = ax.scatter(wphi[wallid]*180/np.pi, wpol[wallid]*180/np.pi, 1, 'red')

    ax.set_xlabel(r'Tor [deg]')
    ax.set_ylabel(r'Pol [deg]')

    ax.set_xlim([0, 360])
    ax.set_ylim([-180, 180])

    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([-180, -90, 0, 90, 180])

    return h1

def plotianglehist(ax, iangle, pload):
    y2, x = np.histogram(iangle*180/np.pi, bins=np.linspace(0,90,179), density=True)
    y1, x = np.histogram(iangle*180/np.pi, bins=np.linspace(0,90,179), density=True, weights=pload)
    xc = (x[:-1] + x[1:])/2

    h1,=s1.plot(xc, y2, color='gray')
    h2,=s1.plot(xc, y1, color='black')

    ax.legend([h2,h1], ["Weighted with heat load", "Unweighted"], frameon=False)

    ax.set_xlabel("Angle of incidence [deg]")
    ax.set_xlim(0,90)
    ax.set_xticks([0, 30, 60, 90])
    ax.set_ylim(bottom=0)
    ax.get_yaxis().set_visible(False)

def plotpattern(ax, wphi, wd, bins, pload=None, eload=None):
    import matplotlib as mpl
    
    dens0 = np.histogram2d(wphi[wallid], wd[wallid], bins=bins)[0]
    if pload is not None:
        dens,xg,yg = np.histogram2d(wphi[wallid], wd[wallid], bins=bins, weights=pload)
        title = r"prt/m$^2$"
    else:
        dens,xg,yg = np.histogram2d(wphi[wallid], wd[wallid], bins=bins, weights=eload)
        title = r"J/m$^2$"

    xg = xg[:-1] + (xg[1] - xg[0]) / 2
    yg = yg[:-1] + (yg[1] - yg[0]) / 2
    
    cmap = mpl.colormaps["Reds"].copy()
    cmap.set_bad(color=[1.0,1.0,1.0])

    a = dens/dens0
    h1 = ax.pcolormesh(xg*180/np.pi, yg, a.T, cmap=cmap,
                       norm=mpl.colors.LogNorm(vmin=np.nanmin(a), vmax=np.nanmax(a)))

    ax.set_xlabel("Toroidal angle [deg]")
    ax.set_ylabel("Distance along wall [m]")
    cax = plt.colorbar(h1, ax=ax, location='top')
    cax.set_label("Mean load " + title)

    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    
    return (h1, cax)

def plothist(ax, warea, eload):
    idx = np.argsort(-eload)
    y = (eload[idx])
    x = np.cumsum(warea[idx])

    ax.fill_between(x, 0, y, color='red')
    ax.plot(x, y, linewidth=2, color='black')
    ax.set_xlim(left=0)
    ax.set_yscale('log')

    ax.set_xlabel(r"Area [m$^2$]")
    ax.set_ylabel(r"Minimum load [J/m$^2$]")
    

def plotmesh(wallmesh, r0, z0, phicam=0, load='eload', bounds=None, ax=None):

    # Removes bottom mesh so it is not in the way of camera
    if bounds is not None:
        wallmesh = wallmesh.clip_box(bounds)

    # This makes non-loaded wall elements appear grey instead red-tinted
    cmap = mpl.colormaps["Reds"].copy()
    cmap.set_bad(color=[0.9,0.9,0.9])

    if live:
        p = pv.Plotter()
    else:
        p = pv.Plotter(off_screen=True)
        p.store_image = True

    a = wallmesh.cell_data[load]
    cmin = np.amin(a[a>0])
    cmax = np.amax(a[a>0])
    p.add_mesh(wallmesh, scalars=load, cmap=cmap, clim=[cmin,cmax])
    if not live:
        p.remove_scalar_bar()

    p.camera.position = (r0*np.cos(phicam), r0*np.sin(phicam), z0-5)
    p.camera.focal_point = (0.9*r0*np.cos(phicam), 0.9*r0*np.sin(phicam), z0+10)
    p.show()

    if ax is not None:
        ax.imshow(p.image)

        norm = mpl.colors.Normalize(vmin=cmin,vmax=cmax)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cax = plt.colorbar(sm,ax=ax,location='top')
        cax.set_label(r"Energy load J/m$^2$")

        ax.set_xticks([])
        ax.set_yticks([])


if __name__ == "__main__":
    import pyvista as pv# Uncomment if you don't have access

    ## Settings ##

    # Make plots or just output the 0D quantities
    plot = True

    # Axis (R,z) coordinates (these determine poloidal angle in plots)
    r0 = 6.2
    z0 = 0.1

    # Wall input data
    wallin   = "iterwall_offset20cm.h5"

    # Wall loads
    wallload = "wallload.h5"

    # Particle restart file containing the output (or None)
    partout  = "part_out.h5"
    #partout = None

    # Group to be plotted
    group = "001"

    # Define a contour along the wall for plotting losses along it
    limiter = np.array([[4,4,4.2,4.9,5.7,6.5,6.8],
                        [3,4,4.2,4.6,4.4,3.8,3.5]]).T
    #limiter = None

    # Choose whether to show interactive view using pyvista or final product with pyplot
    live = False

    # Filename if you want to save mesh & loads in vtk format
    fnvtk = None
    
    ## End of settings ##

    try:
        pv
        pyvista = True
    except:
        pyvista = False

    out = read_wallinput(wallin, limiter, pyvista)
    wr = out[0]; wz = out[1]; wphi = out[2]; wpol = out[3]; warea = out[4]

    i = 5
    if limiter is not None:
        wd = out[i]
        i += 1

    wallmesh = None
    if pyvista:
        wallmesh = out[i]

    del out

    [wallid, eload, pload, iangle, quantities] = read_loads(wallload, warea, wallmesh=wallmesh, fnout=fnvtk)

    if partout is not None:
        [r_prt, z_prt, phi_prt, pol_prt, lost, pload_prt, wettedid, pload_prt, quantities_prt] = read_markers(partout, group, warea)
        
        
    # Plot
    if plot:
    
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        cm = 1/2.54
        params = {'legend.fontsize': 12,
                  'axes.labelsize':  12,
                  'axes.titlesize':  12,
                  'xtick.labelsize': 12,
                  'ytick.labelsize': 12,
                  'font.size' : 12,
                  'text.usetex' : True}
        plt.rcParams.update(params)


        fig1 = plt.figure()
        gs = GridSpec(1,2,figure=fig1,width_ratios=[1, 2], right=0.85)
        s1 = fig1.add_subplot(gs[0])
        s2 = fig1.add_subplot(gs[1])

        plot_wallcontour(s1, wr, wz, r0, z0, wallmesh=wallmesh, phi=0.1)
        if partout is None:
            plotrz(s1, wr, wz, wallid)
            plotphipol(s2, wphi, wpol, wallid)
        else:
            plotrz_prt(s1, r_prt, z_prt, lost)
            plotphipol_prt(s2, phi_prt, pol_prt, lost)
        
        s2.annotate("Wetted area: "        + r"${0:s}$".format(as_si(quantities[0],2)) + r" m$^2$",     (10,-90) )
        s2.annotate("Peak particle load: " + r"${0:s}$".format(as_si(quantities[1],2)) + r" prt/m$^2$", (10,-110) )
        s2.annotate("Peak energy load: "   + r"${0:s}$".format(as_si(quantities[2],2)) + r" J/m$^2$",   (10,-130) )

        s2.yaxis.set_label_position("right")
        s2.yaxis.tick_right()

        if limiter is not None:
            s1.plot(limiter[:,0], limiter[:,1])
    
            fig2 = plt.figure()
            s1 = fig2.add_subplot(1,1,1)
            plotpattern(s1, wphi, wd, bins=[90,40], pload=None, eload=eload)

        bounds = [-20, 20, -20, 20, -10, -2] # Comment to remove clipping

        fig3 = plt.figure()
        s1 = fig3.add_subplot(1,1,1)
        plothist(s1, warea[wallid], eload)
        
        ax = None
        if not live:
            fig4 = plt.figure(figsize=(6,6))
            s1 = fig4.add_subplot(1,1,1)
        plotmesh(wallmesh, r0, z0, phicam=330 * np.pi / 180, load='eloadlog10', bounds=bounds, ax=s1)

        fig5 = plt.figure()
        s1 = fig5.add_subplot(1,1,1)
        plotianglehist(s1, iangle, pload)

        plt.tight_layout()

        # Uncomment to save figures
        fig1.savefig("wallloadscatter.png", format="png", dpi=96*2)
        fig2.savefig("wallloadmean.png", format="png", dpi=96*2)
        fig3.savefig("wallloadhist.png", format="png", dpi=96*2)
        fig4.savefig("wallload3d.png", format="png", dpi=96*2)
        fig5.savefig("angleofincidence.png", format="png", dpi=96*2)
        plt.show()
