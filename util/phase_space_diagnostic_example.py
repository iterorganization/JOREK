#%%

# Simple script to both visualize power exchange and integrate it to show energy conservation between kernels and particles
# and also to visualize interactively the 4D distribution function
# Can be run as notebook as well.

import h5py 
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from mpl_toolkits.mplot3d import Axes3D
# %%

b= h5py.File("fourD_dist_00000.h5")
c=h5py.File("power_exchange_00000.h5")
#%%
mvpar = np.array(c['grids/mgrid_1'])
mmu = np.array(c['grids/mgrid_2'])
vpar = np.array(c['grids/grid_1'])
mu = np.array(c['grids/grid_2'])
pex = np.array(c['values'])
plt.pcolormesh(mvpar,mmu,pex)
p_ex_tot = integrate.trapezoid(integrate.trapezoid(pex,vpar,axis=0),mu,axis=0)
print("Power exchanged:" ,p_ex_tot)
#%%
R1D = np.array(b['grids/grid_1'])
Z1D = np.array(b['grids/grid_2'])
P1D = np.array(b['grids/grid_3'])
E1D = np.array(b['grids/grid_4'])
mR = np.array(b['grids/mgrid_1'])
mR_red = mR[:,:,0,0]
mZ = np.array(b['grids/mgrid_2'])
mZ_red = mZ[:,:,0,0]
mP = np.array(b['grids/mgrid_3'])
mE = np.array(b['grids/mgrid_4'])
val = np.array(b['values'])
RZproj = integrate.trapezoid(integrate.trapezoid(val,P1D,axis=2),E1D,axis=2)
total = integrate.trapezoid(integrate.trapezoid(RZproj,R1D,axis=0),Z1D,axis=0)
print("Total amount of particles",total)
# %%

# %%
F0= val/np.max(val)
class ModifiableF:
    def __init__(self,F0,R1D,Z1D, mP,mE,mR_red,mZ_red,RZproj):
        fig = plt.figure()
        ax1=fig.add_subplot(121)
        ax2=fig.add_subplot(122)
        self.pcRZ= ax1.pcolormesh(mR_red,mZ_red,RZproj,shading="gouraud") 
        self.pcolormesh = ax2.pcolormesh(mP[0,0,:,:],mE[0,0,:,:],F0[0,0,:,:],shading="gouraud",vmin = np.min(F0),vmax=np.max(F0))
        self.pl1, = ax1.plot([R1D[0]],[Z1D[0]],'ro')
        self.F0 = F0
        self.R1D = R1D
        self.Z1D = Z1D
        self.cid = self.pcolormesh.figure.canvas.mpl_connect('button_press_event', self)
        ax1.set_aspect('equal')
    def __call__(self,event):
        if event.inaxes!=self.pcRZ.axes: return
        Rindex=np.argmin(np.abs(R1D-event.xdata))
        Zindex=np.argmin(np.abs(Z1D-event.ydata))
        print(Rindex,Zindex)
        self.pcolormesh.set_array(self.F0[Rindex,Zindex,:,:].ravel())
        self.pl1.set_data([event.xdata],[event.ydata])
        self.pl1.figure.canvas.draw()
Fmod = ModifiableF(val,R1D,Z1D,mP,mE,mR_red,mZ_red,RZproj)



plt.show()
# %%
