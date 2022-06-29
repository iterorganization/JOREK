#!/usr/bin/env python
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse



class C(object):
    pass
arg=C()

parser = argparse.ArgumentParser(description=__doc__);
parser.add_argument('-dir', type=str, help="Absolute path to directory.") 
parser.parse_args(namespace=arg)

scan = [[],[]]
cummulative = []

# result from standard gmres: 31,101,78,78

for p in Path(arg.dir).glob('**/*.dat'):
  fname = f"{p.name}"
  n_pars = fname.replace("_tmp_newton_niterates.dat","").replace('d','E').split('_')
  try:
   n_pars = [float(n_pars_i) for n_pars_i in n_pars]
  except:
   print(n_pars)
   exit()
  #if n_pars[-1] == 1e-3: continue
  #if n_pars[-1] == 1e-4: continue
  tng = np.loadtxt(p,unpack=True)
  tng.astype(int)
  #if max(tng[0])<139: continue  # throw away stopped simulations
  #if max(tng[2])>199: continue  # throw away non convergent simulations
  #if max(tng[1])> 49: continue
  scan[0].append(n_pars) 
  scan[1].append(tng)
  cummulative.append( [*n_pars, sum(tng[1]), int(sum(tng[2])-231)])#sum([31,101,78,78]))] )
  
  
cummulative = np.array(cummulative).T
#print(cummulative)
mins1 = [sum(tng[1]) for tng in scan[1]]
mins2 = [sum(tng[2]) for tng in scan[1]]

for imin in np.flatnonzero(np.array(mins2)==min(mins2)):
#  if 200 in scan[1][imin][2]: continue
#  if max(scan[1][imin][2]-100)>0: continue
  
  # create dict of tstep entries
  n_pars= scan[0][imin]
  tstep = scan[1][imin][0]
  ni    = scan[1][imin][1]
  giter = scan[1][imin][2] 
 
  print("NEWTON PARAMETERS newton_start={}, newton_alpha={}, newton_gamma={}, newton_eps_0={}".format(*n_pars))
  for t in np.unique(tstep):
    ni_t    = ni[tstep==t]
    giter_t = giter[tstep==t]
    print("At tstep={0:4d}, newton_i={1:3d}, total_gmres_iter={2:4d}".format(int(t),int(max(ni_t)),int(sum(giter_t))))
    

# find minima in cummulative

# select newton_start=0.6, newton_eps_gmres=1e-6
for newton_start in np.unique(cummulative[0]):#[0.55, 0.6, 0.65]:
 for newton_eps_0 in np.unique(cummulative[3]): #[1e-6, 1e-5, 1e-4, 1e-3]:
  #print(cummulative[0]==newton_start)
  cummulative_t = np.array([
  cummulative[i][
                 (cummulative[0]==newton_start)&(cummulative[3]==newton_eps_0)#&(cummulative[1]<=0.6)&(cummulative[1]>=0.325)&(cummulative[2]<=1.3)&(cummulative[2]>=0.3)#&(cummulative[-1]>450)
                 ] for i in range(len(cummulative))])
  x,y,z=cummulative_t[1],cummulative_t[2],cummulative_t[-1]
  data=pd.DataFrame(data={'newton_alpha':x,'newton_gamma':y,'additional_gmres_iter':z}).drop_duplicates()
  data.additional_gmres_iter = data.additional_gmres_iter.astype(int)
  #print(data.to_string())
  data=data.pivot(index='newton_alpha',columns='newton_gamma',
                  values='additional_gmres_iter')
  #data = data.astype(pd.Int32Dtype())
  #print(data.astype(pd.Int32Dtype()))
  ax=plt.axes()
  sns.heatmap(data,cmap="YlGnBu",annot=True,fmt='0.0f',cbar_kws={'label':'additional_gmres_iter'}, ax=ax)
  ax.set_title(r'newton_start={0:.2f}, newton_eps_0={1:.0E}'.format(newton_start,cummulative_t[3][0]))
  plt.show()
  del data

