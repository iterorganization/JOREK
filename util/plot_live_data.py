#!/usr/bin/env python
""" Reads and plot data from macroscopic_variables.dat
    https://www.jorek.eu/wiki/doku.php?id=plot_live_data.py
    iholod@ipp.mpg.de
"""

import os.path
import argparse
import numpy as np

from matplotlib import rcParams
import itertools

import matplotlib.pyplot as plt
if 'classic' in plt.style.available: plt.style.use('classic')
import matplotlib.colors as colors
import matplotlib.cm as cmx

def plot1d(x,y,xlbl,ylbl,**kwargs):
    """ 1D plot
        takes x,y,xlbl,ylbl as input
        extra arguments: xmin,xmax,ymin,ymax,title,fname
        extra arguments x1,y1
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    
    plt.rc('font', size=10)
    #plt.rcParams.update({"text.usetex": True})
    f, ax = plt.subplots(figsize=(8, 6), dpi=120)
    lbl0 = None
    lbl1 = None
   
    legend=[]
    
    if type(x)!=list:
        x=[x]
        y=[y]
    if "legend" in kwargs: 
        legend=kwargs["legend"]
        if (legend==None):
            printLegend = False
        else:
            printLegend = True
            if type(legend)!=list:
                legend=[legend]
            if len(legend)!=len(x):
                print("error with legend")
                return()
    else:
        printLegend=False
    
    for ip in range(len(x)):
        label = None
        if (printLegend): label = legend[ip]
        ax.plot(x[ip],y[ip],linewidth=1,marker=".",markersize=5,label=label)

#     ax.yaxis.set_major_formatter(FormatStrFormatter('%1.2e'))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 3))
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    if "title" in kwargs:
        ax.set_title(kwargs["title"])
    # ax1.legend(loc='best', shadow=True)
    if "xmin" in kwargs:
        ax.set_xlim([kwargs["xmin"], ax.get_xlim()[-1]])
    if "xmax" in kwargs:
        ax.set_xlim([ax.get_xlim()[0],kwargs["xmax"]])
    if "ymin" in kwargs:
        ax.set_ylim([kwargs["ymin"], ax.get_ylim()[-1]])
    if "ymax" in kwargs:
        ax.set_ylim([ax.get_ylim()[0],kwargs["ymax"]])
    if "logy" in kwargs:
        if (kwargs["logy"]): ax.set_yscale("log")

    if printLegend: ax.legend(loc='best', shadow=True)
    
    if (("fname" in kwargs) and (kwargs["fname"]!=None)):
        fname = kwargs["fname"]
        f.savefig(fname+".png",dpi=200,bbox_inches="tight")
        f.savefig(fname+".eps",bbox_inches="tight")
        
    return(f,ax)

dirname=os.getcwd()
#dirname ='C:\\Users\\iholod\\tmp'

class C(object):
    pass
arg=C()

parser = argparse.ArgumentParser(description=__doc__);
parser.add_argument('-legend', action='store_false', help="Print legend");
parser.add_argument('-q', type=str, help="Plot given quantity");
parser.add_argument('-f', type=str, help="File name");
parser.add_argument('-l', action='store_true', help="List plottable quantities");
parser.add_argument('-ps', action='store_true', help="Save plot into eps file");
parser.add_argument('-nology', action='store_true', help="Linear y axis");
 
parser.parse_args(namespace=arg)

addLegend = arg.legend

logy = True
logy = not arg.nology

if (arg.q!=None):
    varname = arg.q
else:
    varname = "energies"

if (arg.f!=None):
    fname = arg.f
else:
    fname="macroscopic_vars.dat"
    
fname = os.path.join(dirname, fname)
fid = open(fname,'r',encoding = "ISO-8859-1")

if (arg.l):
    for line in fid:
        line = ' '.join(line.split()).split()
        if (len(line)==0): continue
        if ("@plottable:"==line[0]):
            plottable = str(' '.join(line[1:]))
            break
    print(plottable)
    exit()
    
t = []
dat = []

for line in fid:
    line = ' '.join(line.split()).split()
    if (len(line)==0): continue
    if ("@"+varname+"_xlabel:"==line[0]):
        xlbl = str(' '.join(line[1:]))
    if ("@"+varname+"_ylabel:"==line[0]):
        ylbl = str(' '.join(line[1:]))        
    if ("@n_"+varname+":"==line[0]):
        nvar = int(line[-1])
    if ("@"+varname+":"==line[0]):
        if ("%" in line[1]):
            lbl = line
        else:
            t.append(float(line[1]))
            dat.append([float(line[i]) for i in range(2,len(line))])

fid.close()

dat = np.array(dat)
tt = [t for i in range(nvar)]
vars = [dat[:,i] for i in range(nvar)]

if (not addLegend):
    legend = None
else:
    legend=lbl[2:]
    
if (not arg.ps):
    fname = None
else:
    fname = os.path.join(dirname, varname+"_plot")
    
f,ax = plot1d(tt,vars,xlbl,ylbl,logy=logy,legend=legend,fname=fname)

plt.show()

exit()
