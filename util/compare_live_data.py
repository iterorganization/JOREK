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


def latex_float(f):
 float_str="{0:.2g}".format(f)
 if "e" in float_str:
  base, exp = float_str.split("e")
  return r"{0}\times 10^{{{1}}}".format(base,int(exp))
 else:
  return float_str

def merge_2plots(f1,ax1,f2,ax2):
 fig, (ax2,ax4) = plt.subplots(2,1,sharex=True)
 for l in ax1.lines:
  x=l.get_xdata()
  y=l.get_ydata()
  label=l.get_label()
  marker=l.get_marker()
  markersize=l.get_markersize()
  ax3.plot(x,y,marker=marker,markersize=markersize,label=label)
 for l in ax2.lines:
  x=l.get_xdata()
  y=l.get_ydata()
  label=l.get_label()
  marker=l.get_marker()
  markersize=l.get_markersize()
  ax4.plot(x,y,marker=marker,markersize=markersize,label=label)
# set ylim, xlim
 ax3.set_ylim(ax1.get_ylim());ax3.set_xlim(ax1.get_xlim())
 ax4.set_ylim(ax2.get_ylim());ax4.set_xlim(ax2.get_xlim())
 ax3.set_yscale(ax1.get_yscale())
 ax4.set_yscale(ax2.get_yscale())
 ax3.set_ylabel(ax1.get_ylabel())
 ax4.set_ylabel(ax2.get_ylabel())
 ax4.set_xlabel(ax2.get_xlabel())
 ax3.legend(loc='best',shadow=True)
 ax4.legend(loc='best',shadow=True)
 plt.close(f1)
 plt.close(f2)
 return fig, (ax3,ax4)

def plot1d(x,y,xlbl,ylbl,**kwargs):
    """ 1D plot
        takes x,y,xlbl,ylbl as input
        extra arguments: xmin,xmax,ymin,ymax,title,fname
        extra arguments x1,y1
    """
    import matplotlib.pyplot as plt
    import itertools
    from matplotlib.ticker import FormatStrFormatter
    cm = plt.get_cmap('gist_rainbow')
    plt.rc('font', size=8)
    #plt.rcParams.update({"text.usetex": True})
    f, ax = plt.subplots(figsize=(8, 6), dpi=120)
    NUM_COLOURS = len(x)
    marker=itertools.cycle((',','+','.','o','*','^'))
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
        if (printLegend): label = r"${}$".format(legend[ip].replace('"',''))
        ax.plot(x[ip],y[ip],linewidth=1,marker=next(marker),
                markersize=5,label=label,color=cm(ip/NUM_COLOURS))

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


def read_data(dirname):
 varname="energies"
 fname="macroscopic_vars.dat"
 fname = os.path.join(dirname, fname)
 fid = open(fname,'r',encoding = "ISO-8859-1")

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
 try:
  logfile=os.path.join(dirname,"logfile")
  fid=open(logfile,'r')
  for line in fid:
   if ' tstep ' in line: 
    tstep = float(line.replace('\n','').split("=")[-1])
    tstep=[r'$\delta_t={}$'.format(latex_float(tstep))]; 
    break
  fid.close()
 except:
  pass
 return t, dat, xlbl, ylbl, nvar, lbl, tstep

dirname=os.getcwd()
#dirname ='C:\\Users\\iholod\\tmp'

class C(object):
    pass
arg=C()

parser = argparse.ArgumentParser(description=__doc__);
parser.add_argument('-legend', action='store_false', help="Print legend");
parser.add_argument('-nology', action='store_true', help="Linear y axis");
parser.add_argument('-compare', type=str, nargs='+') 
parser.parse_args(namespace=arg)

addLegend = arg.legend

logy = True
logy = not arg.nology

ttt, data, legends, variables, tsteps = [], [], [], [], []   
for f in arg.compare: 
 
 t, dat, xlbl, ylbl, nvar, lbl, tstep= read_data(f)
 dat    = np.array(dat)
 tt     = [t for i in range(nvar)]
 vars   = [dat[:,i] for i in range(nvar)]
 legend = lbl[2:]


 ttt.extend(tt)
 data.extend(dat)
 legends.extend(legend)
 variables.extend(vars)
 tsteps.extend(tstep)

f,ax = plot1d(ttt,variables,xlbl,ylbl,logy=logy,legend=legends) 
#if len(arg.compare)==1:
# h, l = ax.get_legend_handles_labels()
# ph=[plt.plot([],marker="",ls="")[0]]
# handles=ph[:1]+h[:nvar]
# labels=tsteps[:1]+l[:nvar]
# print(handles, labels)
# 
# leg = plt.legend(handles,labels,ncol=len(arg.compare),loc='best')
# for vpack in leg._legend_handle_box.get_children():
#  for hpack in vpack.get_children()[:1]:
#   hpack.get_children()[0].set_width(0)



if len(arg.compare)>0:
 h, l = ax.get_legend_handles_labels()
 ph=[plt.plot([],marker="",ls="")[0]]*len(arg.compare)
 #handles=ph[:1]+h[:nvar]+ph[1:]+h[nvar:]
 #labels=tsteps[:1]+l[:nvar]+tsteps[1:]+l[nvar:]
 handles2=np.array([ph[i:i+1]+h[i*nvar:(i+1)*nvar] for i in range(len(arg.compare))]).flatten()
 labels2 =np.array([tsteps[i:i+1]+l[i*nvar:(i+1)*nvar] for i in range(len(arg.compare))]).flatten()
# print(labels,labels2)

 leg = plt.legend(handles2,labels2,ncol=len(arg.compare),loc='best')
 for vpack in leg._legend_handle_box.get_children():
  for hpack in vpack.get_children()[:1]:
   hpack.get_children()[0].set_width(0)

#if len(arg.compare)>2:
# h, l = ax.get_legend_handles_labels()
# ph=[plt.plot([],marker="",ls="")[0]]*len(arg.compare)
#
# handles=[ph[i-1:i]+h[(i-1)*nvar:nvar*i] for i in range(len(arg.compare))]
# 
# labels=tsteps[:1]+l[:nvar]+tsteps[1:]+l[nvar:]
#
# leg = plt.legend(handles,labels,ncol=len(arg.compare),loc='best')
# for vpack in leg._legend_handle_box.get_children():
#  for hpack in vpack.get_children()[:1]:
#   hpack.get_children()[0].set_width(0)




plt.show()
 
exit()
