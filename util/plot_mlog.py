#!/usr/bin/env python
""" Reads memory usage files and optionally 
    JOREK timeline file. Plot memory usage vs time with 
    optional time line marks
    https://www.jorek.eu/wiki/doku.php?id=plot_mlog.py
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
    
    plt.rc('font', size=12) 
    f, ax = plt.subplots(figsize=(6,4), dpi=160)
    lbl0 = None
    lbl1 = None
   
    legend=[]
    
    if type(x)!=list:
        x=[x]
        y=[y]
    if "legend" in kwargs: 
        legend=kwargs["legend"]
        if type(legend)!=list:
            legend=[legend]
        if len(legend)!=len(x): 
            return()

    for ip in range(len(x)):
        if "legend" in kwargs: 
            label = legend[ip]
        else:
            label=None
            
        ax.plot(x[ip],y[ip],linewidth= 1,label=label)

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
        ax.set_yscale("log")

    if len(legend)>0: ax.legend(loc='best', shadow=True)
    if "fname" in kwargs:
        fname = kwargs["fname"]+".png"
        f.savefig(fname,dpi=200,bbox_inches="tight")
    # plt.show()
    return(f,ax)

dirname=os.getcwd()

class C(object):
    pass
arg=C()

parser = argparse.ArgumentParser(description=__doc__);
parser.add_argument('-N', type=int, help="number of nodes");
parser.add_argument('-tline', type=str, help="print time line");
parser.add_argument('-legend', action='store_true', help="print legend");
 
parser.parse_args(namespace=arg)

nodes = True
if arg.N==None:
    nodes = False
    nN = 0
else:
    nN = arg.N
if arg.tline==None:
    tline = False
else:
    tline = True
    tfile = arg.tline
addLegend = arg.legend

fnames=[];
lst_dir=os.listdir(dirname)
pid=[]
mem_rank = True
for fid in lst_dir:
    if fid.startswith('mem_rank') and fid.endswith('.log'):
        dum = int(fid[9:-4])
        pid.append(dum)
    elif fid.startswith('mem_pid') and fid.endswith('.log'):
        mem_rank = False
        dum = int(fid[8:-4])
        pid.append(dum)


tstamp = []
event = []
if tline:
    if tfile in lst_dir:
        infile=open(tfile)
        for line in infile:
            line = line.split()
            if (line[0]=="#"): continue
            if len(line)>1:
                tstamp.append(float(line[0]))
                event.append(str(line[1]))

nP = len(pid)
tt = []
mm = []
legend=[]

Nlist = {}
for i in range(nP):
    tt.append([])
    mm.append([])
    legend.append("PID" + str(i))
    if (mem_rank):
        fid = "mem_rank_" + str(pid[i]) + ".log"
    else:
        fid = "mem_pid_" + str(pid[i]) + ".log"

    fid = os.path.join(dirname,fid)
    with open(fid) as infile:
        for line in infile:
            lsp = line.split()
            if (len(lsp)==0): continue
            if (lsp[0] == "#"):
                if (len(lsp)==5):
                    nodeid=int(lsp[2])
                    if nodeid not in Nlist.keys():
                        Nlist[nodeid]=[]
                    Nlist[nodeid].append(i)
                    continue
                else: continue
            tt[i].append(float(lsp[0]))
            mm[i].append(1e-6*float(lsp[1]))

t0 = []
for i in range(nP):
    tt[i]=np.array(tt[i])
    mm[i]=np.array(mm[i])
    t0.append(min(tt[i]))

t0 = min(t0)
for i in range(nP):
    tt[i] = (tt[i] - t0)/1000

if (False):
    # absolute
    f,ax=plot1d(tt,mm,xlbl="time (s)",ylbl="memory (GB)",legend=legend)
else:
    # cumulative
    mmc = []; ttc =[]
    mmc.append(mm[0]); ttc.append(tt[0])
    for i in range(1, nP):
        x = np.array(tt[0])     
        mmc.append(np.interp(x,tt[i],mm[i]) + mmc[i-1])
        ttc.append(x)
    if (addLegend):
        f,ax=plot1d(ttc,mmc,xlbl="time (s)",ylbl="memory (GB)",legend=legend)
    else:
        #f,ax=plot1d(ttc,mmc,xlbl="time (s)",ylbl="memory (GB)")
	#ax.plot(ttc[-1],mmc[-1],'k',linewidth= 2)
        
        f,ax=plot1d(ttc[-1],mmc[-1],xlbl="time (s)",ylbl="memory (GB)")
        ax.properties()['children'][0].set_color('black')
        ax.properties()['children'][0].set(color='black',linewidth=2)
        print("Max value = ",max(mmc[-1]))
#    ax.set_ylim([0,400])
#    ax.set_ylim([0,2000])


ymax = ax.get_ylim()[-1]
for i in range(len(tstamp)):
    x=np.ones(10)*(tstamp[i]-t0)/1000
    y=np.linspace(0,ymax,10)
    ax.plot(x,y,"k--",linewidth=1)
    ax.text(x[0],0.98*y[-1], event[i], fontsize=10, rotation=90,verticalalignment='top')

f.savefig("memhist.png",dpi=200,bbox_inches="tight")
f.savefig("memhist.eps",bbox_inches="tight")

# memory per node
if nodes:
    # distributes pid between nodes
    nR = int(nP/nN)
    for i in range(nN):
        Nlist[i]=[]
        for j in range(i*nR,(i+1)*nR):
            Nlist[i].append(j)
else:
    nN = len(Nlist.keys())

print(Nlist)

if nN>1:
    mN = []; tN=[]; lN=[]
    for i in range(nN):
        mmax= np.zeros(len(tt[0]))
        mN.append(np.zeros(len(tt[0])))
        tN.append(tt[0])
        lN.append("Node {}".format(i))
        for j in Nlist[i]:
            # mmax = np.maximum(mmax,np.interp(tt[0],tt[j],mm[j]))
            mN[i][:] += np.interp(tt[0],tt[j],mm[j])
        # mN[i][:] = mmax*len(Nlist[i])

    if (addLegend):
        f1,ax1=plot1d(tN,mN,xlbl="time (s)",ylbl="memory (GB)",legend=lN)
    else:
        f1,ax1=plot1d(tN,mN,xlbl="time (s)",ylbl="memory (GB)")
    #ax1.set_ylim([0,140])
    f1.savefig("mpernode.png",dpi=200,bbox_inches="tight")

plt.show()

exit()
