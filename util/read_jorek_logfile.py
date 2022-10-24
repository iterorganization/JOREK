#!/usr/bin/env python3
"""find variable in JOREK logfile"""
import numpy as np
import time
import os
import argparse

mydir=os.getcwd()

class C(object):
    pass
arg=C()

parser = argparse.ArgumentParser(description=__doc__);
parser.add_argument('-fname', type=str, help="input file name");
parser.add_argument('-text', type=str, help="string to find");
parser.add_argument('-sum', action='store_true', help="Sum of entries");
parser.add_argument('-n', type=int, help="max number of entries");
parser.add_argument('-e0', type=int, help="first entry");
parser.add_argument('-e1', type=int, help="last entry");
 
parser.parse_args(namespace=arg)

if arg.fname!=None:
	fid = os.path.join(mydir,arg.fname)
else:
	fid=os.path.join(mydir,'logfile.out')
	
if arg.text!=None:
	sstr = arg.text
else:
	sstr= 'gmres/solve'

if arg.n!=None:
	nmax = arg.n
else:
	nmax = 1

if arg.e0!=None:
        e0 = max(arg.e0,0)
else:
        e0 = 0

if arg.e1!=None:
        e1 = max(arg.e1,0)
else:
        e1 = -1        
        

print(fid)
print(sstr)

dat1 = []
with open(fid,encoding = "ISO-8859-1") as infile:
	for line in infile:
		if sstr in line:
			line = line.strip().split()
			dat1.append(float(line[-1]))
			print(line)

dat1 = np.array(dat1)
print("Number of etries found: {}".format(len(dat1)))
if (nmax==1) :
	nmax = len(dat1)
else:
	nmax = min(nmax,len(dat1))
if (e1>0):
    nmax = e1

print(dat1[:nmax])
if (arg.sum):
        print("Sum = {}".format(np.sum(dat1[e0:nmax])))
else:
        print("Average = {}".format(np.mean(dat1[e0:nmax])))
        print("Maximum = {}".format(np.max(dat1[e0:nmax])))
        print("Minimum = {}".format(np.min(dat1[e0:nmax])))

exit()
