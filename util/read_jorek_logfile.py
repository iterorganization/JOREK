#!/usr/bin/env python3
""" Search JOREK log file for lines containing specific string and 
    get average/sum for the values found at the end of such lines
    https://www.jorek.eu/wiki/doku.php?id=read_jorek_logfile.py
    iholod@ipp.mpg.de
"""

import numpy as np
import pylab as plt
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

print(dat1[:nmax])
if (arg.sum):
	print("Sum = {}".format(np.sum(dat1[:nmax])))
else:
	print("Average = {}".format(np.mean(dat1[:nmax])))

exit()
