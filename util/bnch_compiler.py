#!/usr/bin/env python3
import os
import sys
import argparse

# usage example
# ./util/bnch_compiler.py -input <path_to_input> -restart <path_to_restart>


def generate_job_script(dirname, nN, nR, nT, tag, jorek_name):
    fname = os.path.join(dirname,"jscipt.job")
    with open(fname,'w') as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=jrk_{}\n".format(tag))
        f.write("#SBATCH --nodes={:d}\n".format(nN))
        f.write("#SBATCH --ntasks={:d}\n".format(nR))
        f.write("#SBATCH --cpus-per-task={:d}\n".format(nT))
        f.write("#SBATCH --time=00:10:00\n")
        f.write("#SBATCH --account=FUSIO_HLST_1\n")
        f.write("#SBATCH --partition=skl_fua_prod\n")
        f.write("\n")
        f.write("export OMP_NUM_THREADS={:d}\n".format(nT))
        f.write("\n")
        f.write("mpirun ./{} < jorek.in | tee logfile.out \n".format(jorek_name))

    f.close()

def generate_log_file(dirname,message):
    fname = os.path.join(dirname,"test.log")
    with open(fname,'w') as f:
        f.write(message)
    f.close()    

def p2f(a):
    """ convert to fortran data type string """
    if type(a) is bool:
        if a:
            s = ".t."
        else:  s = ".f."
    else:
        s = str(a)
    return(s)

def generate_input_file(dirname,**kwargs):
    fname = os.path.join(dirname,"jorek.in")
    with open(fname,'w') as f:
        f.write("&in1\n")

        for var in kwargs.keys():
            f.write(f"{var} = {p2f(kwargs[var])}\n")
            
        f.write("/\n")
    f.close()
    
    
####### start ######
os.system("which python3")

jorek_dir = os.getcwd()

class C(object):
    pass
arg=C()

parser = argparse.ArgumentParser(description=__doc__);
parser.add_argument('-input', type=str, help="input file name");
parser.add_argument('-restart', type=str, help="restart file name");
parser.add_argument('-check', dest='check', action='store_true', help="option for retrieving results")
parser.add_argument('-testing', dest='testing', action='store_true', help="option for dry run")


parser.parse_args(namespace=arg)

testing = arg.testing
check = arg.check
if check: testing = True


if arg.input!=None:
    fid = os.path.join(jorek_dir,arg.input)
else:
    fid=os.path.join(jorek_dir,'input')
    
if arg.restart!=None:
    frestart = os.path.join(jorek_dir,arg.restart)
    if (os.path.isfile(frestart)):
        restart = True
    else:
        print("Restart file not found")
        exit()
else:
    restart = False
    
print("input file: {}".format(fid))

# dictionary for jorek inputs
jorek_input = {}


with open(fid,encoding = "ISO-8859-1") as infile:
    for line in infile:
        x=line.find("!")
        if (x!=-1):
            # print("Replacing {} wtih {}".format(line,line[0:x]))
            line = line[0:x]
            
        line = line.strip().split()
        if len(line)>1:
            jorek_input[line[0]] = line[-1]    

jorek_defaults = {**jorek_input}

kwargsJorek = {**jorek_defaults}

# basic compiler flags
flags0  = [" FLAGS=\"-qopenmp "]
fflags0 = [" FFLAGS=\"-diag-disable 8889 -warn all -warn nointerfaces -warn nounused -fpp -r8 -module .mod -DMPI_VERSION=3 "]

build_options = [\
                 flags0[0] + "-O2 -march=skylake" + "\"" + fflags0[0] + "-vecabi=compat -mcmodel=medium -align " + "\"" \
                # , " USE_PASTIX=1 USE_STRUMPACK=0 ",\
                # , " USE_PASTIX=0 USE_STRUMPACK=1 ",\
                # , flags0[0] + "-O2 -xCORE-AVX512" + "\"" + fflags0[0] + "-vecabi=compat -mcmodel=medium -align " + "\"" \
                ]

# set [nN,nR,nT]
run_options = [[1,2,24],[1,4,12],[1,8,6]]

# dictionary for input arguments
input_options = {"use_strumpack":True, "use_pastix":True}

# execute configuration script
jorek_name = "jorek_model303"
cmd = "./util/config.sh model=303 n_tor=3 n_period=1 n_plane=4"
os.system(cmd) if (not testing)  else print(cmd)


benchmark_dir = os.path.join(jorek_dir,"bnch_compiler")
if not os.path.isdir(benchmark_dir):
    os.makedirs(benchmark_dir)


# go through build_options   
for build_id in range(len(build_options)):
# compile jorek       
    os.chdir(jorek_dir)
    cmd = "make cleanall; make -j" + build_options[build_id] #+ " printsettings"
    os.system(cmd) if (not testing)  else print(cmd)
    
    # go through run_options
    input_id = 0
    for input in input_options:
        for run_id in range(len(run_options)):
            
            tag = "build_" + str(build_id) + "_input_" + str(input_id) + "_run_" + str(run_id)
            
            test_dir = os.path.join(benchmark_dir, tag)
    
            if os.path.isdir(test_dir):
                cmd = "rm -rf " + test_dir
                os.system(cmd)
            if not os.path.isdir(test_dir):
                os.makedirs(test_dir)
                
            cmd = "cp " + os.path.join(jorek_dir,jorek_name) + " " + test_dir
            os.system(cmd) if (not testing)  else print(cmd)
            
            if restart:
                cmd = "cp " + os.path.join(jorek_dir,frestart) + " " + os.path.join(test_dir,"jorek_restart.h5")
                os.system(cmd) if (not testing)  else print(cmd)
                kwargsJorek["restart"] = True
            else:
                kwargsJorek["restart"] = False
        
            if ("USE_PASTIX=0 USE_STRUMPACK=1" in build_options[build_id]):
                kwargsJorek["use_pastix"] = False
                kwargsJorek["use_strumpack"] = True
                kwargsJorek["centralize_harm_mat"] = False
            if ("USE_PASTIX=1 USE_STRUMPACK=0" in build_options[build_id]):
                kwargsJorek["use_pastix"] = True
                kwargsJorek["use_strumpack"] = False
            
            kwargsJorek[input] = input_options[input]
                
            [nN,nR,nT] = run_options[run_id]
            
            generate_job_script(test_dir,nN,nR,nT,tag,jorek_name)
            
            generate_log_file(test_dir,str(build_options[build_id]))
            
            generate_input_file(test_dir,**kwargsJorek)
            
            os.chdir(test_dir)
            cmd = "sbatch jscipt.job"
            os.system(cmd) if (not testing)  else print(cmd)
            
            if check:
                cmd = os.path.join(jorek_dir,"util","read_jorek_logfile.py") + " -fname logfile.out -text \"ITERATION\""
                print(cmd)
                # os.system(cmd)
            
        input_id += 1
