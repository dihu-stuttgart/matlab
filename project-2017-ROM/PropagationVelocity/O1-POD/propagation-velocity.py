import math
import numpy as np
import sys
import os

#folder
folder="./"

#get the gridsize and time step of the two files for which the propagation velocity should be computed

if len(sys.argv) >1:
   num_files=len(sys.argv)-1
   n_elem=sys.argv[1]
   timestep=sys.argv[2] 
else:
   # if no command line arguments are given
   n_elem=2048
   timestep=0.0005
###############################################################################
def read_outVm_file(filename):
    with open (filename,'r') as f:
      lines=f.read().splitlines()
      Vm=np.array(map(float,lines))
      #print "Vm: ", Vm
    return Vm 
###############################################################################
def node_indx(filename):
    if "_node_" in filename:
       pos_s=filename.find("_node_")+6
       #print "pos_s: ",pos_s
    if "_dt_" in filename:
       pos_e=filename.find("_dt_")
       #print "pos_e: ",pos_e
    indx=int(filename[pos_s:pos_e])
    print "indx: ",indx
    return indx
###############################################################################
def PropagationVelocity(Vm_1,Vm_2,timestep,dL):
    Vm_1_max=np.amax(Vm_1) 
    print "Vm_1_max:", Vm_1_max
    time_1=np.argmax(Vm_1)*timestep
    print "time_1:", time_1
    Vm_2_max=np.amax(Vm_2) 
    print "Vm_2_max:", Vm_2_max 
    time_2=np.argmax(Vm_2)*timestep
    print "time_2 :", time_2
    
    vel=dL*0.01/(abs(time_1-time_2)*0.001) #velocity in (m/s)
    print "velocity: ",vel
    return vel
###############################################################################
ls=os.listdir(folder)
print ls
phrase1="out_Vm"
phrase2="_dt_"+str(timestep)
print phrase2
phrase3="_n_"+str(n_elem)
print phrase3 
condition=lambda name: phrase1 in name and phrase2 in name and phrase3 in name and ".txt" in name
ls_files=list(np.extract(map(condition,ls),ls))
num_files=len(ls_files)
print "num_files: ",num_files
print "ls_files",ls_files
if (num_files > 2):
   print "more than 2 files are given to extract the propagation velocity."
else:
   Vm_1=read_outVm_file(ls_files[0])
   indx_1=node_indx(ls_files[0])
   Vm_2=read_outVm_file(ls_files[1])
   indx_2=node_indx(ls_files[1])
   L=1
   dL=(abs(indx_1-indx_2)-1)/float(n_elem)*L
   print "dL: ",dL
   PropagationVelocity(Vm_1,Vm_2,float(timestep),dL)
