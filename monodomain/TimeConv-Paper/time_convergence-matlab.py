import math
import numpy as np
from numpy import linalg as LA
import matplotlib
#matplotlib.use('Agg') for saving figure in eps format
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import sys
import os 
from scipy import stats

# name of the folder
folder="./"

# get the grid size for which the convergence rate should be computed
if len(sys.argv) >1:
   gridsize=sys.argv[1]
else:
   gridsize=2048

print "grid size: ", gridsize

###############################################################################
# extract time step from the name of files
def get_timestep(filename):
  #print filename
  j_end=len(filename)
  #print j_end
  for j in range(j_end):
    if filename[j]=='d':
      j_s=j+3
      #print "j_s",j_s
    if filename[j]=='n':
      j_e=j-1
      #print "j_e:",j_e
  
  #print "filename[j_s:j_e]",filename[j_s:j_e]
  timestep=float(filename[j_s:j_e])
  print "timestep",timestep
  return timestep
###############################################################################
# get the data for the time convergence study from the current folder 
def get_data_TimeConv(order):

  global folder
   
  ls = os.listdir(folder)

  phrase="O"+str(order)
  condition= lambda name: phrase in name and str(gridsize) in name
  ls_selected=list(np.extract(map(condition,ls),ls))
  ls_conv=sorted(ls_selected)
  print "order: ", order
  print "ls_conv:", ls_conv

  num_runs=len(ls_conv)
  print "No. of runs: ", num_runs
  
  #get the spatial coordinates
  #with open("Out_x_n_"+str(gridsize)+".txt") as f:  
   # lines=f.readlines()
   # x=np.array(lines)
    #print "x: ",x

  Vm=[0.0 for i in range(num_runs)]
  timesteps=[0.0 for i in range(num_runs)]

  for i in range(num_runs):
    #print "i: ",i
    #print"file name: ", ls_conv[i]
    timesteps[i]=get_timestep(ls_conv[i])
    with open(ls_conv[i],'r') as f:
      lines=f.read().splitlines()  
      Vm[i]=np.array(map(float,lines))
      #print "i, Vm: ",i,Vm[i]
    
  return Vm,timesteps
###############################################################################
def makeplot_TimeConv(x,y,labels,order):

   global folder
   imagename="{}".format(folder+"image-Time-Conv-n"+str(gridsize)+"-O"+str(order)+".png")
        
   fig=plt.figure(1)
   for i in range(len(y)):
      plt.plot(x,y[i],label=str(labels[i]))
   plt.title("Vm")
   plt.legend()
   plt.show()    
   fig.savefig(imagename)
   print "Figure saved to {}".format(imagename)
###############################################################################
# The smallest time step is considered as the reference
def TimeConv_err_L2_rel(Vm,Vm_ref):
  num_runs=len(Vm)  
  err=[0.0 for i in range(num_runs)]
  
  Vm_ref_L2=LA.norm(Vm_ref,2)
  #print Vm_ref_L2
  for i in range(num_runs):
    #print i
    err[i]=LA.norm(Vm[i]-Vm_ref,2)/Vm_ref_L2
    #print err
  return err
###############################################################################
def makeplot_TimeConv_err(x,y,labels):

  global folder
  #imagename="{}".format(folder+"errL2-Time-Conv-n"+str(gridsize)+".eps")
  imagename="{}".format(folder+"errL2-Time-Conv-n"+str(gridsize)+".png")    
  fig=plt.figure(1,figsize=[8,7])

  symbols=['ro','bo']
  sl=[1,2]
  sl_labels=["slope 1","slope 2"]
  sl_lines=['r:','b:']

  #intercepts=[0.8,2]
  intercepts=[1,2.35]
  for i in range(len(y)):
    plt.loglog(x[i],y[i],symbols[i],label=labels[i])
    #slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x[i]),np.log10(y[i]))
    #print "i, slope: ", i,slope
    #slope, intercept=np.polyfit(np.log10(x[i]),np.log10(y[i]),1)
    #print "i,slope,intercept: ", i,slope,intercept
    #abline_values=[pow(10,slope*j+intercept) for j  in np.log10(x[i])]
    #plt.loglog(x[i],abline_values,'k:')
    abline_values=[pow(10,sl[i]*j+intercepts[i]) for j  in np.log10(x[i])]
    plt.loglog(x[i],abline_values,sl_lines[i],label=sl_labels[i],linewidth=3)
  
  #plt.title("Vm")
  plt.xlabel('Time step, $dt$',{'fontsize':20})
  plt.xticks(size=18)
  plt.xlim(left=0.0001,right=1)
  plt.ylabel('Relative L2-norm error of $V_m$',{'fontsize':20})
  plt.yticks(size=18)
  plt.legend()
  plt.show()    
  #fig.savefig(imagename,format='eps')
  fig.savefig(imagename)
  print "Figure saved to {}".format(imagename)

################################################################################
Vm_1,timesteps_1=get_data_TimeConv(1) # read the data for the 1st order method

Vm,ts=get_data_TimeConv(2) # read the data for the 2nd order method
Vm_ref=Vm[0] # selecting the smallest timestep as the reference
num_run=len(Vm)
num_run_2=len(Vm)-1 
Vm_2=Vm[1:num_run]
timesteps_2=ts[1:num_run]
#print "time steps: ", timesteps_2

#makeplot_TimeConv(x,Vm_1,timesteps_1)
#makeplot_TimeConv(x,Vm,ts) # plot including the reference data

timesteps=[1.0 for i in range(2)]
timesteps[0]=timesteps_1
timesteps[1]=timesteps_2

err_L2=[1.0 for i in range(2)]
err_L2[0]=TimeConv_err_L2_rel(Vm_1,Vm_ref)
err_L2[1]=TimeConv_err_L2_rel(Vm_2,Vm_ref)
#print "err_L2", err_L2[1]

labels=[None for i in range(2)]
labels[0]="1st order"
labels[1]="2nd order"

makeplot_TimeConv_err(timesteps,err_L2,labels)




