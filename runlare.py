#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for running multiple instances of lare2d, varying certain parameters

"""


import os
import shutil
import numpy as np

nruns = 1   #number of runs
res = 512   #resolution in the x direction (y direction half of this so should be even)
nplots = 500  #number of outputs
output_directory = './Data/'    #this needs altering in control.f90 as well
tmax = 250 #maximum time in days
ncores = 1   #number of MPI cores. Too many cores at low resolution causes the SDF module to go a bit screwy sometimes. This appears to be a bug with lare, not my stuff...


#INITIAL CONDITIONS
energy = 1.0   #initial internal energy
vout = 1.0  #initial outflow speed
grav = 1.0   #gravitational acceleration
bfact = 1.0  #initial magnetic field strength factor
density = 1.0   #initial fluid density

shearfact = 0.2   #maximum shearing velocity
vmin = 5e-3
vmax = 1e-2
eta0 = np.geomspace(vmin, vmax, nruns)*shearfact   #supergranular diffusion rate (variable)

var_array = np.ones((nruns, 11))
var_array[:,0] = np.linspace(0,nruns-1,nruns)
var_array[:,1] = res
var_array[:,3] = nplots
var_array[:,4] = energy
var_array[:,5] = vout
var_array[:,6] = bfact
var_array[:,7] = density
var_array[:,8] = shearfact
var_array[:,9] = eta0
var_array[:,10] = grav

var_array[:,2] = tmax/(24.2*var_array[:,7])

shutil.rmtree(output_directory)
os.mkdir(output_directory)

print('Initial conditions set up, compiling and running...')
for run in range(0,nruns):
    if os.path.isdir('%s%d/' % (output_directory,run)):
        shutil.rmtree('%s%d/' % (output_directory,run))

    os.mkdir('%s%d/' % (output_directory,run))

    var = var_array[run]

    if run == 0:  #only compile the first time
        os.system('make COMPILER=gfortran')
    np.savetxt('variables.txt', var)
    os.system('/usr/lib64/openmpi/bin/mpiexec -np %d ./bin/lare2d' % ncores)
    os.system('rm variables.txt')
