#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27
"""

import os
import shutil
import numpy as np
import sdf
import matplotlib.pyplot as plt
import time

nsnaps = 100
run_number = 0

if os.path.isdir('./plots'):
    for i in range(nsnaps):
        if os.path.isfile('./plots/%04d.png' % i):
            os.remove('./plots/%04d.png' % i)
else:
    os.mkdir('./plots')

go = True
allvmax = [-1e6,-1e6,-1e6,-1e6]
allvmin = [1e6,1e6,1e6,1e6]

allvs = []

allofs = []
allens = []
alljs = []
allbs = []
ts = []

#determine maxima and minima for plotting
allvmaxs = [[],[],[],[],[],[]]
allvmins = [[],[],[],[],[],[]]

data_directory = './Data/'

for i in range(0, nsnaps-1, 1):   #Establish bounds for the colourmaps
    if os.path.isfile('%s%d/%04d.sdf' % (data_directory, run_number, i+1)) :
        dat = sdf.read('%s%d/%04d.sdf' % (data_directory, run_number, i))

        t = dat.__dict__['Header']['time']

        # Read in coordinates:
        grid = dat.__dict__['Grid_Grid']
        x1 = grid.data[0]
        y1 = grid.data[1]

        # Grid spacing:
        dx = np.mean(x1[1:] - x1[:-1])
        dy = np.mean(y1[1:] - y1[:-1])

        # Grid sizes:
        nx = np.size(x1, 0)
        ny = np.size(y1, 0)


        # Generate cell-centre coordinate arrays, including ghost cells:
        xc = np.linspace(x1[0]-0.5*dx, x1[-1]+0.5*dx, nx+1)
        yc = np.linspace(y1[0]-0.5*dy, y1[-1]+0.5*dy, ny+1)

        xs = np.linspace(x1[0], x1[-1], nx)
        ys = np.linspace(y1[0], y1[-1], ny)


        # Read density:
        rho = (dat.__dict__['Fluid_Rho']).data
        bx = (dat.__dict__['Magnetic_Field_Bx']).data
        by = (dat.__dict__['Magnetic_Field_By']).data
        bz = (dat.__dict__['Magnetic_Field_Bz']).data
        vy = (dat.__dict__['Velocity_Vy']).data
        en = (dat.__dict__['Fluid_Energy']).data
        jz = (dat.__dict__['Current_Jz']).data

        def magfield(bx, by, bz):
            bx1 = 0.5*(bx[1:,:] + bx[:-1,:])
            by1 = 0.5*(by[:,1:] + by[:,:-1])
            bz1 = bz
            return (bx1**2 + by1**2+ bz1**2)


        toplot = [vy, magfield(bx,by,bz), np.abs(jz), np.abs(bz), rho, en ]

        for k in range(len(toplot)):
            allvmaxs[k].append(np.percentile(toplot[k], 95))
            allvmins[k].append(np.percentile(toplot[k], 5))

allvmaxs = np.array(allvmaxs)
allvmins = np.array(allvmins)

vmaxs = []
vmins = []
for k in range(len(allvmaxs)):
    vmaxs.append(np.percentile(allvmaxs[k], 95))
    vmins.append(np.percentile(allvmaxs[k], 5))

vmaxs[0] = max(vmaxs[0], -vmins[0])
vmins[0] = -vmaxs[0]

vmins[1] = 0.0; vmins[2] = 0.0; vmins[3] = 0.0; vmins[4] = 0.0; vmins[5] = 0.0

while go:
    for i in range(0,nsnaps-1,1):
        if os.path.isfile('%s%d/%04d.sdf' % (data_directory, run_number, i+1)) :
            if not os.path.isfile('./plots/%04d.png' % i):
                print('Making plot', i, 'from run number', run_number)
                if i == nsnaps - 2:
                    go = False
                fig_width = 20
                dat = sdf.read('%s%d/%04d.sdf' % (data_directory, run_number, i))

                t = dat.__dict__['Header']['time']
                
                # Read in coordinates:
                grid = dat.__dict__['Grid_Grid']
                x1 = grid.data[0]
                y1 = grid.data[1]
            
                # Grid spacing:
                dx = np.mean(x1[1:] - x1[:-1])
                dy = np.mean(y1[1:] - y1[:-1])

                # Grid sizes:
                nx = np.size(x1, 0)
                ny = np.size(y1, 0)
                

                # Generate cell-centre coordinate arrays, including ghost cells:
                xc = np.linspace(x1[0]-0.5*dx, x1[-1]+0.5*dx, nx+1)
                yc = np.linspace(y1[0]-0.5*dy, y1[-1]+0.5*dy, ny+1)
            
                xs = np.linspace(x1[0], x1[-1], nx)
                ys = np.linspace(y1[0], y1[-1], ny)
            

                # Read density:
                rho = (dat.__dict__['Fluid_Rho']).data
                bx = (dat.__dict__['Magnetic_Field_Bx']).data
                by = (dat.__dict__['Magnetic_Field_By']).data
                bz = (dat.__dict__['Magnetic_Field_Bz']).data
                vy = (dat.__dict__['Velocity_Vy']).data
                en = (dat.__dict__['Fluid_Energy']).data
                jz = (dat.__dict__['Current_Jz']).data

                def divergence(bx, by):
                    div = np.zeros((nx-1, ny-1))
                    div = bx[1:,:]*dy - bx[:-1,:]*dy + by[:,1:]*dx - by[:,:-1]*dx
                    return div
                div = divergence(bx, by)



                def openflux(by):
                    return np.sum(np.abs(by[-1]))*dx

                def energy(en):
                    return np.sum(np.abs(en[:]))*dx*dy

                def current(jz):
                    return np.sum(np.abs(jz[:]))*dx*dy
                
                def magnetic_energy(bx, by, bz):
                    bx1 = 0.5*(bx[1:,:] + bx[:-1,:])
                    by1 = 0.5*(by[:,1:] + by[:,:-1])
                    bz1 = bz
                    #print(np.shape(bz1))

                    #print(np.sum(bx1**2))
                    #print(np.sum(by1**2))
                    #print(np.sum(bz1**2))
                    return 0.5*(np.sum(bx1**2 + by1**2+bz1**2))*dx*dy

                def magfield(bx, by, bz):
                    bx1 = 0.5*(bx[1:,:] + bx[:-1,:])
                    by1 = 0.5*(by[:,1:] + by[:,:-1])
                    bz1 = bz
                    return (bx1**2 + by1**2+ bz1**2)

                allofs.append(openflux(by))
                #allens.append(energy(en))
                allens.append(magnetic_energy(bx, by, bz))
                alljs.append(current(jz))
                allbs.append(np.max(np.abs(bz)))
                ts.append(t)
                if True:
                    fig = plt.figure(figsize = (fig_width,0.6*fig_width*(yc[-1] - yc[0])/(xc[-1] - xc[0])))
                    gs = fig.add_gridspec(4,3)

                    allvs.append(vy[0])
                    toplot = [vy, magfield(bx,by,bz), np.abs(jz), np.abs(bz)]
                    allvmin[0] = 0.0; allvmax[0] = 1.5
                    allvmin[1] = 0.0; allvmax[1] = 0.00025
                    allvmin[2] = 0.0; allvmax[2] = 80.0

                    for k in range(4):
                        allvmin[k] = min(allvmin[k], np.min(toplot[k]))
                        allvmax[k] = max(allvmax[k], np.max(toplot[k]))
                    titles = ['Vy', 'Magnetic Energy', 'J_z']
                    for k in range(3):
                        ax = plt.subplot(gs[2*(k//2):2*(k//2)+2,k%2])
                        if np.max(toplot[k]) > allvmax[k]:
                            allvmax[k] = np.max(toplot[k])
                        if np.min(toplot[k]) < allvmin[k]:
                            allvmin[k] = np.min(toplot[k])
                        #allvmin[2] = 0.; allvmin[3] = 0.
                        #allvmax[2] = 1.; allvmax[3] = 1.

                        im = ax.pcolormesh(xs, ys, toplot[k].T, cmap = 'plasma', rasterized=True, vmin = vmins[k], vmax = vmaxs[k])
                        fig.colorbar(im)
                        ax.set_title(titles[k])



                    def find_a(bx, by):   #find the vector potential a from the (hopefully) divergence-free magnetic fields
                        a = np.zeros((nx, ny))
                        for j in range(1,ny):
                            a[0,j] = a[0,j-1] -dy*bx[0,j-1]

                        for i in range(1,nx):
                            a[i,0] = a[i-1,0] + dx*by[i-1,0]
                            for j in range(1,ny):
                                a[i,j] = a[i,j-1] -dy*bx[i,j-1]

                        return a
                    
                    
                    a = find_a(bx, by)
                    
                    if True:
                        ax = plt.subplot(gs[2:4,1])
                        im = ax.pcolormesh(xs, ys, bz.T, cmap = 'plasma',rasterized=True,vmin = vmins[3], vmax = vmaxs[3])
                        #im = ax.pcolormesh(xs, ys, a.T, cmap = 'seismic',rasterized=True)#,vmin = allvmin[k], vmax = allvmax[k])
                        if np.max(a) > 1e-6:
                            
                            ax.contour(xs, ys, a.T, 64, colors = 'black', linewidths = 0.5)
                        fig.colorbar(im)
                        ax.set_title('Magnetic Field Into Plane')
                    
                    if False:
                        ax = plt.subplot(gs[2:4,1])
                        #im = ax.pcolormesh(xs, ys, .T, cmap = 'Reds',rasterized=True)#,vmin = allvmin[k], vmax = allvmax[k])
                        im = ax.pcolormesh(xs, ys, a.T, cmap = 'Reds',rasterized=True)#,vmin = allvmin[k], vmax = allvmax[k])
                        if np.max(a) > 1e-6:
                            ax.contour(xs, ys, a.T, 64, colors = 'black', linewidths = 0.5)
                        fig.colorbar(im)
                        ax.set_title('Magnetic Field Vector Potential')
                    
                    if False:
                        ax = plt.subplot(gs[2:4,1])
                        #im = ax.pcolormesh(xs, ys, .T, cmap = 'Reds',rasterized=True)#,vmin = allvmin[k], vmax = allvmax[k])
                        im = ax.pcolormesh(xc[1:-1], yc[1:-1], div.T, cmap = 'Reds',rasterized=True)#,vmin = allvmin[k], vmax = allvmax[k])

                        fig.colorbar(im)
                        ax.set_title('Magnetic Field Divergence')



                    ax = plt.subplot(gs[:2,2])
                    im = ax.pcolormesh(xs, ys, rho.T, cmap = 'plasma',rasterized=True,vmin = vmins[4], vmax = vmaxs[4])
                    fig.colorbar(im)
                    ax.set_title('Density')

                    ax = plt.subplot(gs[2:4,2])
                    im = ax.pcolormesh(xs, ys, en.T, cmap = 'plasma',rasterized=True,vmin = vmins[5], vmax = vmaxs[5])
                    fig.colorbar(im)
                    ax.set_title('Internal Energy')

                    
                    plt.tight_layout()
                    plt.savefig('./plots/%04d.png' % i, bbox_inches='tight')
                    #plt.show()
                    plt.close()

    time.sleep(1.)
