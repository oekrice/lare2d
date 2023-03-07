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
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
import random
import shutil
import pfield_2d as pfield

class diagnostics():
    def __init__(self, run, data_root):
        #Can probably do something with this and do it properly. I do have many hours.

        self.diag_titles = ['Time', 'Isrope', 'Poloidal Rope Flux', 'Rope Height', 'Magnetic Energy', 'Internal Energy', 'Maximum Sheared Field', 'Open Flux', 'Average Current', 'Maximum Current', 'Axial Rope Flux', 'Axial Rope Current', 'Helicity', 'Rope Area', 'Current-Carrying Helicity 1', 'Relative Helicity 1', 'Free Energy 1', 'Current-Carrying Helicity 2', 'Relative Helicity 2', 'Free Energy 2']
        self.diags = []
        #Diagnostic titles are extended to include some of the stuff that only applies to full MHD
        self.data_root = data_root
        self.shear_fact = 0.2
        self.do_potential = True
        self.A_prev = 0.0
        self.b_prev = 0.0
        i = 0
        self.xbasis = []
        self.ybasis = []
        self.m2 = []
        for i in range(1000):
        #while i < 1000:
            if os.path.isfile('%s%d/%04d.sdf' % (self.data_root, run, i + 1)) :
                print(i)
                dat = sdf.read('%s%d/%04d.sdf' % (self.data_root, run, i))

                self.t = dat.__dict__['Header']['time']*(24.2*self.shear_fact)

                # Read in coordinates:
                self.grid = dat.__dict__['Grid_Grid']
                self.x1 = self.grid.data[0]
                self.y1 = self.grid.data[1]

                # Grid spacing:
                self.dx = np.mean(self.x1[1:] - self.x1[:-1])
                self.dy = np.mean(self.y1[1:] - self.y1[:-1])

                # Grid sizes:
                self.nx = np.size(self.x1, 0)
                self.ny = np.size(self.y1, 0)
                # Generate cell-centre coordinate arrays, including ghost cells:
                self.xc = np.linspace(self.x1[0]-0.5*self.dx, self.x1[-1]+0.5*self.dx, self.nx+1)
                self.yc = np.linspace(self.y1[0]-0.5*self.dy, self.y1[-1]+0.5*self.dy, self.ny+1)

                self.xs = np.linspace(self.x1[0], self.x1[-1], self.nx)
                self.ys = np.linspace(self.y1[0], self.y1[-1], self.ny)

                # Read density:
                self.rho = (dat.__dict__['Fluid_Rho']).data
                self.en = (dat.__dict__['Fluid_Energy']).data

                self.bx = (dat.__dict__['Magnetic_Field_Bx']).data
                self.by = (dat.__dict__['Magnetic_Field_By']).data
                self.bz = (dat.__dict__['Magnetic_Field_Bz']).data

                self.ax = self.find_ax(self.bz)
                self.ay = np.zeros((self.nx, self.ny-1))
                self.az = self.find_az(self.bx, self.by)

                self.jx, self.jy, self.jz = self.calculate_current()

                self.rope_height = self.findrope()

                if self.rope_height == 0 or self.rope_height == self.ys[-1]:
                    self.isrope = 0
                else:
                    self.isrope = 1
                self.check_divergence()   #checks the domain to be divergence-free. If it's not, the simulations are nonsense
                self.check_magfield()   #checks that the potential field integration has been fine

                if self.do_potential:
                    self.bxp, self.byp, self.bzp = self.find_potential_field()
                else:
                    self.bxp, self.byp, self.bzp = 0.0*self.bx, 0.0*self.by, 0.0*self.bz

                self.axp = self.find_ax(self.bzp)
                self.ayp = np.zeros((self.nx, self.ny-1))
                self.azp = self.find_az(self.bxp, self.byp)

                self.helicity_calculations()

                self.ropemask = self.az < self.az[self.nx//2,0]

                self.rope_area = np.sum(self.ropemask)*self.dx*self.dy
                self.update_diagnostic_data()
                #i = i + 1
            #else:
                #time.sleep(1.0)
            if i%10 == 0:
                np.save('diags/diags%d.npy' % run, self.diags, allow_pickle = True)
                np.save('diags/diagtitles.npy', self.diag_titles)
        #plt.plot(self.diags[0], self.diags[3])
        #plt.show()

    def find_potential_field(self):
        lbound = self.by[:,0]
        ubound = self.by[:,-1]
        #side boundaries are zero (hopefully) - will check this is definitely true
        bxp = 0*self.bx
        byp = 0*self.by
        bzp = 0*self.bz

        pfield_data = pfield.compute_pfield(self.grid, ubound, lbound, xbasis = self.xbasis, ybasis = self.ybasis, m2 = self.m2, max_mode = 1e6)
        bxp = pfield_data.bxp
        byp = pfield_data.byp
        bzp = 0.0*self.bz
        bzp[:,:] = (np.sum(self.bz)/np.size(self.bz))

        if len(self.xbasis) == 0:
            self.xbasis = pfield_data.xbasis
            self.ybasis = pfield_data.ybasis
            self.m2 = pfield_data.m2

        return bxp, byp, bzp

    def helicity_calculations(self):
        #I think it's best to do all the helicity calclulations in here. B and A require averaging first
        def avga(ax, ay, az):
            ax0 = 0.5*(ax[:,1:] + ax[:,:-1])
            ay0 = 0.5*(ay[1:,:] + ay[:-1,:])
            az0 = 0.25*(az[1:,1:] + az[:-1,1:] + az[1:,:-1] + az[:-1,:-1])
            return ax0, ay0, az0

        def avgb(bx, by, bz):
            bx0 = 0.5*(bx[1:,:] + bx[:-1,:])
            by0 = 0.5*(by[:,1:] + by[:,:-1])
            bz0 = bz[:,:]
            return bx0, by0, bz0

        ax1, ay1, az1 = avga(self.ax, self.ay, self.az)
        axp1, ayp1, azp1 = avga(self.axp, self.ayp, self.azp)

        bx1, by1, bz1 = avgb(self.bx, self.by, self.bz)
        bxp1, byp1, bzp1 = avgb(self.bxp, self.byp, self.bzp)

        self.helicity = -np.sum(ax1*bx1 + ay1*by1 + az1*bz1)*self.dx*self.dy
        #Let helicity 1 be the definition from the original 2D parameter study. helicity 2 be the one with no bz component, as in the 3D (probably?!) and the axisymmetric polar runs
        self.hr_1 = -np.sum((ax1 + axp1)*(bx1 - bxp1) + (ay1 + ayp1)*(by1 - byp1) + (az1 + azp1)*(bz1 - bzp1))*self.dx*self.dy
        self.hj_1 = -np.sum((ax1 - axp1)*(bx1 - bxp1) + (ay1 - ayp1)*(by1 - byp1) + (az1 - azp1)*(bz1 - bzp1))*self.dx*self.dy

        self.hr_2 = -np.sum((ax1 + axp1)*(bx1 - bxp1) + (ay1 + ayp1)*(by1 - byp1) + (az1 + azp1)*(bz1))*self.dx*self.dy
        self.hj_2 = -np.sum((ax1 - axp1)*(bx1 - bxp1) + (ay1 - ayp1)*(by1 - byp1) + (az1 - azp1)*(bz1))*self.dx*self.dy

        self.free_energy_1 = 0.5*(np.sum(bx1**2 + by1**2+bz1**2) - np.sum(bxp1**2 + byp1**2+bzp1**2))*self.dx*self.dy
        self.free_energy_2 = 0.5*(np.sum(bx1**2 + by1**2+bz1**2) - np.sum(bxp1**2 + byp1**2))*self.dx*self.dy


        plt.show()

    def calculate_current(self):
        jx = 0*self.ax
        jy = 0*self.ay
        jz = 0*self.az
        jx[:,1:-1] = (self.bz[:,1:] - self.bz[:,:-1])/self.dy
        jy[1:-1,:] = -(self.bz[1:,:] - self.bz[:-1,:])/self.dx
        jz[1:-1,1:-1] = (self.by[1:,1:-1] - self.by[:-1,1:-1])/self.dx - (self.bx[1:-1,1:] - self.bx[1:-1,:-1])/self.dy
        return jx, jy, jz

    def update_diagnostic_data(self):
        #Collates all the data that has been collected and appends to the big alldiags file (in the same format as the axisymmetric one)
        if len(self.diags) == 0:
            ndiags = len(self.diag_titles)
            for i in range(ndiags):
                self.diags.append([])
        for n in range(len(self.diags)):
            if self.diag_titles[n] == 'Time':
                self.diags[n].append(self.t)
            if self.diag_titles[n] == 'Isrope':
                self.diags[n].append(self.isrope)
            if self.diag_titles[n] == 'Poloidal Rope Flux':
                self.diags[n].append(self.ropeflux())
            if self.diag_titles[n] == 'Rope Height':
                self.diags[n].append(self.findrope())
            if self.diag_titles[n] == 'Magnetic Energy':
                self.diags[n].append(self.magnetic_energy())
            if self.diag_titles[n] == 'Internal Energy':
                self.diags[n].append(self.internal_energy())
            if self.diag_titles[n] == 'Maximum Sheared Field':
                self.diags[n].append(np.max(np.abs(self.bz)))
            if self.diag_titles[n] == 'Open Flux':
                self.diags[n].append(self.find_openflux())
            if self.diag_titles[n] == 'Average Current':
                self.diags[n].append(self.avgcurrent())
            if self.diag_titles[n] == 'Maximum Current':
                self.diags[n].append(self.maxcurrent())
            if self.diag_titles[n] == 'Axial Rope Flux':
                self.diags[n].append(self.axial_flux())
            if self.diag_titles[n] == 'Axial Rope Current':
                self.diags[n].append(self.axial_current())
            if self.diag_titles[n] == 'Helicity':
                self.diags[n].append(self.helicity)
            if self.diag_titles[n] == 'Free Energy 1':
                self.diags[n].append(self.free_energy_1)
            if self.diag_titles[n] == 'Free Energy 2':
                self.diags[n].append(self.free_energy_2)
            if self.diag_titles[n] == 'Relative Helicity 1':
                self.diags[n].append(self.hr_1)
            if self.diag_titles[n] == 'Relative Helicity 2':
                self.diags[n].append(self.hr_2)
            if self.diag_titles[n] == 'Current-Carrying Helicity 1':
                self.diags[n].append(self.hj_1)
            if self.diag_titles[n] == 'Current-Carrying Helicity 2':
                self.diags[n].append(self.hj_2)
            if self.diag_titles[n] == 'Rope Area':
                self.diags[n].append(self.rope_area)

        if False:
            for i in range(len(self.diags)):
                print(self.diag_titles[i], self.diags[i][-1])

    def check_divergence(self):
        div = np.zeros((self.nx-1, self.ny-1))
        div = self.bx[1:,:]*self.dy - self.bx[:-1,:]*self.dy + self.by[:,1:]*self.dx - self.by[:,:-1]*self.dx
        if np.max(np.abs(div)) > 1e-14:
            print('Not divergence-free')
        return div

    def find_az(self, bx, by):   #find the vector potential a from the (hopefully) divergence-free magnetic fields
        az = np.zeros((self.nx, self.ny))
        for j in range(1,self.ny):
            az[0,j] = az[0,j-1] + self.dy*bx[0,j-1]

        for i in range(1,self.nx):
            az[i,0] = az[i-1,0] - self.dx*by[i-1,0]
            for j in range(1,self.ny):
                az[i,j] = az[i,j-1] + self.dy*bx[i,j-1]
        return az

    def find_ax(self, bz):
        ax = np.zeros((self.nx-1, self.ny))
        for j in range(1,self.ny):
            ax[:,j] = ax[:,j-1] - self.dy*bz[:,j-1]
        return ax

    def check_magfield(self):
        #checks that the integrated vector potentials are actually accurate
        bxc = 0*self.bx
        byc = 0*self.by
        bzc = 0*self.bz
        bxc[:,:] = (self.az[:,1:] - self.az[:,:-1])/self.dy
        byc[:,:] = -(self.az[1:,:] - self.az[:-1,:])/self.dx
        bzc[:,:] = (self.ay[1:,:] - self.ay[:-1,:])/self.dx - (self.ax[:,1:] - self.ax[:,:-1])/self.dy

        if np.max(np.abs(self.bx - bxc)) > 1e-10:   #this is quite a high tolerance but I think it's close enough
            print('Integration wrong', np.max(np.abs(self.bx - bxc)))
        if np.max(np.abs(self.by - byc)) > 1e-10:
            print('Integration wrong', np.max(np.abs(self.by - byc)))
        if np.max(np.abs(self.by - byc)) > 1e-10:
            print('Integration wrong', np.max(np.abs(self.by - byc)))

    def findrope(self):

        aslice = self.az[self.nx//2,:]   #want to find highest point of this

        max_index = np.where(aslice == min(aslice))

        if len(max_index[0]) == 0:
            return 0.0
        else:
            max_index = max_index[0][0]
            if max_index == 0:
                return 0
            elif max_index == len(aslice) - 1:
                return self.ys[-1]
            else:
                pts = aslice[max_index-1:max_index+2]
                poly = np.polyfit(self.ys[max_index-1:max_index+2], aslice[max_index-1:max_index+2],deg = 2)
                xfit = np.linspace(self.ys[max_index-1], self.ys[max_index+1], 100)
                avals = np.polyval(poly, xfit)
                return xfit[np.where(avals == min(avals))[0]][0]

    def find_openflux(self):
        return np.sum(np.abs(self.by[-1]))*self.dx

    def internal_energy(self):
        return np.sum(np.abs(self.en[:]))*self.dx*self.dy

    def ropeflux(self):
        if self.isrope:
            return max(0, max(-self.az[self.nx//2,:] + self.az[self.nx//2,0]))
        else:
            return 0.0

    def magnetic_energy(self):
        bx1 = 0.5*(self.bx[1:,:] + self.bx[:-1,:])
        by1 = 0.5*(self.by[:,1:] + self.by[:,:-1])
        bz1 = self.bz
        return 0.5*(np.sum(bx1**2 + by1**2+bz1**2))*self.dx*self.dy

    def avgcurrent(self):
        return np.sum(np.abs(self.jz[:]))*self.dx*self.dy/(2.0)

    def maxcurrent(self):
        return np.max(np.abs(self.jz[:]))


    def axial_flux(self):
        bz1 = 0.25*(self.bz[1:,1:] + self.bz[:-1,1:] + self.bz[1:,:-1] + self.bz[:-1, :-1])
        ropeb = bz1[self.ropemask[1:-1,1:-1]]
        return np.sum(np.abs(ropeb))*self.dx*self.dy

    def axial_current(self):
        ropej = self.jz[self.ropemask]
        return np.sum(np.abs(ropej))*self.dx*self.dy
