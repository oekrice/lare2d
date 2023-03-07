import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.linalg import eigh_tridiagonal

nx = 256   #number of CELLS, not grid points
ny = 128


ubound = np.zeros(nx)
lbound = np.zeros(nx)


class compute_pfield():
    def __init__(self, grid, ubound, lbound, xbasis = [], ybasis = [], m2 = 0.0, max_mode = 1e6):
        #use the same grid notation as in the proper diagnostic calculator, even though it can be bodged for now
        self.xs = grid.data[0]
        self.ys = grid.data[1]

        self.nx = np.size(self.xs, 0)
        self.ny = np.size(self.ys, 0)

        self.xc = np.zeros(self.nx + 1)
        self.yc = np.zeros(self.ny + 1)

        self.xc[1:-1] = 0.5*(self.xs[1:] + self.xs[:-1])
        self.yc[1:-1] = 0.5*(self.ys[1:] + self.ys[:-1])
        self.xc[0] = self.xc[1] - (self.xc[2] - self.xc[1])
        self.yc[0] = self.yc[0] - (self.yc[2] - self.yc[2])

        self.xc[-1] = self.xc[-2] + (self.xc[-2] - self.xc[-3])
        self.yc[-1] = self.yc[-2] + (self.yc[-2] - self.yc[-3])

        self.dx = np.sum(self.xs[1:] - self.xs[:-1])/len(self.xs[1:])
        self.dy = np.sum(self.ys[1:] - self.ys[:-1])/len(self.ys[1:])

        self.ubound = ubound
        self.lbound = lbound

        self.max_mode = min(max_mode, self.nx - 1)
        self.ubound_transform = 0.0*ubound
        self.lbound_transform = 0.0*lbound

        if len(xbasis) == 0.0 or len(ybasis) == 0.0:   #Calculate new basis functions here. Not necessary if it's already been done for this grid
            self.xbasis = np.zeros((self.nx-1, self.nx+1))
            #find eigenvalues and basis vectors (eigenvalues are m^2 and numbered k)
            self.m2, self.xbasis[:,1:-1] = self.find_eigenstuff()
            self.xbasis[:,0] = self.xbasis[:,1]; self.xbasis[:,-1] = self.xbasis[:,-2]

            self.ybasis = self.find_ys(self.m2)

        else:
            self.xbasis = xbasis; self.ybasis = ybasis; self.m2 = m2
        for k in range(self.nx-1):   #basis fourier transform on the two boundaries. This step needs doing twice
            self.ubound_transform[k] = self.coeff(self.ubound, k)
            self.lbound_transform[k] = self.coeff(self.lbound, k)

        self.phi = self.make_phi()

        #self.test_phi()
        self.bxp = (self.phi[1:,1:-1] - self.phi[:-1,1:-1])/self.dx
        self.byp = (self.phi[1:-1,1:] - self.phi[1:-1,:-1])/self.dy
        #plt.pcolormesh(self.xs, self.ys, self.phi[1:-1,1:-1].T)
        #plt.pcolormesh(self.xs, self.yc, self.byp.T)
        #plt.show()

    def find_ys(self, eigenvalues):
        #finds suitable (approximately hyperbolic) functions in the y direction
        #CAN SPEED THIS UP WITH ARRAYS!
        ybasis = np.zeros((self.nx-1, self.ny+1))  #number of x modes plus dimension (with ghosts) in the y direction
        for k in range(1,self.max_mode):  #run through the modes
            ybasis[k][0] = 1e-10
            ybasis[k][1] = 1e-10  #ensure the lower boundary derivative is zero
            for i in range(1, self.ny):
                ybasis[k][i+1] = eigenvalues[k]*ybasis[k][i]*self.dy**2
                ybasis[k][i+1] += 2*ybasis[k][i] - ybasis[k][i-1]
            dfact = (ybasis[k][-1] - ybasis[k][-2])/self.dy
            ybasis[k] = ybasis[k]/dfact
        return ybasis

    def make_phi(self):
        phi = np.zeros((self.nx+1, self.ny+1))
        for k in range(self.nx-1):
            phi = phi + self.ubound_transform[k]*self.xbasis[k,:][:,np.newaxis]*self.ybasis[k,:][np.newaxis,:]   #upper bound
            phi = phi - self.lbound_transform[k]*self.xbasis[k,:][:,np.newaxis]*np.flip(self.ybasis[k,:])[np.newaxis,:]   #lower bound
        return phi

    def find_eigenstuff(self):
        # Uses scipy tridiagonal solver to find the numerical approximations to the sine functions that have the desired properties.
        #Generates a matrix etc. then solves. Should only need to do this once for a given resolution. Doesn't depend on boundary conditions etc.
        d = 2*np.ones(self.nx-1)
        e = -1*np.ones(self.nx-2)
        d[0] = 1.0; d[-1] = 1.0
        w, v = eigh_tridiagonal(d, e)
        m2 = w/(self.dx**2)
        return m2, v.T

    def second_derivative(self, array, d):
        return (array[:-2] - 2*array[1:-1] + array[2:])/d**2

    def fcheck(self, bound_trans):
        bcheck = 0.0*ubound
        for k in range(self.nx -1):
            bcheck = bcheck + bound_trans[k]*self.xbasis[k,1:-1]
        return bcheck

    def mode(self, m):
        return np.sin(0.5*m*np.pi*self.xc/self.xs[-1])

    def coeff(self, bound, k):
        return np.sum(bound*self.xbasis[k,1:-1])

    def test_phi(self):
        d2x = (self.phi[:-2,1:-1] - 2*self.phi[1:-1,1:-1] + self.phi[2:,1:-1])/self.dx**2
        d2y = (self.phi[1:-1,:-2] - 2*self.phi[1:-1,1:-1] + self.phi[1:-1,2:])/self.dy**2
        print('Max laplacian', np.max(np.abs(d2x + d2y)))
        ubound_test = (self.phi[:,-1] - self.phi[:,-2])/self.dy
        print('Max ubound error', np.max(np.abs(self.ubound - ubound_test[1:-1])))
        lbound_test = (self.phi[:,1] - self.phi[:,0])/self.dy
        print('Max lbound error', np.max(np.abs(self.lbound - lbound_test[1:-1])))
        print('Max left side flux', np.max(np.abs(self.phi[1,1:-1] - self.phi[0,1:-1])))
        print('Max right side flux', np.max(np.abs(self.phi[-1,1:-1] - self.phi[-2,1:-1])))


class make_test_grid():
    def __init__(self, nx, ny):
        self.nx = nx + 1
        self.ny = ny + 1   #number of grid points in the main block
        self.data = [[],[]]
        self.data[0] = np.linspace(-1.,1.,self.nx)
        self.data[1] = np.linspace(0.,1.,self.ny)

