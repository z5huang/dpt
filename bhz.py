import sys
import os
if'../lib/' not in sys.path:
    sys.path.append('../lib/')
import numpy as np
import math
#from numpy import linalg as la
from scipy import linalg as la
#from npext import *
import npext
import cmath
import itertools as it
import kpath

def chern(u, loop_kx = True, loop_ky = True):
    """Compute Chern numbers of band wave functions u, where u[nx,ny] is

    the unitary matrix that diagonalizes H(kx,ky)

    u must have shape (nx,ny,nband,nband)

    If loop_kx == False, then do not include cells connection the last
    kx (presumably 2pi) to the first kx (0). Same for loop_ky
    """

    # 2 - 1 (x,y)
    # |   |
    # 3 - 4

    # vectorized version. In ipython, use '%time <statement>' to time a
    # statement. E.g., a 9-band system with 200x200 k mesh takes about 110ms
    # 'user time'
    YY,XX,ROW,COL = range(4)
    u12 = np.sum(u.conj() * np.roll(u, 1, axis=XX), axis=ROW)
    u41 = np.sum(np.roll(u, -1, axis=YY).conj() * u, axis=ROW)
    u34 = np.roll(u12, -1, axis=YY).conj()
    u23 = np.roll(u41, 1, axis=XX).conj()
    berry_curvature = np.angle(u12 * u23 * u34 * u41)

    x0,y0 = 1-int(loop_kx), 1-int(loop_ky)
    res = np.sum(berry_curvature[x0:,y0:], axis=(XX,YY))
    
    # non-vectorized version. Significantly slower but makes it cleaner. A
    # 9-bandsystem with 200x200 k mesh takes about 4.6s 'user time'
    ## nx,ny,nband = np.shape(u)[:3]
    ## x0,y0 = 1-int(loop_kx), 1-int(loop_ky)
    ## res = np.zeros(nband)
    ## for x in range(x0,nx):
    ##     for y in range(y0,ny):
    ##         # 2 - 1 (x,y)
    ##         # |   |
    ##         # 3 - 4
    ##         
    ##         u1 = u[x,y]
    ##         u2 = u[x-1,y]
    ##         u3 = u[x-1,y-1]
    ##         u4 = u[x,y-1]
    ## 
    ##         # u12 = [ <a1|a2>, <b1|b2>, <c1|c2>, ... ], where a,b,c... are band
    ##         # indices
    ##         u12 = np.sum(u1.conj()*u2, axis=0)
    ##         u23 = np.sum(u2.conj()*u3, axis=0)
    ##         u34 = np.sum(u3.conj()*u4, axis=0)
    ##         u41 = np.sum(u4.conj()*u1, axis=0)
    ##         res += np.angle(u12 * u23 * u34 * u41)

    return res/(2*np.pi)


def patch_chern_test(u):
    
    nx,ny,nrow,ncol = u.shape
    flux = np.zeros(u.shape[-1],dtype=np.complex)
    for x in range(nx-1):
        for y in range(ny-1):
            # 4 - 3   /|\
            # |   |    |
            # 1 - 2  (x,y)-->
            u1 = u[x,y]
            u2 = u[x+1,y]
            u3 = u[x+1,y+1]
            u4 = u[x,y+1]

            u12 = np.sum(u1.conj() * u2, axis=0)
            u23 = np.sum(u2.conj() * u3, axis=0)
            u34 = np.sum(u3.conj() * u4, axis=0)
            u41 = np.sum(u4.conj() * u1, axis=0)

            flux += np.angle(u12 * u23 * u34 * u41)

    # wf on 4 zore edges
    b_edge = u[:,0]
    r_edge = u[-1]
    t_edge = u[:,-1]
    l_edge = u[0]

    # berry factor over the edges
    berry = 1
    for x in range(nx-1):
        # bottom
        berry *= np.sum(b_edge[x].conj() * b_edge[x+1], axis=0)
        # top
        berry *= np.sum(t_edge[x] * (t_edge[x+1].conj()), axis=0)
    for y in range(ny-1):
        # right
        berry *= np.sum(r_edge[y].conj() * r_edge[y+1], axis=0)
        # left
        berry *= np.sum(l_edge[y] * (l_edge[y+1].conj()), axis=0)
    berry_phase = np.angle(berry)
    return np.array((flux - berry_phase, flux, berry_phase))/(np.pi * 2)

def patch_chern(u):
    # 4 - 3   /|\
    # |   |    |
    # 1 - 2  (x,y)-->

    u1 = u[:-1,:-1]
    u2 = u[1:, :-1]
    u3 = u[1:, 1: ]
    u4 = u[:-1, 1:]

    u12 = np.sum(u1.conj() * u2, axis=2)
    u23 = np.sum(u2.conj() * u3, axis=2)
    u34 = np.sum(u3.conj() * u4, axis=2)
    u41 = np.sum(u4.conj() * u1, axis=2)

    curvature = np.angle(u12 * u23 * u34 * u41)
    flux = np.sum(curvature, axis=(0,1))

    # now to the boundary
    b_edge = u12[:,0]
    r_edge = u23[-1]
    t_edge = u34[:,-1]
    l_edge = u41[0]

    t_prod = np.prod(t_edge, axis=0)
    b_prod = np.prod(b_edge, axis=0)
    l_prod = np.prod(l_edge, axis=0)
    r_prod = np.prod(r_edge, axis=0)
    
    boundary_phase = np.angle(t_prod * b_prod * l_prod * r_prod)
    return np.array((flux - boundary_phase, flux, boundary_phase))/(2*np.pi)

def patch_chern_wrong(u):
    """ Compute the Chern numbers of band wavefunctions u (shape: nx,ny,nrow, ncol) over a patch of BZ, i.e., not periodic over either boundaries. """

    # 2 - 1 (x,y)
    # |   |
    # 3 - 4
    YY,XX,ROW,COL = range(4)
    u12 = np.sum(u.conj() * np.roll(u, 1, axis=XX), axis=ROW)
    u41 = np.sum(np.roll(u, -1, axis=YY).conj() * u, axis=ROW)
    u34 = np.roll(u12, -1, axis=YY).conj()
    u23 = np.roll(u41, 1, axis=XX).conj()
    berry_curvature = np.angle(u12 * u23 * u34 * u41)[1:,1:]

    berry_flux = np.sum(berry_curvature, axis=(XX,YY))

    # berry phase factors along the four zone edges
    t_edge = np.prod(u12[0,1:], axis=0)
    b_edge = np.prod(u12[-1,1:], axis=0).conj()
    l_edge = np.prod(u41[:-1,0], axis=0).conj()
    r_edge = np.prod(u41[:-1,-1], axis=0)
    boundary_berry_phase = np.angle(t_edge * b_edge * l_edge * r_edge)
    return np.array([berry_flux - boundary_berry_phase, berry_flux, boundary_berry_phase])/(2*np.pi)

    
    
    
    


def hk(kx,ky,p,q,t=0,swap_xy=False):
    """k-space Hofstadter Hamiltonian on square lattice, flux per plaquette is 2pi*p/q

    t: next nearest hopping, i.e. diagonals of the square plaquette
    swap_xy: if true, swap kx and ky before proceeding. Useful for calculating winding number in the other direction, say
    """
    
    if swap_xy:
        kx,ky = ky,kx

    res = np.zeros((q,q),dtype=np.complex)
    i,j = np.indices((q,q))
    phi = 2*np.pi*p/q
    res[i==j] = 2*np.cos(np.arange(1,q+1) * phi + kx)

    sup_diag = 1 + 2*t*np.cos(np.arange(1.5,q+0.5) * phi + kx)
    res[i==j-1] = sup_diag
    res[i==j+1] = sup_diag

    corner = 1 + 2*t*np.cos(0.5*phi + kx)
    res[0,-1]+= np.exp(1j * ky) * corner
    res[-1,0]+= np.exp(-1j * ky) * corner
    return res

def memo_u(nx,ny,p,q,t,swap_xy=False,h=hk,kxfrac=1,kyfrac=1,kx0=0,ky0=0,basis=None):
    """ memoize the wavefunctions on the momentum mesh

        if basis != None, then store wavefunctions in the basis of basis[x,y]
    """
    res = np.zeros((nx+1,ny+1,q,q),dtype=np.complex)
    dkx = np.pi*2/(nx*kxfrac)
    dky = np.pi*2/(ny*kyfrac)
    for x in range(nx+1):
        kx = dkx * x + kx0
        for y in range(ny+1):
            ky = dky * y + ky0
            hh = h(kx,ky,p,q,t,swap_xy)
            eig,u = la.eigh(hh)
            if basis == None:
                res[x,y] = u
            else:
                res[x,y] = np.dot(npext.dagger(basis[x,y]),u)
    return res

def chern_hof(p,q,t=0,swap_xy=False,nx=50,ny=50,h=hk,kxfrac=1,kyfrac=1):
    """ compute the Chern numbers of the Hofstadter model, cf. def hk() """
    uu = memo_u(nx,ny,p,q,t,swap_xy,h,kxfrac,kyfrac)
    return chern(uu)
