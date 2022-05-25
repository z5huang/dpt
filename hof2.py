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
#from scipy.special import binom as binom
import kpath
#from chern import chern

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
    c,f,b = patch_chern(uu)
    print('Chern = %s\nTotal flux = %s\nContour phase = %s\n'%(c,f,b))

def parallelize(u,uref):
    """ adjust the phases in the column vectors u, so that they are parallel with the corresponding columns in uref """
    ncol = u.shape[1]

    # <u(1)|ref(1)>, <u(2)|ref(2)>, ..., <u(ncol) | ref(ncol)>
    overlaps = np.sum(np.conj(u) * uref, 0)

    phase_factors = overlaps / np.abs(overlaps)
    return u*phase_factors

def memo_u_quench(nx,ny,q,pi,pf,ti,tf,swap_xy=False,h=hk,kxfrac=1,kyfrac=1,kx0=0, ky0=0):
    """ memoize the relative wavefunctions before and after the quench. I.e., the initial wavefunctions in the basis of the quenched states """

    ###
    # Parallel gauge for both wavefunctions. i.e., for any kx, the ky line is
    # parallel gauged. For ky = 0, the kx line is parallel gauged.
    res = np.zeros((nx+1,ny+1,q,q),dtype=np.complex)
    dkx = np.pi*2/(nx*kxfrac)
    dky = np.pi*2/(ny*kyfrac)
    # wavefunctions at the same kx and previous ky
    ui_prev_ky = np.diag(np.ones(q,dtype=np.complex))
    uf_prev_ky = np.diag(np.ones(q,dtype=np.complex))
    # wavefunctions at ky=0 and the previous kx
    ui_prev_kx = np.diag(np.ones(q,dtype=np.complex))
    uf_prev_kx = np.diag(np.ones(q,dtype=np.complex))
    for x in range(nx+1):
        kx = dkx * x + kx0
        for y in range(ny+1):
            ky = dky * y + ky0
            hi = h(kx,ky,pi,q,ti,swap_xy)
            eig,ui = la.eigh(hi)
            hf = h(kx,ky,pf,q,tf,swap_xy)
            eig,uf = la.eigh(hf)
            if y == 0:
                ui_para = parallelize(ui,ui_prev_kx)
                uf_para = parallelize(uf,uf_prev_kx)
                ui_prev_kx = ui_para
                uf_prev_kx = uf_para
            else:
                ui_para = parallelize(ui,ui_prev_ky)
                uf_para = parallelize(uf,uf_prev_ky)
            ui_prev_ky = ui_para
            uf_prev_ky = uf_para
            res[x,y] = np.dot(npext.dagger(uf_para), ui_para)
    return res

def patch_chern_hof_parallel(q,pi,ti,tf,pf=None,nx=50,ny=None,swap_xy=False,h=hk,kx0=0,ky0=0,kxfrac=1,kyfrac=1):
    """ compute the Chern numbers for the relative wavefunctions, using a parallel-gauged final basis"""
    if pf == None:
        pf = pi
    if ny == None:
        ny = nx
    print('pi = %g, pf = %g, ti = %g, tf = %g'%(pi,pf,ti,tf))
    uu = memo_u_quench(nx,ny,q,pi,pf,ti,tf,swap_xy,h,kxfrac,kyfrac,kx0,ky0)
    return patch_chern(uu)
