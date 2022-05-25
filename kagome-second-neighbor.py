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

def patch_chern(u, return_full=True):
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
    
    # the gauge-independent boundary phase mod 2pi
    boundary_phase_mod = np.angle(t_prod * b_prod * l_prod * r_prod)
    # the gauge-dependent boundary phase
    boundary_phase = 0
    for edge in (b_edge,r_edge,t_edge,l_edge):
        boundary_phase += np.sum(np.angle(edge), axis=0)
    
    if return_full:
        # return 'C, flux, boundary_phase_mod, boundary_phase, curvature, connection_x, connection_y'
        # NB: connection is meaningful only when a smooth gauge is used
        nkx,nky,nband = curvature.shape
        connection_x = np.angle(u12) * nkx
        connection_y = -np.angle(u41) * nky
        return np.array((flux - boundary_phase, flux, boundary_phase_mod, boundary_phase, curvature*nkx*nky, connection_x, connection_y))/(2*np.pi)
    else:
        # just return C
        return (flux-boundary_phase)/(2*np.pi)

def fix_gauge(u, elt=0):
    """fix the gauge of an eigenbasis u, such that the elt'th element of

    each eigenstate is real positive. I.e., the elt'th row of u is
    adjusted to be real. NB: numerically nothing is exactly zero so this
    is always doable
    """
    phase_factors = u[elt] / np.abs(u[elt])
    return u / phase_factors

def chern(u, leave_out_kx_bound = False, leave_out_ky_bound = False):
    """Compute Chern numbers of band wave functions u, where u[nx,ny] is
    the unitary matrix that diagonalizes H(kx,ky)

    u must have shape (nx,ny,nband,nband)

    If leave_out_kx_bound, then do not include cells connecting between kx = 0 with 2pi.
    Same for leave_out_ky_bound
    """
    nx,ny,nband = np.shape(u)[:3]

    res = np.zeros(nband)

    if leave_out_kx_bound == True:
        x0 = 1
    else:
        x0 = 0
    if leave_out_ky_bound == True:
        y0 = 1
    else:
        y0 = 0

    for x in range(x0,nx):
        for y in range(y0,ny):
            # 2 - 1 (x,y)
            # |   |
            # 3 - 4
            
            u1 = u[x,y]
            u2 = u[x-1,y]
            u3 = u[x-1,y-1]
            u4 = u[x,y-1]

            for b in np.arange(nband):
                u12 = np.dot(u1[:,b].conj(), u2[:,b])
                u23 = np.dot(u2[:,b].conj(), u3[:,b])
                u34 = np.dot(u3[:,b].conj(), u4[:,b])
                u41 = np.dot(u4[:,b].conj(), u1[:,b])
                res[b] += np.angle(u12 * u23 * u34 * u41)
    return res/(2*np.pi)

class Kagome: 
    def __init__(self,phi,t=0,m=0):
        """ t: second neighbor hopping strength """
        self.phi = phi
        self.t = t
        self.m = m
        # mb = mb or ma
        # mc = mc or 2*mb
        # self.ma = ma
        # self.mb = mb
        # self.mc = mc
    def hk(self,k1,k2):
        # k1 = k dot a1, k2 = k dot a2. Note that a3 = -(a1 + a2). See onenote
        res = np.zeros((3,3),dtype=np.complex)

        eiphi = np.exp(1j*self.phi/3)
        eimk = np.exp(-1j * np.array((k1,k2,-k1-k2)))
        res[0,1] = eiphi * (1 + eimk[0]) + self.t * eiphi.conj() * (eimk[1] + eimk[2])
        res[1,2] = eiphi * (1 + eimk[1]) + self.t * eiphi.conj() * (eimk[2] + eimk[0])
        res[2,0] = eiphi * (1 + eimk[2]) + self.t * eiphi.conj() * (eimk[0] + eimk[1])
        res += npext.dagger(res)
        res[0,0] = 0
        res[1,1] = self.m
        res[2,2] = self.m * 2

        return res
    def export_hk(self,nx=50,ny=50,fn='kagomi-erg'):
        with open(fn + '.par', 'w') as par:
            txt = """
            nx = %d
            ny = %d
            m = %g
            phi = %g
            t = %g
            """%(nx,ny,self.m,self.phi,self.t)
            print(txt,file=par)

        labels = 'k1 k2 kx ky erg1 erg2 ...'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])

        with open(fn + '.dat', 'w') as dat:
            print(header, file=dat)
            for k1 in np.linspace(0,2*np.pi,nx):
                for k2 in np.linspace(0,2*np.pi,ny):
                    kx = k2
                    ky = -(2*k1 + k2)/np.sqrt(3)
                    print('\r k1/pi = %g    k2/pi = %g           '%(k1/np.pi,k2/np.pi), end='')
                    hh = self.hk(k1,k2)
                    eig = la.eigvalsh(hh)
                    print('%g\t%g\t%g\t%g\t%s'%(k1,k2,kx,ky,'\t'.join([ '%g'%erg for erg in eig ])), file=dat)
                dat.write('\n')

def memo_u(nkx,nky,model,kxr=(0,2*np.pi), kyr=(0,2*np.pi),swap_xy=False,fix_gauge_elt=0,basis=None):
    """
    nkx,nky: k-space mesh size. Note that for kxr from 0 to 2pi, e.g., real space lattice size is actually nkx-1, because in k space, both 0 and 2pi are included.
    kxr: range of kx, inclusive on both limits
    kyr: range of ky, inclusive on both limits
    model: something that supports model.hk(kx,ky)
    fix_gauge_elt: if < 0, then do not fix, otherwise, pass it as elt to fix_gauge_elt()
    basis: if not None, then use basis[x,y] as the basis of the eigenstates
    """
    hk00 = model.hk(0,0)
    nband = hk00.shape[0]
    res = np.zeros((nkx,nky,nband,nband),dtype=np.complex)

    for x,kx in enumerate(np.linspace(*kxr, nkx, endpoint=True)):
        for y,ky in enumerate(np.linspace(*kyr, nky, endpoint=True)):
            print('\rmemo_u: x, y = %d, %d\t        '%(x,y), end='')
            if swap_xy == True:
                h = model.hk(ky,kx)
            else:
                h = model.hk(kx,ky)
            eig,u = la.eigh(h)
            if fix_gauge_elt >= 0:
                u = fix_gauge(u, fix_gauge_elt)
            if basis != None:
                u = np.dot(npext.dagger(basis[x,y]), u)
            res[x,y] = u
    return res

def chernkag(phi,t,m,nkx,nky,fix_gauge_elt=1):
    kag = Kagome(phi=phi,t=t,m=m)
    u = memo_u(nkx,nky,kag, fix_gauge_elt=fix_gauge_elt)
    return chern(u)
