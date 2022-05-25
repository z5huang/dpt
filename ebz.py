###
# extended BZ chern number calculation for relative wave functions
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

one2 = np.identity(2, dtype=np.complex)
one4 = np.identity(4, dtype=np.complex)
sx = np.array([[0,1],[1,0]], dtype=np.complex)
sy = np.array([[0,-1j],[1j,0]], dtype=np.complex)
sz = np.array([[1,0],[0,-1]], dtype=np.complex)

def chern(u):
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
    return flux/(2*np.pi)

class Haldane:
    def __init__(self,m=0.5,phi=0.3*np.pi,t1=0.3,t2=0.3,t3=0.3):
        self.m,self.phi = m,phi
        self.t1,self.t2,self.t3 = t1,t2,t3
    def hk(self,kx,ky):
        """ haldane model """
        bx = -1 - np.cos(kx) - np.cos(ky)
        by = -np.sin(kx) - np.sin(ky)
        bz = self.m + 2*np.sin(self.phi)*(
            self.t1 * np.sin(kx)
            - self.t2 * np.sin(ky)
            + self.t3 * np.sin(ky - kx))
        omega = -2*np.cos(self.phi)*(
            self.t1 * np.cos(kx)
            + self.t2 * np.cos(ky)
            + self.t3 * np.cos(ky - kx))
        return omega*one2 + bx*sx + by*sy + bz*sz
class Hof:
    def __init__(self,p=1,q=3,t=0):
        self.p,self.q,self.t = p,q,t
    def hk(self,kx,ky):
        res = np.zeros((self.q,self.q),dtype=np.complex)
        i,j = np.indices((self.q,self.q))
        phi = 2*np.pi*self.p/self.q
        res[i==j] = 2*np.cos(np.arange(1,self.q+1) * phi + kx)

        sup_diag = 1 + 2*self.t*np.cos(np.arange(1.5,self.q+0.5) * phi + kx)
        res[i==j-1] = sup_diag
        res[i==j+1] = sup_diag

        corner = 1 + 2*self.t*np.cos(0.5*phi + kx)
        res[0,-1]+= np.exp(1j * ky) * corner
        res[-1,0]+= np.exp(-1j * ky) * corner
        return res

def fix_gauge(u, elt=0):
    """fix the gauge of an eigenbasis u, such that the elt'th element of

    each eigenstate is real positive. I.e., the elt'th row of u is
    adjusted to be real. NB: numerically nothing is exactly zero so this
    is always doable
    """
    phase_factors = u[elt] / np.abs(u[elt])
    return u / phase_factors

def memo_u(nkx,nky,model,kxr=(0,2*np.pi), kyr=(0,2*np.pi),swap_xy=False,fix_gauge_elt=0,basis=None):
    """
    nkx,nky: k-space mesh size. Note that for kxr from 0 to 2pi, e.g., real space lattice size is actually nkx-1, because in k space, both 0 and 2pi are included.
    kxr: range of kx, inclusive on both limits
    kyr: range of ky, inclusive on both limits
    model: something that supports model.hk(kx,ky)
    fix_gauge_elt: if < 0, then do not fix, otherwise, pass it as elt to fix_gauge_elt()
    NB: gauge is fixed before basis transformation
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

def parallelize(u,uref):
    """ adjust the phases in the column vectors u, so that they are parallel with the corresponding columns in uref """
    ncol = u.shape[1]

    # <u(1)|ref(1)>, <u(2)|ref(2)>, ..., <u(ncol) | ref(ncol)>
    overlaps = np.sum(np.conj(u) * uref, 0)

    phase_factors = overlaps / np.abs(overlaps)
    return u*phase_factors

def parallelize_path(all_u):
    """ parallelization of a path of u's """
    ntot = all_u.shape[0]
    res = np.zeros(all_u.shape,dtype=np.complex)
    res[0] = all_u[0]
    for n in range(ntot-1):
        res[n+1] = parallelize(all_u[n+1], res[n])
    return res

def smooth_path(all_u, ref_berry_phase = None):
    """smoothify a path of u's. After 'smoothification', the first and

    last elements of all_u will be parallel

    if ref_berry_phase == None, then just use the phases of
    berry_factors as the corresponding berry phases. Otherwise, the
    Berry phases are such that the phase change from the ref_berry_phase
    be small (i.e. in the branch (-pi,pi])

    """
    u_para = parallelize_path(all_u)
    u1,uN = u_para[[0,-1]]
    ntot,nband = all_u.shape[0:2]

    nseg = ntot - 1 # number of segments in the path

    berry_factors = (u1.conj() * uN).sum(axis=-2)
    if ref_berry_phase == None:
        berry_phase = np.angle( berry_factors )
    else:
        phase_change = np.angle(berry_factors * np.exp(-1j*ref_berry_phase))
        berry_phase = ref_berry_phase + phase_change
    twist = np.exp(-1j * np.outer(np.arange(0,ntot), berry_phase)/nseg).reshape(ntot,1,nband)
    return u_para*twist,berry_phase

def memo_u_smooth(nkx,nky,model,swap_xy=False):
    """memoize u with smooth gauge, that is, for all kx, the gauge is

    smooth and periodic along ky, and for ky=0, the gauge is smooth and
    periodic along kx.

    nkx,nky: k space mesh size, including both 0 and 2pi

    """
    hk00 = model.hk(0,0)
    nband = hk00.shape[0]
    u = memo_u(nkx,nky,model,kxr=(0, np.pi*2),kyr=(0, np.pi*2),swap_xy=swap_xy, fix_gauge_elt=-1)
    
    u[:,0] = smooth_path(u[:,0])[0]

    ref_berry_phase = np.zeros(nband)
    for x,u_kx in enumerate(u):
        u[x],ref_berry_phase = smooth_path(u_kx, ref_berry_phase)
    return u

def advance_basis(u, cvec=None, nbz=1):
    """advance a smooth basis, u, by a number of BZs, nbz. cvec: vector

    of Chern numbers

    if cvec == None, then compute it from the u basis
    """
    nkx,nky,nband = u.shape[0:3]

    if cvec == None:
        cvec = chern(u)

    # the vector of nbz * ky for all ky
    kyvec = np.linspace(0, nbz*2*np.pi, nky, endpoint=True)
    phases = np.outer(kyvec, cvec).reshape(1,nky,1,nband)
    return np.exp(1j * phases) * u

def advance_relU(relU, cvec, nbz=1):
    """advance wf in a smooth basis, relU, by a number of BZs, nbz. cvec: vector

    of Chern numbers

    """
    nkx,nky,nband = relU.shape[0:3]
    # the vector of nbz * ky for all ky
    kyvec = np.linspace(0, nbz*2*np.pi, nky, endpoint=True)
    phases = np.outer(kyvec, cvec).reshape(1,nky,nband,1)
    return np.exp(-1j * phases) * relU
    
def quench_chern_ebz(mi,mf,nkx,nky,nbz=None,swap_xy=False,fix_gauge_elt=0):

    uf = memo_u_smooth(nkx,nky,mf,swap_xy)
    cvec = chern(uf)
    ui = memo_u(nkx,nky,mi,kxr=(0,np.pi*2),kyr=(0,np.pi*2),swap_xy=swap_xy,fix_gauge_elt=fix_gauge_elt,basis=uf)

    if nbz == None:
        nbz = nky - 1 # number of BZs in the extended BZ

    nband = uf.shape[-1]
    flux = np.zeros(nband)

    for bz in range(nbz):
        ru = advance_relU(ui, cvec, bz)
        flux += chern(ru)
    return (flux, ui, ru)

################################################################
## test codes
def phase_adv(upath):
    """ compute phase advancement for each step along a path """
    ntot,nband = upath.shape[0:2]
    for n in range(ntot-1):
        arg = np.angle((upath[n].conj() * upath[n+1]).sum(axis=-2))
        print('angle between ', n, 'and ', n+1, ': ', arg)

def phase_paths(upath1, upath2):
    """ compute the relative phase of two paths at each step """
    ntot,nband = upath1.shape[0:2]
    for n in range(ntot):
        arg = np.angle((upath1[n].conj() * upath2[n]).sum(axis=-2))
        print('phase at step ', n, ': ', arg)
        
