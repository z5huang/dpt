###
# Some code to look at relation between wavefunction zeros (i.e. nodes) and Berry flux (a gauge-independent property)

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
def export_u_and_curvature(nkx,nky,model,
                           kxr=(0,2*np.pi),kyr=(0,2*np.pi),
                           swap_xy=False,fix_gauge_elt=0,
                           fn_wf = 'wf', fn_curv = 'curv', fn_conn = 'conn'):
    """
    export wavefunction magnitude, corresponding curvature, and connection vector, of all bands of 'model'

    nkx and nky are k-space mesh size, not related to real space lattice size. If kx includes both 0 and 2pi, then lattice size is instead nkx-1
    """
    u = memo_u(nkx,nky,model,kxr,kyr,swap_xy,fix_gauge_elt)
    cc,ff,bbm,bb,curv,conn_x, conn_y = patch_chern(u,return_full=True)

    labels_wf = 'kx ky <1|psi1>_abs <1|psi1>_phase <2|psi1>_abs <2|psi1>_phase ... <1|psi2>_abs <1|psi2>_phase <2|psi2>_abs ...'
    header_wf = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels_wf.split())])    

    labels_curv = 'kx ky curvature1/2pi curvature2/2pi ...'
    header_curv = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels_curv.split())])

    labels_conn = 'kx ky conn1_x conn1_y conn2_x conn2_y ...'
    header_conn = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels_conn.split())])

    with open(fn_wf + '.dat', 'w') as dat:
        print(header_wf, file=dat)
        for x,kx in enumerate(np.linspace(*kxr, nkx, endpoint=True)):
            for y,ky in enumerate(np.linspace(*kyr, nky, endpoint=True)):
                print('\r nkx = %d / %d    nky = %d / %d    '%(x,nkx,y,nky), end='')
                print('%g\t%g\t%s'%(kx,ky,  '\t'.join([ '%g\t%g'%(np.abs(wf),np.angle(wf)) for wf in u[x,y].transpose().flatten() ])), file=dat)
            dat.write('\n')
                
    with open(fn_curv + '.dat', 'w') as dat:
        print(header_curv, file=dat)
        for x,kx in enumerate(np.linspace(*kxr, nkx-1, endpoint=False)):
            for y,ky in enumerate(np.linspace(*kyr, nky-1, endpoint=False)):
                print('\r nkx = %d / %d    nky = %d / %d    '%(x,nkx-1,y,nky-1), end='')
                print('%g\t%g\t%s'%(kx,ky,'\t'.join([ '%g'%curv_k for curv_k in curv[x,y].flatten() ])), file=dat)
            dat.write('\n')

    with open(fn_conn + '.dat', 'w') as dat:
        print(header_conn, file=dat)
        for x,kx in enumerate(np.linspace(*kxr, nkx-1, endpoint=False)):
            for y,ky in enumerate(np.linspace(*kyr, nky-1, endpoint=False)):
                print('\r nkx = %d / %d    nky = %d / %d    '%(x,nkx-1,y,nky-1), end='')
                print('%g\t%g\t'%(kx,ky), file=dat, end='')
                print('\t'.join( [ '%g\t%g\t'%(cx,cy) for cx,cy in zip(conn_x[x,y],conn_y[x,y]) ]), file=dat)
            dat.write('\n')

    print('\nChern:\t %s\nTotal Flux:\t %s\nBoundary flux (mod 2pi):\t %s\nBoundary flux (in current gauge):\t %s'%(cc,ff,bbm,bb))
                
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

def memo_u_parallel(nkx,nky,model,kxr=(0, 2*np.pi), kyr=(0, 2*np.pi), swap_xy=False):
    """ memoize u with parallel gauge: at fixed kx, the gauge is parallel along ky, with bounds given by kyr. For ky = kyr[0], i.e. the lower bound, the gauge is parallel along kx

    nkx and nky are k-space mesh size, not related to real space lattice size. If kx includes both 0 and 2pi, then lattice size is instead nkx-1
    """
    hk00 = model.hk(0,0)
    nband = hk00.shape[0]
    #res = np.zeros((nkx,nky,nband,nband),dtype=np.complex)

    u = memo_u(nkx,nky,model,kxr,kyr,swap_xy,fix_gauge_elt=-1)
    u[:,0] = parallelize_path(u[:,0])
    res = np.asarray([ parallelize_path(ux) for ux in u ])
    return res

def patch_chern_quench(mi,mf,nkx,nky,kxr=(0,2*np.pi),kyr=(0,2*np.pi),swap_xy=False,fix_gauge_elt=0):
    """ calculate the 'quench chern number'

    mi: initial model
    mf: final model

    nkx and nky are k-space mesh size, not related to real space lattice size. If kx includes both 0 and 2pi, then lattice size is instead nkx-1    
    """
    uf = memo_u_parallel(nkx,nky,mf,kxr,kyr,swap_xy)
    ui = memo_u(nkx,nky,mi,kxr,kyr,swap_xy,fix_gauge_elt,basis=uf)
    return (patch_chern(ui), uf,ui)

def export_rel_u_and_curvature(nkx,nky,mi,mf,
                           kxr=(0,2*np.pi),kyr=(0,2*np.pi),
                           swap_xy=False,fix_gauge_elt=0,
                           fn_wf = 'wf', fn_curv = 'curv', fn_conn = 'conn'):
    """
    export relative wavefunction magnitude, corresponding curvature, and connection vector, of all bands of 'model'

    nkx and nky are k-space mesh size, not related to real space lattice size. If kx includes both 0 and 2pi, then lattice size is instead nkx-1
    """

    pc,uf,u = patch_chern_quench(mi,mf,nkx,nky,kxr,kyr,swap_xy,fix_gauge_elt)
    cc,ff,bbm,bb,curv,conn_x, conn_y = pc

    labels_wf = 'kx ky <1|psi1>_abs <1|psi1>_phase <2|psi1>_abs <2|psi1>_phase ... <1|psi2>_abs <1|psi2>_phase <2|psi2>_abs ...'
    header_wf = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels_wf.split())])    

    labels_curv = 'kx ky curvature1/2pi curvature2/2pi ...'
    header_curv = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels_curv.split())])

    labels_conn = 'kx ky conn1_x conn1_y conn2_x conn2_y ...'
    header_conn = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels_conn.split())])

    with open(fn_wf + '.dat', 'w') as dat:
        print(header_wf, file=dat)
        for x,kx in enumerate(np.linspace(*kxr, nkx, endpoint=True)):
            for y,ky in enumerate(np.linspace(*kyr, nky, endpoint=True)):
                print('\r nkx = %d / %d    nky = %d / %d    '%(x,nkx,y,nky), end='')
                print('%g\t%g\t%s'%(kx,ky,  '\t'.join([ '%g\t%g'%(np.abs(wf),np.angle(wf)) for wf in u[x,y].transpose().flatten() ])), file=dat)
            dat.write('\n')
                
    with open(fn_curv + '.dat', 'w') as dat:
        print(header_curv, file=dat)
        for x,kx in enumerate(np.linspace(*kxr, nkx-1, endpoint=False)):
            for y,ky in enumerate(np.linspace(*kyr, nky-1, endpoint=False)):
                print('\r nkx = %d / %d    nky = %d / %d    '%(x,nkx-1,y,nky-1), end='')
                print('%g\t%g\t%s'%(kx,ky,'\t'.join([ '%g'%curv_k for curv_k in curv[x,y].flatten() ])), file=dat)
            dat.write('\n')

    with open(fn_conn + '.dat', 'w') as dat:
        print(header_conn, file=dat)
        for x,kx in enumerate(np.linspace(*kxr, nkx-1, endpoint=False)):
            for y,ky in enumerate(np.linspace(*kyr, nky-1, endpoint=False)):
                print('\r nkx = %d / %d    nky = %d / %d    '%(x,nkx-1,y,nky-1), end='')
                print('%g\t%g\t'%(kx,ky), file=dat, end='')
                print('\t'.join( [ '%g\t%g\t'%(cx,cy) for cx,cy in zip(conn_x[x,y],conn_y[x,y]) ]), file=dat)
            dat.write('\n')

    print('\nChern:\t %s\nTotal Flux:\t %s\nBoundary flux (mod 2pi):\t %s\nBoundary flux (in current gauge):\t %s'%(cc,ff,bbm,bb))

def relU_next(relU, cvec):
    """ given a set of relative wavefunctions, relU, compute its value in the next BZ, i.e., kx -> kx + 2pi, where cvec is a vector of Chern numbers

    relU.shape = (nkx,nky, nband, nband), where nkx,nky are k-space mesh size -- implying that real space mesh size is nkx-1,nky-1, assuming both kx and ky are from 0 to 2pi, *inclusive*
    cvec.shape = (nband)
    """

    nkx,nky,nband = relU.shape[:3]

    kyvec = np.linspace(0,np.pi*2,nky,endpoint=True)

    advance = np.exp(-1j * np.outer(kyvec, cvec)).reshape((1,nky,nband,1))

    return advance * relU
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
        
        

def quench_chern_extended_BZ(mi,mf,nkx,nky,nbz = None, swap_xy=False):
    """
    calculate the Chern number in the extended BZ, [0, 2pi*nbz] x [0, 2pi]

    mi: initial model
    mf: final model

    nkx,nky: k-space mesh size
    nbz: number of consecutive BZ's in the kx direction
    """
    uf = memo_u_smooth(nkx,nky,mf,swap_xy)
    cvec = patch_chern(uf, return_full=True)[1]
    print('cvec = ',cvec)
    
    ui = memo_u(nkx,nky,mi,kxr=(0,np.pi*2),kyr=(0,np.pi*2),swap_xy=swap_xy,fix_gauge_elt=1,basis=uf)
    uii = ui # remember initial ui

    if nbz == None:
        nbz = nky - 1 # number of BZs in the extended BZ

    nband = uf.shape[-1]
    flux,bp_mod,bp = np.zeros((3,nband))
    for bz in range(nbz):
        pc = patch_chern(ui,return_full=True)
        flux += pc[1]
        bp_mod += pc[2]
        bp += pc[3]
        ui = relU_next(ui,cvec)
    return (flux,bp,bp_mod,ui,uii,cvec)
