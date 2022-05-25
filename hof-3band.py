# 3-band Hofstadter model with onsite mass to open up trivial gaps
import sys
import os
if'../lib/' not in sys.path:
    sys.path.append('../lib/')
import numpy as np
import math
#from numpy import linalg as la
from scipy import linalg as la
#from scipy.optimize import minimize as spmin
import scipy.optimize as so
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

def brent_premin(f,xrange,nx_premin):
    """minimize a function using brent(), with the bracket calculated on

    a finer mesh given by nx_premin
    """
    xvals = np.linspace(*xrange, nx_premin)
    dx = xvals[1] - xvals[0]
    fvals = [ f(x) for x in xvals ]
    imin = np.argmin(fvals)
    #print('imin = %d'%imin)
    xmin = xvals[imin]
    xleft = xmin - dx
    xright = xmin + dx
    while f(xleft) == f(xmin):
        xleft -= dx
        print('bracket failed')
    while f(xright) == f(xmin):
        xright += dx
    bracket = (xleft, xmin, xright)
    #bracket = (xmin - dx, xmin, xmin + dx)
    #print('fvals = %g, %g, %g'%(f(xmin - dx), f(xmin), f(xmin + dx)))
    return so.brent(f,brack=bracket, full_output=True)

def patch_chern(u, return_full=True):
    """CAUTION: in populating u, one should slightly shift the k-mesh,

    say by 1e-6 in each direction, to avoid numerically sick cases in
    calculating the boundary_phase below. Otherwise, even with a smooth
    gauge, the boundary_phase may not vanish.
    """
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

    # now to the boundary. brtl: bottom, right, top, left
    b_edge = u12[:,0]    # i.e., for all kx, but for ky = 0
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
class Hof:
    def __init__(self,p=1,q=3,t=0,m=0.1):
        self.p,self.q,self.t = p,q,t
        self.m = m
    def hk(self,kx,ky):
        res = np.zeros((self.q,self.q),dtype=np.complex)
        i,j = np.indices((self.q,self.q))
        phi = 2*np.pi*self.p/self.q
        res[i==j] = 2*np.cos(np.arange(1,self.q+1) * phi + kx) + np.arange(self.q)*self.m

        sup_diag = 1 + 2*self.t*np.cos(np.arange(1.5,self.q+0.5) * phi + kx)
        res[i==j-1] = sup_diag
        res[i==j+1] = sup_diag

        corner = 1 + 2*self.t*np.cos(0.5*phi + kx)
        res[0,-1]+= np.exp(1j * ky) * corner
        res[-1,0]+= np.exp(-1j * ky) * corner
        return res

    def export_hk(self,nx=50,ny=50,fn='hof-erg'):
        with open(fn + '.par', 'w') as par:
            txt = """
            nx = %d
            ny = %d
            p = %d
            q = %d
            m = %g
            """%(nx,ny,self.p,self.q,self.m)
            print(txt,file=par)

        labels = 'kx ky erg1 erg2 ...'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])

        with open(fn + '.dat', 'w') as dat:
            print(header, file=dat)
            for k1 in np.linspace(0,2*np.pi,nx):
                for k2 in np.linspace(0,2*np.pi,ny):
                    print('\r k1/pi = %g    k2/pi = %g           '%(k1/np.pi,k2/np.pi), end='')
                    hh = self.hk(k1,k2)
                    eig = la.eigvalsh(hh)
                    print('%g\t%g\t%s'%(k1,k2,'\t'.join([ '%g'%erg for erg in eig ])), file=dat)
                dat.write('\n')
    def cubic_discriminant(self,kx,ky):
        """ return discriminant of the cubic equation, whose zeros correspond to the energy eigenvalues """

        phi = 2*np.pi*self.p/self.q
        dd = 2*np.cos(np.arange(1,self.q+1) * phi + kx) + np.arange(self.q)*self.m
        vv = 1 + 2*self.t*np.cos(np.arange(1.5,self.q+1.5) * phi + kx)

        a = -1
        b = np.sum(dd)
        c = - np.sum( dd * np.roll(dd,1) ) + np.sum(vv*vv)
        d = np.prod(dd) - np.sum(dd * np.roll(vv*vv, -1)) + 2*np.prod(vv) * np.cos(ky)
        res = 18 * a*b*c*d - 4*b**3*d + b**2*c**2 - 4*a*c**3 - 27*a**2*d**2
        return res

    def min_cubic_discriminant(self,nx,ny):
        kgrid = npext.n2kgrid([nx,ny])
        all_disc = [ self.cubic_discriminant(kx,ky) for kx,ky in kgrid ]
        return np.min(all_disc)

    def export_cubic_discriminant(self,nx=100,ny=2,tmin=-4,tmax=4,dt=0.01,mmin=-1,mmax=1,dm=0.1,fn='cubic_discriminant'):
        # it is good to keep ny=2, because it seems like gap collapsing always occurs at ky = 0 or pi
        with open(fn + '.par', 'w') as par:
            txt = """
            p = %d
            q = %d
            nx = %d
            ny = %d
            tmin = %g
            tmax = %g
            dt = %g
            mmin = %g
            mmax = %g
            dm = %g
            """%(self.p,self.q,nx,ny,tmin,tmax,dt,mmin,mmax,dm)
            print(txt,file=par)

        labels = 't m Del'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        
        with open(fn + '.dat', 'w') as dat:
            print(header, file=dat)
            for m in np.arange(mmin,mmax,dm):
                for t in np.arange(tmin,tmax,dt):
                    print('\rt = %g  m = %g        '%(t,m), end='')
                    self.t = t
                    self.m = m
                    disc = self.min_cubic_discriminant(nx,ny)
                    print('%g\t%g\t%g'%(t,m,disc), file=dat)
                dat.write('\n')
                dat.flush()
        
    def gap_n(self,kx,ky,n):
        eig = la.eigvalsh(self.hk(kx,ky))
        return eig[n+1] - eig[n]
    def gapn_gen(self,ky,n):
        def res(kx):
            eig = la.eigvalsh(self.hk(kx,ky))
            return eig[n+1] - eig[n]
        return res
    def min_gap_n_ky(self,ngap,ky,nx_premin=50):
        """ use brent() to minimize the nth gap, with bracket computed using pre_min() """
        gapn = self.gapn_gen(ky,ngap)
        return brent_premin(gapn,(0,2*np.pi),nx_premin)[:-1]
    def min_gap_n(self,ngap,ny,nx_premin=50):
        kgrid = np.linspace(0,np.pi*2,ny,endpoint=False)
        res = [ self.min_gap_n_ky(ngap,ky,nx_premin) for ky in kgrid ]
        kxgap = np.array(res)[:,1]

        idx = np.argmin(kxgap)
        return res[idx]
        
    def export_min_gap(self,nx_premin=50,ny=2,tmin=0,tmax=4,dt=0.01,mmin=0,mmax=6,dm=0.1,fn='min_gap'):
        """ for each ky, look for minimal gap using the min_gap_ky_recursive """
        # it is good to keep ny=2, because it seems like gap collapsing always occurs at ky = 0 or pi
        with open(fn + '.par', 'w') as par:
            txt = """
            p = %d
            q = %d
            nx_premin = %d
            ny = %d
            tmin = %g
            tmax = %g
            dt = %g
            mmin = %g
            mmax = %g
            dm = %g
            """%(self.p,self.q,nx_premin,ny,tmin,tmax,dt,mmin,mmax,dm)
            print(txt,file=par)

        labels = 't m gap1 kx_itr1 gap2 kx_itr2 ...'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        
        with open(fn + '.dat', 'w') as dat:
            print(header, file=dat)
            for m in np.arange(mmin,mmax,dm):
                for t in np.arange(tmin,tmax,dt):
                    print('\rt = %g  m = %g        '%(t,m), end='')
                    self.t = t
                    self.m = m
                    print('%g\t%g\t'%(t,m),end='',file=dat)
                    for gap in np.arange(self.q-1):
                        kxmin,mgap,kxitr = self.min_gap_n(gap,ny,nx_premin)
                        print('%g\t%d\t'%(mgap,kxitr), end='', file=dat)
                    dat.write('\n')
                dat.write('\n')
                dat.flush()

###    def min_gap_n_ky_recursive(self,n,nx,ky,tolerance=1e-6,zone=1):
###        """ find the minimal of the nth gap of given ky """
###        k_low, k_high = 0, 2*np.pi
###        kgrid = np.linspace(k_low, k_high, nx, endpoint=True)
###        dk = (k_high - k_low) / (nx - 1)
###        
###        all_kgap = [ self.gap_n(kx,ky,n) for kx in kgrid ]
###        idx = np.argmin(all_kgap)
###
###        itr = 1
###        mgap = all_kgap[idx]
###        kmin = kgrid[idx]
###        k_low = kmin - zone*dk
###        k_high = kmin + zone*dk
###        dk = (k_high - k_low)/(nx -1)
###
###        mgap_diff = mgap
###        while mgap_diff > tolerance:
###            itr += 1
###            kgrid = np.linspace(k_low, k_high, nx, endpoint=True)
###            all_kgap = [ self.min_gap_k(kx,ky) for kx in kgrid ]
###            idx = np.argmin(all_kgap)
###
###            mgap_old = mgap
###            mgap = all_kgap[idx]
###            mgap_diff = np.abs(mgap_old - mgap)
###            kmin = kgrid[idx]
###            k_low = kmin - zone*dk
###            k_high = kmin + zone*dk
###            dk = (k_high - k_low)/(nx-1)
###
###        return (mgap,kmin,itr)
###    def min_gap_n_recursive(self,n,nx,ny,tolerance = 1e-6, zone=1):
###        kgrid = np.linspace(0,np.pi*2,ny,endpoint=False)
###        res = [ self.min_gap_n_ky_recursive(n,nx,ky,tolerance,zone) for ky in kgrid ]
###        kxgap = np.array(res)[:,0]
###
###        idx = np.argmin(kxgap)
###        return res[idx]
###    def min_gap_recursive_2(self,nx,ny,tolerance=1e-6,zone=1):
###        """ calculate the minimal value of each gap, then return the minimal among all gaps """
###        kgrid = np.linspace(0,np.pi*2,ny,endpoint=False)
###        res = [ self.min_gap_n_recursive(n, nx,ny,tolerance,zone) for n in range(self.q-1) ]
###        kxgap = np.array(res)[:,0]
###
###        idx = np.argmin(kxgap)
###        return res[idx]
###
###    def min_gap_k(self,kx,ky):
###        """ find the minimal gap at a k point """
###        eigs = la.eigvalsh(self.hk(kx,ky))
###        gaps = [ eigs[i+1] - eigs[i] for i in range(len(eigs) - 1) ]
###        return np.min(gaps)
###    def min_gap(self,nx,ny):
###        kgrid = npext.n2kgrid([nx,ny])
###        all_kgap = [ self.min_gap_k(kx,ky) for kx,ky in kgrid ]
###        return np.min(all_kgap)
###    def min_gap_ky_recursive(self,nx,ky,tolerance = 1e-6, zone=10):
###        """find minimal gap recursively in the kx direction, until
###
###        convergence reaches tolerance
###
###        At each stage, search for the minimal gap on a nx x ny lattice
###        within recursively smaller boundary
###
###        1: 0 to 2pi -> kmin
###        2: kmin - dk*zone to kmin + dk*zone -> kmin'
###        etc.
###        """
###        k_low, k_high = 0, 2*np.pi
###        kgrid = np.linspace(k_low, k_high, nx, endpoint=True)
###        dk = (k_high - k_low) / (nx - 1)
###        
###        all_kgap = [ self.min_gap_k(kx,ky) for kx in kgrid ]
###        idx = np.argmin(all_kgap)
###        #if (idx == 0) or (idx == (nx-1)):
###        #    k_low = kgrid[-2] - 2*np.pi   # wrap around 0 and 2pi
###        #    k_high = kgrid[1]
###        #else:
###        #    k_low = kgrid[idx - 1]
###        #    k_high = kgrid[idx + 1]
###
###
###        itr = 1
###        mgap = all_kgap[idx]
###        kmin = kgrid[idx]
###        k_low = kmin - zone*dk
###        k_high = kmin + zone*dk
###        dk = (k_high - k_low)/(nx -1)
###
###        mgap_diff = mgap
###        while mgap_diff > tolerance:
###            itr += 1
###            kgrid = np.linspace(k_low, k_high, nx, endpoint=True)
###            all_kgap = [ self.min_gap_k(kx,ky) for kx in kgrid ]
###            idx = np.argmin(all_kgap)
###
###            mgap_old = mgap
###            mgap = all_kgap[idx]
###            mgap_diff = np.abs(mgap_old - mgap)
###            kmin = kgrid[idx]
###            k_low = kmin - zone*dk
###            k_high = kmin + zone*dk
###            dk = (k_high - k_low)/(nx-1)
###
###            ## # if minimal occurs at the lower bound, then keep the lower bound
###            ## # as the new lower bound. Otherwise, use 1 grid prior as the lower
###            ## # bound
###            ## idx_low = idx - 1
###            ## k_low = kgrid[ 0 if idx == 0 else idx - 1 ]
###            ## # if minimal occurs at the upper bound, then keep the upper bound
###            ## # as the new upper bound. Otherwise use 1 grid next as the upper
###            ## # bound
###            ## k_high = kgrid[ idx if idx == nx-1 else idx+1 ]
###        return (mgap,kmin,itr)
###        
###    def min_gap_recursive(self,nx,ny,tolerance = 1e-6, zone=10):
###        kgrid = np.linspace(0,np.pi*2,ny,endpoint=False)
###        res = [ self.min_gap_ky_recursive(nx,ky,tolerance,zone) for ky in kgrid ]
###        kxgap = np.array(res)[:,0]
###
###        idx = np.argmin(kxgap)
###        return res[idx]
###        
###    def export_min_gap(self,nx=100,ny=2,tmin=-4,tmax=4,dt=0.01,mmin=-1,mmax=1,dm=0.1,fn='min_gap'):
###        # it is good to keep ny=2, because it seems like gap collapsing always occurs at ky = 0 or pi
###        with open(fn + '.par', 'w') as par:
###            txt = """
###            p = %d
###            q = %d
###            nx = %d
###            ny = %d
###            tmin = %g
###            tmax = %g
###            dt = %g
###            mmin = %g
###            mmax = %g
###            dm = %g
###            """%(self.p,self.q,nx,ny,tmin,tmax,dt,mmin,mmax,dm)
###            print(txt,file=par)
###
###        labels = 't m min_gap'
###        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
###        
###        with open(fn + '.dat', 'w') as dat:
###            print(header, file=dat)
###            for m in np.arange(mmin,mmax,dm):
###                for t in np.arange(tmin,tmax,dt):
###                    print('\rt = %g  m = %g        '%(t,m), end='')
###                    self.t = t
###                    self.m = m
###                    mgap = self.min_gap(nx,ny)
###                    print('%g\t%g\t%g'%(t,m,mgap), file=dat)
###                dat.write('\n')
###                dat.flush()

###    def export_min_gap_recursive(self,nx=50,ny=2,tmin=0,tmax=4,dt=0.01,mmin=0,mmax=6,dm=0.1,fn='min_gap',tolerance=1e-6, zone=1):
###        """ for each ky, look for minimal gap using the min_gap_ky_recursive """
###        # it is good to keep ny=2, because it seems like gap collapsing always occurs at ky = 0 or pi
###        with open(fn + '.par', 'w') as par:
###            txt = """
###            p = %d
###            q = %d
###            nx = %d
###            ny = %d
###            tmin = %g
###            tmax = %g
###            dt = %g
###            mmin = %g
###            mmax = %g
###            dm = %g
###            tolerance = %g
###            zone = %d
###            """%(self.p,self.q,nx,ny,tmin,tmax,dt,mmin,mmax,dm,tolerance, zone)
###            print(txt,file=par)
###
###        labels = 't m gap1 kx_itr1 gap2 kx_itr2 ...'
###        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
###        
###        with open(fn + '.dat', 'w') as dat:
###            print(header, file=dat)
###            for m in np.arange(mmin,mmax,dm):
###                for t in np.arange(tmin,tmax,dt):
###                    print('\rt = %g  m = %g        '%(t,m), end='')
###                    self.t = t
###                    self.m = m
###                    print('%g\t%g\t'%(t,m),end='',file=dat)
###                    for gap in np.arange(self.q-1):
###                        mgap,kxmin,kxitr = self.min_gap_n_recursive(gap,nx,ny,tolerance,zone)
###                        print('%g\t%d\t'%(mgap,kxitr), end='', file=dat)
###                    dat.write('\n')
###                dat.write('\n')
###                dat.flush()
        
###    def export_min_gap_n_recursive(self,n,nx=50,ny=2,tmin=0,tmax=4,dt=0.01,mmin=0,mmax=6,dm=0.1,fn='min_gap',tolerance=1e-6, zone=1):
###        """ for each ky, look for minimal nth gap using the min_gap_n_ky_recursive """
###        # it is good to keep ny=2, because it seems like gap collapsing always occurs at ky = 0 or pi
###        with open(fn + '.par', 'w') as par:
###            txt = """
###            p = %d
###            q = %d
###            nx = %d
###            ny = %d
###            gap = %d
###            tmin = %g
###            tmax = %g
###            dt = %g
###            mmin = %g
###            mmax = %g
###            dm = %g
###            """%(self.p,self.q,nx,ny,n,tmin,tmax,dt,mmin,mmax,dm)
###            print(txt,file=par)
###
###        labels = 't m min_gap kx_itr'
###        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
###        
###        with open(fn + '.dat', 'w') as dat:
###            print(header, file=dat)
###            for m in np.arange(mmin,mmax,dm):
###                for t in np.arange(tmin,tmax,dt):
###                    print('\rt = %g  m = %g        '%(t,m), end='')
###                    self.t = t
###                    self.m = m
###                    mgap,kxmin,kxitr = self.min_gap_n_recursive(n,nx,ny,tolerance,zone)
###                    print('%g\t%g\t%g\t%d'%(t,m,mgap,kxitr), file=dat)
###                dat.write('\n')
###                dat.flush()

def memo_u(nkx,nky,model,kxr=(0,2*np.pi), kyr=(0,2*np.pi),swap_xy=False,fix_gauge_elt=0,basis=None, epsilon_k=1e-6):
    """nkx,nky: k-space mesh size. Note that for kxr from 0 to 2pi, e.g., real space lattice size is actually nkx-1, because in k space, both 0 and 2pi are included.

    kxr: range of kx, inclusive on both limits
    kyr: range of ky, inclusive on both limits
    model: something that supports model.hk(kx,ky)
    fix_gauge_elt: if < 0, then do not fix, otherwise, pass it as elt to fix_gauge_elt()
    basis: if not None, then use basis[x,y] as the basis of the eigenstates

    epsilon_k: shift of k-mesh. In very rare cases, something numerically sick may happen at k=0, pi, 2pi, etc.. This is resolved by shifting the k-mesh by a tiny amount

    Example of usefulness of epsilon_k: Consider Hofstadter with
    p/q=1/3, t=0.5, m=0. Let uu = memo_u(50,50, epsilon_k=0). Then
    patch_chern(uu)[0] -> [0,-1,-1], whereas the correct one should be
    [2,-1,-1]. Setting epsilon_k = 1e-10, say, resolves this
    problem. The problem is because numerically, exp(ix2pi) is not
    exactly 1, so in the algorithm of patch_chern, the top and bottom
    edges are not exactly conj to each other at two spots, where the
    wavefunction overlap analytically should be -1, but differ in the
    top and bottom edge, causing one phase to be pi, the other to be
    -pi, hence they don't cancel each other.
    """
    hk00 = model.hk(0,0)
    nband = hk00.shape[0]
    res = np.zeros((nkx,nky,nband,nband),dtype=np.complex)

    for x,kx in enumerate(np.linspace(*kxr, nkx, endpoint=True) + epsilon_k):
        for y,ky in enumerate(np.linspace(*kyr, nky, endpoint=True) + epsilon_k):
            print('\rmemo_u: x, y = %d, %d\t        '%(x,y), end='')
            if swap_xy == True:
                h = model.hk(ky,kx)
            else:
                h = model.hk(kx,ky)
            eig,u = la.eigh(h)
            if basis != None:
                u = np.dot(npext.dagger(basis[x,y]), u)
            if fix_gauge_elt >= 0: # fix the gauge in the relative basis
                u = fix_gauge(u, fix_gauge_elt)
            res[x,y] = u
    return res
def hofc(t,m,nkx=50,nky=50):
    hof = Hof(1,3,t,m)
    uu = memo_u(nkx,nky,hof)
    print('')
    res = patch_chern(uu)[0]
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

def memo_u_smooth(nkx,nky,model,swap_xy=False,epsilon_k=1e-6,kxr=(0,2*np.pi), kyr=(0, 2*np.pi)):
    """memoize u with smooth gauge, that is, for all kx, the gauge is

    smooth and periodic along ky, and for ky=0, the gauge is smooth and
    periodic along kx.

    nkx,nky: k space mesh size, including both 0 and 2pi

    """
    hk00 = model.hk(0,0)
    nband = hk00.shape[0]
    u = memo_u(nkx,nky,model,kxr=kxr,kyr=kyr,swap_xy=swap_xy, fix_gauge_elt=-1, epsilon_k=epsilon_k)
    
    u[:,0] = smooth_path(u[:,0])[0]

    ref_berry_phase = np.zeros(nband)
    for x,u_kx in enumerate(u):
        u[x],ref_berry_phase = smooth_path(u_kx, ref_berry_phase)
    return u
        
def patch_chern_quench(imod,fmod,nkx,nky,fix_gauge_elt,swap_xy=False,kxr=(0, np.pi*2), kyr=(0, np.pi*2), epsilon_k=1e-6):
    """ calculate the 'quench chern number'

    imod: initial model
    fmod: final model

    nkx and nky are k-space mesh size, not related to real space lattice size. If kx includes both 0 and 2pi, then lattice size is instead nkx-1    
    """
    uf = memo_u_smooth(nkx=nkx,nky=nky,model=fmod,swap_xy=swap_xy, kxr=kxr, kyr=kyr, epsilon_k=epsilon_k)
    ui = memo_u(nkx=nkx,nky=nky,model=imod,swap_xy=swap_xy,fix_gauge_elt=fix_gauge_elt,basis=uf, kxr=kxr, kyr=kyr, epsilon_k=epsilon_k)
    return (patch_chern(ui), uf,ui)

def oneband_quench_coef(kx,ky,imod,fmod,bidx):
    """
    compute the single band wf coefficients after quench.
    
    kx,ky: momenta
    imod,fmod: initial and final model
    bidx: band index to be filled initially
    """
    hi = imod.hk(kx,ky)
    hf = fmod.hk(kx,ky)

    ui = la.eigh(hi)[1][:,bidx]
    uf = la.eigh(hf)[1]

    return np.dot(npext.dagger(uf),ui)

def export_oneband_quench_coef(imod,fmod,bidx=0,nx=50,ny=50,kxfrac=1,kyfrac=1, fn='qcoef'):
    """
    cf. oneband_quench_coef
    nx,ny: lattice size
    kxfrac: fraction of the kx axis to compute. e.g.: 3 means kx ranges from [0,2pi]/3, etc
    kyfrac: similar to kxfrac
    """

    if bidx > imod.q:
        print('Ensure bidx < q')
        return -1
    
    with open(fn + '.par', 'w') as par:
        txt = """
        nx = %d
        ny = %d
        kxfrac = %g
        kyfrac = %g
        q  = %d
        p_ii = %d
        t_ii = %g
        p_f = %d
        t_f = %g
        m_i = %g
        m_f = %g
        bidx = %d
        """%(nx,ny,kxfrac,kyfrac,imod.q,imod.p,imod.t,fmod.p,fmod.t,imod.m,fmod.m,bidx)
        print(txt, file=par)


    labels = 'kx ky max_abs abs(coef1) angle(coef1) abs(coef2) angle(coef2) ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
    
    dkx = np.pi * 2/(nx*kxfrac)
    dky = np.pi * 2/(ny*kyfrac)
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        for x in np.arange(nx):
            kx = x*dkx
            for y in np.arange(ny):
                ky = y*dky
                print('\r nx = %d / %d    ny = %d / %d     '%(x,nx,y,ny),end='')
                coef = oneband_quench_coef(kx,ky,imod,fmod,bidx)
                amp = np.abs(coef)
                ang = np.angle(coef)
                maxamp = np.max(amp)
                print('%g\t%g\t%g\t%s'%(kx,ky,maxamp,'\t'.join([ '%g\t%g'%(rr,tt) for rr,tt in zip(amp,ang)])), file=dat)
            dat.write('\n')
    return None

def export_rel_u_and_curvature(nkx,nky,imod,fmod,
                           swap_xy=False,fix_gauge_elt=0,
                           fn_wf = 'wf', fn_curv = 'curv', fn_conn = 'conn', kxr=(0,2*np.pi), kyr=(0,2*np.pi), epsilon_k=1e-6):
    """
    export relative wavefunction magnitude, corresponding curvature, and connection vector, of all bands of 'model'

    nkx and nky are k-space mesh size, not related to real space lattice size. If kx includes both 0 and 2pi, then lattice size is instead nkx-1
    """

    pc,uf,u = patch_chern_quench(imod,fmod,nkx,nky,fix_gauge_elt,swap_xy,kxr,kyr,epsilon_k)
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

def memo_proj_kx(model,kx,ny=50,ky0=0,swap_xy=False):
    dky = np.pi * 2 / ny
    q=model.q
    all_p = np.zeros((q,ny,q,q),dtype=np.complex)
    for y in range(ny):
        ky = dky * y + ky0
        if swap_xy:
            hh = model.hk(ky,kx)
        else:
            hh = model.hk(kx,ky)
        eig,u = la.eigh(hh)
        for b in range(q):
            psi = u[:,b]
            all_p[b,y] = np.outer(psi,psi.conj())
    return all_p
def wilson_loop(p_path):
    res = np.diag(np.ones(p_path.shape[1], dtype=np.complex))
    for p in p_path:
        res = np.dot(res, p)
    return res
    
def export_berry_phase(model,nx=50,ny=50,ky0=0,fn='berry',swap_xy=False):
    with open(fn + '.par', 'w') as par:
        txt = """
        nx = %d
        ny = %d
        p = %d
        q = %d
        t = %g
        m = %g
        ky0 = %g
        swap_xy = %d
        """%(nx,ny,model.p,model.q,model.t,model.m,ky0, swap_xy)
        print(txt, file=par)
    labels = 'kx abs1 ang1 abs2 ang2 ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
    
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        dkx = 2*np.pi/nx
        for x in np.arange(nx):
            kx = dkx * x
            print('%g\t'%kx, file=dat, end='')
            print('\rkx = %g            '%kx, end='')

            all_p = memo_proj_kx(model,kx,ny,ky0,swap_xy)
            for p_band in all_p:
                ww = wilson_loop(p_band)
                eig = np.trace(ww)
                r,theta = cmath.polar(eig)
                print('%g\t%g\t'%(r,theta), file=dat, end='')
            dat.write('\n')
    return None

################################################################
## non-generic stuff
################################################################
def patch_chern_quench_mat(ti,mi,tf,mf,nkx=50,nky=None):
    if nky == None:
        nky = nkx
    ihof = Hof(p=1,q=3,t=ti,m=mi)
    fhof = Hof(p=1,q=3,t=tf,m=mf)

    res = np.zeros((3,3),dtype=np.double)
    for row in range(3):
        res[row] = patch_chern_quench(ihof,fhof,nkx,nky,row)[0][0]
    return res
def wf_gen(ihof,fhof,ky,iband,fbasis):
    """ generates a function that will return the wavefunction overlap between two Hamiltonians, as a function of kx """
    def res(kx):
        ui = la.eigh(ihof.hk(kx,ky))[-1][:,iband]
        uf = la.eigh(fhof.hk(kx,ky))[-1][:,fbasis]
        return np.abs(np.dot(uf.conj(), ui))
    return res

def wf_maxmin(ihof,fhof,iband,all_ky=(0,np.pi), nx_premin=50):
    """compute the minimal for each component of the relative wavefunction(initial state on the final basis), and return the max of these minima.

    If the max is zero, then all components touch zero at some point,
    hence there is protected DPT


    nx_premin: mesh size for computing the brent() bracket. [0,2pi] is divided into nx_premin segments, and the bracket is the one around the minimal
    """

    res = np.zeros(ihof.q, dtype=np.double)
    for fbasis in range(ihof.q):
        wf_all_ky = [ wf_gen(ihof,fhof,ky,iband,fbasis) for ky in all_ky ]
        min_all_ky = [ brent_premin(wf, (0,2*np.pi), nx_premin)[1] for wf in wf_all_ky ]
        res[fbasis] = np.min(min_all_ky)
    return np.max(res)

def export_wf_maxmin(ihof=Hof(t=0.25,m=1),
                  mf=1,tf_from=0,tf_to=4,dtf=0.1,
                  all_ky=(0,np.pi),
                  nx_premin=50,
                  fn='wf_maxmin'):
    with open(fn + '.par', 'w') as par:
        txt = """
        p = %d
        q = %d

        ti = %g
        mi = %g

        mf = %g
        tf_from = %g
        tf_to = %g
        dtf = %g
        """%(ihof.p, ihof.q, ihof.t, ihof.m, mf, tf_from, tf_to, dtf)
        print(txt,file=par)
    labels = 't iband1 iband2 ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])

    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)

        fhof = Hof(p=ihof.p, q=ihof.q, m=mf)
        for t in np.arange(tf_from,tf_to,dtf):
            fhof.t = t
            dat.write('%g\t'%t)
            for iband in np.arange(ihof.q):
                print('\rtf = %g\t\tiband = %d            '%(t,iband),end='')
                res = wf_maxmin(ihof,fhof,iband=iband,all_ky = all_ky, nx_premin=nx_premin)
                dat.write('%g\t'%res)
            dat.write('\n')
    
                
    
