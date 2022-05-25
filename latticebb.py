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
from scipy.special import binom as binom
import kpath

one2 = np.identity(2,dtype=np.complex)
one4 = np.identity(4, dtype=np.complex)
sx = np.array([[0,1],[1,0]],dtype=np.complex)
sy = np.array([[0,-1j],[1j,0]],dtype=np.complex)
sz = np.array([[1,0],[0,1]],dtype=np.complex)

gg1 = np.kron(sz,sy)
gg2 = -np.kron(sz,sx)
gg4 = np.kron(sx,one2)
gg12 = 1j * np.dot(gg1,gg2)

################################################################
## Chern number stuff. Should move to external files later
def chern(u):
    """Compute Chern numbers of band wave functions u, where u[nx,ny] is
    the unitary matrix that diagonalizes H(kx,ky)

    u must have shape (nx,ny,nband,nband)
    """
    nx,ny,nband = np.shape(u)[:3]

    res = np.zeros(nband)

    for x in range(nx):
        for y in range(ny):
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
################################################################

def hk(kx,ky,ds,eta):
    """ k-space latticebb model with Del_D = 0, making it effectively 2D. See latticebb.pdf, p1.

    """
    ssx,ssy = np.sin([kx,ky])
    ccx,ccy = np.cos([kx,ky])
    res = ssx * gg1 + ssy * gg2 + ds * (3 - ccx - ccy) * gg4 + eta * gg12
    return res
    
def export_hk(nx,ny,ds,eta,h=hk,fn='hk'):
    with open(fn + '.par', 'w') as par:
        txt = """
        nx = %d
        ny = %d
        ds = %g
        eta = %g
        """%(nx,ny,ds,eta)
        print(txt, file=par)

    labels = 'kx ky erg1 erg2 ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])


    dkx = np.pi * 2 / nx
    dky = np.pi * 2 / ny
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        for x in np.arange(nx):
            kx = x*dkx
            for y in np.arange(ny):
                ky = y*dky
                print('\r nx = %d / %d    ny = %d / %d     '%(x,nx,y,ny),end='')
                hh = h(kx,ky,ds,eta)
                eig = la.eigvalsh(hh)
                print('%g\t%g\t%s'%(kx,ky,'\t'.join([ '%g'%erg for erg in eig ])), file=dat)
            dat.write('\n')
    print('\n')
    print(chern_bb(ds,eta,nx,ny))
#def export_hk_k(kx,ky,p,q,tlist=np.arange(-5,5,0.005),fn='hk_k'):
#    with open(fn + '.par', 'w') as par:
#        txt = """
#        kx = %g
#        ky = %g
#        p = %d
#        q = %d
#        """%(kx,ky,p,q)
#        print(txt, file=par)
#
#    labels = 't erg1 erg2 ...'
#    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
#
#
#    with open(fn + '.dat', 'w') as dat:
#        print(header, file=dat)
#        for t in tlist:
#            print('\r t = %g      '%t,end='')
#            hh = hk(kx,ky,p,q,t)
#            eig = la.eigvalsh(hh)
#            print('%g\t%s'%(t,'\t'.join([ '%g'%erg for erg in eig ])), file=dat)
def export_erg_vs_t(p,q,tmax=2,dt=0.005,h=hk,fn='erg_vs_t'):
    with open(fn + '.par', 'w') as par:
        txt = """
        p = %d
        q = %d
        """%(p,q)
        print(txt, file=par)

    labels = 'kx ky t erg1 erg2 ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])


    tlist = np.arange(-tmax,tmax,dt)
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        for kx in (0, np.pi/q):
            for ky in (0,np.pi):
                for t in tlist:
                    print('\rkx, ky = %g,%g\t t = %g      '%(kx,ky,t),end='')
                    hh = h(kx,ky,p,q,t)
                    eig = la.eigvalsh(hh)
                    print('%g\t%g\t%g\t%s'%(kx,ky,t,'\t'.join([ '%g'%erg for erg in eig ])), file=dat)
                dat.write('\n')

def memo_u(nx,ny,ds,eta,h):
    res = np.zeros((nx,ny,4,4),dtype=np.complex)
    dkx = np.pi*2/nx
    dky = np.pi*2/ny
    for x in range(nx):
        kx = dkx * x
        for y in range(ny):
            ky = dky * y
            hh = h(kx,ky,ds,eta)
            eig,u = la.eigh(hh)
            res[x,y] = u
    return res
def chern_bb(ds,eta,nx=50,ny=50,h=hk):
    """ compute the Chern numbers of the latticebb model, cf. def hk() """
    uu = memo_u(nx,ny,ds,eta,h)
    return chern(uu)

def oneband_quench_coef(kx,ky,q,pi,ti,pf,tf,bidx):
    """
    compute the single band wf coefficients after quench.
    
    kx,ky: momenta
    q: for both H before and after quench, in order for the number of bands to remain the same
    pi,ti: p,q,t before quench
    pf,tf: p,q,t after quench
    bidx: band index to be filled initially
    """
    hi = hk(kx,ky,pi,q,ti)
    hf = hk(kx,ky,pf,q,tf)

    ui = la.eigh(hi)[1][:,bidx]
    uf = la.eigh(hf)[1]

    return np.dot(npext.dagger(uf),ui)

def export_oneband_quench_coef(nx,ny,kxfrac,kyfrac,q,pi,pf,ti,tf,bidx, fn='qcoef'):
    """
    cf. oneband_quench_coef
    nx,ny: lattice size
    kxfrac: fraction of the kx axis to compute. e.g.: 3 means kx ranges from [0,2pi]/3, etc
    kyfrac: similar to kxfrac
    """

    if bidx > q:
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
        bidx = %d
        """%(nx,ny,kxfrac,kyfrac,q,pi,ti,pf,tf,bidx)
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
                coef = oneband_quench_coef(kx,ky,q,pi,ti,pf,tf,bidx)
                amp = np.abs(coef)
                ang = np.angle(coef)
                maxamp = np.max(amp)
                print('%g\t%g\t%g\t%s'%(kx,ky,maxamp,'\t'.join([ '%g\t%g'%(rr,tt) for rr,tt in zip(amp,ang)])), file=dat)
            dat.write('\n')
    return None

def export_oneband_quench_coef_kline(q,pi,pf,ti,tf,bidx,knodes=[[0,0],[0,1],[1,1],[0,0]],dk=0.01, fn='qcoef_kline'):
    """
    cf. oneband_quench_coef
    knodes, dk: for kx, in unit of pi/q; for ky, in unit of pi
    """

    if bidx > q:
        print('Ensure bidx < q')
        return -1
    
    kk = kpath.kpath(knodes,dk)
    def ticname(kpoint):
        kx,ky = kpoint
        res = ''
        if kx == 0:
            res = """'$(0\ , \ %g)$'"""%(ky)
        else:
            res = """'$(\\\\frac{%g}{%d}\ , \ %g)$'"""%(kx,q,ky)
        return res
    tics,info,nk = kpath.kpath_aux(knodes,dk,ticname)
    with open(fn + '.par', 'w') as par:
        txt = """
        tics = "%s"
        info = "%s"
        nk = %d
        q  = %d
        p_ii = %d
        t_ii = %g
        p_f = %d
        t_f = %g
        bidx = %d
        """%(tics,info,nk,q,pi,ti,pf,tf,bidx)
        print(txt, file=par)


    labels = 'kx ky max_abs abs(coef1) angle(coef1) abs(coef2) angle(coef2) ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
    
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        for xx,yy in kk:
            kx = xx*np.pi/q
            ky = yy*np.pi
            print('\r kx = %g    ky = %g     '%(kx,ky),end='')
            coef = oneband_quench_coef(kx,ky,q,pi,ti,pf,tf,bidx)
            amp = np.abs(coef)
            ang = np.angle(coef)
            maxamp = np.max(amp)
            print('%g\t%g\t%g\t%s'%(kx,ky,maxamp,'\t'.join([ '%g\t%g'%(rr,tt) for rr,tt in zip(amp,ang)])), file=dat)
    return None


    
def max_quench_coef(ui,uf):
    """ compute the maximal quench coefficient (wavefunction mod squared)

    ui: initial wave functions. shape(ui) = (ntotal, nfilled) for nfilled bands
    uf: final (quenched) wavefunctions. shape(uf) = (ntotal,ntotal)
    """
    ntotal,nfilled = np.shape(ui)

    # projection onto initially filled bands
    proji = np.dot(ui, npext.dagger(ui))
    # proji in the quenched basis
    uu = np.dot(npext.dagger(uf),np.dot(proji, uf))

    
    res = -1
    max_idx = None
    for final_bands in it.combinations(range(ntotal),nfilled):
        idx = np.ix_(final_bands,final_bands)
        coef = np.abs(la.det(uu[idx]))
        if coef > res:
            res = coef
            max_idx = final_bands
    return(res,max_idx)
def export_max_quench_coef(nx,ny,kxfrac,kyfrac,q,pi,pf,ti,tf,bands,fn='qcoef-max'):
    """
    nx,ny: k-space lattice size
    kxfrac, kyfrac: fraction in k axis to be computed. e.g., kxfrac=3 means kx runs over [0,2pi]/3
    q: muc size
    pi,pf: initial and final p
    ti,tf: initial and final t
    bands: indices of initially filled bands
    """

    with open(fn + '.par', 'w') as par:
        txt = """
        nx = %d
        ny = %d
        q  = %d
        pi = %d
        ti = %g
        pf = %d
        tf = %g
        bands = %s
        """%(nx,ny,q,pi,ti,pf,tf,bands)
        print(txt, file=par)

    labels = 'kx ky max_coef target_bands'
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
                hi = hk(kx,ky,pi,q,ti)
                hf = hk(kx,ky,pf,q,tf)
                ui = la.eigh(hi)[1][:,bands]
                uf = la.eigh(hf)[1]

                max_coef,target_bands = max_quench_coef(ui,uf)
                print('%g\t%g\t%g\t%s'%(kx,ky,max_coef,target_bands),file=dat)
            dat.write('\n')
    return None
