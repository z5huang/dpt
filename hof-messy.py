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
from chern import chern

################################################################
## Chern number stuff. Should move to external files later
#def chern(u):
#    """Compute Chern numbers of band wave functions u, where u[nx,ny] is
#    the unitary matrix that diagonalizes H(kx,ky)
#
#    u must have shape (nx,ny,nband,nband)
#    """
#    nx,ny,nband = np.shape(u)[:3]
#
#    res = np.zeros(nband)
#
#    for x in range(nx):
#        for y in range(ny):
#            # 2 - 1 (x,y)
#            # |   |
#            # 3 - 4
#            
#            u1 = u[x,y]
#            u2 = u[x-1,y]
#            u3 = u[x-1,y-1]
#            u4 = u[x,y-1]
#
#            for b in np.arange(nband):
#                u12 = np.dot(u1[:,b].conj(), u2[:,b])
#                u23 = np.dot(u2[:,b].conj(), u3[:,b])
#                u34 = np.dot(u3[:,b].conj(), u4[:,b])
#                u41 = np.dot(u4[:,b].conj(), u1[:,b])
#                res[b] += np.angle(u12 * u23 * u34 * u41)
#    return res/(2*np.pi)
################################################################

def hk(kx,ky,p,q,t=0):
    """ k-space Hofstadter Hamiltonian on square lattice, flux per plaquette is 2pi*p/q

    t: next nearest hopping, i.e. diagonals of the square plaquette

    """
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
def hk_tri(kx,ky,p,q,t=0):
    """ k-space Hofstadter Hamiltonian on square lattice with triangular hopping, flux per plaquette is 2pi*p/q

    t: triangular next nearest hopping, i.e. one diagonal of each square plaquette

    """
    res = np.zeros((q,q),dtype=np.complex)
    i,j = np.indices((q,q))
    phi = 2*np.pi*p/q
    res[i==j] = 2*np.cos(np.arange(1,q+1) * phi + kx)

    sup_diag = 1 + t*np.exp(1j * (np.arange(1.5,q+0.5) * phi + kx))
    res[i==j-1] = np.conjugate(sup_diag)
    res[i==j+1] = sup_diag

    corner = 1 + t*np.exp(1j * (0.5*phi + kx))
    res[0,-1]+= np.exp(1j * ky) * corner
    res[-1,0]+= np.exp(-1j * ky) * np.conjugate(corner)
    return res
def hk_no_diag_flux(kx,ky,p,q,t=0):
    """ k-space Hofstadter Hamiltonian on square lattice, flux per plaquette is 2pi*p/q

    t: next nearest hopping, i.e. diagonals of the square plaquette

    """
    res = np.zeros((q,q),dtype=np.complex)
    i,j = np.indices((q,q))
    phi = 2*np.pi*p/q
    res[i==j] = 2*np.cos(np.arange(1,q+1) * phi + kx)

    sup_diag = 1 + 2*t*np.cos(kx)
    res[i==j-1] = sup_diag
    res[i==j+1] = sup_diag

    corner = 1 + 2*t*np.cos(kx)
    res[0,-1]+= np.exp(1j * ky) * corner
    res[-1,0]+= np.exp(-1j * ky) * corner
    return res

def export_hk(nx,ny,p,q,t,h=hk,fn='hk', export_wf = None):
    with open(fn + '.par', 'w') as par:
        txt = """
        nx = %d
        ny = %d
        p = %d
        q = %d
        t = %g
        """%(nx,ny,p,q,t)
        print(txt, file=par)

    labels = 'kx ky erg1 erg2 ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])


    dkx = np.pi * 2 / nx
    dky = np.pi * 2 / ny
    if export_wf == True:
        dat_wf = open(fn + '_wf.dat', 'w')
        labels_wf = 'kx ky <1|psi1>_abs <1|psi1>_phase <2|psi1>_abs <2|psi1>_phase ... <1|psi2>_abs <1|psi2>_phase <2|psi2>_abs ...'
        header_wf = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels_wf.split())])
        print(header_wf, file=dat_wf)
    
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        for x in np.arange(nx):
            kx = x*dkx
            for y in np.arange(ny):
                ky = y*dky
                print('\r nx = %d / %d    ny = %d / %d     '%(x,nx,y,ny),end='')
                hh = h(kx,ky,p,q,t)
                if export_wf == True:
                    eig,u = la.eigh(hh)
                    print('%g\t%g\t%s'%(kx,ky,  '\t'.join([ '%g\t%g'%(np.abs(wf),np.angle(wf)) for wf in u.transpose().flatten() ])), file=dat_wf)
                else:
                    eig = la.eigvalsh(hh)
                print('%g\t%g\t%s'%(kx,ky,'\t'.join([ '%g'%erg for erg in eig ])), file=dat)
            dat.write('\n')
            if export_wf == True:
                dat_wf.write('\n')
    if export_wf == True:
        dat_wf.close()
    print('\n')
    print(chern_hof(p,q,t,nx,ny))
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

def memo_u(nx,ny,p,q,t,h,kxfrac=1,kyfrac=1):
    res = np.zeros((nx,ny,q,q),dtype=np.complex)
    dkx = np.pi*2/(nx*kxfrac)
    dky = np.pi*2/(ny*kyfrac)
    for x in range(nx):
        kx = dkx * x
        for y in range(ny):
            ky = dky * y
            hh = h(kx,ky,p,q,t)
            eig,u = la.eigh(hh)
            res[x,y] = u
    return res
def chern_hof(p,q,t=0,nx=50,ny=50,h=hk,kxfrac=1,kyfrac=1):
    """ compute the Chern numbers of the Hofstadter model, cf. def hk() """
    uu = memo_u(nx,ny,p,q,t,h,kxfrac,kyfrac)
    return chern(uu)


def parallelize(u,uref):
    """ adjust the phases in the column vectors u, so that they are parallel with the corresponding columns in uref """
    ncol = u.shape[1]

    # <u(1)|ref(1)>, <u(2)|ref(2)>, ..., <u(ncol) | ref(ncol)>
    overlaps = np.sum(np.conj(u) * uref, 0)

    phase_factors = overlaps / np.abs(overlaps)
    return u*phase_factors


def memo_u_quench(nx,ny,q,pi,pf,ti,tf,h,kxfrac=1,kyfrac=1):
    """ memoize the relative wavefunctions before and after the quench. I.e., the initial wavefunctions in the basis of the quenched states """

    ###
    # Parallel gauge for both wavefunctions. i.e., for any kx, the ky line is
    # parallel gauged. For ky = 0, the kx line is parallel gauged.
    res = np.zeros((nx,ny,q,q),dtype=np.complex)
    dkx = np.pi*2/(nx*kxfrac)
    dky = np.pi*2/(ny*kyfrac)
    # wavefunctions at the same kx and previous ky
    ui_prev_ky = np.diag(np.ones(q,dtype=np.complex))
    uf_prev_ky = np.diag(np.ones(q,dtype=np.complex))
    # wavefunctions at ky=0 and the previous kx
    ui_prev_kx = np.diag(np.ones(q,dtype=np.complex))
    uf_prev_kx = np.diag(np.ones(q,dtype=np.complex))
    for x in range(nx):
        kx = dkx * x
        for y in range(ny):
            ky = dky * y
            hi = h(kx,ky,pi,q,ti)
            eig,ui = la.eigh(hi)
            hf = h(kx,ky,pf,q,tf)
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
def chern_hof_quench(q,pi,ti,tf,pf=None,nx=50,ny=None,h=hk,kxfrac=1,kyfrac=1):
    """ compute the Chern numbers for the relative wavefunctions, i.e., of the initial states in the basis of the final states """
    if pf == None:
        pf = pi
    if ny == None:
        ny = nx
    print('pi = %g, pf = %g, ti = %g, tf = %g'%(pi,pf,ti,tf))
    uu = memo_u_quench(nx,ny,q,pi,pf,ti,tf,h,kxfrac,kyfrac)
    return chern(uu)
    


#def parallelize_all(all_u,all_uref):
#    """ adjust the phases in the column vectors u, so that they are parallel with the corresponding columns in uref
#
#    all_u = [u1, u2, ...]
#    all_uref = [uref1, uref2, ...]
#    """
#    ncol = all_u.shape[-1]
#
#    # <u(1)|ref(1)>, <u(2)|ref(2)>, ..., <u(ncol) | ref(ncol)>
#    overlaps = np.sum(np.conj(all_u) * all_uref, -2) # sum over the row index
#
#    phase_factors = (overlaps / np.abs(overlaps)).reshape(all_u.shape[0], 1, ncol)
#    return all_u*phase_factors


def parallelize_path(all_u):
    """ parallelization of a path of u's """
    ntot = all_u.shape[0]
    res = np.zeros(all_u.shape,dtype=np.complex)
    res[0] = all_u[0]
    for n in range(ntot-1):
        res[n+1] = parallelize(all_u[n+1], res[n])
    return res
        
def gen_connection(u_path, smooth_twist = True, end_with_u1 = True, ref_berry_phase = None):
    """generate the connection matrix from a path of u, u_path =

    [u1,u2,..., uN] returns [u12, u23, u34, ..., u(N,1)]. See OneNote
    for meaning of smooth_twist and end_with_u1
    """


    npoints,nbands = u_path.shape[:2]

    u_para = parallelize_path(u_path)
    u1,uN = u_para[[0,-1]]
    berry_factors = (u1.conj() * uN).sum(axis = -2)
    if ref_berry_phase == None:
        berry_phase = np.angle( berry_factors )
    else:
        phase_change = np.angle(berry_factors * np.exp(-1j*ref_berry_phase))
        berry_phase = ref_berry_phase + phase_change
    if smooth_twist:
        twist_factor = np.exp(1j * berry_phase/npoints)
    else:
        twist_factor = 1

    res = np.zeros(u_path.shape,dtype=np.complex)
    for n in np.arange(npoints-1):
        ua,ub = u_para[[n,n+1]]
        res[n] = np.dot(ua * twist_factor, npext.dagger(ub))
    if smooth_twist:
        res[-1] = np.dot(uN * (twist_factor * np.exp(-1j * berry_phase)), npext.dagger(u1))
    else:
        uend = u1
        if end_with_u1 == False:
            uend = u1 * np.exp(1j * berry_phase)
        res[-1] = np.dot(uN, npext.dagger(uend))
    return (res, berry_phase)
def rel_wilson_loop(p_path, connection):
    """ generate a wilson loop from a path of projectors p_path, and the corresponding connection
    
    p_path = [p1, p2, ..., pN],
    connection = [u12, u23, ..., u(N,1)]
    """
    res = np.diag(np.ones(p_path.shape[1], dtype=np.complex))
    for p,u in zip(p_path,connection):
        res = np.dot(res, np.dot(p,u))
    return res
        


def memo_kx(q,pi,pf,ti,tf,kx, ny=50,ky0=0,h=hk):
    "memoize initial projectors and final wavefunctions at fixed kx"
    dky = np.pi * 2 / ny
    all_p = np.zeros((q,ny,q,q),dtype=np.complex)
    all_u = np.zeros((ny,q,q), dtype=np.complex)

    for y in range(ny):
        ky = dky * y + ky0
        hi = h(kx,ky,pi,q,ti)
        eigi,ui = la.eigh(hi)
        for b in range(q):
            psi = ui[:,b]
            all_p[b,y] = np.outer(psi,psi.conj())
        hf = h(kx,ky,pf,q,tf)
        eigf,uf = la.eigh(hf)
        all_u[y] = uf
    parallelize_path(all_u)
    return (all_p,all_u)
        
def export_rel_wilson_loop(q,pi,pf,ti,tf, nx=50,ny=50, ky0=0, smooth_twist=True, end_with_u1=True, fn='wilson',h=hk):

    with open(fn + '.par', 'w') as par:
        txt = """
        nx = %d
        ny = %d
        pi = %d
        pf = %d
        q = %d
        ti = %d
        tf = %d
        ky0 = %g
        smooth_twist = %d
        end_with_u1 = %d
        """%(nx,ny,pi,pf,q,ti,tf,ky0,smooth_twist,end_with_u1)
        print(txt, file=par)

    labels = 'kx abs1 ang1 abs2 ang2 ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
    
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        dkx = 2*np.pi/nx

        ref_berry_phase = None
        for x in np.arange(nx):
            kx = dkx * x
            print('%g\t'%kx, file=dat, end='')
            print('\rkx = %g                 '%kx, end='')
            
            all_p,all_u = memo_kx(q,pi,pf,ti,tf,kx,ny,ky0,h)
            #connection, ref_berry_phase = gen_connection(all_u, smooth_twist, end_with_u1, ref_berry_phase)
            connection, ref_berry_phase = gen_connection(all_u, smooth_twist, end_with_u1)
            for p_band in all_p:
                ww = rel_wilson_loop(p_band, connection)
                #ww = np.dot(ww, p_band[0])
                eig = np.trace(ww)
                r,theta = np.abs(eig),np.angle(eig)
                print('%g\t%g\t'%(r,theta), file=dat, end='')
            dat.write('\n')
    return None
        
def memo_proj_kx(p,q,t,kx,ny=50,ky0=0,h=hk):
    dky = np.pi * 2 / ny
    all_p = np.zeros((q,ny,q,q),dtype=np.complex)
    for y in range(ny):
        ky = dky * y + ky0
        hh = h(kx,ky,p,q,t)
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
    
    
def export_berry_phase(p,q,t,nx=50,ny=50,ky0=0,fn='berry',h=hk):
    with open(fn + '.par', 'w') as par:
        txt = """
        nx = %d
        ny = %d
        p = %d
        q = %d
        t = %g
        ky0 = %g
        """%(nx,ny,p,q,t,ky0)
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

            all_p = memo_proj_kx(p,q,t,kx,ny,ky0,h)
            for p_band in all_p:
                ww = wilson_loop(p_band)
                eig = np.trace(ww)
                r,theta = cmath.polar(eig)
                print('%g\t%g\t'%(r,theta), file=dat, end='')
            dat.write('\n')
    return None
                
            

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
