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

################################################################
## Chern number stuff. Should move to external files later
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
################################################################

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

def export_hk(nx,ny,p,q,t,swap_xy=False,h=hk,fn='hk', export_wf = None):
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
                hh = h(kx,ky,p,q,t,swap_xy)
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

def export_erg_vs_t(p,q,tmax=2,dt=0.005,h=hk,fn='erg_vs_t'):
    """plot energy at the high symmetry momenta (0,0), (0,pi), (pi/q, 0)

    and (pi/q, pi) as a function of the diagonal hopping t
    """
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
                    hh = h(kx,ky,p,q,t,swap_xy=False)
                    eig = la.eigvalsh(hh)
                    print('%g\t%g\t%g\t%s'%(kx,ky,t,'\t'.join([ '%g'%erg for erg in eig ])), file=dat)
                dat.write('\n')

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
        kx = dkx * x + kx0
        for y in range(ny):
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
def chern_hof_quench(q,pi,ti,tf,pf=None,nx=50,ny=None,swap_xy=False,h=hk,kxfrac=1,kyfrac=1, kx0=0, ky0=0):
    """ compute the Chern numbers for the relative wavefunctions, i.e., of the initial states in the basis of the final states """
    if pf == None:
        pf = pi
    if ny == None:
        ny = nx
    print('pi = %g, pf = %g, ti = %g, tf = %g'%(pi,pf,ti,tf))
    uu = memo_u_quench(nx,ny,q,pi,pf,ti,tf,swap_xy,h,kxfrac,kyfrac,kx0,ky0)
    return chern(uu)

def parallelize_path(all_u):
    """ parallelization of a path of u's """
    ntot = all_u.shape[0]
    res = np.zeros(all_u.shape,dtype=np.complex)
    res[0] = all_u[0]
    for n in range(ntot-1):
        res[n+1] = parallelize(all_u[n+1], res[n])
    return res

def smooth_path(all_u, ref_berry_phase = None, u_inclusive=True):
    """smoothify a path of u's

    if ref_berry_phase == None, then just use the phases of
    berry_factors as the corresponding berry phases. Otherwise, the
    Berry phases are such that the phase change from the ref_berry_phase
    be small (i.e. in the branch (-pi,pi])

    if u_inclusive == True, it means the first and last member of all_u
    are equivalent. This would influence the calculation of the phase
    adjustment
    """
    u_para = parallelize_path(all_u)
    u1,uN = u_para[[0,-1]]
    ntot,nband = all_u.shape[0:2]

    # nk is the actual number of inequivalent k points
    if u_inclusive == True:
        nk = ntot - 1
    else:
        nk = ntot
    berry_factors = (u1.conj() * uN).sum(axis=-2)
    if ref_berry_phase == None:
        berry_phase = np.angle( berry_factors )
    else:
        phase_change = np.angle(berry_factors * np.exp(-1j*ref_berry_phase))
        berry_phase = ref_berry_phase + phase_change
    twist = np.exp(-1j * np.outer(np.arange(0,ntot), berry_phase)/nk).reshape(ntot,1,nband)
    return u_para*twist,berry_phase

def memo_u_smooth(nx,ny,p,q,t,swap_xy=False,h=hk,kx0=0,ky0=0,kxfrac=1,kyfrac=1):
    """memoize a smooth gauged eigenbasis. The gauge is smooth along the
    kx direction for all ky values, with initial point at kx0. For kx =
    kx0, the gauge is smoothed along the ky direction, with initial
    point at ky0"""

    res = np.zeros((nx+1,ny+1,q,q),dtype=np.complex)
    dkx = np.pi*2/(nx*kxfrac)
    dky = np.pi*2/(ny*kyfrac)

    # smooth the kx0 line first
    # NB: use ny + 1, so that ky = [0, dky, 2dky, ..., 2pi], i.e., inclusive of both 0 and 2pi
    for y in range(ny+1):                 
        ky = dky * y + ky0
        hh = h(kx0,ky,p,q,t,swap_xy)
        eig,u = la.eigh(hh)
        res[0,y] = u
    res[0] = smooth_path(res[0], ref_berry_phase=None)[0]

    ref_berry_phase = None
    for y in range(ny+1):
        ky = dky * y + ky0
        # NB: use nx + 1, so that kx = [0, dkx, 2dkx, ..., 2pi], i.e., inclusive of both 0 and 2pi
        for x in range(1,nx+1):             # skip kx=kx0 so starting from x=1
            kx = dkx * x + kx0
            
            hh = h(kx,ky,p,q,t,swap_xy)
            eig,u = la.eigh(hh)
            res[x,y] = u
        # smooth the ky line
        res[:,y] , ref_berry_phase = smooth_path(res[:,y], ref_berry_phase)
    return res
def chern_hof_smooth(q,pi,ti,tf,pf=None,nx=50,ny=None,swap_xy=True,h=hk,kx0=0,ky0=0,leave_out_kx_bound=True, leave_out_ky_bound=True):
    """ compute the Chern numbers for the relative wavefunctions, using a smooth-gauged final basis"""
    if pf == None:
        pf = pi
    if ny == None:
        ny = nx
    print('pi = %g, pf = %g, ti = %g, tf = %g'%(pi,pf,ti,tf))
    u_smooth = memo_u_smooth(nx,ny,pf,q,tf,swap_xy,h,kx0,ky0)
    uu = memo_u(nx,ny,pi,q,ti,swap_xy=swap_xy,h=h,kx0=kx0,ky0=ky0,basis=u_smooth)
    return chern(uu,leave_out_kx_bound, leave_out_ky_bound)
def patch_chern_hof_smooth(q,pi,ti,tf,pf=None,nx=50,ny=None,swap_xy=True,h=hk,kx0=0,ky0=0,kxfrac=1,kyfrac=1):
    """ compute the Chern numbers for the relative wavefunctions, using a smooth-gauged final basis"""
    if pf == None:
        pf = pi
    if ny == None:
        ny = nx
    print('pi = %g, pf = %g, ti = %g, tf = %g'%(pi,pf,ti,tf))
    u_smooth = memo_u_smooth(nx,ny,pf,q,tf,swap_xy,h,kx0,ky0,kxfrac,kyfrac)
    uu = memo_u(nx,ny,pi,q,ti,swap_xy=swap_xy,h=h,kx0=kx0,ky0=ky0,basis=u_smooth)
    return patch_chern(uu)
    
def gen_connection(u_path, smooth_twist = False, end_with_u1 = True, ref_berry_phase = None):
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
        


def memo_kx(q,pi,pf,ti,tf,kx, ny=50,ky0=0,swap_xy=False,h=hk):
    "memoize initial projectors and final wavefunctions at fixed kx"
    dky = np.pi * 2 / ny
    all_p = np.zeros((q,ny,q,q),dtype=np.complex)
    all_u = np.zeros((ny,q,q), dtype=np.complex)

    for y in range(ny):
        ky = dky * y + ky0
        hi = h(kx,ky,pi,q,ti,swap_xy)
        eigi,ui = la.eigh(hi)
        for b in range(q):
            psi = ui[:,b]
            all_p[b,y] = np.outer(psi,psi.conj())
        hf = h(kx,ky,pf,q,tf,swap_xy)
        eigf,uf = la.eigh(hf)
        all_u[y] = uf
    parallelize_path(all_u)
    return (all_p,all_u)
        
def export_rel_wilson_loop(q,pi,pf,ti,tf, nx=50,ny=50, ny0=None, ky00=0, smooth_twist=False, end_with_u1=True, fn='wilson',swap_xy=False, h=hk, basis_winding=False):
    """ export eigenvalues of the relative wilson loop on a mesh of kx and ky0, where ky0 is the initial point of the loop (or: the end point corresponding to the link which is not parallel transported)


    nx: number of kx points
    ny: number of ky points used for the Wilson loop
    ny0: number of ky0 mesh
    ky00: the starting point of ky0. Useful if we only want to export a single ky0 by setting ny0 = 1
    basis_winding: allow the Berry phase of the final eigenstates to also wind as a function of ky, so that it is smooth when crossing pi
    
    """
    if ny0 == None:
        ny0 = ny
    with open(fn + '.par', 'w') as par:
        txt = """
        nx = %d
        ny = %d
        ny0 = %d
        pi = %d
        pf = %d
        q = %d
        ti = %d
        tf = %d
        ky00 = %g
        smooth_twist = %d
        end_with_u1 = %d
        """%(nx,ny,ny0,pi,pf,q,ti,tf,ky00,smooth_twist,end_with_u1)
        print(txt, file=par)

    labels = 'kx ky0 abs1 ang1 abs2 ang2 ...'
    header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
    
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        dkx = 2*np.pi/nx
        dky0 = 2*np.pi/ny0

        ref_berry_phase = None
        for x in np.arange(nx):
            kx = dkx * x
            for y in np.arange(ny0):
                ky0 = dky0 * (y + ky00)
                print('%g\t%g\t'%(kx,ky0), file=dat, end='')
                print('\rkx = %g\tky0 = %g\t                 '%(kx,ky0), end='')
                
                all_p,all_u = memo_kx(q,pi,pf,ti,tf,kx,ny,ky0,swap_xy,h)
                if basis_winding:
                    connection, ref_berry_phase = gen_connection(all_u, smooth_twist, end_with_u1, ref_berry_phase)
                else:
                    connection, ref_berry_phase = gen_connection(all_u, smooth_twist, end_with_u1,ref_berry_phase=None)
                for p_band in all_p:
                    ww = rel_wilson_loop(p_band, connection)
                    #ww = np.dot(ww, p_band[0])
                    eig = np.trace(ww)
                    r,theta = np.abs(eig),np.angle(eig)
                    print('%g\t%g\t'%(r,theta), file=dat, end='')
                dat.write('\n')
            if ny0 != 1:
                dat.write('\n')
    return None
        
def memo_proj_kx(p,q,t,kx,ny=50,ky0=0,swap_xy=False,h=hk):
    dky = np.pi * 2 / ny
    all_p = np.zeros((q,ny,q,q),dtype=np.complex)
    for y in range(ny):
        ky = dky * y + ky0
        hh = h(kx,ky,p,q,t,swap_xy)
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
    
    
def export_berry_phase(p,q,t,nx=50,ny=50,ky0=0,fn='berry',swap_xy=False,h=hk):
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

            all_p = memo_proj_kx(p,q,t,kx,ny,ky0,swap_xy,h)
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

def export_oneband_quench_coef(q,pi,pf,ti,tf,bidx=0,nx=50,ny=50,kxfrac=1,kyfrac=1, fn='qcoef'):
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
