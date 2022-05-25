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

def hk(k,t,m,d):
    res = np.zeros((3,3),dtype=np.double)
    ckt = np.cos(k)+t
    sk = np.sin(k)
    res[0,1]=res[1,0] = ckt + d
    res[1,2]=res[2,1] = ckt - d

    res[0,0] = sk + m
    res[2,2] = sk - m
    res[1,1] = -sk
    return res

def hk2(k,t,m=0,d=0):
    # m and d are dummy
    res = np.zeros((3,3), dtype=np.complex)
    eik = np.exp(1j*k)
    a = 1 + t*eik
    b = t + 1*eik
    res[0,1] = res[2,1] = a;
    res[1,0] = res[1,2] = np.conjugate(a);
    #res[0,1] = res[1,2] = a
    #res[1,0] = res[2,1] = np.conjugate(a)
    return res

def export_hk(t=0.8, m=0.1, d=0, nk=100, fn='ssh3-erg'):
    with open(fn + '.dat', 'w') as dat:
        labels = 'k erg1 erg2 erg3'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])

        print(header, file=dat)
        for k in np.linspace(0,2*np.pi,nk):
            hh = hk(k,t,m,d)
            eig=la.eigvalsh(hh)
            print('%g\t%s'%(k, '\t'.join([ '%g'%erg for erg in eig ])), file=dat)
            dat.write('\n')
    print(np.angle(berry(t,m,d,nk)))
def berry(t=0.8,m=0.1,d=0,nk=100):
    id3 = np.identity(3, dtype=np.complex)
    wilson = [id3,id3,id3]
    for k in np.linspace(0,2*np.pi, nk, endpoint=False):
        hh = hk(k,t,m,d)
        eig,u = la.eigh(hh)
        for band in range(3):
            wilson[band] =  np.dot(npext.v2proj(u[:,band]), wilson[band])
    res = [ np.trace(w) for w in wilson ]
    return res
        
def export_overlap(iband = 0, ti=0.8, mi=0.5, di=0, tf=100, mf=0,df=0, nk=100, fn='ssh3-overlap'):
    with open(fn + '.par', 'w') as par:
        txt = """
        iband = %d
        ti = %g
        mi = %g
        di = %g
        tf = %g
        mf = %g
        df = %g
        nk = %d
        """%(iband, ti,mi,di,tf,mf,df,nk)
        print(txt, file=par)

    with open(fn + '.dat', 'w') as dat:
        labels = 'k <phi_1|psi>_square <phi_2|psi>_square <phi_3|psi>_square'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])

        print(header, file=dat)
        for k in np.linspace(0, 2*np.pi,nk):
            hi = hk(k,ti,mi,di)
            eig, ui = la.eigh(hi)

            hf = hk(k,tf,mf,df)
            eig, uf = la.eigh(hf)

            overlap = np.dot(npext.dagger(uf), ui[:,iband])
            overlap2 = overlap * np.conjugate(overlap)
            print('%g\t%s'%(k, '\t'.join([ '%g'%oo for oo in overlap2 ])), file=dat)

