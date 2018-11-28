import numpy as np
import random
import datetime
import h5py

import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 

from matplotlib import rcParams

from wetzel_utils.utility.utility_catalog import indices_tree
from wetzel_utils.utility.utility_array import elements #returns index numbers of elements of an array that fall within given limts
from subhalo_io_hack import TreepmClass

from progressbar import ProgressBar

zis=np.arange(35)
treepm=TreepmClass()
subcat,hostcat=treepm.read(zis=zis,catalog_kind='both')

M0wid=0.5

def get_hostis0(M0,M0wid):
    hostis=elements(hostcat[0]['m.fof'],lim=[M0-M0wid/2.,M0+M0wid/2.])
    return hostis

def getsubs(hosti,zi):
    subis=np.where(subcat[zi]['halo.i']==hosti)[0]
    return subis

def getinfall(subi,zi):
    zi_inf=subcat[zi]['inf.last.zi'][subi]
    subi_inf=subcat[zi]['inf.last.i'][subi]
    return subi_inf,zi_inf

def mM0s_add_f(hosti0,excl=False):
    mM0s=[]
    M0=hostcat[0]['m.fof'][hosti0]
    subis_prev=[]
    for zi in zis:
        #print'\n\nsnapshot %i'%zi
        hosti=indices_tree(hostcat,0,zi,hosti0)
        #print'\nhost idx: %i'%hosti 
        if hosti<1:
            return mM0s
        subis=getsubs(hosti,zi)
        isnew=np.array([not subi in subis_prev for subi in subis],
                       dtype=bool)
        subis_inf,zis_inf=getinfall(subis,zi)
        infhaps=subis_inf>0 #If the subhalo is part of the host for the whole simulation, the indices will be -1.
        subis_inf,zis_inf=subis_inf[infhaps],zis_inf[infhaps]
        ms=[subcat[zi_inf][msubtype][subi_inf] 
            for zi_inf,subi_inf 
            in zip(zis_inf,subis_inf)] 
        mM0s_add=ms-M0
        mM0s+=list(mM0s_add)
        subis_prev=subis
        
        if not excl:
            nhosti=hostcat[zi]['par.n.i'][hosti]
            while nhosti>0:
                #print'next host idx: %i'%nhosti
                m=hostcat[zi]['m.fof'][nhosti]
                mM0_add=m-M0
                #print mM0s[-1]
                mM0s+=list(mM0_add.flatten())
                #print mM0s[-2:]
                nhosti=hostcat[zi]['par.n.i'][nhosti]
    return mM0s

def accr_rates(M0,N=None,excl=False):  
    mM0s=[]
    hostis0=get_hostis0(M0,M0wid)
    #print hostis0
    print'%i host halos in mass range.'%len(hostis0)
    if not N is None:
        random.seed(1)
        hostis0=random.sample(hostis0,N)
    print'evaluating %i'%len(hostis0)
    #print hostis0
    Nhostis0=len(hostis0)
    pbar=ProgressBar()
    for hosti0 in pbar(hostis0):
        #print'\n\n\nhost0 idx: %i'%hosti0
        mM0s+=mM0s_add_f(hosti0,excl)
   
    timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
    #filename='./dat/true_m_M0_{0:0.0f}_{2}_{1}.h5'.format(mmid,timestmp,self.mkind)
    filename='./dat/mM0accr_{1}_{0}.h5'.format(timestmp,M0)
    f = h5py.File(filename, 'w')
    f.create_dataset('mM0s', data=mM0s)
    f.create_dataset('Nhostis0',data=Nhostis0)
    f.close()

    return
