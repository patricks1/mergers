import numpy as np
import random
import datetime
import h5py

'''
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 

from matplotlib import rcParams
'''

from wetzel_utils.utility.utility_catalog import indices_tree
from wetzel_utils.utility.utility_array import elements #returns index numbers of elements of an array that fall within given limts
from subhalo_io_hack import TreepmClass

from progressbar import ProgressBar

allzis=np.arange(35)
treepm=TreepmClass()
#subcat,hostcat=treepm.read(zis=allzis,catalog_kind='both')
hostcat=treepm.read(zis=allzis,catalog_kind='halo')

msubtype='m.bound'
Mwid=0.5

def get_hostis(M,Mwid,zi):
    hostis=elements(hostcat[zi]['m.fof'],lim=[M-Mwid/2.,M+Mwid/2.])
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
    for zi in allzis:
        #print'\n\nsnapshot %i'%zi
        hosti=indices_tree(hostcat,0,zi,hosti0)
        #print'\nhost idx: %i'%hosti 
        if hosti<1:
            return mM0s
        '''
        subis=getsubs(hosti,zi)
        isnew=np.array([not subi in subis_prev for subi in subis],
                       dtype=bool)
        subis_inf,zis_inf=getinfall(subis,zi)
        infhaps=subis_inf>0 #If the subhalo is part of the host for the whole simulation, the indices will be -1.
        #subis_inf,zis_inf=subis_inf[infhaps],zis_inf[infhaps]
        subis_b4inf,zis_b4inf=subis_inf[infhaps],zis_inf[infhaps]
        zis_inf=zis_b4inf-1
        subis_inf=[indices_tree(subcat,zi_b4inf,zi_inf,subi_b4inf) 
                   for zi_b4inf,zi_inf,subi_b4inf 
                   in zip(zis_b4inf,zis_inf,subis_b4inf)]
        ms=[subcat[zi_inf][msubtype][subi_inf] 
            for zi_inf,subi_inf 
            in zip(zis_inf,subis_inf)] 
        mM0s_add=ms-M0
        mM0s+=list(mM0s_add)
        subis_prev=subis
        '''

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

def accr_ratesz(M,zibeg,ziend,N=None):
    '''
    timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
    filename='./dat/mMaccr_{1:d}_{0}.h5'.format(timestmp,M)
    f = h5py.File(filename, 'w')
    '''

    zis=np.arange(zibeg,ziend+1)
    print'zis:'
    print zis
    mMs=[]
    #I might not need to do the following---make it so the program searches for host halos that match M at all zis being evaluated. Put those hostis into a dictionary corresponding to their zis.
    #hostis_mat={}
    #hostis_mat[0]=get_hostis(M,Mwid,zis[0])
    ''' 
    hostis=get_hostis(M,Mwid,zis[0]) 
    if not N is None:
        random.seed(1)
        hostis=random.sample(hostis,N)
    '''
    #subis_counted={} #making this is a dictionary so I can have a subis_counted[-1] check when evaluating zi=0
    #subis_counted[-1]=[]
    Nhost=0
    for zi in zis[1:]: #The [1:] effectively makes it so we're correctly finding mergers that happen BETWEEN the zis specified.
        print'snapshot %i'%zi
        hostis=get_hostis(M,Mwid,zi)
        print'%i halos in range.'%len(hostis)
        if not N is None:
            random.seed(1)
            hostis=random.sample(hostis,N)
        print'evaluating %i'%len(hostis) 
        #print'host indices:'
        #print hostis
        Nhost+=len(hostis)
        pbar=ProgressBar()
        for hosti in hostis:
        #for hosti in pbar(hostis):
            #subis_counted[zi]=[]
            M=hostcat[zi]['m.fof'][hosti]
            '''
            #print'prev counted:'
            #print subis_counted[zi-1]
            subis=getsubs(hosti,zi)
            #uncounted=np.array([not subi in subis_counted[zi-1] for subi in subis],
            #                   dtype=bool)
            #print'not yet counted:'
            #print subis[uncounted]
            subis_b4inf,zis_b4inf=getinfall(subis,zi)
            infhaps=zis_b4inf>=0
            subis=subis_b4inf[infhaps]
            zis_b4inf=zis_b4inf[infhaps]
            zis_inf=zis_b4inf-1
            subis_inf=[indices_tree(subcat,zi_b4inf,zi_inf,subi_b4inf) 
                       for zi_b4inf,zi_inf,subi_b4inf 
                       in zip(zis_b4inf,zis_inf,subis_b4inf)]
            infnow=zis_inf==zi
            
            ##DEBUGGING
            return subis_inf,subis,infnow
            ##END DEBUGGING
            
            #print infnow
            print subis
            print subis_inf
            print zis_inf
            print'test whether subis_inf=subis[infnow]:'
            print subis_inf[infnow]==subis[infnow]
            print ''
            ms=subcat[zi][msubtype][subis[infnow]]
            mMs_add=ms-M
            mMs+=list(mMs_add)
            #subis_counted[zi]+=list(subis[infhap])
            '''
            nhosti=hostcat[zi]['par.n.i'][hosti]
            while nhosti>0:
                #print'next host idx: %i'%nhosti
                m=hostcat[zi]['m.fof'][nhosti]
                mM_add=m-M
                #print'mM: %f'%mM_add
                #print mM0s[-1]
                mMs+=list(mM_add.flatten())
                #print mM0s[-2:]
                nhosti=hostcat[zi]['par.n.i'][nhosti]
    zbeg=hostcat.snap[zis[0]][1]
    zend=hostcat.snap[zis[-1]][1]
    return mMs,Nhost,zbeg,zend
    #return mMs,Nhost

def accr_ratesz_new(M,zibeg,ziend,N=None):
    mtype='m.fof'
    zis=np.arange(zibeg,ziend+1)
    zbeg=hostcat.snap[zibeg][1]
    zend=hostcat.snap[ziend][1]
    print'zis:'
    print zis
    mMs=[]
    Nhost=0
    pbar=ProgressBar()
    for zi in zis[:-1]: 
        print'snapshot %i'%zi
        i_primaries=elements(hostcat[zi][mtype],lim=[M-Mwid/2.,M+Mwid/2.])
        print'%i halos in range.'%len(i_primaries)
        if not N is None:
            random.seed(1)
            i_primaries=random.sample(i_primaries,N)
        print'evaluating %i'%len(i_primaries)
        Nhost+=len(i_primaries)
        
        for i_primary in i_primaries:
            M=hostcat[zi][mtype][i_primary]
            is_prog=hostcat[zi+1]['chi.i']==i_primary
            is_prog[hostcat[zi]['par.i'][i_primary]]=False
            ms=hostcat[zi+1][mtype][is_prog]
            mMs_add=ms-M
            mMs+=list(mMs_add)
    return mMs,Nhost,zbeg,zend


def accr_rates0(M0,N=None,excl=False):  
    timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
    filename='./dat/mM0accr_{1:d}_{0}_test.h5'.format(timestmp,M0)
    f = h5py.File(filename, 'w')

    mM0s=[]
    hostis0=get_hostis(M0,Mwid,zi=0)
    #print hostis0
    print'%i host halos in mass range.'%len(hostis0)
    if not N is None:
        random.seed(1)
        hostis0=random.sample(hostis0,N)
    print'evaluating %i'%len(hostis0)
    #print hostis0
    Nhost0=len(hostis0)
    pbar=ProgressBar()
    for hosti0 in pbar(hostis0):
        #print'\n\n\nhost0 idx: %i'%hosti0
        mM0s+=mM0s_add_f(hosti0,excl=excl)
   
    f.create_dataset('mM0s', data=mM0s)
    f.create_dataset('Nhost0',data=Nhost0)
    f.close()

    return
