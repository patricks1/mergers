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

def get_hostis(M,Mwid,zi,mtype='m.fof'):
    hostis=elements(hostcat[zi][mtype],lim=[M-Mwid/2.,M+Mwid/2.])
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
        
        #The following commented block of code uses the getinfall method to count when subhalos entered host halos. It may come in handy later, so I'm not deleting it.
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

def mM0s_add_f_new(hosti0):
    mtype='m.fof'
    mM0s=[]
    M0=hostcat[0][mtype][hosti0]
    for zi in allzis[:-1]: 
        i_primary=indices_tree(hostcat,0,zi,hosti0)
        if i_primary<1:
            return mM0s
        is_prog=hostcat[zi+1]['chi.i']==i_primary
        is_prog[hostcat[zi]['par.i'][i_primary]]=False #is_prog is a bool array. This line sets the elements corresponding to the main progentior to false, so we're not incorrectly counting a merger of the main progenitor with itself.
        ms=hostcat[zi+1][mtype][is_prog]
        mM0s_add=ms-M0
        mM0s+=list(mM0s_add)
    return mM0s

def accr_ratesz_partner_meth(M,zibeg,ziend,N=None,mtype='m.200c'): #This uses the par.n.i method. 
    '''
    timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
    filename='./dat/mMaccr_{1:d}_{0}.h5'.format(timestmp,M)
    f = h5py.File(filename, 'w')
    '''

    zis=np.arange(zibeg,ziend+1)
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
        hostis=get_hostis(M,Mwid,zi-1,mtype=mtype) #get child indices whose masses are in range AFTER the merger
        hostis=hostcat[zi-1]['par.i'][hostis] #get the host indices for the current snapshot
        print'%i halos in range.'%len(hostis)
        if not N is None:
            #random.seed(1)
            hostis=random.sample(hostis,N)
        print'evaluating %i'%len(hostis) 
        #print'host indices:'
        #print hostis
        Nhost+=len(hostis)
        #print'Nhost=%i'%Nhost
        pbar=ProgressBar()
        hostis_used=[]
        #for hosti in hostis:
        for hosti in pbar(hostis):
            #subis_counted[zi]=[]
            if hosti in hostis_used:
                #print'duplicate'
                Nhost-=1
                continue
            chii=hostcat[zi]['chi.i'][hosti]
            M=hostcat[zi-1][mtype][chii]
            nhosti=hostcat[zi]['par.n.i'][hosti]
            while nhosti>-1:
                #print'next host idx: %i'%nhosti
                m=hostcat[zi][mtype][nhosti]
                mM_add=m-M
                #print'mM: %f'%mM_add
                #print mM0s[-1]
                mMs+=list(mM_add.flatten())
                #print mM0s[-2:]
                nhosti=hostcat[zi]['par.n.i'][nhosti]
                #print nhosti
                #print type(nhosti)
                hostis_used+=[nhosti]
        #print'Nhost=%i'%Nhost
    zbeg=hostcat.snap[zis[0]][1]
    zend=hostcat.snap[zis[-1]][1]
    #print'ending Nhost=%i'%Nhost
    return mMs,Nhost,zbeg,zend
    #return mMs,Nhost

def accr_ratesz_parent_meth(M,zibeg,ziend,N=None,mtype='m.200c'):
    zis=np.arange(zibeg,ziend+1)
    zbeg=hostcat.snap[zibeg][1]
    zend=hostcat.snap[ziend][1]
    #print'zis:'
    #print zis
    mMs=[]
    Nhost=0
    pbar=ProgressBar()
    for zi in zis[:-1]: 
        print'\nsnapshot %i'%zi
        i_primaries=elements(hostcat[zi][mtype],lim=[M-Mwid/2.,M+Mwid/2.])
        print'%i halos in range.'%len(i_primaries)
        if not N is None:
            random.seed(1)
            i_primaries=random.sample(i_primaries,N)
        print'evaluating %i'%len(i_primaries)
        Nhost+=len(i_primaries)
        
        pbar=ProgressBar()
        for i_primary in pbar(i_primaries):
            M=hostcat[zi][mtype][i_primary]
            is_prog=hostcat[zi+1]['chi.i']==i_primary
            is_prog[hostcat[zi]['par.i'][i_primary]]=False
            ms=hostcat[zi+1][mtype][is_prog]
            mMs_add=ms-M
            mMs+=list(mMs_add)
    return mMs,Nhost,zbeg,zend

'''
def accr_ratesz_subs(M,zibeg,ziend,N=None):
    mtype='m.max'
    zis=np.arange(zibeg,ziend+1)
    zbeg=hostcat.snap[zibeg][1]
    zend=hostcat.snap[ziend][1]
    #print'zis:'
    #print zis
    mMs=[]
    Nhost=0
    pbar=ProgressBar()
    for zi in zis[:-1]: 
        print'\nsnapshot %i'%zi
        i_primaries=elements(hostcat[zi][mtype],lim=[M-Mwid/2.,M+Mwid/2.])
        print'%i halos in range.'%len(i_primaries)
        if not N is None:
            random.seed(1)
            i_primaries=random.sample(i_primaries,N)
        print'evaluating %i'%len(i_primaries)
        Nhost+=len(i_primaries)
        
        pbar=ProgressBar()
        for i_primary in pbar(i_primaries):
            M=hostcat[zi][mtype][i_primary]
            is_prog=hostcat[zi+1]['chi.i']==i_primary
            is_prog[hostcat[zi]['par.i'][i_primary]]=False
            ms=hostcat[zi+1][mtype][is_prog]
            mMs_add=ms-M
            mMs+=list(mMs_add)
    return mMs,Nhost,zbeg,zend
'''

def accr_rates0(M0,N=None,excl=False):  
    timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
    filename='./dat/mM0accr_{1:d}_{0}.h5'.format(timestmp,M0)
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
        mM0s+=mM0s_add_f_new(hosti0)
   
    f.create_dataset('mM0s', data=mM0s)
    f.create_dataset('Nhost0',data=Nhost0)
    f.close()

    return

def compare(M,zibeg,ziend,N=None):
    mtype='m.fof'
    zis=np.arange(zibeg,ziend+1)
    hostis=elements(hostcat[zibeg][mtype],lim=[M-Mwid/2.,M+Mwid/2.])
    if not N is None:
        hostis=random.sample(hostis,N)
    for zi in zis:
        print'snapshot: %i'%zi 
        pbar=ProgressBar()
        for hosti in pbar(hostis):
            print'host index: %i'%hosti
            isprog_parent=hostcat[zi+1]['chi.i']==hosti
            #print type(isprog_parent)
            #print isprog_parent
            pari=hostcat[zi]['par.i'][hosti]
            isprog_parent[pari]=False

            
            parni=hostcat[zi+1]['par.n.i'][pari]
            parnis=[]
            while parni>-1:
                parnis+=[parni]
                parni=hostcat[zi+1]['par.n.i'][parni]
                #hostis_used+=[nhosti]
            isprog_partner=np.repeat(False,
                                     len(hostcat[zi+1][mtype]))
            isprog_partner[parnis]=True
            if not np.array_equal(isprog_partner,isprog_parent):
                print('problem with host index {0:d} parents at snapshot {1:d}'.format(hosti,
                                                                                       zi))
    return isprog_parent,isprog_partner 
