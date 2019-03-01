import numpy as np
import random
import datetime
import time
import h5py

from wetzel_utils.utility.utility_catalog import indices_tree
from wetzel_utils.utility.utility_array import elements #returns index numbers of elements of an array that fall within given limts
from subhalo_io_hack import TreepmClass
from progressbar import ProgressBar
import my_sham_hack as sham
#from centralms import sham_hack as sham
#from treepm import sham

import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 
#from matplotlib import rcParams

class shamedTreepmClass(TreepmClass):
    def __init__(self,scat=0.,dis_mf=0.,source='li-march',
                 Mwid=0.5,catkind='subhalo',seed=None):
        TreepmClass.__init__(self)
        #super(shamedTreepmClass,self).__init__()
        self.smtype='m.max'
        self.hmtype='m.200c'
        self.gmtype='m.star'
        self.allzis=np.arange(35)
        self.shamzis=self.allzis[:-1]
        self.scat=scat
        self.dis_mf=dis_mf
        self.source=source
        self.Mwid=Mwid
        self.seed=seed
        self.subisread=False
        self.hostisread=False
        self.mtreebuilt=False
        if catkind=='subhalo':
            self.readsub()
        elif catkind=='halo':
            self.readhost()
        elif catkind=='both':
            self.readboth()

    def write_galmMs(self,M0cond,condtype,zibeg,ziend,Mtime,N=None,app=''):
        Mtime=str(Mtime)
        tstmp='{:%Y%m%d}'.format(datetime.datetime.now())
        if (self.seed is None) or (not self.scat):
            seedstr=''
        else:
            seedstr='_'+str(self.seed)+'seed'
        fname='./dat/gal_mM{3}_M{4}0_{0:0.2f}_scat{2:0.2f}{6}_{1}{5}.h5'\
              .format(M0cond,tstmp,self.scat,Mtime,condtype,app,seedstr)
        with h5py.File(fname,'w') as f:
            out=self.gal_mMs(M0cond,condtype,zibeg,ziend,Mtime,N)
            f.create_dataset('mMs',data=out[0])
            f.create_dataset('ms',data=out[1])
            f.create_dataset('Ms',data=out[2])
            f.create_dataset('iprim0s_merg',data=out[3])
            f.create_dataset('iprim0s_keys',data=out[4])
            f.create_dataset('Nprim0',data=out[5])
            f.create_dataset('zbeg',data=out[6])
            f.create_dataset('zend',data=out[7])
            f.create_dataset('zs',data=out[8])

    def merg_rates_z(self,M,zibeg,ziend,N=None):
        if not self.hostisread:
            self.readhost()
        if not self.subisread:
            self.readsub()
        zis=np.arange(zibeg,ziend+1)
        zbeg=self.hostcat.snap[zibeg][1]
        zend=self.hostcat.snap[ziend][1]
        mMs=[]
        Nhost=0
        pbar=ProgressBar()
        for zi in zis[:-1]:
            print'\nsnapshot %i'%zi
            hostis=elements(self.hostcat[zi][self.hmtype],lim=[M-self.Mwid/2.,M+self.Mwid/2.])
            cenis=self.hostcat[zi]['cen.i'][hostis]
            print'%i halos in range.'%len(hostis)
            if not N is None:
                random.seed(1)
                cenis=random.sample(cenis,N)
            print'evaluating %i'%len(cenis)
            Nhost+=len(cenis)

            pbar=ProgressBar()
            for ceni,hosti in pbar(zip(cenis,hostis)):
                M=self.hostcat[zi][self.hmtype][hosti]
                is_prog=self.subcat[zi+1]['chi.i']==ceni
                is_prog[self.subcat[zi]['par.i'][ceni]]=False
                ms=self.subcat[zi+1][self.smtype][is_prog]
                mMs_add=ms-M
                mMs+=list(mMs_add)
        return mMs,Nhost,zbeg,zend

    def merg_rates_z_gal(self,M,zibeg,ziend,N=None):
        #Get list of m_{*,z}/M_{*,z} for galaxies whose M_{host hal,z}
        #match a given mass at z.
        #
        #Things to note:
        #    -This is most useful for calculating dn/dz. The matching
        #     process checks the M condition at EACH redshift, NOT just once
        #     at z=0, because the
        #     assumption is that dz is very small.
        #    -Usually, ziend=zibeg+1.
        #    -z(m)!=z(M). This is a functional choice because the merger happens
        #     between snapshots, so we can't ever know the simulaneous mass of
        #     the accreted galaxy and the main progenitor galaxy. zi(m) will 
        #     always be the zi before the merger because there is no data on the
        #     accreted galaxy post-merger. The script
        #     could have taken zi(M) to also be the zi prior to the merger, in
        #     which case zi(m)==zi(M).
        #     However, I thought it better to keep zi(M)==zi(m)-1 so
        #     that when ziend-zibeg=1 we get what we would expect for M(z);
        #     M(z)=M(zibeg), where zis is monotonically increasing i.e. each
        #     step moves backward in time.
        #    -The mass matching process looks for HOST galaxies whose 
        #     M_{host hal,z} matches M. This, in effect, isolates CENTRAL
        #     subhalos as
        #     opposed to all subhalos. 
        #    -This method subsequently takes the M_z in m/M_z
        #     to be the galaxy mass sham'd to the HOST halo.
        #    -m.max of the central subhalo usually but not always comes close
        #     to/matches m.200c of the host halo.
        #    -Other methods in this script 
        #     match M_z or M_0 directly to the mass of subhalos (m.max usually)
        #     and therefore encompass all mergers, not just those with centrals.
        #
        #Also note:
        #    It is probably possible to combine merg_rates_z_gal and 
        #    merg_rates_z_hal with a simple switch in the function parameters.
        #    I just first
        #    need to make sure everything but the hmtype/gmtype inputs
        #    throughout the two methods are the same.
        if not self.hostisread:
            self.readhost()
        if not self.subisread:
            self.readsub()
        zis=np.arange(zibeg,ziend+1)
        zbeg=self.hostcat.snap[zibeg][1]
        zend=self.hostcat.snap[ziend][1]
        mMs=[]
        Nhost=0
        pbar=ProgressBar()
        for zi in zis[:-1]:
            print'\nsnapshot %i'%zi
            #Get host indices whose masses at z match M:
            hostis=elements(self.hostcat[zi][self.hmtype],
                            lim=[M-self.Mwid/2.,M+self.Mwid/2.])
            #Get the corresponding central subhalo indices from the subhalo 
            #cat:
            cenis=self.hostcat[zi]['cen.i'][hostis]
            print'%i halos in range.'%len(hostis)
            if not N is None:
                random.seed(1)
                cenis=random.sample(cenis,N)
            print'evaluating %i'%len(cenis)
            Nhost+=len(cenis)

            pbar=ProgressBar()
            for ceni,hosti in pbar(zip(cenis,hostis)):
                #Take M_z to be the galaxy mass sham'd to the HOST halo at z:
                M=self.hostcat[zi][self.gmtype][hosti]
                is_prog=self.subcat[zi+1]['chi.i']==ceni
                is_prog[self.subcat[zi]['par.i'][ceni]]=False
                #Take m_z to be the mass of the accreted galaxy at zi+1: 
                ms=self.subcat[zi+1][self.gmtype][is_prog]
                mMs_add=ms-M
                mMs+=list(mMs_add)
        return mMs,Nhost,zbeg,zend

    def merg_tree(self,zibeg,ziend):
        #Get list of m_*/M_* 
        #M0cond can be based on either galaxy or halo mass at z=0.
        #t0=time.time()
        if not self.subisread:
            self.readsub()
        zis=np.arange(zibeg,ziend+1)
        #print zis[:0:-1]
        zbeg=self.subcat.snap[zibeg][1]
        zend=self.subcat.snap[ziend][1]
        #i_primaries0=self.get_primaries0(condtype,M0cond)
        #i0s=np.arange(len(self.subcat[zis[-1]][self.gmtype]))

        '''
        print'%i galaxies'%len(i0s)
        if not N is None:
            random.seed(1)
            i0s=random.sample(i0s,N)
        N0=len(i0s)
        print'evaluating %i'%N0
        ''' 
        
        #pbar=ProgressBar()
        #Evaluate zis in reverse order from highest redshift through second to
        #lowest redshift:
        #the second lowest:
        for zi in zis[:0:-1]:
        #for zi in pbar(zis[:0:-1]):
            self.subcat[zi-1]['par.tree']={}
            his=np.arange(len(self.subcat[zi][self.gmtype]))
            #i_primaries=indices_tree(self.subcat,zis[0],zi,i0s)

            #isolate primaries that exist and have mass at zi:
            #i_primaries=i_primaries[i_primaries>=0.]
            chiis=indices_tree(self.subcat,zi,zi-1,his)
            his=his[chiis>=0]
            chiis=chiis[chiis>=0]
            m_chis=self.subcat[zi-1][self.gmtype][chiis]
            ms=self.subcat[zi][self.smtype][his]
            his=his[(ms>0.) & (m_chis>0.)]
            
            print'{0:d} galaxies in snapshot {1:d}'.format(len(his),zi)

            for hi,chii in zip(his,chiis):
                #print hi
                '''
                This shouldn't be necessary given that I'm removing primaries
                where hi<0, so I'm commenting it. I'll delete it after testing.
                #if the parent doesn't exist at zi, go to the next parent:
                if hi<0: 
                    continue
                '''
                #chii=self.subcat[zi]['chi.i'][hi]
                #prim_chii=indices_tree(self.subcat,zi,zi-1,hi)
                #isolate the accreted galaxies:
                #notprim=chii!=prim_chii
                #chii=chii[notprim]
                
                '''
                if hi==6877:
                    print zi
                    print hi
                    print chii
                    print self.subcat[zi-1]['par.i'][chii]
                    print hi==self.subcat[zi-1]['par.i'][chii]
                '''
                if hi==self.subcat[zi-1]['par.i'][chii]: 
                    continue
                #print'it worked'
                #if len(chii)==0:
                #    continue
                if chii not in self.subcat[zi-1]['par.tree']:
                    self.subcat[zi-1]['par.tree'][chii]=list([hi])
                else:
                    branch=self.subcat[zi-1]['par.tree'][chii]
                    branch.extend(hi)
        self.mtreebuilt=True

    def gal_mMs_fromtree(self,M0cond,condtype,zibeg,ziend,Mtime,N=None):
        if not self.mtreebuilt:
            self.merg_tree(zibeg,ziend)

    def gal_mMs(self,M0cond,condtype,zibeg,ziend,Mtime,N=None):
        #Get list of m_*/M_* 
        #M0cond can be based on either galaxy or halo mass at z=0.
        #t0=time.time()
        if not self.subisread:
            self.readsub()
        zis=np.arange(zibeg,ziend+1)
        zbeg=self.subcat.snap[zibeg][1]
        zend=self.subcat.snap[ziend][1]
        mMs=[]
        ms=[]
        Ms=[]
        zs=[]
        iprim0s_merg=[]
        
        i_primaries0=self.get_primaries0(condtype,M0cond)
        
        print'%i galaxies in range'%len(i_primaries0)
        if not N is None:
            random.seed(1)
            i_primaries0=random.sample(i_primaries0,N)
        Nprim0=len(i_primaries0)
        print'evaluating %i'%Nprim0

        pbar=ProgressBar()
        for i_primary0 in pbar(i_primaries0):
            out=self.mMs_add_gal_f(i_primary0,zis,Mtime)
            #print out
            mMs+=out[0]
            ms+=out[1]
            Ms+=out[2]
            zs+=out[3]
            iprim0s_merg+=list(np.repeat(i_primary0,len(out[0])))
        #t1=time.time()
        #print'run time: {0}'.format(t1-t0)
        return mMs,ms,Ms,iprim0s_merg,i_primaries0,Nprim0,zbeg,zend,zs

    '''
    def quench_dat(mus,mMs,iprim0_merg,iprim0_keys):
        lmuwid=0.01 #log width of the merger ratio
        lmus=np.log10(np.array(mus))
        iscen=self.subcat[0][ilk]==1
        M0s=self.subcat[0][gmtype][iscen]
        minM0=np.min(M0s)
        print'min: %f'%minM0
        maxM0=np.max(M0s)
        print'max: %f'%maxM0
        Mbins=np.arange(minM0,maxM0,)
        for mu in mus:
            Ns=[]
            for M in Mbins:
    '''

    def get_primaries0(self,condtype,M0cond):
        if condtype=='gal':
            i_primaries0=self.get_is(M0cond,condtype)
        elif condtype=='subhal':
            i_primaries0=self.get_is(M0cond,condtype)
        elif condtype in ['censubhal','cengal']:
            i_primaries0=self.get_cenis(M0cond,condtype)
        else:
            raise ValueError('Condition type is "cengal", "gal", "subhal", or '
                             '"censubhal".')
        return i_primaries0

    def get_is(self,M,condtype,zi=0):
        if condtype=='subhal':
            mtype=self.smtype
        elif condtype=='gal':
            mtype=self.gmtype
        else:
            raise ValueError('Unexpected condition type error')

        il=elements(self.subcat[zi][mtype],
                       lim=[M-self.Mwid/2.,M+self.Mwid/2.])
        return il

    def get_cenis(self,M,condtype,zi=0):
        if condtype=='censubhal':
            mtype=self.smtype
        elif condtype=='cengal':
            mtype=self.gmtype
        else:
            raise ValueError('Unexpected condition type error')
        
        iscen=self.subcat[zi]['ilk']==1
        ms=self.subcat[zi][mtype]
        inrange=(ms>M-self.Mwid/2.)&(ms<M+self.Mwid/2.)
        Nsubs=len(self.subcat[zi][mtype])
        cenis=np.arange(Nsubs)[iscen & inrange]
        return cenis

    def mMs_add_gal_f(self,i_primary0,zis,Mtime):
        #Get a list of m_gal/M_gal ratios spanning zis for a single z=0 subhalo
        #M_gal can either be M0 or Mz depending on Mtime.
        if not Mtime in ['0','z']:
            raise ValueError('Mtime is either "0" or "z".')
        mMs=[]
        ms=[]
        Ms=[]
        zs=[]
        zi0=zis[0]
        if Mtime=='0':
            M=self.subcat[zi0][self.gmtype][i_primary0]
        for zi in zis[:-1]:
            #index of the main progenitor at zi:
            i_primary=indices_tree(self.subcat,zi0,zi,i_primary0)
            if i_primary<0:
                return mMs,ms,Ms,zs
            #index of the main progenitor at zi+1:
            pari=self.subcat[zi]['par.i'][i_primary] 
            
            Mz=self.subcat[zi+1][self.gmtype][pari]
            #For some reason, it sometimes happens that Mz==0 while
            #secondary progenitors exist at z that also have 0 mass. This would
            #lead to m==M, in which case m-M and m/M equal 0. (Remember m and M 
            #are actually log masses.) Therefore, as soon as Mz==0, the method
            #writes no more data and returns what it has. (This didn't show up
            #before we were doing m/Mz, because with m/M0, M0 is never 0, and
            #0-M0 shows up so far on the left of the plot that we never noticed
            #it.)
            if Mz==0:
                return mMs,ms,Ms,zs
            
            if Mtime=='z':
                M=Mz #When we're doing m/Mz, set the denominator to Mz.
            is_prog=self.subcat[zi+1]['chi.i']==i_primary
            #Remove the main progenitor from ms in the m/M calculation:
            is_prog[pari]=False 
            ms_add=self.subcat[zi+1][self.gmtype][is_prog]
            mMs_add=ms_add-M
            Ms_add=np.repeat(M,len(ms_add))
            z=self.subcat[zi+1].snap['z']
            zs_add=np.repeat(z,len(ms_add))
            ms+=list(ms_add)
            Ms+=list(Ms_add)
            mMs+=list(mMs_add)
            zs+=list(zs_add)
            #print zs
        #mMs=list(-np.abs(np.array(mMs)))
        return mMs,ms,Ms,zs

    def readhost(self):
        self.hostcat=self.read(zis=self.allzis,catalog_kind='halo')
        sham.assign(self.hostcat,scat=self.scat,dis_mf=self.dis_mf,
                    source=self.source,
                    sham_prop=self.hmtype,zis=self.shamzis,seed=self.seed)
        self.hostisread=True

    def readsub(self):
        self.subcat=self.read(zis=self.allzis,catalog_kind='subhalo')
        sham.assign(self.subcat,scat=self.scat,dis_mf=self.dis_mf,
                    source=self.source,
                    #sham_prop=self.smtype,zis=self.shamzis)
                    sham_prop=self.smtype,zis=self.shamzis,seed=self.seed)
        self.subisread=True

    def readboth(self):
        self.subcat,self.hostcat=self.read(zis=self.allzis,
                                           catalog_kind='both')
        sham.assign(self.subcat,scat=self.scat,dis_mf=self.dis_mf,
                    source=self.source,
                    sham_prop=self.smtype,zis=self.shamzis,seed=self.seed)
        sham.assign(self.hostcat,scat=self.scat,dis_mf=self.dis_mf,
                    source=self.source,
                    sham_prop=self.hmtype,zis=self.shamzis,seed=self.seed)
        self.subisread=True
        self.hostisread=True

    def dndz(self,mMs,Nhost,dz):
        bins=np.linspace(-5,0,100)
        fig=plt.figure(figsize=(12,10))
        ax=fig.add_subplot(111)
        n,bins,patches=ax.hist(mMs,bins=bins,
                               weights=np.repeat(1./float(Nhost)/dz,
                                                 len(mMs)),
                               cumulative=-1,ec='k')
        ax.set_yscale('log')
        midbins=(bins[:-1]+bins[1:])/2.
        plt.clf()
        return n,midbins

    def count(self,mMs,Nhost):
        bins=np.linspace(-5,0,100)
        fig=plt.figure(figsize=(12,10))
        ax=fig.add_subplot(111)
        n,bins,patches=ax.hist(mMs,bins=bins,
                               weights=np.repeat(1./float(Nhost),
                                                 len(mMs)),
                               cumulative=-1,ec='k')
        ax.set_yscale('log')
        midbins=(bins[:-1]+bins[1:])/2.
        plt.clf()
        return n,midbins

###############################################################################

def read_galmMs(Mcond,condtype,Mtime,scat,tstmp,seed=None,apnd=''):
    Mtime=str(Mtime)
    if (seed is None) or (not scat):
        seedstr=''
    else:
        seedstr='_'+str(seed)+'seed'
    fname='/home/users/staudt/projects/mergers/dat/' \
          'gal_mM{4}_M{3}0_{0:0.2f}_scat{2:0.2f}{5}_{1}{6}.h5' \
          .format(Mcond,tstmp,scat,condtype,Mtime,seedstr,apnd)
    print fname
    with h5py.File(fname,'r') as f:
        mMs=np.array(f['mMs'])
        ms=np.array(f['ms'])
        Ms=np.array(f['Ms'])
        iprim0s_merg=np.array(f['iprim0s_merg'])
        iprim0s_keys=np.array(f['iprim0s_keys'])
        Nprim0=(f['Nprim0']).value
        zbeg=(f['zbeg']).value
        zend=(f['zend']).value
        #Using this for now because of a typo in the mMswrite code. Change
        #later:
        zs=np.array(f['za'])
        #zs=np.array(f['zs'])
    return mMs,ms,Ms,iprim0s_merg,iprim0s_keys,Nprim0,zbeg,zend,zs

def write_gal_hgram_dat(Mconds,condtype,Mtime,rngs,zmax,
                        scat,mMs_tstmp,seed=None,
                        readapnd='',apnd='',tstmp0=None):
    Mtime=str(Mtime)
    if (seed is None) or (not scat):
        seedstr=''
    else:
        seedstr='_'+str(seed)+'seed'
    if tstmp0 is None:
        tstmt0=mMs_tstmp
    tstmp='{:%Y%m%d}'.format(datetime.datetime.now())
    fname='./dat/gal_hgram_dat{2}_{0}{1}.h5'.format(tstmp,apnd,seedstr)
    with h5py.File(fname,'w') as f: 
        dat=[]
        Nprim0s=[]
        for Mcond in Mconds:
            mMsdat_0scat=read_galmMs(Mcond,condtype,Mtime,0,tstmp0,seed)
            mMsdat_scat=read_galmMs(Mcond,condtype,Mtime,scat,mMs_tstmp,seed,
                                    readapnd)
            Mcond_dat=[]
            pbar=ProgressBar()
            for rng in pbar(rngs):
                Ns_0scat=gal_hgram_dat_noz(rng,mMsdat_0scat[0],
                                       mMsdat_0scat[3],
                                       mMsdat_0scat[4])
                Ns_scat=gal_hgram_dat(rng,mMsdat_scat[0],
                                      mMsdat_scat[3],
                                      mMsdat_scat[4],
                                      zmax,mMsdat_scat[8])
                dname0='Ns_{0:0.1f}M_{1:0.1f}-{2:0.1f}_0scat'\
                                 .format(Mcond,rng[0],rng[1])
                dname='Ns_{0:0.1f}M_{1:0.1f}-{2:0.1f}_scat'\
                                 .format(Mcond,rng[0],rng[1])
                f.create_dataset(dname0,data=Ns_0scat)
                f.create_dataset(dname,data=Ns_scat)
            #number of primary ending galaxies for the 0scat and scat data:
            Nprim0s.append([mMsdat_0scat[5],mMsdat_scat[5]])
        f.create_dataset('Nprim0s',data=Nprim0s)
        f.create_dataset('Mconds',data=Mconds)
        f.create_dataset('rngs',data=rngs)

def write_quench_dat(condtype,Mtime,mus,scat,mMs_tstmp):
    #unfinished function. Might delete later.
    lmuwid=0.01 #log width of the merger ratio
    lmus=np.log10(np.array(mus))

def gal_hgram_dat(rng,mMs,iprim0s_merg,iprim0s_keys,zmax,zs):
    #Generate a list of the number of times a merger happens within the 
    #given range. Each datapoint represents the merger count for a single
    #main progenitor.
    mMlo=rng[0]
    mMhi=rng[1]
    inrange=(mMs>mMlo)&(mMs<mMhi)&(zs<=zmax)
    Ns=np.array([np.sum(ismp & inrange)
                 for ismp in
                 (iprim0s_merg==key for key in iprim0s_keys)])
    return Ns

def gal_hgram_dat_noz(rng,mMs,iprim0s_merg,iprim0s_keys):
    #Generate a list of the number of times a merger happens within the 
    #given range. Each datapoint represents the merger count for a single
    #main progenitor.
    mMlo=rng[0]
    mMhi=rng[1]
    inrange=(mMs>mMlo)&(mMs<mMhi)
    Ns=np.array([np.sum(ismp & inrange)
                 for ismp in
                 (iprim0s_merg==key for key in iprim0s_keys)])
    return Ns
   
