import numpy as np
import random
import datetime
import h5py
import pandas as pd
import scipy
import operator

#from wetzel_utils.utility.utility_catalog import indices_tree
from wtreepm3.utility.utility_catalog import indices_tree
#returns index numbers of elements of an array that fall within given limts:
#from wetzel_utils.utility.utility_array import elements 
from wtreepm3.utility.utility_array import elements 
#from subhalo_io_hack import TreepmClass
from wtreepm3.subhalo_io import TreepmClass
from progressbar import ProgressBar
import my_sham_hack3 as sham
#from centralms import sham_hack as sham
#from treepm import sham
import staudt_utils3 as sutils

import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 
#from matplotlib import rcParams

class GalTreepmClass(object):
    def __init__(self,haltpm,scat=0.,dis_mf=0.,source='rand',Mwid=0.5,
                 haltyp='subhalo',shamzibeg=0,shamziend=34,seed=None,mmin=5.,
                 conscat=True,sham_prop='m.max'):
        if haltyp=='subhalo':
            haltpm.readsub()
            halcat=haltpm.subcat
        elif haltyp=='host':
            haltpm.readhost()
            halcat=haltpm.hostcat

        self.galcat={}
        self.shamzis=np.arange(shamzibeg,shamziend+1)
        self.scat=scat
        self.gmtype='m.star'
        self.dis_mf=dis_mf
        self.source=source
        self.Mwid=Mwid
        self.mmin=mmin
        self.conscat=conscat
        self.seed=seed
        self.mtree_gal_blt=False
        self.mptree_gal_blt=False

        self.run_sham(halcat,sham_prop)

    def hgram_dat_ft(self,mu_rng,Mcond,zibeg=1,ziend=34,Mwid=None,
                     indices=False,cutoff=-np.inf):
        #Calculate the cumulative probabilities for a galaxy with M(zi=0)=Mcond 
        #to experience a certain number of mergers of raio m/M within mu_rng 
        
        #indices=True includes hi0s_merg in the output.
        #mu_rng is the range of log(mu)'s to include.
        #Given the habit so far of keeping both masses and mass ratios in log,
        #cutuff is the log(M/M0) after which we will stop counting mergers in
        #the given main progenitor.
        #IT IS EXTREMELY IMPORTANT THAT THE USER INPUT CUTOFF WITH THE CORRECT
        #SIGN. E.G. M/M0=1/10 implies cutoff=-1. INCLUDING THE NEGATIVE IS 
        #CRITICAL.

        zis=np.arange(zibeg,ziend+1)
        cat=self.galcat
        mtype=self.gmtype
        mp_brstr='gal.mp.branch'
        m_brstr='gal.merg.branch'
        M0s=cat[0][mtype]
        if Mwid is None:
            Mwid=self.Mwid
        inrange=(M0s<Mcond+Mwid/2.) & (M0s>Mcond-Mwid/2.)
        hi0s=np.arange(len(M0s))[inrange]
        M0s=M0s[inrange]

        #For the histogram later. Need the arrays for the if stmt now.
        bins=np.arange(-0.5,9.5,1)
        count_ax=(bins[1:]+bins[:-1])/2. #x-axis of merger counts
        #print(len(hi0s))
        #print(hi0s)
        if len(hi0s)==0:
            return (np.repeat(np.nan,len(count_ax)), 
                    np.repeat(np.nan,len(count_ax)))

        Ns=[] #initialize the number of mergers in each count bin
        hi0s_merg=[] #z=0 halos that have any mergers in their history

        def Ns_add_f(hi0,M0):
            N=0 #starting off the number of hi0's mergers in the bin at 0
            for zi in zis:
                #2 parents will share a child, so get the mp branch for zi-1
                mpbranch=cat[zi-1][mp_brstr]
                mergbranch=cat[zi][m_brstr]
                if hi0 in mpbranch:  
                    chii=mpbranch[hi0]
                    if chii in mergbranch:
                        merg_is=mergbranch[chii]
                    else:
                        continue
                else:
                    continue
                primi=cat[zi][mp_brstr][hi0]
                M=cat[zi][mtype][primi]
                Mfrac=M-M0
                #If Mfrac falls outside the acceptable range, stop evaluating
                #hi0, and return its N.
                if cutoff and Mfrac <= cutoff:
                    return N
                else:
                    ms=cat[zi][mtype][merg_is]
                    #print(hi0)
                    #print(ms)
                    #print(M)
                    #print('')
                    mus=-np.abs(ms-M)
                    in_mu_rng=(mus>mu_rng[0])&(mus<mu_rng[1])
                    N+=np.sum(in_mu_rng) 
            return N 

        for hi0,M0 in zip(hi0s,M0s):
            Ns_add=Ns_add_f(hi0,M0)
            Ns+=[Ns_add]
            if Ns_add>0:
                hi0s_merg+=[hi0]

        fig=plt.figure()
        ax=fig.add_subplot(111)
        Ps=ax.hist(Ns,bins=bins,cumulative=-1,
                   weights=np.repeat(1./float(len(hi0s)),len(Ns)))[0]
        plt.close()
        
        if indices==True:
            output=(count_ax,Ps,hi0s_merg)
        else:
            output=(count_ax,Ps)
        return output

    def quench_frac_ft(self,min_mu,zibeg=1,ziend=34,
                       M0s=np.linspace(9.7,12.,7.),cutoff=-np.inf):
        #Run hgram_dat_ft for each Mcond in M0s to find the fraction of
        #galaxies that experience at least one mergers of ratio m/M>min_mu as a
        #function of M0

        #M0s=np.linspace(9.5,11.5,17)
        #Mwid=np.average(M0s[1:]-M0s[:-1])
        Mwid=0.1
        fracs=[]
        print('running quench data')
        for M0 in M0s: 
            #print(M0)
            hgram_dat=self.hgram_dat_ft(mu_rng=[min_mu,0.],Mcond=M0,
                                   zibeg=zibeg,ziend=ziend,Mwid=Mwid,
                                   indices=False,cutoff=cutoff)
            #print(hgram_dat[0])
            #print(hgram_dat[1])
            #pull the fraction corresponding to the "at least one" count
            frac=hgram_dat[1][1] 
            fracs+=[frac]
        return M0s,fracs
    
    def run_sham(self,halcat,sham_prop):
        smf_mtrx=pd.read_csv('/home/users/staudt/projects/mergers/data/smf.csv')  
        print('running SHAM')
        if isinstance(self.source,(list,np.ndarray,pd.core.series.Series)):
            pbar=ProgressBar()
            for zi in pbar(self.shamzis):
                source=self.source[zi]
                sham.assign(halcat,scat=self.scat,dis_mf=self.dis_mf,
                            source=source,       
                            sham_prop=sham_prop, 
                            zis=[zi],            
                            seed=self.seed,mmin=self.mmin,const=self.conscat,
                            galcat=self.galcat)
        elif self.source=='rand':
            self.source={}
            pbar=ProgressBar()
            for zi in pbar(self.shamzis):
                #Choose random smf from those allowed at this z.
                allowed=np.array(smf_mtrx.loc[zi,:]).astype('bool')
                all_smf=np.arange(len(allowed))
                choices=smf_mtrx.columns.values[allowed]
                source=random.choice(choices) 
                self.source[zi]=source
                sham.assign(halcat,scat=self.scat,dis_mf=self.dis_mf,
                            source=source,       
                            sham_prop=sham_prop, 
                            zis=[zi],            
                            seed=self.seed,mmin=self.mmin,const=self.conscat,
                            galcat=self.galcat)
        else:
            sham.assign(halcat,scat=self.scat,dis_mf=self.dis_mf,
                        source=self.source,  
                        sham_prop=sham_prop, 
                        zis=self.shamzis,
                        seed=self.seed,mmin=self.mmin,const=self.conscat,
                        galcat=self.galcat)
        return

class HalTreepmClass(TreepmClass):
    def __init__(self,dis_mf=0.,
                 Mwid=0.5,catkind='subhalo',shamzibeg=0,shamziend=34):
    #def __init__(self,scat=0.,dis_mf=0.,source='rand',
    #             Mwid=0.5,catkind='subhalo',shamzibeg=0,shamziend=34,
    #             seed=None,mmin=5.,conscat=True):
        TreepmClass.__init__(self)
        self.smtype='m.max'
        self.hmtype='m.200c'
        self.allzis=np.arange(35)
        self.Mwid=Mwid
        self.mtree_host_blt=False
        self.mtree_sub_blt=False
        self.mptree_host_blt=False
        self.mptree_sub_blt=False
        self.subisread=False
        self.hostisread=False
        if catkind=='subhalo':
            self.readsub()
        elif catkind=='halo':
            self.readhost()
        elif catkind=='both':
            self.readboth()

    def write_mMs(self,M0cond,condtype,zibeg,ziend,Mtime,mutype='gal',
                     N=None,app=''):
        Mtime=str(Mtime)
        tstmp='{:%Y%m%d}'.format(datetime.datetime.now())
        if (self.seed is None) or (not self.scat):
            seedstr=''
        else:
            seedstr='_'+str(self.seed)+'seed'
        if mutype=='hal':
            scatstr=''
        else:
            scatstr='_scat{0:0.2f}'.format(self.scat)
        fname='./dat/{7}_mM{3}_M{4}0_{0:0.2f}{2}{6}_{1}{5}.h5'\
              .format(M0cond,tstmp,scatstr,Mtime,condtype,app,seedstr,mutype)
        with h5py.File(fname,'w') as f:
            if self.subisread:
                f.create_dataset('snap',data=self.subcat.snap)
            elif self.hostisread:
                f.create_dataset('snap',data=self.hostcat.snap)
            out=self.mMs_f(M0cond,condtype,zibeg,ziend,Mtime,mutype,N)
            f.create_dataset('mMs',data=out[0])
            f.create_dataset('ms',data=out[1])
            f.create_dataset('Ms',data=out[2])
            f.create_dataset('iprim0s_merg',data=out[3])
            f.create_dataset('iprim0s_keys',data=out[4])
            f.create_dataset('Nprim0',data=out[5])
            f.create_dataset('zbeg',data=out[6])
            f.create_dataset('zend',data=out[7])
            f.create_dataset('zs',data=out[8])

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
            print('\nsnapshot %i'%zi)
            #Get host indices whose masses at z match M:
            hostis=elements(self.hostcat[zi][self.hmtype],
                            lim=[M-self.Mwid/2.,M+self.Mwid/2.])
            #Get the corresponding central subhalo indices from the subhalo 
            #cat:
            cenis=self.hostcat[zi]['cen.i'][hostis]
            print('%i halos in range.'%len(hostis))
            if not N is None:
                random.seed(1)
                cenis=random.sample(cenis,N)
            print('evaluating %i'%len(cenis))
            Nhost+=len(cenis)

            pbar=ProgressBar()
            for ceni,hosti in pbar(list(zip(cenis,hostis))):
                #Take M_z to be the galaxy mass sham'd to the HOST halo at z:
                M=self.hostcat[zi][self.gmtype][hosti]
                is_prog=self.subcat[zi+1]['chi.i']==ceni
                is_prog[self.subcat[zi]['par.i'][ceni]]=False
                #Take m_z to be the mass of the accreted galaxy at zi+1: 
                ms=self.subcat[zi+1][self.gmtype][is_prog]
                mMs_add=ms-M
                mMs+=list(mMs_add)
        return mMs,Nhost,zbeg,zend

    def merg_tree(self,zibeg=0,ziend=34,typ='gal',gal_tpm=None,prnt=False):
        if typ not in ['gal','subhal','host']:
            raise ValueError('typ must be "gal", "subhal" or "host".')
        if typ in ['gal','subhal']:
            if not self.subisread:
                self.readsub()
            cat=self.subcat
            if typ=='gal':
                if gal_tpm is None:
                    raise ValueError('type is "gal", but gal_tpm not provided')
                mtype=gal_tpm.gmtype
                #Make a pointer to the cat into which we'll put the merger tree
                objcat=gal_tpm.galcat 
            elif typ=='subhal':
                mtype=self.smtype
                #Make a pointer to the cat into which we'll put the merger tree
                objcat=self.subcat
        elif typ=='host':
            if not self.hostisread:
                self.readhost()
            cat=self.hostcat
            #Make a pointer to the cat into which we'll put the merger tree
            objcat=self.hostcat
            mtype=self.hmtype
        
        zis=np.arange(zibeg,ziend+1)
        zbeg=cat.snap[zibeg][1]
        zend=cat.snap[ziend][1]
        hi0s=np.arange(len(objcat[0][mtype]))

        print('building merger tree')
        #Evaluate zis in reverse order from highest redshift through second to
        #lowest redshift:
        for zi in zis[:0:-1]:
            #initialize merg.branch and par.branch branches at zi as a 
            #dictionaries:
            if typ=='gal':
                brstr='gal.merg.branch'
            elif typ=='subhal':
                brstr='sub.merg.branch'
            elif typ=='host':
                brstr='merg.branch'
            objcat[zi][brstr]={}
            mergbranch=objcat[zi][brstr]

            his=np.arange(len(objcat[zi][mtype]))
            #Take each halo at zi and find out either what it merged into
            #or what it's new index is at zi-1:
            chiis=indices_tree(cat,zi,zi-1,his)
            
            if prnt:
                if typ in ['host','subhal']:
                    print('{0:d} halos in snapshot {1:d}'.format(len(his),zi))
                elif typ=='gal': 
                    print('{0:d} galaxies in snapshot {1:d}'.format(len(his),zi))

            allms_chi=objcat[zi-1][mtype]
            allms=objcat[zi][mtype]
            if (len(allms_chi)==0)or(len(allms)==0):
                #At high z, there may be no mass array.
                continue

            #Only evaluate where the halo's child exists:
            mask=chiis>=0
            #Don't count primary parent as a merger 
            paris=indices_tree(cat,zi-1,zi,chiis[mask])
            mask[paris]=False
            his=his[mask]
            chiis=chiis[mask]
            
            #Only evaluate where both the primary and secondary object have
            #positive masses
            primis=cat[zi-1]['par.i'][chiis]
            prim_ms=objcat[zi][mtype][primis]
            sec_ms=objcat[zi][mtype][his]
            notnan=~(np.isnan(sec_ms) | np.isnan(prim_ms))
            #Get a mask for notnan objects that have both prim and sec masses:
            hasmasses=prim_ms[notnan]*sec_ms[notnan]>0.
            #Apply the notnan mask to get the notnan arrays and then apply the
            #hasmasses mask.
            his=his[notnan][hasmasses]
            chiis=chiis[notnan][hasmasses]

            '''
            #Only evaluate where both halo and its child have positive masses
            ms_chi=allms_chi[chiis]
            ms=allms[his]
            notnan=~(np.isnan(ms_chi) | np.isnan(ms))
            hasmasses=(ms[notnan]>0.) & (ms_chi[notnan]>0.)
            chiis=chiis[notnan][hasmasses]
            his=his[notnan][hasmasses]
            '''
            
            if prnt:
                if typ in ['host','subhal']:
                    print('{0:d} merge with larger halos'.format(len(his),zi))
                elif typ=='gal':
                    print('{0:d} merge with larger galaxies'.format(len(his),zi))

            for hi,chii in zip(his,chiis):
                if chii not in mergbranch: 
                    mergbranch[chii]=[hi]
                else:
                    mergbranch[chii]+=[hi]

            '''
            #Rename branch keys with 0 indices of main progenitor
            #Removing this because some mergers are not mergers into the main
            #progenitor but mergers into a galaxy that will itself merge into
            #the main progenitor
            i1s=indices_tree(cat,zi-1,0,mergbranch.keys())
            for chii,i0 in zip(mergbranch.keys(),i0s):
                mergbranch[i0]=mergbranch.pop(chii)
            '''
        if typ=='gal':
            gal_tpm.mtree_gal_blt=True
        elif typ=='subhal':
            self.mtree_sub_blt=True
        elif typ=='host':
            self.mtree_host_blt=True

    def mp_tree(self,typ='gal',zibeg=0,ziend=34,gal_tpm=None):
        if typ not in ['gal','subhal','host']:
            raise ValueError('typ must be "gal", "subhal" or "host".')
        if typ in ['gal','subhal']:
            if not self.subisread:
                self.readsub()
            #Given that all the index information is stored in the halo
            #catalog, I'm keeping mp_tree() and merg_tree() in 
            #HalTreepmClass(), which contains the halo catalog. This then
            #requires me to add a gal_tpm argument so, if we're building a 
            #galaxy mp_tree, the program knows the tpm
            #into which it should put the mp_tree.
            cat=self.subcat
            if typ=='gal':
                if gal_tpm is None:
                    raise ValueError('type is "gal", but gal_tpm not provided')
                mtype=gal_tpm.gmtype
                #Make a pointer to the cat into which we'll put the main
                #progenitor tree:
                objcat=gal_tpm.galcat
            elif typ=='subhal':
                mtype=self.smtype
                #Make a pointer to the cat into which we'll put the main
                #progenitor tree:
                objcat=self.subcat
        elif typ=='host':
            if not self.hostisread:
                self.readhost()
            cat=self.hostcat
            #Make a pointer to the at into which we'll put the main
            #progenitor tree:
            objcat=self.hostcat
            mtype=self.hmtype

        zis=np.arange(zibeg,ziend+1)
        zbeg=cat.snap[zibeg][1]
        zend=cat.snap[ziend][1]

        hi0s=np.arange(len(objcat[0][mtype]))

        #Evaluate zis in ascending order, tracing back z=0 main progenitors so
        #we have a list of main progenitors at every snapshot.
        print('building main progenitor tree:')
        pbar=ProgressBar()
        for zi in pbar(zis):
            #initialize mp.branch
            if typ=='gal':
                brstr='gal.mp.branch'
            elif typ=='subhal':
                brstr='sub.mp.branch'
            elif typ=='host':
                brstr='mp.branch'
            objcat[zi][brstr]={}
            mpbranch=objcat[zi][brstr]

            #fill this zi branch of the main-progenitor tree
            #Only add paris for mps that have positive indices and mass at z.
            #Make a copy of hi0s. "exatz" stands for "exists at z":
            hi0s_exatz=hi0s.copy() 
            paris=indices_tree(cat,0,zi,hi0s_exatz)

            exists=paris>=0
            paris=paris[exists]
            hi0s_exatz=hi0s_exatz[exists]

            ms=objcat[zi][mtype]
            if len(ms)==0:
                #At high z, there may be no mass array.
                continue
            ms=ms[paris]

            #Sometimes, especially (and maybe only) because of SHAM, masses are
            #nan. Get rid of the galaxies/halos corresponding to those.
            notnan=~np.isnan(ms)
            paris=paris[notnan]
            hi0s_exatz=hi0s_exatz[notnan]
            ms=ms[notnan]

            hasmass=ms>0.
            paris=paris[hasmass]
            hi0s_exatz=hi0s_exatz[hasmass]

            for hi0,pari in zip(hi0s_exatz,paris):
                mpbranch[hi0]=pari
        if typ=='gal':
            gal_tpm.mptree_gal_blt=True
        elif typ=='subhal':
            self.mptree_sub_blt=True
        elif typ=='host':
            self.mptree_host_blt=True
        return

    def gal_mMs_fromtree(self,M0cond,condtype,zibeg,ziend,Mtime,N=None):
        if not self.mtreebuilt:
            self.merg_tree(zibeg,ziend)
        zis=np.arange(zibeg,ziend+1)
        zbeg=self.subcat.snap[zibeg][1]
        zend=self.subcat.snap[ziend][1]
        mMs=[]
        ms=[]
        Ms=[]
        zs=[]
        iprim0s_merg=[]
        
        i_primaries0=self.get_primaries0(condtype,M0cond)
        
        for zi in zis[1:]:
            i_primary=indices_tree(self.subcat,0,zi,i_primaries0)
            Mz=self.subcat[zi]['m.star'][i_primary]
            if Mz==0.:
                continue
            ms_add=self.subcat[zi]['par.tree'][i_primary]
        return

    def dNdx_ofz(self,Mcond,mu_cond,typ,dx='dz',zibeg=1,ziend=34,
                 through=False,zi_Mcond=0,galtpm=None):
        #returns dNdz(z) or dNdt(z) around specific mu=m/M and M0 values.
        if not typ in ['gal','cengal','host','subhal','censubhal']:             
            raise ValueError('typ must be "gal", "cengal", "subhal", '          
                             '"censubhal", or "host"')                
        if typ in ['gal','cengal']:                                             
            if galtpm is None:
                raise ValueError('Galaxy type specified but no galaxy catalog'
                                 'provided')
            mtype=galtpm.gmtype                                                   
            #set the merger-branch and main-progenitor-branch strings           
            m_brstr='gal.merg.branch'                                           
            mp_brstr='gal.mp.branch'                                            
            self.readsub()                                                      
            halcat=self.subcat                                                     
            mcat=galtpm.galcat
            if not galtpm.mptree_gal_blt:                                         
                galtpm.mp_tree(typ='gal')                               
            if not galtpm.mtree_gal_blt:                                          
                galtpm.merg_tree(typ='gal')                             
        elif typ in ['subhal','censubhal']:                                     
            if galtpm is not None:
                raise ValueError('Halo type specified, but the user provided a'
                                 'galaxy catalog.')
            mtype=self.smtype                                                   
            #set the merger-branch and main-progenitor-branch strings           
            m_brstr='sub.merg.branch'                                           
            mp_brstr='sub.mp.branch'                                            
            self.readsub()                                                      
            halcat=self.subcat                                                     
            mcat=halcat
            if not self.mptree_sub_blt:                                         
                self.mp_tree(typ='subhal')                               
            if not self.mtree_sub_blt:                                          
                self.merg_tree(typ='subhal')                             
        elif typ=='host':                                                       
            if galtpm is not None:
                raise ValueError('Halo type specified, but the user provided a'
                                 'galaxy catalog.')
            mtype=self.hmtype                                                   
            #set the merger-branch and main-progenitor-branch strings           
            m_brstr='merg.branch'                                               
            mp_brstr='mp.branch'                                                
            self.readhost()                             
            halcat=self.hostcat                                                    
            mcat=halcat
            if not self.mptree_host_blt:                                        
                self.mp_tree(typ='host')                               
            if not self.mtree_host_blt:                                         
                self.merg_tree(typ='host')                  
        
        zis=np.arange(zibeg,ziend+1)   
        allM0s=mcat[zi_Mcond][mtype]
        inrange=(allM0s<Mcond+self.Mwid/2.) & (allM0s>Mcond-self.Mwid/2.)
        if typ in ['cengal','subhal']:
            iscen=halcat[zi_Mcond]['ilk']==1
            hi0s=np.arange(len(allM0s))[inrange & iscen]
        else:
            hi0s=np.arange(len(allM0s))[inrange]    
        if through:
            hi0s=through_f(mcat,hi0s,mtype,ziend)
        #M0s=mcat[0][mtype][hi0s] #Do I need this?
        
        if zi_Mcond!=0:
            #In order to implement the ability to set an M condition at a zi other
            #than 0 while still using the mp_tree indexed to 0, the code finds hi0s
            #based on the conditions at zi_Mcond and then uses indices_tree to find 
            #out what the indices are at zi=0.
            #
            #I might run into an issue where a halo/gal exists at zi_Mcond but not
            #at zi=0. I'll come back and fix this later if that issue arrises.
            hi0s=indices_tree(halcat,zi_Mcond,0,hi0s)

        Nprim0=len(hi0s)
        dNdxs=[]
        zs=[halcat[zi].snap['z'] for zi in zis]
        ts=[halcat[zi].snap['t'] for zi in zis] #not currently used

        for zi in zis:
            dN=0.
            mergbranch=mcat[zi][m_brstr]
            #Set the mpbranch one snapshot behind, because the keys         
            #of the mergbranch are the children into which two or more          
            #hals/gals merg.     
            mpbranch=mcat[zi-1][mp_brstr]
           
            for hi0 in hi0s:
                if hi0 in mpbranch:                
                    chii=mpbranch[hi0]
                    if chii in mergbranch:
                        merg_is=mergbranch[chii]
                    else:
                        continue
                else:
                    continue
                    
                primi=mcat[zi][mp_brstr][hi0]
                if primi<0:
                    print('hi0 {0:d} has a negative index at zi='
                          '{1:d}'.format(hi0,zi))
                    continue
                M=mcat[zi][mtype][primi]
                ms=mcat[zi][mtype][merg_is]
                
                mus=-np.abs(ms-M)
                in_mu_range=mus>=mu_cond
                dN+=np.sum(in_mu_range)
                
            if dx=='dz':
                dx_val=halcat[zi].snap['z']-halcat[zi-1].snap['z']
                #Can't use zs[zi] because zs aren't necessarily indexed to zis.
            elif dx=='dt':
                dx_val=halcat[zi-1].snap['t']-halcat[zi].snap['t']
            dNdx=float(dN)/dx_val
            #print(zi)
            #print(halcat[zi].snap['z'])
            #print('dN: {0}'.format(dN))
            #print('dx: {0:0.2f}'.format(dx))
            #print('dNdx: {0:0.4f}'.format(dNdx))
            #print('')
            dNdxs+=[dNdx]
        dNdxs=np.array(dNdxs)
        dNdxs/=float(Nprim0)
        return zs,dNdxs

    def gal_mMs(self,M0cond,condtype,zibeg,ziend,Mtime,N=None):
        #***CANDIDATE FOR REMOVAL BECAUSE I THINK GAL_MmS_FROMTREE() COVERS
        #THIS***

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
        
        i_primaries0=self.get_primaries0(condtype,M0cond,throughout=True)
        
        print('%i galaxies in range'%len(i_primaries0))
        if not N is None:
            random.seed(1)
            i_primaries0=random.sample(i_primaries0,N)
        Nprim0=len(i_primaries0)
        print('evaluating %i'%Nprim0)

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

    def mMs_f(self,M0cond,condtype,zibeg,ziend,Mtime,mutype='gal',N=None):
        #***CANDIDATE FOR REMOVAL ONCE I WRITE AN MmS_FROMTREE THAT CAN HANDLE
        #EITHER GALAXIES OR HALOS. RIGHT NOW GAL_MmS_FROMTREE IS GOOD, BUT FOR
        #SOME REASON I WROTE IT SO IT ONLY HANDLES GALAXIES.
        
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
        
        i_primaries0=self.get_primaries0(condtype,M0cond,throughout=True)
        
        print('%i galaxies in range'%len(i_primaries0))
        if not N is None:
            random.seed(1)
            i_primaries0=random.sample(i_primaries0,N)
        Nprim0=len(i_primaries0)
        print('evaluating %i'%Nprim0)

        pbar=ProgressBar()
        for i_primary0 in pbar(i_primaries0):
            out=self.mMs_add_f(i_primary0,zis,Mtime,mutype)
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
    
    def dndm(self,snap=0,paramed=False,typ='gal',galtpm=None,num=50):

        #Get dn/dlogm at given snapshot.
        #num sets number of mass bins.
        #This is different from dndm in the sham module because this checks the
        #dndm that actually resulted within the mock, as opposed to what dndm
        #should be in theory.
        if typ=='gal':
            if galtpm is None:
                raise ValueError('Type is "gal", but no galtpm is provided.')
            galcat=galtpm.galcat
            allms=galcat[snap][galtpm.gmtype] #Get all masses at snap
        elif typ=='hal':
            allms=self.subcat[snap][self.smtype] #Get all masses at snap
        mask1=allms!=0.
        mask2=~np.isnan(allms)
        mask=mask1*mask2 #Get halo index for every halo with nonzero mass and ID not equal to -1.
        ms=allms[mask] #Get the corresponding masses                            
        #self.mstest=ms
        lower=np.min(ms) #Lower histogram bound                                 
        upper=np.max(ms) #Upper histogram bound                                 
        bins=np.linspace(lower,upper,num) #Mass bins in log                     
        binw=np.average(bins[1:]-bins[:-1]) #Bin width in dex                   
        midbins=(bins[:-1]+bins[1:])/2.
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        N,binsout,patches=ax.hist(ms,bins)                                     
        #plt.clf()                                                               
        plt.close()                                                               
        
        #number density i.e. number per volume per bin size:
        if paramed:
            n=N/(250.)**3./binw 
        else:
            n=N/(250./0.7)**3./binw
        return midbins,n                                 

    def shmr_sham(self,zi=0,galtpm=None,source=None,scat=None,
                  mmin=5.):
        #shmr found by integrating the simulated dndm for galaxies and halos.
        if galtpm is None:
            if source is None:
                raise ValueError('galtpm not provided but ' 
                                 'no SMF name'
                                 'is provided in the source argument.')
            if scat is None:
                raise ValueError('galtpm not provided but '
                                 'no scatter value is'
                                 ' specified.')
            #Creating a temporary tpm instance into which to put galaxy masses
            #for the dndm calcualation
            galtpm=TreepmClass()
            galtpm.gmtype='m.star'
            galtpm.galcat={}

            sham.assign(self.subcat,scat=scat,source=source,zis=[zi],seed=None,
                        mmin=mmin,const=False,galcat=galtpm.galcat)
            ms_gal,dndms_gal=self.dndm(zi,typ='gal',paramed=False,
                                       galtpm=galtpm,num=500)
        else:
            if source or scat:
                raise ValueError('You provided a galtpm but also specified '
                                 'attributes meant for when galtpm is not '
                                 'provided.')
            ms_gal,dndms_gal=self.dndm(zi,typ='gal',paramed=False,
                                       galtpm=galtpm)
            
        ms_hal,dndms_hal=self.dndm(zi,typ='hal',paramed=False,num=500)
        
        #reverse the arrays so dndms are in ascending order:
        ms_gal=ms_gal[::-1]
        dndms_gal=dndms_gal[::-1]
        ms_hal=ms_hal[::-1]
        dndms_hal=dndms_hal[::-1]
        
        dm_gal=(ms_gal.max()-ms_gal.min())/float(len(ms_gal))
        dm_hal=(ms_hal.max()-ms_hal.min())/float(len(ms_hal))

        #Need a negative sign in front because masses are in descending order.
        numdens_gal=-scipy.integrate.cumtrapz(dndms_gal,ms_gal,initial=0)
        numdens_hal=-scipy.integrate.cumtrapz(dndms_hal,ms_hal,initial=0)

        #Switch to lognumdens.
        numdens_gal=np.log10(numdens_gal)
        numdens_hal=np.log10(numdens_hal)

        #maximum difference between a target numden and the numden that lookup
        #finds in the lookup array in order for a result to register:
        step=np.average(numdens_gal[1:]-numdens_gal[:-1])
        #Initialize ms_hal corresponding to ms_gal.
        ms_hal_shmr=np.zeros(len(ms_gal))
        for i in range(len(ms_gal)):
            numden=numdens_gal[i]
            ms_hal_shmr[i]=sutils.lookup(numden,numdens_hal,ms_hal,step)
        
        return ms_hal_shmr,ms_gal#,numdens_gal

    def shmr_avg(self,ms_hal,galtpm=None,zi=0,Mwid=0.004,source=None,scat=None,
                 mmin=5.):
        #SHMR based on the average resulting galaxy mass for givin halo-mass
        #buckets
        if galtpm is None:
            if source is None:
                raise ValueError('galtpm not provided but ' 
                                 'no SMF name'
                                 'is provided in the source argument.')
            if scat is None:
                raise ValueError('galtpm not provided but '
                                 'no scatter value is'
                                 ' specified.')
            #Creating a temporary tpm instance into which to put galaxy masses
            #for the avg mass calcualations
            galtpm=TreepmClass()
            galtpm.gmtype='m.star'
            galtpm.galcat={}
            sham.assign(self.subcat,scat=scat,source=source,zis=[zi],seed=None,
                        mmin=mmin,const=False,galcat=galtpm.galcat)
        else:
            if source or scat:
                raise ValueError('You provided a galtpm but also specified '
                                 'attributes meant for when galtpm is not '
                                 'provided.')
        ms_gal=np.zeros(len(ms_hal))
        halcatz=self.subcat[zi]
        galcatz=galtpm.galcat[zi]
        all_hal_ms=halcatz[self.smtype]
        all_gal_ms=galcatz[galtpm.gmtype]
        notzero=all_gal_ms>0
        notnan=~np.isnan(all_gal_ms)
        all_hal_ms=all_hal_ms[notzero & notnan]
        all_gal_ms=all_gal_ms[notzero & notnan]
        for i,m_hal in enumerate(ms_hal):
            matches_m_hal=((all_hal_ms>m_hal-Mwid/2.) 
                           & (all_hal_ms<m_hal+Mwid/2.))
            ms_gal[i]=np.average(all_gal_ms[matches_m_hal])
        return ms_gal

    def get_primaries0(self,condtype,M0cond,throughout=False):
        if condtype in ['gal','subhal']:
            i_primaries0=self.get_is(M0cond,condtype,throughout=throughout)
        elif condtype in ['censubhal','cengal']:
            i_primaries0=self.get_cenis(M0cond,condtype,throughout=throughout)
        elif condtype=='hal':
            i_primaries0=self.get_halis(M0cond,throughout=throughout)
        else:
            raise ValueError('Condition type is "cengal", "gal", "subhal", '
                             '"censubhal", or "hal".')
        return i_primaries0

    def get_is(self,M,condtype,zi=0,throughout=False):
        if condtype=='subhal':
            mtype=self.smtype
        elif condtype=='gal':
            mtype=self.gmtype
        else:
            raise ValueError('Unexpected condition type error')
        
        il=elements(self.subcat[zi][mtype],
                       lim=[M-self.Mwid/2.,M+self.Mwid/2.])
        if throughout:
            ziend=self.allzis[-1]
            il_fin=indices_tree(self.subcat,zi,ziend,il)
            hasfin=il_fin>=0
            il=il[hasfin]
        return il

    def get_cenis(self,M,condtype,zi=0,throughout=False):
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
        if throughout:
            ziend=self.allzis[-1]
            cenis_fin=indices_tree(self.subcat,zi,ziend,cenis)
            hasfin=cenis_fin>=0
            cenis=cenis[hasfin]
        return cenis

    def get_halis(self,M,zi=0,throughout=False):
        mtype=self.hmtype
        il=elements(self.hostcat[zi][mtype],
                    lim=[M-self.Mwid/2.,M+self.Mwid/2.])
        if throughout:
            ziend=self.allzis[-1]
            il_fin=indices_tree(self.subcat,zi,ziend,il)
            hasfin=il_fin>=0
            il=il[hasfin]
        return il

    def mMs_add_f(self,i_primary0,zis,Mtime,mutype):
        #CANDIDATE FOR REMOVAL ONCE I REMOVE ALL FUNCTIONS THAT USE THIS
        #FUNCTION, BECAUSE THIS DOES NOT USE THE TREE.

        #Get a list of m_gal/M_gal ratios spanning zis for a single z=0 subhalo
        #M_gal can either be M0 or Mz depending on Mtime.
        if not Mtime in ['0','z']:
            raise ValueError('Mtime is either "0" or "z".')
        if not mutype in ['gal','hal']:
            raise ValueError('mutype is either "gal" or "hal".')
        mMs=[]
        ms=[]
        Ms=[]
        zs=[]
        zi0=zis[0]
        if Mtime=='0':
            if mutype=='gal':
                M=self.subcat[zi0][self.gmtype][i_primary0]
            elif mutype=='hal':
                M=self.subcat[zi0][self.smtype][i_primary0]
        for zi in zis[:-1]:
            #index of the main progenitor at zi:
            i_primary=indices_tree(self.subcat,zi0,zi,i_primary0)
            if i_primary<0:
                return mMs,ms,Ms,zs
            #index of the main progenitor at zi+1:
            pari=self.subcat[zi]['par.i'][i_primary] 
            
            if mutype=='gal':
                Mz=self.subcat[zi+1][self.gmtype][pari]
            elif mutype=='hal':
                Mz=self.subcat[zi+1][self.smtype][pari]
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
            if mutype=='gal':
                ms_add=self.subcat[zi+1][self.gmtype][is_prog]
            elif mutype=='hal':
                ms_add=self.subcat[zi+1][self.smtype][is_prog]
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
        if not self.hostisread:
            self.hostcat=self.read(zis=self.allzis,catalog_kind='halo')
            self.hostisread=True
        return

    def readsub(self):
        if not self.subisread:
            self.subcat=self.read(zis=self.allzis,catalog_kind='subhalo')
            self.subisread=True
        return

    def readboth(self):
        self.readhost()
        self.readsub()
        return

###############################################################################

def read_mMs(mutype,Mcond,condtype,Mtime,scat,tstmp,seed=None,apnd=''):
    Mtime=str(Mtime)
    if (seed is None) or (not scat):
        seedstr=''
    else:
        seedstr='_'+str(seed)+'seed'
    if mutype=='hal':
        scatstr=''
    else:
        scatstr='_scat{0:0.2f}'.format(scat)
    fname='/home/users/staudt/projects/mergers/dat/' \
          '{7}_mM{4}_M{3}0_{0:0.2f}{2}{5}_{1}{6}.h5' \
          .format(Mcond,tstmp,scatstr,condtype,Mtime,seedstr,apnd,mutype)
    print(fname)

    with h5py.File(fname,'r') as f:
        mMs=np.array(f['mMs'])
        ms=np.array(f['ms'])
        Ms=np.array(f['Ms'])
        iprim0s_merg=np.array(f['iprim0s_merg'])
        iprim0s_keys=np.array(f['iprim0s_keys'])
        Nprim0=(f['Nprim0']).value
        zbeg=(f['zbeg']).value
        zend=(f['zend']).value
        zs=np.array(f['zs'])
        snap=np.array(f['snap'])
    return mMs,ms,Ms,iprim0s_merg,iprim0s_keys,Nprim0,zbeg,zend,zs,snap

def write_gal_hgram_dat(Mconds,condtype,Mtime,rngs,zmax,zmin,
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
                                      zmax,zmin,mMsdat_scat[8])
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

def gal_hgram_dat(rng,mMs,iprim0s_merg,iprim0s_keys,zmax,zmin,zs):
    #Generate a list of the number of times a merger happens within the 
    #given range. Each datapoint represents the merger count for a single
    #main progenitor.
    mMlo=rng[0]
    mMhi=rng[1]
    inrange=(mMs>mMlo)&(mMs<mMhi)&(zs<=zmax)&(zs>=zmin)
    Ns=np.array([np.sum(ismp & inrange)
                 for ismp in
                 (iprim0s_merg==key for key in iprim0s_keys)])
    return Ns

def dndz_ofz(mMs,zs,mu_min,zkeys,Nprim0,
             zibeg,ziend):
    #The user should feed this function the zibeg and ziend corresponding to 
    #the mergers' LABEL. Keep in mind that mergers are labeled with the z
    #immediately before the merger occurs. Therefore, a merger that occurs
    #between zi=1 and zi=0 will be labeled with z[1].
    #**1 is the lowest zibeg the user should feed into this function.
    dNdzs=[]
    zaxis=(zkeys[zibeg:ziend+1]+zkeys[zibeg-1:ziend])/2.
    #print zaxis
    zis=np.arange(zibeg,ziend+1)
    in_murange=mMs>=mu_min
    for zi in zis:
        zmin=zkeys[zi-1]
        zmax=zkeys[zi]
        #print'zi: %d'%zi
        #print'zmin: %f'%zmin
        #print'zmax: %f'%zmax
        dz=zmax-zmin
        #Mergers are labeled with the z of the snapshot right before the 
        #merger, discreetly, so we only need to search for equality, not a 
        #range:
        atz=zs==zmax
        #print'atz: %d'%np.sum(atz)
        inrange=in_murange & atz
        #print'inrange: %d'%np.sum(inrange)
        #print np.sum(inrange)
        dNdz=[np.sum(inrange)/dz/float(Nprim0)]
        #print dNdz
        dNdzs+=dNdz
    return dNdzs,zaxis

def count(mMs,Nhost):
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

def through_f(cat,hi0s,mtype,ziend,mcat=None):                                                
    if mcat is None:
        mcat=cat
    #take only main progenitors that exist at ziend                             
    hiends=indices_tree(cat,0,ziend,hi0s)                                       
    goesthrough=hiends>=0                                                       
    hi0s=hi0s[goesthrough]                                                      
    hiends=hiends[goesthrough]                                                  
    #take only main progenitors whose mass is not nan at                        
    #ziend                                                                      
    msend=mcat[ziend][mtype][hiends]                                             
    notnan=~np.isnan(msend)                                                     
    hi0s=hi0s[notnan]                                                           
    hiends=hiends[notnan]                                                       
    msend=msend[notnan]                                                         
    #take only main progenitors that have non-zero mass at                      
    #ziend                                                                      
    pos_msend=msend>0.                                                          
    hi0s=hi0s[pos_msend]                                                        
    return hi0s

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

def dNdx_ofmu(self,Mcond,typ,Mtime='z',dx='dz',forcem200=False,zibeg=0,
        ziend=1,through=False):
    if not typ in ['gal','cengal','host','subhal','censubhal']:             
        raise ValueError('typ must be "gal", "cengal", "subhal", '          
                         '"censubhal", or "host"')                
    if typ in ['gal','cengal']:                                             
        mtype=self.gmtype                                                   
        #set the merger-branch and main-progenitor-branch strings           
        m_brstr='gal.merg.branch'                                           
        mp_brstr='gal.mp.branch'                                            
        self.readsub()                                                      
        cat=self.subcat                                                     
        if not self.mptree_gal_blt:                                         
            self.mp_tree(typ='gal')                               
        if not self.mtree_gal_blt:                                          
            self.merg_tree(typ='gal')                             
    elif typ in ['subhal','censubhal']:                                     
        mtype=self.smtype                                                   
        #set the merger-branch and main-progenitor-branch strings           
        m_brstr='sub.merg.branch'                                           
        mp_brstr='sub.mp.branch'                                            
        self.readsub()                                                      
        cat=self.subcat                                                     
        if not self.mptree_sub_blt:                                         
            self.mp_tree(typ='subhal')                               
        if not self.mtree_sub_blt:                                          
            self.merg_tree(typ='subhal')                             
    elif typ=='host':                                                       
        mtype=self.hmtype                                                   
        #set the merger-branch and main-progenitor-branch strings           
        m_brstr='merg.branch'                                               
        mp_brstr='mp.branch'                                                
        self.readhost()                             
        cat=self.hostcat                                                    
        if not self.mptree_host_blt:                                        
            self.mp_tree(typ='host')                               
        if not self.mtree_host_blt:                                         
            self.merg_tree(typ='host')                  
    mus,dNs=N_mu_ft(self,Mcond,typ,Mtime,forcem200,zibeg,ziend,through)
    if dx=='dz':
        zbeg=cat[zibeg].snap['z']
        zend=cat[ziend].snap['z']
        dx=zend-zbeg
    if dx=='dt':
        tbeg=cat[zibeg].snap['t']
        tend=cat[ziend].snap['t']
        dx=tbeg-tend
    dNdxs=dNs/dx
    return mus,dNdxs 

def N_mu_ft(self,M0,typ,Mtime='z',forcem200=False,zibeg=0,ziend=34,
            through=False):
    #Plots the cumulative number of mergers per                             
    #halo/galaxy, N(>mu), where mu is m/M.                                  
    #"ft" stands for "function using tree"                                  
    if not typ in ['gal','cengal','host','subhal','censubhal']:             
        raise ValueError('typ must be "gal", "cengal", "subhal", '          
                         '"censubhal", or "host"')                          
    if typ in ['gal','cengal']:                                             
        if forcem200:
            raise ValueError('forcem200 and galaxy analysis are incompatible')
        mtype=self.gmtype                                                   
        #set the merger-branch and main-progenitor-branch strings           
        m_brstr='gal.merg.branch'                                           
        mp_brstr='gal.mp.branch'                                            
        self.readsub()                                                      
        cat=self.subcat                                                     
        if not self.mptree_gal_blt:                                         
            self.mp_tree(typ='gal')                               
        if not self.mtree_gal_blt:                                          
            self.merg_tree(typ='gal')                             
    elif typ in ['subhal','censubhal']:                                     
        mtype=self.smtype                                                   
        #set the merger-branch and main-progenitor-branch strings           
        m_brstr='sub.merg.branch'                                           
        mp_brstr='sub.mp.branch'                                            
        self.readsub()                                                      
        cat=self.subcat                                                     
        if not self.mptree_sub_blt:                                         
            self.mp_tree(typ='subhal')                               
        if not self.mtree_sub_blt:                                          
            self.merg_tree(typ='subhal')                             
    elif typ=='host':                                                       
        mtype=self.hmtype                                                   
        #set the merger-branch and main-progenitor-branch strings           
        m_brstr='merg.branch'                                               
        mp_brstr='mp.branch'                                                
        self.readhost()                             
        cat=self.hostcat                                                    
        if not self.mptree_host_blt:                                        
            self.mp_tree(typ=typ)                               
        if not self.mtree_host_blt:                                         
            self.merg_tree(typ=typ)                             

    #Filter for gals/hals that meet the mass condition:
    if forcem200 & (typ!='host'):
        self.readhost()
        hosti0s=self.subcat[0]['halo.i']
        allM0s=self.hostcat[0]['m.200c'][hosti0s]
    else:
        allM0s=cat[0][mtype]                                                    
    inrange=(allM0s<M0+self.Mwid/2.) & (allM0s>M0-self.Mwid/2.)             
    hi0s=np.arange(len(allM0s))[inrange]

    if typ in ['cengal','censubhal']:                                       
        iscen=cat[0]['ilk'][hi0s]==1                                        
        hi0s=hi0s[iscen]                                                    
        #print'%d halos match the M0 condition and are centrals'%hi0s.size
    
    if through:                                                             
        #Filter for gals/hals that exist for all zi up to ziend
        hi0s=through_f(cat,hi0s,mtype,ziend)

    allzis=np.arange(zibeg,ziend+1)                                   
    mzMs=[]                                                                 

    print('\nrunning N(>mu)')                                                
    
    self.checkdic={}
    self.mcheckdic={}
    self.Mcheckdic={}
    #Check the main progenitor tree at each redshift to see if hi0 existed
    #then. If it does, check for mergers in the merger tree.
    pbar=ProgressBar()
    for hi0 in pbar(hi0s):                                                        
        for zi in allzis[1:]:                                               
            mergbranch=cat[zi][m_brstr]                                     
            mpbranch=cat[zi-1][mp_brstr]                                    
            if hi0 in mpbranch:                                             
                chii=mpbranch[hi0]                                          
                if chii in mergbranch:                                      
                    #If the mp exists at z and it had mergers there, 
                    #get the indices
                    #of the merged secondary progenitors.
                    merg_is=mergbranch[chii]                                
                else:                                                       
                    #If no mergers at z, move to next z.
                    continue                                                
            else:                                                           
                #If the main progenitor doesn't exist at z, move to the next z.
                #I would assume that if the mp doesn't exist at z, it won't
                #exist at any of the following zs, but I'm letting it run
                #through to be safe.
                continue                                                    
            #If we are forcing masses to be m200
            if forcem200 & (typ!='host'):
                self.readhost()
                inf_zis=cat[zi]['inf.last.zi'][merg_is]
                inf_is=cat[zi]['inf.last.i'][merg_is]
                inf_hostis=np.array([cat[inf_zi]['halo.i'][inf_i] 
                                     for inf_zi,inf_i in zip(inf_zis,inf_is)])
                ms=np.array([self.hostcat[inf_zi]['m.200c'][inf_hosti]
                            for inf_zi,inf_hosti in zip(inf_zis,inf_hostis)])
                if -1 in inf_zis:
                    #In some cases, subhalos form inside a host halo and don't
                    #have an infall time, per se. Get rid of those:
                    is_imac=inf_zis==-1
                    inf_zis=inf_zis[~is_imac]
                    if len(inf_zis)==0:
                        continue
                    inf_is=inf_is[~is_imac]
                    inf_hostis=inf_hostis[~is_imac]
                    ms=ms[~is_imac]
                if Mtime=='z':
                    primi=cat[zi][mp_brstr][hi0]                                
                    hosti=self.subcat[zi]['halo.i'][primi]
                    #M=self.hostcat[zi]['m.200c'][hosti]

                    #Instead of having a single mass M,
                    #Get an array corresponding to the main progenitor mass at
                    #the infall times for each halo that merged at z, and
                    #name that array M:
                    inf_mpis=np.array([indices_tree(self.hostcat,zi,
                                                    inf_zi,hosti)
                                       for inf_zi in inf_zis])
                    isneg=inf_mpis<0
                    if np.sum(isneg)>0:
                        #In some cases, for some reason likely beyond human
                        #comprehensibility, subhalos are not sitting in a host 
                        #at the time of infall but
                        #are not considered hosts themselves. Get rid of those:
                        inf_zis=inf_zis[~isneg]
                        if len(inf_zis)==0:
                            continue
                        inf_is=inf_is[~isneg]
                        inf_hostis=inf_hostis[~isneg]
                        ms=ms[~isneg]
                        inf_mpis=inf_mpis[~isneg]
                    M=np.array([self.hostcat[inf_zi]['m.200c'][inf_mpi]
                                for inf_zi,inf_mpi 
                                in zip(inf_zis,inf_mpis)])
                elif Mtime=='0':
                    hosti0=self.subcat[0]['halo.i'][hi0]
                    M=self.hostcat[0]['m.200c'][hosti0]
            else:
                ms=cat[zi][mtype][merg_is]                                      
                if Mtime=='z':                                                  
                    primi=cat[zi][mp_brstr][hi0] 
                    M=cat[zi][mtype][primi]                                     
                elif Mtime=='0':                                                      
                    M0=cat[0][mtype][hi0]                                               
                    M=M0                                                            
            mzMs_add=list(ms-M)                                             
            mzMs+=mzMs_add                                                  
            ###THIS IS FOR TESTING. WE CAN REMOVE IT LATER.###
            if hi0 in self.checkdic:
                self.checkdic[hi0]+=mzMs_add
                self.mcheckdic[hi0]+=list(ms)
            else:
                self.checkdic[hi0]=mzMs_add
                self.Mcheckdic[hi0]=M
                self.mcheckdic[hi0]=list(ms)
            ##################################################
        '''
        for zi in allzis[1:]:                                               
            mergbranch=cat[zi][m_brstr]                                     
            mpbranch=cat[zi-1][mp_brstr]                                    
            if hi0 in mpbranch:                                             
                chii=mpbranch[hi0]                                          
                if chii in mergbranch:                                      
                    merg_is=mergbranch[chii]                                
                else:                                                       
                    continue                                                
            else:                                                           
                continue                                                    
            ms=cat[zi][mtype][merg_is]
            primi=cat[zi-1]['par.i'][chii]
            hosti=cat[zi]['halo.i'][primi]
            M=self.hostcat[zi]['m.200c'][hosti]
            mzMs_add=list(ms-M)
            mzMs+=mzMs_add
            if hi0 in self.checkdic:
                self.checkdic[hi0]+=mzMs_add
                self.mcheckdic[hi0]+=list(ms)
            else:
                self.checkdic[hi0]=mzMs_add
                self.Mcheckdic[hi0]=M
                self.mcheckdic[hi0]=list(ms)
        '''
    mzMs=-np.abs(mzMs)             
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    bins=np.linspace(-5,0,200)
    Ns,muaxis=ax.hist(mzMs,bins=bins,cumulative=-1,               
                      weights=np.repeat(1./float(len(hi0s)),len(mzMs)))[:2]                    
    #Ns,muaxis=ax.hist(mzMs,bins=250,cumulative=-1,               
    #                  weights=np.repeat(1./float(len(hi0s)),len(mzMs))
    #                 )[:2]                    
    plt.close()
    muaxis=(muaxis[1:]+muaxis[:-1])/2.
    return muaxis,Ns

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
    ###TESTING
    self.mrz_checkdic={}
    self.mrz_mcheckdic={}
    self.mrz_Mcheckdic={}
    ###TESTING
    for zi in zis[:-1]:
        print('\nsnapshot %i'%zi)
        hostis=elements(self.hostcat[zi][self.hmtype],
                        lim=[M-self.Mwid/2.,M+self.Mwid/2.])
        cenis=self.hostcat[zi]['cen.i'][hostis]
        randis=np.arange(cenis.size)
        print('%i halos in range.'%len(hostis))
        if not N is None:
            #random.seed(1)
            randis=random.sample(randis,N)
        cenis=cenis[randis]
        hostis=hostis[randis]
        print('evaluating %i'%len(cenis))
        Nhost+=len(cenis)

        pbar=ProgressBar()
        for ceni,hosti in pbar(list(zip(cenis,hostis))):
            M=self.hostcat[zi][self.hmtype][hosti]
            if ceni==89483:
                print('zi: %d'%zi)
                print('hosti: %d'%hosti)
            is_prog=self.subcat[zi+1]['chi.i']==ceni
            is_prog[self.subcat[zi]['par.i'][ceni]]=False
            ms=self.subcat[zi+1][self.smtype][is_prog]
            ###TESTING###
            mergis_frmeth=np.arange(len(self.subcat[zi+1][self.smtype]))\
                          [is_prog]
            ispos=ms>0.
            ms=ms[ispos]
            mergis_frmeth=mergis_frmeth[ispos]
            if (np.sum(is_prog)>0)&(np.sum(ispos)>0):
                mergis_frtree=np.sort(self.subcat[zi+1]['sub.merg.branch']\
                                      [ceni])
                mergis_frmeth_oth=np.sort(np.arange((self.subcat[zi+1]\
                                                    [self.smtype]).size)\
                                          [is_prog][ispos])
                mMs_add=ms-M
                if ceni in self.mrz_checkdic:
                    self.mrz_checkdic[ceni]+=list(mMs_add)
                    self.mrz_mcheckdic[ceni]+=list(ms)
                else:
                    self.mrz_checkdic[ceni]=list(mMs_add)
                    self.mrz_Mcheckdic[ceni]=M
                    self.mrz_mcheckdic[ceni]=list(ms)
                if not np.array_equal(mergis_frtree,mergis_frmeth):
                    print('Arrays are not equal.')
                    print('%d: central index at z=0'%ceni)
                    print('merger list from mrz method:')
                    print(mergis_frmeth)
                    print(mergis_frmeth_oth)
                    print('merger list from tree method:')
                    print(mergis_frtree)
                    #return mergis_frmeth,mergis_frtree
            ###END TESTING###
            mMs_add=ms-M
            mMs+=list(mMs_add)
    return mMs,Nhost,zbeg,zend

def sv_tpm(self):
    tstmp='{:%Y%m%d}'.format(datetime.datetime.now())
    fname='tpm_{0}'.format(tstmp)
    with h5py.File(fname,'w') as f:
        f.create_dataset('tpm',data=self.subcat[1]['gal.merg.branch'])
    return

def bld_smf_compo(main):
    smf_ar=pd.read_csv('/home/users/staudt/projects/mergers/data/smf_composites.csv')[main]
    return smf_ar
