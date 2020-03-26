import numpy as np
import random
import datetime
import h5py

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['font.family']='serif'
rcParams['mathtext.fontset'] = 'dejavuserif'

from wetzel_utils.utility.utility_catalog import indices_tree
from wetzel_utils.utility.utility_array import elements #returns index numbers of elements of an array that fall within given limts
from treepm.subhalo_io import TreepmClass
from subhalo_io_hack import TreepmClass as TreepmClass_h

from progressbar import ProgressBar

class Mergerdata:
    def __init__(self,mhal0,mwidth,catalog_kind='subhalo',n=None,zis=np.arange(35)): #inputs <mass of interest>,<acceptable width around that mass>,<size of random subsample to take>,<snapshot indices>
        self.mhal0=mhal0
        self.mwidth=mwidth
        self.snapshotindex=zis 
        self.catalog_kind=catalog_kind
        treepm=TreepmClass_h() #create a treepm instance from Wetzel's code
        #treepm=TreepmClass() #create a treepm instance from Wetzel's code
        self.cat=treepm.read(zis=self.snapshotindex,catalog_kind=self.catalog_kind) #populate the halo catalog using Wetzel's code
        self.mkind=self.cat.info['m.kind'] #Get the dictionary identifier for the halo mass list. Subhalos use 'm.max'. Halos use 'm.fof'. 
        self.zarray=np.array([self.cat.snap[snapshotnum][1] for snapshotnum in self.snapshotindex]) #look in column 1 of the snapshot table to get the redshifts
        
        self.halis={} #create an empty dictionary in which to put the index numbers of halos at the various snapshots
        self.halis[0]=elements(self.cat[0][self.mkind],
                               lim=[self.mhal0-1./2.*self.mwidth,self.mhal0+1./2.*self.mwidth]) #fill the halindexnums dictionary at snapshot 0 with the index numbers of halos that fall within the desired mass range
        #self.cat.info['m.kind'] is there because the halo masses are stored either as 'm.max' if catalog_kind='subhalo' or 'm.fof' if catalog_kind='halo'
        if(self.catalog_kind=='subhalo'):
            is_cen=self.cat[0]['ilk'][self.halis[0]]==1
            self.halis[0]=self.halis[0][is_cen]

        random.seed(1) 
        if n is None:
            n=len(self.halis[0])
            print('\n{0:0.0f} halos in range. Evaluating all.'.
                  format(len(self.halis[0])))
        else:
            print('\n{0:0.0f} halos in range. Evaluating a sample of {1:0.0f}.'.
                  format(len(self.halis[0]),n))
            self.halis[0]=random.sample(self.halis[0],n) #Take a radom sample if the sample size is specified.
    
    def populatemap(self,mkind,zis,n=None):
        t0=self.cat.snap[0][2]
        self.lbt=t0-[self.cat.snap[snap][2] for snap in zis]
        indices_tree(self.cat,zis[0],zis[-1])
        for snapshotnum in zis[1:]: #fill the halis dictionary at all snapshots after snapshot 0 with the index numbers of the parents of the snapshot 0 halos
            self.halis[snapshotnum]=indices_tree(self.cat,0,snapshotnum,self.halis[0])
        massmap=np.zeros([len(self.halis[0]),len(zis)])
        pbar=ProgressBar()
        for snapshotnum in pbar(zis):    
            for massmaprow, halindexnum in enumerate(self.halis[snapshotnum]): #massmap rows are halos. massmap columns are snapshots (each at a specific redshift). This loop populates the massmap by putting each halo's mass at the given redshift in the corresponding row and column.
                if halindexnum<0: #If the parent ID is negative, the halo doesn't exist yet, so keep its mass at 0.
                    continue
                else:
                    massmap[massmaprow,snapshotnum]=self.cat[snapshotnum][mkind][halindexnum]
        #avgratio=np.average(massmap[:,snap]/massmap[:,0] for snap in zis)
        self.avgratios=[np.average(10.**massmap[:,snap]/10.**massmap[:,0]) for snap in zis]
        self.avg_gal_m=[np.average(massmap[:,snap]) for snap in zis]
        self.massmap=massmap
    
    def gendistributionnew(self):
        #self.mratios=np.array([])
        #self.ms=np.array([])
        mratios=[]
        ms=[]

        pbar=ProgressBar()
        for finhali in pbar(self.halis[0]):
            m0=self.cat[0][self.mkind][finhali]
            mainprogi_prevsnap=finhali
            for snap in self.snapshotindex[1:]:
                mainprogi_snap=self.cat[snap-1]['par.i'][mainprogi_prevsnap] #Parent index of main progenitor at snap-1 is the main progenitor index at snap.
                if mainprogi_snap==-1:
                    break #If the main progenitor index reads -1 at snap, it means the main progenitor formed at snap-1. Therefore break the loop and move on to the next finhali. 
                progenis_snap=np.where(self.cat[snap]['chi.i']==mainprogi_prevsnap)[0] #Progenitors of the main progenitor are identified by finding at snap halos whose child index is the main progenitor index at snap-1 
                mask_mainprog=progenis_snap!=mainprogi_snap
                mask_zeromass=self.cat[snap][self.mkind][progenis_snap]!=0.
                mask=mask_mainprog*mask_zeromass
                progenis_snap=progenis_snap[mask] #Remove the main progenitor and zero-mass halos from calculations
                ms_add=self.cat[snap][self.mkind][progenis_snap]
                #print type(ms_add)
                mratios_add=10.**ms_add/10.**m0 
                ms+=list(ms_add)
                mratios+=list(mratios_add)
                #self.ms=np.append(self.ms,ms_add)
                #self.mratios=np.append(self.mratios,mratios_add)
                mainprogi_prevsnap=mainprogi_snap #Prepare the loop for the next iteration by moving forward mainprogi_prevsnap.

        self.ms=np.array(ms)
        self.mratios=np.array(mratios)
        
        num=50.
        lowerbound=min(self.mratios)
        upperbound=max(self.mratios)
        bins=np.logspace(np.log10(lowerbound),np.log10(upperbound),num)
        midbins=(bins[:-1]+bins[1:])/2
        n,binsout,patches=plt.hist(self.mratios,bins)
        plt.clf()
        n=n/len(self.halis[0])
       
        ngeqmratio=np.array([n[i:].sum() for i in np.arange(len(n))]) #n(>=m/m0) --- avg number of mergers that occured with ratios greater than or equal to the given m/m0
        fig=plt.figure(figsize=(8,6))
        ax=fig.add_subplot(1,1,1)
        #ax.plot(midbins,ngeqmratio,'s',mec='b',mfc='None')
        #ax.semilogy(midbins,ngeqmratio,'s',mec='b',mfc='None')
        ax.loglog(midbins,ngeqmratio,'s',mec='b',mfc='None')
        ax.set_title('$M_0=10^{{{0:0.0f}}}$ halo mergers'.
                     format(self.mhal0),fontsize=16)
        ax.set_xlabel('$m/M_0$',fontsize=14)
        ax.set_ylabel('$N\left(\geq\\frac{m}{M_0}\\right)$',fontsize=14)
        ax.tick_params(axis='both',labelsize=12)
        plt.show()

    def gendistribution(self):
        self.fillfamilytrees()
        self.mratios=np.array([])
        self.ms=np.array([])
        for finhali in self.familytrees.keys():
            m0=self.cat[0][self.mkind][finhali] #mass of the halo at z=0
            for snap in self.snapshotindex[1:]:
                if self.familytrees[finhali].get(snap) is None:
                    break 
                else:
                    #print 'snapshot: %d' %snap
                    for hali in self.familytrees[finhali][snap]:
                        #print 'hali: %d' %hali
                        chii=self.cat[snap]['chi.i'][hali]
                        #print 'chii: %d' %chii
                        mainprogi=self.cat[snap-1]['par.i'][chii]
                        '''
                        print('main progenitor index of child: {0}'.
                              format(mainprogi))
                        '''
                        m=self.cat[snap][self.mkind][hali]
                        if (hali!=mainprogi) and (m!=0.): #if the given halo is not the main progenitor and the given halo's mass is non-zero
                            self.ms=np.append(self.ms,m)
                            mratioadd=10.**m/10.**m0
                            self.mratios=np.append(self.mratios,mratioadd)
                            #print('m/m0: {0:0.4f}'.format(mratioadd))
        num=50.
        lowerbound=min(self.mratios)
        upperbound=max(self.mratios)
        bins=np.logspace(np.log10(lowerbound),np.log10(upperbound),num)
        midbins=(bins[:-1]+bins[1:])/2
        n,binsout,patches=plt.hist(self.mratios,bins)
        plt.clf()
        n=n/len(self.halis[0])
       
        ngeqmratio=np.array([n[i:].sum() for i in np.arange(len(n))]) #n(>=m/m0) --- avg number of mergers that occured with ratios greater than or equal to the given m/m0
        fig=plt.figure(figsize=(8,6))
        ax=fig.add_subplot(1,1,1)
        #ax.plot(midbins,ngeqmratio,'s',mec='b',mfc='None')
        #ax.semilogy(midbins,ngeqmratio,'s',mec='b',mfc='None')
        ax.loglog(midbins,ngeqmratio,'s',mec='b',mfc='None')
        ax.set_title('$M_0=10^{{{0:0.0f}}}$ halo mergers'.
                     format(self.mhal0),fontsize=16)
        ax.set_xlabel('$m/M_0$',fontsize=14)
        ax.set_ylabel('$N\left(\geq\\frac{m}{M_0}\\right)$',fontsize=14)
        ax.tick_params(axis='both',labelsize=12)
        plt.show()

    def genaccdist(self): #Generate unqualified accretion distribution
        m1array=np.zeros([len(self.halis[0]),len(self.snapshotindex)])*np.nan #initialize the array containing the accreted  masses of each merger
        m1array[:,1:]=self.massmap[:,:-1]-self.massmap[:,1:]
        m0array=self.massmap[:,0]
        mergerratioarray=[10.**m1array[rownum,snapshotnum]/10.**m0array[rownum]
                         for rownum in np.arange(self.massmap[:,0].size)
                         for snapshotnum in self.snapshotindex[1:]]
        return m1array, m0array, mergerratioarray
    
    def fillfamilytrees(self): #Fill a familytrees dictionary with halo indexes. Also fill an identically shaped mtree dictionary with the masses corresponding to the halos in familytrees.
        self.familytrees={} #Initialize familytrees as a dictionary so we can lookup family trees by final halo index number. (Note that index number is not halo ID. It is the position of the halo in the catalog.)
        #The keys for familytrees will be index numbers of the halos in the objective mass range at z=0 i.e. halindexnumbs[0]. The values for any given key in self.familytrees will be a second-tier dictionary. The keys for that second tier dictionary will be snapshots. The values will be arrays containing the index numbers of the halos that were progenitors at the given snapshot for the given redshift-zero halo. 
        self.mtree={} #Initialize the mass dictionary
        
        def fillfinhali(): #This function will be used later to fill finhali's tree once the snapshot 1 branch has been filled.            
            for snap in self.snapshotindex[2:]: 
                #print('snap: %d' % snap)
                halis_snap=np.array([],dtype=int) #Initialize an empty array into which we will append progenitor halo indices for the given snapshot every time we evaluate each halo in the previous snapshot of the tree.
                ms_snap=np.array([])
                for hali_prevsnap in self.familytrees[finhali][snap-1]:
                    #print'hali_prevsnap: %d' %hali_prevsnap
                    addhalis_snap=np.where(self.cat[snap]
                                           ['chi.i']
                                           ==hali_prevsnap)[0] #indices of the halos at snap that are the parents of the parent at snap-1 currently being evaluated
                    mask=self.cat[snap]['halo.i'][addhalis_snap]!=-1 #When a halo is newly formed at snap-1, it has a parent at snap with halo ID -1. This and the next line get rid of those 'halo.i'=-1 halos at snap.
                    addhalis_snap=addhalis_snap[mask]
                    addms_snap=self.cat[snap][self.mkind][addhalis_snap]
                    halis_snap=np.append(halis_snap,addhalis_snap) #Indices of every halo at snap that was a progenitor for the halos that were found to be progenitors at snap-1.
                    ms_snap=np.append(ms_snap,addms_snap)
                    if len(halis_snap)==0:
                        return #If the mask removes all elements of addhalis_snap, there are no progenitors at snap, so break the whole loop and go to the next finhali.
                    else:
                        self.mtree[finhali][snap]=ms_snap
                        self.familytrees[finhali][snap]=halis_snap #Put halis_snap into familytrees according to finalhaloindexnum and snap.                      
        
        print('generating family trees for {0:0.0f} main progenitors'.format(len(self.halis[0])))
        pbar=ProgressBar()
        for finhali in pbar(self.halis[0]): #for every halo at redshift 0
            #print('final halo index: %d' % finhali)
            #The next two lines initialize dictionaries within familytrees and ms for the given finhali, the keys for which will be snapshot indices.
            self.familytrees[finhali]={}
            self.mtree[finhali]={}
            paris_one=np.where(self.cat[1]['chi.i']==[finhali])[0] #Anywhere the [child index @ snapshot 1] == [the given halo index @ snapshot 0], that halo @ snapshot 1 is a progenitor of the halo @ snapshot 0. 
            if len(paris_one)==0.:
                continue
            self.familytrees[finhali][1]=paris_one #Fill with all the halos at snapshot 1 that were progentiors for finalhalo.
            self.mtree[finhali][1]=self.cat[1][self.mkind][paris_one] #masses of the parents at snapshot one
            fillfinhali()

    def galmassf(self,snap=0): 
        #Get galaxy  mass number density at given snapshot
        gmtype='m.star'
        if self.catalog_kind=='subhalo':
            allhalids=self.cat[snap]['halo.i'] #Get every halo ID at snap
            mask2=allhalids!=-1
        else:
            mask2=self.cat[snap]['m.fof']>0.
        allms=self.cat[snap][gmtype] #Get all masses at snap
        mask1=allms!=0.
        mask=mask1*mask2 #Get halo index for every halo with nonzero mass and ID not equal to -1.
        ms=allms[mask] #Get the corresponding masses
        num=50. #Number of bins
        lower=np.min(ms) #Lower histogram bound
        upper=np.max(ms) #Upper histogram bound
        bins=np.linspace(lower,upper,num) #Mass bins in log
        binw=np.average(bins[1:]-bins[:-1]) #Bin width in dex
        midbins=(bins[:-1]+bins[1:])/2 
        N,binsout,patches=plt.hist(ms,bins)
        plt.clf()
        n=N/(250.)**3./binw #number density i.e. number per volume per bin size
        #n=N/(250./0.7)**3./binw #number density i.e. number per volume per bin size
        return midbins,n
    
    def halmassf(self,snap=0): #Get halo mass number density at given snapshot
        if self.catalog_kind=='subhalo':
            allhalids=self.cat[snap]['halo.i'] #Get every halo ID at snap
            mask2=allhalids!=-1
        else:
            mask2=self.cat[snap]['m.fof']>0.
        allms=self.cat[snap][self.mkind] #Get all masses at snap
        mask1=allms!=0.
        mask=mask1*mask2 #Get halo index for every halo with nonzero mass and ID not equal to -1.
        ms=allms[mask] #Get the corresponding masses
        num=50. #Number of bins
        lower=np.min(ms) #Lower histogram bound
        upper=np.max(ms) #Upper histogram bound
        bins=np.linspace(lower,upper,num) #Mass bins in log
        binw=np.average(bins[1:]-bins[:-1]) #Bin width in dex
        midbins=(bins[:-1]+bins[1:])/2 
        N,binsout,patches=plt.hist(ms,bins)
        plt.clf()
        n=N/(250.)**3./binw #number density i.e. number per volume per bin size
        #n=N/(250./0.7)**3./binw #number density i.e. number per volume per bin size
        return midbins,n
         
    def tracker(self,zis=np.arange(35)):
        mfrackey='m.frac.min'
        thresh=0.1
        his=np.arange(len(self.cat[self.snapshotindex[-1]]['halo.i']))
        
        '''
        #FOR TESTING####
        n=int(1e5)
        his=np.array(random.sample(his,n)) #Take a radom sample for testing
        #FOR TESTING####
        '''
        
        print'starting with %i halos' %len(his)
        m_M0s=[]
        M0s=[]
        iscen=self.cat[self.snapshotindex[-1]]['ilk'][his]==1 #specifying his here in the first step is usually redundant, but I'm including it to cover special case where we start with less than the full catalog of his
        #mfracsprev=np.repeat(1.,sum(~iscen))
        mfracsprev=np.repeat(1.,len(his))
        pbar=ProgressBar()
        #for zi in pbar(self.snapshotindex[::-1]): 
        for zi in zis[::-1]: 
            print'\nsnapshot %i' %zi
            iscen=self.cat[zi]['ilk'][his]==1
            #print'%i non-centrals'%sum(~iscen)
            #his_nc=his[~iscen] #Do I really want to be removing centrals here? I think so, but verify this.
            mfracs=self.cat[zi][mfrackey][his]
            
            smaller=mfracs<mfracsprev
            destd=mfracs<=thresh #qualifies as destroyed
            isev=smaller & destd & ~iscen #qualifies as an event
            his_ev=his[isev] 
            #print'%i newly destroyed non-central halos' %len(his_ev)
            mmaxs_ev=self.cat[zi]['m.max'][his_ev]
            mfracs_ev=mfracs[isev]
            ms_inst=(mmaxs_ev*mfracs_ev) #The program currently does not use this data
            
            chis=self.cat[zi]['halo.i'][his_ev] #central halo indices---Don't confuse with chiis (child halo indices)
            #fchis=np.array(indices_tree(self.cat,zi,0,chis)) #z=0 indices of central halos
            fchis=indices_tree(self.cat,zi,0,chis) #z=0 indices of central halos
            has_fchi=fchis>=0 #if indices_tree returns a negative for the final central halo index, there is no central at z=0.
            fchis=fchis[has_fchi]
            #print'%i final halos for %i newly destroyed halos' %(len(fchis),len(his_ev))
            #M0s=self.cat[0]['m.max'][fchis]
            M0s_add=self.cat[0]['m.max'][fchis]
            #inmbin=(M0s<self.mhal0+self.mwidth/2.)&(M0s>self.mhal0-self.mwidth/2)
            #ms=ms[has_fchi][inmbin] #Cut the ms array down to elements whose central exists at z=0. Then pull elements of that array (the shape of which matches fchis) whose centrals are in the mbin. 
            ms=mmaxs_ev[has_fchi].flatten() #Cut the mmaxs_ev array down to elements whose central exists at z=0. I'm not going to eliminate elements outside the mbin, because I want to only run this once, and it's probably a better idea to just record m_M0 along with what the M0 is and then narrow down the m_M0's later.
            #The reason I need to flatten the array: If mmaxs_ev has only one element, the mmaxs_ev[has_fchi] array returns a nested array, apparently because it evaluates both the single element and the dtype=float32 descriptor. I have no idea why it does this, but flattening the resulting masked array is a work around. 
            #M0s=M0s[inmbin] #M0s was made from fchis, so we don't need to cut anything before taking inmbin mask.
            m_M0s_add=ms-M0s_add
            m_M0s_add=np.minimum(m_M0s_add,-m_M0s_add)
            #print'new m_M0 ratios:'
            #print m_M0s_add
            '''
            print"from m.max's:"
            print ms
            print"from M0's:"
            print M0s_add
            '''
            M0s+=list(M0s_add)
            m_M0s+=list(m_M0s_add)

            #prep for next run:
            his=his[~destd] #Remove destroyed subhalos from next evaluation.
            mfracsprev=mfracs[~destd] #Remove mfracs for destroyed subhalos 
            #print 'removing %i destroyed halos. %i halos remaining'%(sum(destd),len(his))
            chiis=self.cat[zi]['chi.i'][his] #child indices
            his=chiis

        self.M0s=np.array(M0s)
        self.m_M0s=np.array(m_M0s)
        #print''
        #print m_M0s

	timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
	#filename='./dat/true_m_M0_{0:0.0f}_{2}_{1}.h5'.format(mmid,timestmp,self.mkind)
	filename='./dat/true_m_M0_{1}_{0}.h5'.format(timestmp,self.mkind)
	f = h5py.File(filename, 'w')
	f.create_dataset('m_M0s', data=self.m_M0s)
        f.create_dataset('M0s',data=self.M0s)
        f.create_dataset('N_fc',data=len(self.halis[0]))
	f.close()

    def trackernew(self,zis=np.arange(35)): #trackernew is an attempt to add m/M_z functionality
        mfrackey='m.frac.min'
        thresh=0.1
        his=np.arange(len(self.cat[self.snapshotindex[-1]]['halo.i']))
        
        '''
        #FOR TESTING####
        n=int(1e5)
        his=np.array(random.sample(his,n)) #Take a radom sample for testing
        #FOR TESTING####
        '''
        
        print'starting with %i halos' %len(his)
        m_M0s=[]
        M0s=[]
        m_Ms=[]
        Ms=[]
        iscen=self.cat[self.snapshotindex[-1]]['ilk'][his]==1 #specifying his here in the first step is usually redundant, but I'm including it to cover special case where we start with less than the full catalog of his
        #mfracsprev=np.repeat(1.,sum(~iscen))
        mfracsprev=np.repeat(1.,len(his))
        pbar=ProgressBar()
        #for zi in pbar(self.snapshotindex[::-1]): 
        for zi in zis[::-1]: 
            print'\nsnapshot %i' %zi
            iscen=self.cat[zi]['ilk'][his]==1
            #print'%i non-centrals'%sum(~iscen)
            #his_nc=his[~iscen] #Do I really want to be removing centrals here? I think so, but verify this.
            mfracs=self.cat[zi][mfrackey][his]
            
            smaller=mfracs<mfracsprev
            destd=mfracs<=thresh #qualifies as destroyed
            isev=smaller & destd & ~iscen #qualifies as an event
            his_ev=his[isev] 
            #print'%i newly destroyed non-central halos' %len(his_ev)
            mmaxs_ev=self.cat[zi]['m.max'][his_ev]
            mfracs_ev=mfracs[isev]
            ms_inst=(mmaxs_ev*mfracs_ev) #The program currently does not use this data
            
            chis=self.cat[zi]['halo.i'][his_ev] #central halo indices---Don't confuse with chiis (child halo indices)
            #fchis=np.array(indices_tree(self.cat,zi,0,chis)) #z=0 indices of central halos
            fchis=indices_tree(self.cat,zi,0,chis) #z=0 indices of central halos
            has_fchi=fchis>=0 #if indices_tree returns a negative for the final central halo index, there is no central at z=0.
            ##fchis=fchis[has_fchi]
            #print'%i final halos for %i newly destroyed halos' %(len(fchis),len(his_ev))
            #M0s=self.cat[0]['m.max'][fchis]
            M0s_add=self.cat[0]['m.max'][fchis]
            Ms_add=self.cat[zi]['halo.m'][his_ev]
            #inmbin=(M0s<self.mhal0+self.mwidth/2.)&(M0s>self.mhal0-self.mwidth/2)
            #ms=ms[has_fchi][inmbin] #Cut the ms array down to elements whose central exists at z=0. Then pull elements of that array (the shape of which matches fchis) whose centrals are in the mbin. 
            ms=mmaxs_ev##[has_fchi].flatten() #Cut the mmaxs_ev array down to elements whose central exists at z=0. I'm not going to eliminate elements outside the mbin, because I want to only run this once, and it's probably a better idea to just record m_M0 along with what the M0 is and then narrow down the m_M0's later.
            #The reason I need to flatten the array: If mmaxs_ev has only one element, the mmaxs_ev[has_fchi] array returns a nested array, apparently because it evaluates both the single element and the dtype=float32 descriptor. I have no idea why it does this, but flattening the resulting masked array is a work around. 
            #M0s=M0s[inmbin] #M0s was made from fchis, so we don't need to cut anything before taking inmbin mask.
            m_Ms_add=ms-Ms
            m_Ms_add=np.minimum(m_M_add,-m_M_add)
            m_M0s_add=ms-M0s_add
            m_M0s_add=np.minimum(m_M0s_add,-m_M0s_add)
            #print'new m_M0 ratios:'
            #print m_M0s_add
            '''
            print"from m.max's:"
            print ms
            print"from M0's:"
            print M0s_add
            '''
            M0s+=list(M0s_add)
            m_M0s+=list(m_M0s_add)
            m_Ms+=list(m_M_add)
            Ms+=list(Ms_add)

            #prep for next run:
            his=his[~destd] #Remove destroyed subhalos from next evaluation.
            mfracsprev=mfracs[~destd] #Remove mfracs for destroyed subhalos 
            #print 'removing %i destroyed halos. %i halos remaining'%(sum(destd),len(his))
            chiis=self.cat[zi]['chi.i'][his] #child indices
            his=chiis

        self.M0s=np.array(M0s)
        self.m_M0s=np.array(m_M0s)
        self.m_Ms=np.array(m_Ms)
        self.Ms=np.array(Ms)
        #print''
        #print m_M0s

	timestmp='{:%Y%m%d%H}'.format(datetime.datetime.now())
	#filename='./dat/true_m_M0_{0:0.0f}_{2}_{1}.h5'.format(mmid,timestmp,self.mkind)
	filename='./dat/true_m_M0_{1}_{0}.h5'.format(timestmp,self.mkind)
	f = h5py.File(filename, 'w')
	f.create_dataset('m_M0s', data=self.m_M0s)
        f.create_dataset('M0s',data=self.M0s)
        f.create_dataset('Ms',data=self.Ms)
        f.create_dataset('m_Ms',data=self.m_Ms)
        f.create_dataset('N_fc',data=len(self.halis[0]))
	f.close()

    def tracker_with_opt(self,zis=np.arange(35),kind='merger'):
        mfrackey='m.frac.min'
        thresh=0.1
        his=np.arange(len(self.cat[self.snapshotindex[-1]]['halo.i']))
        
        '''
        #FOR TESTING####
        n=int(1e5)
        his=np.array(random.sample(his,n)) #Take a radom sample for testing
        #FOR TESTING####
        '''
        
        print'starting with %i halos' %len(his)
        m_M0s=[]
        M0s=[]
        
        #This next line might be redundant, because it may have been originally intended to feed the line after itself, and that line after is now commented.
        iscen=self.cat[self.snapshotindex[-1]]['ilk'][his]==1 #specifying his here in the first step is usually redundant, but I'm including it to cover special case where we start with less than the full catalog of his
        #mfracsprev=np.repeat(1.,sum(~iscen))
        
        #mfracsprev=np.repeat(1.,len(his))
        mfracsprev=np.zeros(len(his))
        
        pbar=ProgressBar()
        #for zi in pbar(self.snapshotindex[::-1]): 
        for zi in zis[::-1]: 
            print'\nsnapshot %i' %zi
            iscen=self.cat[zi]['ilk'][his]==1
            #print'%i non-centrals'%sum(~iscen)
            #his_nc=his[~iscen] #Do I really want to be removing centrals here? I think so, but verify this.
            mfracs=self.cat[zi][mfrackey][his]
            
            smaller=mfracs<mfracsprev
            destd=mfracs<=thresh #qualifies as destroyed
            isev=smaller & ~iscen #qualifies as an event
            if kind=='merger':
                isev=smaller & destd & ~iscen #qualifies as an event
            else:
                isev=smaller & ~iscen
            his_ev=his[isev] 
            #print'%i newly destroyed non-central halos' %len(his_ev)
            mmaxs_ev=self.cat[zi]['m.max'][his_ev]
            mfracs_ev=mfracs[isev]
            ms_inst=(mmaxs_ev*mfracs_ev) #The program currently does not use this data
            
            chis=self.cat[zi]['halo.i'][his_ev] #central halo indices---Don't confuse with chiis (child halo indices)
            #fchis=np.array(indices_tree(self.cat,zi,0,chis)) #z=0 indices of central halos
            fchis=indices_tree(self.cat,zi,0,chis) #z=0 indices of central halos
            has_fchi=fchis>=0 #if indices_tree returns a negative for the final central halo index, there is no central at z=0.
            fchis=fchis[has_fchi]
            #print'%i final halos for %i newly destroyed halos' %(len(fchis),len(his_ev))
            #M0s=self.cat[0]['m.max'][fchis]
            M0s_add=self.cat[0]['m.max'][fchis]
            #inmbin=(M0s<self.mhal0+self.mwidth/2.)&(M0s>self.mhal0-self.mwidth/2)
            #ms=ms[has_fchi][inmbin] #Cut the ms array down to elements whose central exists at z=0. Then pull elements of that array (the shape of which matches fchis) whose centrals are in the mbin. 
            ms=mmaxs_ev[has_fchi].flatten() #Cut the mmaxs_ev array down to elements whose central exists at z=0. I'm not going to eliminate elements outside the mbin, because I want to only run this once, and it's probably a better idea to just record m_M0 along with what the M0 is and then narrow down the m_M0's later.
            #The reason I need to flatten the array: If mmaxs_ev has only one element, the mmaxs_ev[has_fchi] array returns a nested array, apparently because it evaluates both the single element and the dtype=float32 descriptor. I have no idea why it does this, but flattening the resulting masked array is a work around. 
            #M0s=M0s[inmbin] #M0s was made from fchis, so we don't need to cut anything before taking inmbin mask.
            m_M0s_add=ms-M0s_add
            m_M0s_add=np.minimum(m_M0s_add,-m_M0s_add)
            #print'new m_M0 ratios:'
            #print m_M0s_add
            '''
            print"from m.max's:"
            print ms
            print"from M0's:"
            print M0s_add
            '''
            M0s+=list(M0s_add)
            m_M0s+=list(m_M0s_add)

            #prep for next run:
            his=his[~isev] #Remove affected subhalos from next evaluation.
            mfracsprev=mfracs[~isev] #Remove mfracs for affected subhalos 
            #print 'removing %i destroyed halos. %i halos remaining'%(sum(destd),len(his))
            chiis=self.cat[zi]['chi.i'][his] #child indices
            his=chiis

        self.M0s=np.array(M0s)
        self.m_M0s=np.array(m_M0s)
        #print''
        #print m_M0s

	timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
	#filename='./dat/true_m_M0_{0:0.0f}_{2}_{1}.h5'.format(mmid,timestmp,self.mkind)
	filename='./dat/true_m_M0_{1}_{0}.h5'.format(timestmp,self.mkind)
	f = h5py.File(filename, 'w')
	f.create_dataset('m_M0s', data=self.m_M0s)
        f.create_dataset('M0s',data=self.M0s)
        f.create_dataset('N_fc',data=len(self.halis[0]))
	f.close()

    def plotnewtracker():
        timestmp=20181031
        fname='/home/users/staudt/projects/mergers/dat/true_m_M0_m.max_{}.h5'.format(timestmp)
        f=h5py.File(fname,'r')
        M0s_zis=np.array(f['M0s'])
        m_M0s_zis=np.array(f['m_M0s'])
