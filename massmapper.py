import numpy as np
import random

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['font.family']='serif'
rcParams['mathtext.fontset'] = 'dejavuserif'

from wetzel_utils.utility.utility_catalog import indices_tree
from wetzel_utils.utility.utility_array import elements #returns index numbers of elements of an array that fall within given limts
from treepm.subhalo_io import TreepmClass

from progressbar import ProgressBar

class Mergerdata:
    def __init__(self,mhal0,mwidth,catalog_kind='subhalo',n=None,zis=np.arange(35)): #inputs <mass of interest>,<acceptable width around that mass>,<size of random subsample to take>,<snapshot indices>
        self.mhal0=mhal0
        self.mwidth=mwidth
        self.snapshotindex=zis 
        self.catalog_kind=catalog_kind
        treepm=TreepmClass() #create a treepm instance from Wetzel's code
        self.cat=treepm.read(zis=self.snapshotindex,catalog_kind=self.catalog_kind) #populate the halo catalog using Wetzel's code
        self.mkind=self.cat.info['m.kind'] #Get the dictionary identifier for the halo mass list. Subhalos use 'm.max'. Halos use 'm.fof'. 
        self.zarray=np.array([self.cat.snap[snapshotnum][1] for snapshotnum in self.snapshotindex]) #look in column 1 of the snapshot table to get the redshifts
        
        self.halis={} #create an empty dictionary in which to put the index numbers of halos at the various snapshots
        self.halis[0]=elements(self.cat[0][self.mkind],
                               lim=[self.mhal0-1./2.*self.mwidth,self.mhal0+1./2.*self.mwidth]) #fill the halindexnums dictionary at snapshot 0 with the index numbers of halos that fall within the desired mass range
        #self.cat.info['m.kind'] is there because the halo masses are stored either as 'm.max' if catalog_kind='subhalo' or 'm.fof' if catalog_kind='halo'
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
    
    def populatemap(self,n=None):
        t0=self.cat.snap[0][2]
        self.lbt=t0-[self.cat.snap[snap][2] for snap in self.snapshotindex]
        indices_tree(self.cat,self.snapshotindex[0],self.snapshotindex[-1])
        for snapshotnum in self.snapshotindex[1:]: #fill the halis dictionary at all snapshots after snapshot 0 with the index numbers of the parents of the snapshot 0 halos
            self.halis[snapshotnum]=indices_tree(self.cat,0,snapshotnum,self.halis[0])
        massmap=np.zeros([len(self.halis[0]),len(self.snapshotindex)])
        pbar=ProgressBar()
        for snapshotnum in pbar(self.snapshotindex):    
            for massmaprow, halindexnum in enumerate(self.halis[snapshotnum]): #massmap rows are halos. massmap columns are snapshots (each at a specific redshift). This loop populates the massmap by putting each halo's mass at the given redshift in the corresponding row and column.
                if halindexnum<0: #If the parent ID is negative, the halo doesn't exist yet, so keep its mass at 0.
                    continue
                else:
                    massmap[massmaprow,snapshotnum]=self.cat[snapshotnum][self.mkind][halindexnum]
        #avgratio=np.average(massmap[:,snap]/massmap[:,0] for snap in self.snapshotindex)
        self.avgratios=[np.average(10.**massmap[:,snap]/10.**massmap[:,0]) for snap in self.snapshotindex]
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

    def halmassf(self,snap=0): #Get halo mass number density at given snapshot
        allhalids=self.cat[snap]['halo.i'] #Get every halo ID at snap
        allms=self.cat[snap][self.mkind] #Get all masses at snap
        mask1=allms!=0.
        mask2=allhalids!=-1
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
        n=N/250**3./binw #number density i.e. number per volume per bin size
        return midbins,n
         
    def tracker(self):
        mfrackey='mfrackey'
        thresh=0.1
        halis=len(self.cat[self.snapshotindex[-1]]['halo.i']):
        m_m0s=[]
            for zi in self.snapshotindex[::-1]     
                iscen=self.cat[zi]['ilk'][halis]==1
                halis=halis[~iscen]
                mfracs=self.cat[zi][mfrackey][halis] 
                mmaxs=self.cat[zi]['m.max'][halis]
                
                diddrop=mfracs<mfracsprev
                destd=mfracs<=thresh #qualifies as destroyed
                ms=(mmaxs*mfracs)[diddrop & destd] 
                mcs=self.cat[zi]
                ms+=msadd

                #prep for next run:
                chiis=self.cat[zi]['chi.i'][halis]
                halis=chiis
                mfracsprev=mfracs
                
