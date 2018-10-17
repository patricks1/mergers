#!/bin/python
'''
'''
import matplotlib as mpl
mpl.use('Agg')

import h5py 
import numpy as np
import datetime
import random

from collections import Counter
try: 
    from treepm import subhalo_io
    from utilities import utility as wUT
except ImportError: 
    print'import fail'
    #pass


def merger_track(mmid, dlogm=0.5,mkind='halo.m'): 
    ''' track the accreted halo masses 

    snap5   * * m3
            |/ 
    snap4   * 
            |
    snap3   * * m2
            |/ 
    snap2   * 
            | 
    snap1   * * m1 
            |/ 
    snap0   * M0 
    
    output [m1/M0, m2/M0, m3/M0...]  to file 

    '''
    snapshots = np.arange(35) 

    # read in treepm subhalos 
    sub = subhalo_io.Treepm.read('subhalo', zis=snapshots)
    nsub0 = len(sub[0]['ilk'])

    # central at snapshot 0 
    is_cen = (sub[0]['ilk'] == 1)
    print('%i out of %i subhalos are centrals' % (np.sum(is_cen), len(is_cen)))
    # within [mmid - 0.5 dlogm, mmid + 0.5 dlogm] at snapshot 0 
    in_mbin = (sub[0][mkind] > mmid - 0.5*dlogm) & (sub[0][mkind] < mmid + 0.5*dlogm)

    ihalos = np.arange(nsub0)[is_cen & in_mbin]
    #random.seed(1)
    #ihalos=np.array(random.sample(ihalos,50))

    print('%i subhalos in range' % len(ihalos))
    m_M0 = [] 
    for isnap in range(snapshots[-1]): 
        # indices of the parent with the highest M_max (i.e. primaries) 
        i_primaries = sub[isnap]['par.i'][ihalos]
        # only keep halos that have parents
        # at snapshot=0 ~40 do not have parents
        has_primary = (i_primaries >= 0) 
        ihalos = ihalos[has_primary]
        i_primaries = i_primaries[has_primary]
        
        # make sure indices match across snapshot i and i+1 
        assert np.sum(ihalos - sub[isnap+1]['chi.i'][i_primaries]) == 0
        
        # identify halos with more than just a primary parent
        # the few lines below is just to speed things up a bit 
        counter = Counter(list(sub[isnap+1]['chi.i']))
        ip1_children_dupl = np.array([ii for ii, count in counter.items() if count >1])
        ihs = np.intersect1d(ihalos, np.array(ip1_children_dupl))
    
        for ih in ihs: 
            # loop through halos with parents and store non-primary halo masses 
            is_parent = np.where(sub[isnap+1]['chi.i'] == ih)[0]
            i_primary = sub[isnap]['par.i'][ih]
            notprimary = (is_parent != i_primary) 
           
            ifh=wUT.catalog.indices_tree(sub, isnap, 0,ih)
            M0=sub[0][mkind][ifh]
            m_M0 += list(sub[isnap+1][mkind][is_parent[notprimary]] - M0)
            #m_M0 += list(sub[isnap+1][mkind][is_parent[notprimary]] - sub[isnap][mkind][ih])
    
        print('snapshot %i, %i m_M0s' % (isnap, len(m_M0)))
        # keep searching through the primaries 
        ihalos = i_primaries
    print('%i halos have primaries in the last snapshot' % len(ihalos)) 
    m_M0 = np.array(m_M0)
    
    timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
    filename='./dat/m_M0_{0:0.0f}_{2}_{1}.h5'.format(mmid,timestmp,mkind)
    f = h5py.File(filename, 'w') 
    f.create_dataset('m_M0', data=m_M0) 
    f.close() 
    return None 


def merger_track_throughout(mmid, dlogm=0.5,mkind='halo.m'): 
    ''' track the accreted halo masses 

    snap5   * * m3
            |/ 
    snap4   * 
            |
    snap3   * * m2
            |/ 
    snap2   * 
            | 
    snap1   * * m1 
            |/ 
    snap0   * M0 
    
    output [m1/M0, m2/M0, m3/M0...]  to file 

    '''
    snapshots = np.arange(35) 

    # read in treepm subhalos 
    sub = subhalo_io.Treepm.read('subhalo', zis=snapshots)
    nsub0 = len(sub[0]['ilk'])

    # central at snapshot 0 
    is_cen = (sub[0]['ilk'] == 1)
    print('%i out of %i subhalos are centrals' % (np.sum(is_cen), nsub0))
    # within [mmid - 0.5 dlogm, mmid + 0.5 dlogm] at snapshot 0 
    in_mbin = (sub[0][mkind] > mmid - 0.5*dlogm) & (sub[0][mkind] < mmid + 0.5*dlogm)
    #print('%i out of %i subhalos are centrals' % (np.sum(is_cen & in_mbin), nsub0))
    # has prognitor until the last snapshot 
    prog_final = wUT.catalog.indices_tree(sub, 0, snapshots[-1], np.arange(nsub0))
    #prog_final = wUT.utilities_catalog.indices_tree(sub, 0, snapshot[-1], np.arange(nsub0))
    has_prog = (prog_final >= 0)
    
    print('%i in range' % np.sum(is_cen & in_mbin))
    ihalos = np.arange(nsub0)[is_cen & in_mbin & has_prog]
    print('%i out of %i subhalos are tracked throughout the snapshots' % (len(ihalos), nsub0))
    m_M0 = [] 
    for isnap in range(snapshots[-1]): 
        # indices of the parent with the highest M_max (i.e. primaries) 
        i_primaries = sub[isnap]['par.i'][ihalos]
        # only keep halos that have parents
        # at snapshot=0 ~40 do not have parents
        has_primary = (i_primaries >= 0) 
        ihalos = ihalos[has_primary]
        i_primaries = i_primaries[has_primary]
        
        # make sure indices match across snapshot i and i+1 
        assert np.sum(ihalos - sub[isnap+1]['chi.i'][i_primaries]) == 0
        
        # identify halos with more than just a primary parent
        # the few lines below is just to speed things up a bit 
        counter = Counter(list(sub[isnap+1]['chi.i']))
        ip1_children_dupl = np.array([ii for ii, count in counter.items() if count >1])
        ihs = np.intersect1d(ihalos, np.array(ip1_children_dupl))
    
        for ih in ihs: 
            # loop through halos with parents and store non-primary halo masses 
            is_parent = np.where(sub[isnap+1]['chi.i'] == ih)[0]
            i_primary = sub[isnap]['par.i'][ih]
            notprimary = (is_parent != i_primary) 
            
            ifh=wUT.catalog.indices_tree(sub, isnap, 0,ih)
            M0=sub[0][mkind][ifh]
            m_M0 += list(sub[isnap+1][mkind][is_parent[notprimary]] - M0)
            #m_M0 += list(sub[isnap+1][mkind][is_parent[notprimary]] - sub[isnap][mkind][ih])
    
        print('snapshot %i, %i m_M0s' % (isnap, len(m_M0)))
        # keep searching through the primaries 
        ihalos = i_primaries
    print('%i halos have primaries in the last snapshot' % len(ihalos)) 
    m_M0 = np.array(m_M0)

    timestmp='{:%Y%m%d}'.format(datetime.datetime.now())
    filename='./dat/m_M0_{0:0.0f}_{2}_through_{1}.h5'.format(mmid,timestmp,mkind)
    f = h5py.File(filename, 'w') 
    f.create_dataset('m_M0', data=m_M0) 
    f.close() 
    return None 


def plotMFaccreted(datfile='m_M0.h5'): 
    ''' plot the mass funciton
    ''' 
    
    #import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['axes.xmargin'] = 1
    mpl.rcParams['xtick.labelsize'] = 'x-large'
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['ytick.labelsize'] = 'x-large'
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['legend.frameon'] = False

    f = h5py.File(datfile, 'r')
    m_M0 = f['m_M0'].value 
    
    fig = plt.figure(figsize=(6,6))
    sub = fig.add_subplot(111)
    bins = np.logspace(-3., np.log10(0.3), 15) 
    sub.hist(m_M0, weights=np.repeat(1./7e4, len(m_M0)), bins=bins, histtype='step', cumulative=-1)
    sub.set_xlabel(r'$m/M_0$', fontsize=25)
    sub.set_xscale("log") 
    sub.set_xlim([1e-3, 0.3]) 
    sub.set_ylabel(r'$n(>m)$', fontsize=25)
    sub.set_yscale("log")
    #sub.set_ylim([0.2, 35]) 
    if datfile=='m_M0_through.h5':
        fig.savefig('m_M0_through_new.png', bbox_inches='tight')
    else:
        fig.savefig('m_M0_new.png', bbox_inches='tight')
    return None 


if __name__=="__main__": 
    #merger_track(12, dlogm=0.5)
    merger_track_throughout(12, dlogm=0.5)
    #plotMFaccreted()
