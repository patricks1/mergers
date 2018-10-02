#!/bin/python
'''
'''
import h5py 
import numpy as np 
from collections import Counter
from treepm import subhalo_io


def merger_track(mmid, dlogm=0.5): 
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
    in_mbin = (sub[0]['halo.m'] > mmid - 0.5*dlogm) & (sub[0]['halo.m'] < mmid + 0.5*dlogm)

    ihalos = np.arange(nsub0)[is_cen & in_mbin]
    m_M0 = [] 
    for isnap in range(snapshots[-1]): 
        # indices of the parent with the highest M_max (i.e. primaries) 
        i_primaries = sub[isnap]['par.i'][ihalos]
        has_primary = (i_primaries >= 0) # at snapshot=0 ~40 do not have parents
        
        # only keep halos that have parents
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

            m_M0 += list(sub[isnap+1]['halo.m'][is_parent[notprimary]] - sub[isnap]['halo.m'][ih])
    
        print('snapshot %i, %i m_M0s' % (isnap, len(m_M0)))
        # keep searching through the primaries 
        ihalos = i_primaries
    
    m_M0 = np.array(m_M0)

    f = h5py.File('m_M0.h5', 'w') 
    f.create_dataset('m_M0', data=m_M0) 
    f.close() 
    return None 


if __name__=="__main__": 
    merger_track(12, dlogm=0.5)
