import numpy as np

from wetzel_utils.utility import array as ut_array 
from numpy import log10, Inf

def indices_tree(cat, zi_start, zi_end, cis=None, get_indices=False):
    '''
    Get parent / child index[s] of input cis at zi_end.
    Assign negative value to [sub]halo if cannot track all the way back.
    Import [sub]halo catalog, starting & ending snapshots (forward or backward), index[s].
    '''
    if zi_start == zi_end:
        return cis
    elif zi_end >= len(cat):
        raise ValueError('zi.end = %d not within %s snapshot limit = %d' %
            (zi_end, cat.info['kind'], len(cat) - 1))
    zi_step, tree_kind = get_tree_direction_info(zi_start, zi_end)
    if cis is None:
        cis_tree = ut_array.arange_length(cat[zi_start][tree_kind]) #index numbers
    else:
        cis_tree = ut_array.arrayize(cis, bit_num=32)
    print'cis_tree:'
    print cis_tree
    print''
    for zi in xrange(zi_start, zi_end, zi_step): #for zi in zis
        ciis = ut_array.elements(cis_tree, [0, Inf]) #returns the index numbers of all non-negative elements of cis_tree (which after the first run is all parent IDs). Parent ID is negative when the parent doesn't exit. 
        print('snapshot {0:0.0f}'.format(zi))
        print'halo IDs:'
        print cat[zi]['halo.i']
        print'parent IDs:'
        print cat[zi]['par.i']
        cis_tree[ciis] = cat[zi][tree_kind][cis_tree[ciis]] #ciis is an array of cis_tree index numbers. cis_tree are parent IDs. Therefore, cis_tree[ciis] is an array of the objective parent IDs.
        print'cis_tree:'
        print cis_tree
        objective_halo_index_num=list(cat[zi]['halo.i']).index(2)
        print'parent ID of halo 2:'
        print cat[zi]['par.i'][objective_halo_index_num]
        print'child ID of halo 11969:'
        objective_halo_index_num=list(cat[zi]['halo.i']).index(11969)
        print cat[zi]['chi.i'][objective_halo_index_num]
        print''
        #on the first loop, this replaces cis_tree with parent IDs. Then the first line of the loop, on the next run, gets rid of all negative parent IDs, i.e. non-existent parent IDs. 
        #before the first run of this line, cis_tree[ciis]==ciis. This line takes the elements ciis of cis_tree and changes their values.
    ciis = ut_array.elements(cis_tree, [0, Inf])
    cis_end = np.zeros(cis_tree.size, cis_tree.dtype) - 1 - cat[zi_end][tree_kind].size
    cis_end[ciis] = cis_tree[ciis]
    if get_indices:
        return ut_array.scalarize(cis_end), ut_array.scalarize(ciis)
    else:
        return ut_array.scalarize(cis_end)

def get_tree_direction_info(zi_start=None, zi_end=None, direction_kind=''):
    '''
    Get snapshot step (+1 or -1) & dictionary key corresponding to parent/child.

    Import starting & ending snapshot indices (forward or backward) OR
    direction kind ('parent', 'child').
    '''
    if direction_kind:
        if direction_kind == 'parent':
            zi_step = 1
            family_key = 'par.i'
        elif direction_kind == 'child':
            zi_step = -1
            family_key = 'chi.i'
    else:
        if zi_end < 0:
            raise ValueError('u.end = %d out of bounds' % zi_end)
        elif zi_end == zi_start:
            raise ValueError('u.end = u.start')
        if zi_end > zi_start:
            zi_step = 1
            family_key = 'par.i'
        elif zi_end < zi_start:
            zi_step = -1
            family_key = 'chi.i'
    return zi_step, family_key

