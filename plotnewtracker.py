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

from progressbar import ProgressBar

tmp=20181107
fname='/home/users/staudt/projects/mergers/dat/true_m_M0_m.max_{}.h5'.format(timestmp)
f=h5py.File(fname,'r')
M0s=np.array(f['M0s'])
m_M0s=np.array(f['m_M0s'])
Ms=np.array(f['Ms'])
