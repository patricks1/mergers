import merg_rates
import h5py

import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 

fname='/home/users/staudt/projects/mergers/dat/simruns/simdat.h5'
try f=h5py.File(fname,'a'):
    action='append'
except:
    f=h5py.File(fname,'w')
    action='create'

tpm=merg_rates.shamedTreepmClass(0.2)
tpm.mp_tree()
tpm.merg_tree()

muax_11,Ns11=N_mu_ft(10.5,'cengal')
muax_105,Ns105=N_mu_ft(10.5,'cengal')
muax_95,Ns95=N_mu_ft(9.5,'cengal')

if action=='create':
    with h5py.File(fname,'w') as f:
        f.create_dataset('Ns11',[muax_110,Ns110])
        f.create_dataset('Ns105',[muax_105,Ns105])
        f.create_dataset('Ns95',[muax_095,Ns095])
elif action=='append':
    with h5py.File(fname,'a') as f:
        Ns11_h5=f['Ns11']
        Ns105=f['Ns105']
        Ns95=f['Ns95']
f.close()
