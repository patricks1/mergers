import merg_rates                                                                  
import h5py                    
import numpy as np
from progressbar import ProgressBar
                                                                                   
import matplotlib as mpl                                                           
mpl.use('Agg')                                                                     
import matplotlib.pyplot as plt     
                                                                                   
fname='/home/users/staudt/projects/mergers/dat/simruns/simdat_0.4.h5'                  
   
tpm=merg_rates.shamedTreepmClass(0.4)                                              
tpm.mp_tree()                                                                      
tpm.merg_tree()                                                                    

###Mergers per galaxy stats####################################################
muax_11,Ns11=merg_rates.N_mu_ft(tpm,10.5,'cengal')                                 
muax_105,Ns105=merg_rates.N_mu_ft(tpm,10.5,'cengal')                               
muax_95,Ns95=merg_rates.N_mu_ft(tpm,9.5,'cengal')

with h5py.File(fname,'a') as f:                                                
    matrix=[muax_11,Ns11]
    if 'Ns11' in f.keys():
        d=(f['Ns11']).value
        newd=np.append(d,[matrix],axis=0)
        del f['Ns11']
        f.create_dataset('Ns11',data=newd)                                    
    else:
        f.create_dataset('Ns11',data=[matrix])                                    
    matrix=[muax_105,Ns105]
    if 'Ns105' in f.keys():
        d=(f['Ns105']).value
        newd=np.append(d,[matrix],axis=0)
        del f['Ns105']
        f.create_dataset('Ns105',data=newd)                                    
    else:
        f.create_dataset('Ns105',data=[matrix])                                    
    matrix=[muax_95,Ns95]
    if 'Ns95' in f.keys():
        d=(f['Ns95']).value
        newd=np.append(d,[matrix],axis=0)
        del f['Ns95']
        f.create_dataset('Ns95',data=newd)                                    
    else:
        f.create_dataset('Ns95',data=[matrix])                                    
###############################################################################

###Histogram stat##############################################################
Mconds=[9.5,10.5,11.]
rngs=[[-5.,-1],[-1,np.log10(1./3.)],[np.log10(1./3.),0]]
pbar=ProgressBar()
for Mcond in pbar(Mconds):
    for rng in rngs:
        key='Ps_{0:0.1f}M_{1:0.1f}-{2:0.1f}'\
            .format(Mcond,rng[0],rng[1])
        #Get [count axis,probabilities] in the form of a
        #single matrix from hgram_dat_ft:
        matrix=merg_rates.hgram_dat_ft(tpm,rng,Mcond,1,34)
        with h5py.File(fname,'a') as f:            
            if key in f.keys():
                Ps_h5=(f[key]).value
                new=np.append(Ps_h5,[matrix],axis=0)
                del f[key]
                f.create_dataset(key,data=new)
            else:
                f.create_dataset(key,data=[matrix])
###############################################################################

####Quench stats###############################################################
mus=[np.log10(1./3.),np.log10(1./10.)]
matrix1=merg_rates.quench_frac_ft(tpm,mus[0])
matrix2=merg_rates.quench_frac_ft(tpm,mus[1])
matrix=[matrix1,matrix2]

with h5py.File(fname,'a') as f:                                                
    if 'quench_dat' in f.keys():
        d=(f['quench_dat']).value
        newd=np.append(d,[matrix],axis=0)
        del f['quench_dat']
        f.create_dataset('quench_dat',data=newd)                                    
    else:
        f.create_dataset('quench_dat',data=[matrix])  
###############################################################################

#dNdz_ofmu#####################################################################
with h5py.File(fname,'a') as f:
    M0s=[9.5,10.5,11.]
    for M0 in M0s:
        key='{0:0.1f}dNdz_ofmu'.format(M0)
        new_dat=merg_rates.dNdx_ofmu(tpm,M0,'cengal','z',
                                 forcem200=False,through=False)
        if key in f.keys():
            dat_ld=(f[key]).value
            dat_apnd=np.append(dat_ld,[new_dat],axis=0)
            del f[key]
            f.create_dataset(key,data=dat_apnd)
        else:
            f.create_dataset(key,data=[new_dat])

#dNdz_ofz######################################################################
with h5py.File(fname,'a') as f:
    M0s=[9.5,10.5,11.]
    mu_conds=[np.log10(1./3),np.log10(1./10.)]
    for M0 in M0s:
        dat3=merg_rates.dNdx_ofz(tpm,M0,mu_conds[0],'cengal',
                                 zibeg=1,ziend=34,through=False)
        dat10=merg_rates.dNdx_ofz(tpm,M0,mu_conds[1],'cengal',
                         zibeg=1,ziend=34,through=False)
        new_dat=[dat3,dat10]
        key='{0:0.1f}dNdz_ofz'.format(M0)
        if key in f.keys():
            dat_ld=(f[key]).value
            dat_apnd=np.append(dat_ld,[new_dat],axis=0)
            del f[key]
            f.create_dataset(key,data=dat_apnd)
        else:
            f.create_dataset(key,data=[new_dat])
