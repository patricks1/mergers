import os
from centralms import catalog as cat

class Staudt_subhalos(cat.Subhalos):
    '''
    def __init__(self, sigma_smhm=None, smf_source=None, nsnap0=None,use_super_defaults=True):
        if use_super_defaults:
            super(Testsubhalo,self).__init__(sigma_smhm=0.2, smf_source='li-march',nsnap0=20)
        else:
            super(Testsubhalo,self).__init__(self,sigma_smhm=sigma_smhm,smf_source=smf_source,nsnap0=nsnap0)
    '''
    def File(self):
        file_name=''.join(['/home/users/staudt/projects/mergers/dat/',
            'Subhalos',
            '.SHAMscat',str(self.sigma_smhm),
            '.smf_',self.smf_source,
            'nsnap0_',str(self.nsnap0),
            '.hdf5'])
        return file_name
