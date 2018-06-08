#test
import os
from centralms import catalog as cat
#from centralms import util as UT
class Testsubhalo(cat.Subhalos):
    def __init__(self, sigma_smhm=0.2, smf_source='li-march',nsnap0=20,use_super_defaults=True):
        if use_super_defaults:
            super(self).__init__(self)
        else:
            super(self).__init__(self,sigma_smhm=sigma_smhm,smf_source=smf_source,nsnap0=nsnap0)
    def File(self):
	file_name=''.join(['/home/users/staudt/projects/mergers',
        	'Subhalos',
		'.SHAMscat',str(self.sigma_smh),
		'.smf_',self.smf_source,
		'nsnap0_',str(self.nsnap0),
		'.hdf5'])
	return file_name
