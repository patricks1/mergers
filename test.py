#test
import os
from centralms import catalog as cat
#from centralms import util as UT
class testsubhalo(cat.Subhalos):
	def File(self):
		file_name=''.join(['/home/users/staudt/projects/mergers',
			'Subhalos',
			'.SHAMscat',str(self.sigma_smh),
			'.smf_',self.smf_source,
			'nsnap0_',str(self.nsnap0),
			'.hdf5'])
		return file_name

print __file__
realpath=os.path.realpath(__file__)
print realpath
print os.path.dirname(realpath)
print UT.dat_dir()
