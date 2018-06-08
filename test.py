#test
import os
from centralms import catalog as cat
from centralms import util as UT
class testsubhalo(cat.Subhalos):
	'''def File(self):
		file_name=''.join([
	'''
print __file__
realpath=os.path.realpath(__file__)
print realpath
print os.path.dirname(realpath)
print UT.dat_dir()
