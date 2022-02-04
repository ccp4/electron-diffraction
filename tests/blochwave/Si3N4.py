from utils import*                      #;imp.reload(dsp)
from blochwave import bloch             ;imp.reload(bloch)
from EDutils import utilities as ut     #;imp.reload(ut)
from multislice import mupy_utils as mut#;imp.reload(mut)
plt.close('all')
cif = 'dat/bloch/Si3N4/Si3N4.cif'
uza = [0,0,1]

crys = mut.import_crys(cif)
lat_vec = np.array(crys.lattice_vectors)
u = lat_vec.T.dot(np.array(uza))
b = bloch.Bloch(cif_file=cif,keV=200,u=u,Smax=0.01,Nmax=5,solve=0,opts='s')
