from utils import*                          ;imp.reload(dsp)
from EDutils import viewers as vw           ;imp.reload(vw)
from multislice import mupy_utils as mut    ;imp.reload(mut)
from multislice import pets as pt           ;imp.reload(pt)
import scipy.optimize as opt
plt.close('all')
pets = 'dat/pets/glycine.pts'

pets = pt.Pets('dat/pets/glycine.pts')
# v0 = mut.Viewer(pets=pets,init='cr',rot=203,Nmax=13,Smax=0.025,frame=14)
# v = mut.Viewer(exp_path='dat/pets/tiff/processed',init='cr',rot=203,Nmax=13,Smax=0.025,frame=14)
# v = vw.Viewer(config='dat/config.pkl')

# v1 = vw.Viewer(path='dat',cif_file='dat/alpha_glycine.cif',init_opts='',
#     Nmax=6,Smax=0.025,tag='vw',dyn_opt=1,
#     pets_opts='SVkr',frame=17,thick=24,rotS=70,rot=205,multi=-1)
v = vw.Viewer(path='dat',cif_file='dat/alpha_glycine.cif',init_opts=' ',
    Nmax=9,Smax=0.025,tag='vw',dyn_opt=0,
    pets_opts='PBSkr',frame=14,thick=300,rot=203,multi=-1,rotS=70)
