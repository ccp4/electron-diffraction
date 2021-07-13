from utils import*                  ;imp.reload(dsp)
from EDutils import viewers as vw   ;imp.reload(vw)
plt.close('all')


v = vw.Viewer(path='dat',cif_file='dat/alpha_glycine.cif',init_opts='',Nmax=6,Smax=0.025,tag='vw',
    pets_opts='SVkr',frame=17,thick=24,rotS=70,rot=205,multi=-1)
# v = vw.Viewer(config='dat/config.pkl')
