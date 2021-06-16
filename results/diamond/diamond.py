from utils import*                  #;imp.reload(dsp)
from EDutils import viewers as vw   ;imp.reload(vw)
plt.close('all')

# v = vw.Viewer(path='dat',cif_file='diamond',init_opts='R')
v = vw.Viewer(config='dat/config.pkl')
