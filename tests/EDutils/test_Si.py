from utils import*                  ;imp.reload(dsp)
from EDutils import viewers as vw   ;imp.reload(vw)
plt.close('all')



# v = vw.Viewer(path='dat/diamond',cif_file='diamond',u=[0,0,1],pets_opts='L',
#     init_opts='R')
v = vw.Viewer(config='dat/diamond/config.pkl')
