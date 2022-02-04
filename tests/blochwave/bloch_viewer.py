from utils import*                  ;imp.reload(dsp)
from EDutils import viewers as vw   ;imp.reload(vw)

# plt.close('all')

# v = vw.Viewer(config='dat/config.pkl')
v = vw.Viewer(config='dat/4_2beam_config.pkl')

# v = vw.Viewer(cif_file='diamond',orient=[0]*2+[0.1]*2,
#     Smax=0.025,Nmax=8,thick=250,path='dat',
#     init_opts='R',pets_opts='S',xylims=3)
