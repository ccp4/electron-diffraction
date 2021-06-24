from utils import*                  #;imp.reload(dsp)
from EDutils import viewers as vw   ;imp.reload(vw)
# from EDutils import viewer_config as vw_cfg       ;imp.reload(vw_cfg)

plt.close('all')

# v = vw.Viewer(path='dat',init_opts='',Nmax=5,tag='test',
#     pets_opts='PMS',frame=18,thick=24,rotS=70,rot=205)

v = vw.Viewer(config='dat/config.pkl')
