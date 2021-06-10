from utils import*                  #;imp.reload(dsp)
from EDutils import viewers as vw   #;imp.reload(vw)
plt.close('all')

v = vw.Viewer(path='dat',init_opts='',Nmax=5,
    pets_opts='PMS',frame=18,thick=24,rotS=70,rot=205)
