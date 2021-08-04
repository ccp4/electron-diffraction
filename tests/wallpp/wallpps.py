from utils import*
from wallpp import config as cg      ;imp.reload(cg)
from wallpp import wallpaper as wp   ;imp.reload(wp)
plt.close('all')
fig_path = 'figures/'


pattern = [[0.2,0.1,1],[0.1,0.3,2]]
pp_args = {'a':2,'b':3,'alpha':90,'pattern':pattern,
    'fract':True,'interp':False,'nh':3,'nk':3}

wallpps={}
fz = lambda x :np.log10(np.maximum(abs(x),1e-5))
for i,pp_type in enumerate(cg.pp_types):
    pp = wp.Wallpaper(pp_type,**pp_args)
    fig,(ax1,ax2) = dsp.create_fig('f',rc='12',pad=5)
    pp.plot('uag',ax=ax1,pOpt='tX',labs=['$x$','$y$'],xylims=[0,5,0,5],
        ms=50,opt='')
    pp.show_structure_factor(title=pp_type,fz=fz,
        caxis=[-5,0],ax=ax2,pOpt='tX',xylims=5,
        name=fig_path+'sf_%s_%s.png' %(str(i).zfill(2),pp_type),opt='sc')
    wallpps[pp_type] = pp
