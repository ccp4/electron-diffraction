from utils import*
from wallpp import config as cg      ;imp.reload(cg)
plt.close('all')
fig_path = 'figures/'


a,b,angle = 2,3,70
ndeg,nh,nk = 30,3,2

pp_args = {'a':a,'b':b,'angle':angle,'pattern':pattern,
        'ndeg':ndeg,'nh':nh,'nk':nk,'path':fig_path}
        'pOpt' :'aue sc','gen':True,'pOpt':'VA'}

wallpps=dict()
for pp_type in cg.pp_types:#['cm']:
    pp = pg.Wallpaper(pp_type,**pp_args)
    pp.plot_unit_cells('ua',name=fig_path+'%s_unitcell.png' %pp_type)
    wallpps[pp_type] = pp
