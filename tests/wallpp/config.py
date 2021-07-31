from utils import*
from matplotlib import patches
from wallpp import config as cg      ;imp.reload(cg)



for pp_type in cg.pp_types:
    verts = cg.df_wallpp.loc[pp_type,'asym_cell']
    pp = patches.Polygon(verts,alpha=0.3,edgecolor='k',linewidth=1.5)
    dsp.stddisp(patches=[pp],title=pp_type)
