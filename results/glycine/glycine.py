from utils import*                  #;imp.reload(dsp)
from multislice import pets as pt   ;imp.reload(pt)
from blochwave import bloch_pp as bl;imp.reload(bl)
from EDutils import utilities as ut ;imp.reload(ut)
# plt.close('all')

Nframes = 30


pets = pt.Pets('dat/pets/glycine.pts')
uvw,alpha = pets.uvw0[:Nframes],pets.alpha[:Nframes]
bloch_args = {'cif_file':'dat/pets/alpha_glycine.cif','Smax':0.025,'Nmax':9,
    'solve':1,'opts':'s'}

cond = '(Sw<1e-2) & (Vga>1e-6)'
refl = [[-3,-1,2],[3,1,-4],[4,1,0],[-3,-1,4],[0,0,-2],[0,0,2]]
hkls = [[h] for h in refl]

rock = bl.bloch_rock(tag='glycine',uvw=uvw,bloch_args=bloch_args,omega=alpha,
    cond=cond,hkls=hkls,thick=500,thicks=(10,1000,50),iTs=slice(2,14),
    i=13,ts0=None,rot=205,
    path='dat/bloch/',opts='',opt='p')

fz = lambda x : -np.log10(np.maximum(x,1e-10))
refl = [[-3,-1,2]]
# df  = pets.rpl
# df0 = df.loc[(df.F<30) & (df.I>1)].sort_values('F')[['h','k','l','F']]
# df0['refl'] = [str(h) for h in df0[['h','k','l']].values]
# refl = pd.unique(df0.refl)
cond = ''#  '(Sw<1e-2) & I>1e-2'
Sw = rock.Sw_vs_theta(refl=refl,cond=cond,fz=fz,
    iTs=slice(2,30,1),opts='tf')
