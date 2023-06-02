from utils import*
from EDutils import utilities as ut;imp.reload(ut)
from blochwave import bloch_pp as bl;imp.reload(bl)
from EDutils import rotate_exp ;imp.reload(rotate_exp )
import os

out = 'dat/rock'
new = 0

rock_file=file=out+'/rock_.pkl'
tag = ''
thicks=(0,100,101)
def test_rock():

    bloch_args = {'cif_file':'diamond','keV':200,'thick':200,'thicks':thicks,
        'Nmax':4,'Smax':0.025,'solve':True,'opts':'stz'}
    uvw = ut.get_uvw_rock(e0=[0,0,1],e1=[2,1],deg=3,npts=5,show=0)
    rock = bl.Bloch_cont(path=out,tag=tag,uvw=uvw,Sargs=bloch_args)

if not os.path.exists(rock_file) or new:test_rock()


r = ut.load_pkl(file=rock_file)
# r.set_beams_vs_thickness(thicks)
df = r.integrate(thick=10)
# b0 = r.load(0)
# thicks=['%.1fA' %z for z in b0.z ]
# df = pd.DataFrame(0,index=r.beams.index,columns=thicks)
#
# for i in range(r.n_simus):
#     print(i)
#     b0 = r.load(i)
#     df.loc[b0.df_G.index] += b0.Iz
