from utils import*
from EDutils import utilities as ut;imp.reload(ut)
from blochwave import bloch_pp as bl;imp.reload(bl)
from EDutils import rotate_exp ;imp.reload(rotate_exp )
import os
from subprocess import check_output

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
rock = ut.load_pkl(file=rock_file)

def test_change_path():
    new_path = 'dat/rock_new_path'
    if os.path.exists(new_path):
        o=check_output('rm -rf %s;mkdir %s' %(new_path,new_path),shell=True).decode()
    o=check_output('cp -r %s/* %s' %(out,new_path),shell=True).decode() #;print(o)
    rock = ut.load_pkl(file=rock_file)
    rock.change_path(new_path)
    rock = ut.load_pkl(file='%s/rock_.pkl' %new_path)
    b=rock.load(0)
    b.save()

test_change_path()

if 0:
    df_int = r.integrate()

    df = rock.integrate(thick=100)
    I_exp = df[df.I > 0.001].copy()
    I_exp = I_exp.drop(str((0,0,0)))
    #noise the data
    I_exp.I*=(np.random.rand(I_exp.shape[0])/5+1)
    I_exp.loc[str((-10,2,5))]=1.1

    r.Rfactor(I_exp)
