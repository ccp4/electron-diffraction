from utils import*                  ;imp.reload(dsp)
from EDutils import utilities as ut ;imp.reload(ut)
from blochwave import bloch_pp as bl;imp.reload(bl)
from utils import pytest_util
import pytest,os
plt.close('all')

tag='big' #do not change as it messes up with the cmp_ref
npts={'big':50,'small':5}[tag]

@pytest.mark.slow
def test_rock():
    uvw = ut.get_uvw_rock(e0=[0,0,1],e1=[2,1],deg=3,npts=npts,show=0)
    rock = bl.Bloch_cont(path=out,tag=tag,uvw=uvw,Sargs=bloch_args)

try:
    out,ref,dir = pytest_util.get_path(__file__)
    bloch_args = {'cif_file':'diamond','keV':200,'thick':200,'thicks':(0,200,500),
        'Nmax':9,'Smax':0.025,'solve':True,'opts':'stz'}
    rock_file=file=out+'/rock_%s.pkl' %tag
    if not os.path.exists(rock_file):test_rock()
    rock = ut.load_pkl(file=rock_file)
except Exception as e:
    print(e)
    pass

def test_load():
    rock.load(0)

@pytest.mark.lvl1
@pytest_util.add_link(__file__)
def test_show_tiff():
    rock.convert2tiff(thick=200,Imax=1e6,aperpixel=0.01)
    vw=rock.show_tiff(cutoff=40,frame=2,pargs={'opt':''},
        h=False)
    return vw.fig,vw.ax

@pytest.mark.lvl1
@pytest_util.add_link(__file__)
def test_plot_rocking_cond():
    rock.set_beams_vs_thickness(thicks=(0,100,100))
    cond=lambda dfG:bl.strong_beams(dfG,tol=1e-2,n=5)
    return rock.plot_rocking(cond=cond,opts='',opt='')

@pytest.mark.lvl1
@pytest_util.add_link(__file__)
def test_plot_rocking_quick():
    return rock.plot_rocking(refl=[str((2,2,0))],zs=[50,100],opts='',opt='')

# @pytest.mark.lvl1
@pytest_util.cmp_ref(__file__)
def test_integrate_rocking():
    refl,h=rock.get_beams(cond=lambda dfG:bl.strong_beams(dfG,tol=1e-2,n=5))
    rock._integrate_rocking(refl=refl,new=1)
    return rock.Iz_dyn[refl[0]]

@pytest.mark.lvl1
@pytest_util.add_link(__file__)
def test_plot_integrate():
    refl,h=rock.get_beams(cond=lambda dfG:bl.strong_beams(dfG,tol=1e-2,n=5))
    # refl=refl[2:3]
    return rock.plot_integrated(refl=refl,opt='')

# rock._build_index()
# test_rock()
# fig,ax=test_plot_rocking_quick()
# test_integrate_rocking()
# test_plot_integrate()

# rock.plot_rocking(refl=['(-2, -2, 0)','(-2, 2, 0)','(2, 2, 0)','(2, -2, 0)'])
# hkls=bl.strong_beams(rock.load(20).df_G,tol=1e-3,n=5,opt='F')
# print(hkls)
# hkls = ['(-2, -2, 0)','(-2, 2, 0)','(2, 2, 0)','(2, -2, 0)']
# i=0
# df = rock.get_frames(hkls[i])
# z,I = rock.get_rocking(refl=hkls,zs=np.arange(25,100,25))#iZs=slice(1,None,50))#zs=np.arange(25,200,50))
# rock.plot_rocking(refl=hkls[i:i+1],zs=np.arange(25,100,25))
# rock.plot_rocking(refl=[hkls[i]],iZs=[-1])#zs=np.arange(25,100,25))
# rock.plot_rocking(refl=hkls,iZs=[-1],x='Sw')


# cond = '(Sw<2e-2) & (Vga>1e-6) & (I>1e-3)'
# cond = '(Sw<1e-4) & (I>1e-4) '
# fz = lambda x : -np.log10(np.maximum(x,1e-10))
# cond = ''#'(Sw<1e-2) & (Vga>1e-6)'
# hkls = [[-8,-2,2],[-5,1,1],[-1,7,-1],[3,7,-1],[-2,6,0]]
# refl = [[-8,-2,2],[-5,1,1],[-1,-7,1],[3,7,-1],[-2,6,0],[7,-3,-1],[5,1,-1]]
# hkls = [[refl[1]]]
# rock = bl.bloch_rock(tag='diamond_rand',uvw=uvw,
#     omega=omega,bloch_args=bloch_args,
#     thicks=(0,800,400),
#     ts0=1.0,refl=refl,cond=cond,
#     zs=np.arange(1,6)*50,hkls=hkls,
#     Z_args={},#{'new':1},
#     path=path,opts=opts,opt=opt)


# rock = ut.load_pkl('dat/bloch/rock_diamond_rand.pkl')
# b = rock.load(0)
# b.get_beam(cond=cond)
