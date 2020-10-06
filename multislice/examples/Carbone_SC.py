import importlib as imp
from utils import*
import multislice.multislice as mupy ;imp.reload(mupy)
import multislice.postprocess as pp
import matplotlib.animation as animation

name ='../dat/Carbone_SC/'
Nz=10
def run(*kwargs):
    return mupy.Multislice(name,tail='test',data='Carbone_SC001A.dat',
        mulslice=1,keV=200,TDS=False,
        NxNy=1024,slice_thick=1.0,repeat=[8,8,Nz],
        Nhk=10,hk_sym=1,**kwargs)
def run_tilts(**kwargs):
    tilts = np.vstack([np.linspace(0,1,30),0*np.ones((30))]).T
    mupy.sweep_var(name,'tilt',tilts,do_prev=0,
        mulslice=1,keV=200,TDS=False,
        NxNy=1024,slice_thick=1.0,repeat=[8,8,Nz],
        Nhk=10,hk_sym=1,**kwargs)

# multi = run(opt='dsrp',fopt='f',ppopt='wB')
# multi = run_tilts(opt='dsr',fopt='f',ppopt='',ssh='tarik-CCP4home',v=1)
# df = pp.update_df_info(name+'df.pkl',files=['beams'])
# df = pd.read_pickle(name+'df.pkl')
multi = pp.load_multi_obj(name+'Carbone_SC_tilt016_mulslice.pkl')
multi.run(ssh_alias='tarik-CCP4home')
# multi.print_log()
# multi.beam_vs_thickness(cm='seismic')
# multi.get_beams(cm='seismic')
# multi.gif_beams()
