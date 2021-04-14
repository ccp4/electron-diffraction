import importlib as imp
from utils import*                      ;imp.reload(dsp)
import multislice.pymultislice as pyms  ;imp.reload(pyms)
import multislice.multi_2D as MS2D      ;imp.reload(MS2D)
import wallpp.lattice as lat            ;imp.reload(lat)
import wallpp.plane_group as pg         ;imp.reload(pg)
import scipy.fftpack as fft

plt.close('all')
path='../../multislice/docs_fig/multi2D/TDS/'
opts='0'

ax,bz = 8,10
pattern = [[2,3,2],[4,6,4],[1.5,2.0,1],[2.5,1.5,1]]
ndeg = 2**10
keV  = 200


Nx,nTDS = 32,32
file = 'dat/TDS_Nx%d_nTDS%d.pkl' %(Nx,nTDS)
if 'r' in opts:
    mp = MS2D.Multi2D(pattern,ax,bz,keV=keV,
            Nx=Nx,dz=2,nz=10,eps=0.25,nx=ndeg,
            TDS=True,nTDS=nTDS,wobble=np.array([0.002,0.025,0.025,0.025,0.05])*2, #[0.01,1:0.025,2:0.025,4:0.05},
            ppopt='',opts='q',v=1)
    mp.save(file)
else :
    mp = pyms.load(file)

mp.Za_show()
mp.Qz_show(iZs=[-1],title='with TDS',xylims=['x',-2,2])

# mp.Bz_show(iBs=np.arange(0,4),sym_opt=1,lw=2)

if '0' in opts:
    # mp = MS2D.Multi2D(pattern,ax,bz,keV=keV,
    #         Nx=Nx,dz=2,nz=10,eps=0.25,nx=ndeg,
    #         TDS=True,nTDS=1,wobble=np.array([0.002,0.025,0.025,0.025,0.05])*0, #[0.01,1:0.025,2:0.025,4:0.05},
    #         ppopt='',opts='q',v=1)

    p1 = pg.Wallpaper('p1',ax,bz,90,pattern,ndeg=ndeg,gen=True,nh=3,nk=3)
    # p1.plot_unit_cells(opts='AV',nh=2,nk=2)
    potential = p1.get_potential_grid_p1()
    #
    # # Loop over 1 unit cell periodic
    mp = MS2D.Multi2D(potential,ax,bz,keV=keV,
            Nx=Nx,dz=2,nz=10,eps=0.25,
            ppopt='',opts='q',v=1)

    # mp.Ewald_show()
    # mp.Tz_show(Vopt='V',lw=2,cmaps='jet')
    mp.Za_show()
    mp.Qz_show(iZs=[-1],xylims=['x',-2,2],title='no TDS')
    # mp.Bz_show(iBs=np.arange(0,4)*mp.Nx,sym_opt=1,lw=2)
    # mp.Bz_show(iBs=np.arange(6),lw=2)
