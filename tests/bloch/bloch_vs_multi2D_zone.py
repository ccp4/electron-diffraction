from utils import*
from blochwave import bloch2D as bl          ;imp.reload(bl)
from multislice import mupy_utils as mut     ;imp.reload(mut)
import multislice.multi_2D as ms             ;imp.reload(ms)
plt.close('all')
opts = 'BM' #(Bloch) (Multislice)

# K = 5.024937810560445
# keV = cst.lam2keV(1/K)
keV=200

fig,axB=dsp.create_fig()
cmap  ='jet'
thick = 100
eps   = 0.1

#### Bloch
if 'B' in opts:
    file = 'dat/p1.py'
    # p1 = mut.import_wallpp(file)
    # p1.plot_unit_cells(opts='uAa',nh=3,nk=3)
    b0  = bl.Bloch2D(file,keV=keV,u=[0,1],thick=thick,eps=eps,
        Nmax=5,Smax=0.1,solve=1)
    # b0.show_beams(opts='B',fz=np.abs)
    # b0.show_ewald()
    # b0.show_beams(opts='SVB')
    b0.show_beams_vs_thickness(thicks=(0,thick,1000),strong=['I'],m={'I':1e4},
        linespec='--',marker='',cm=cmap,ax=axB,opt='')
    # b0.G('Sw')
    # idx = b0.get_Istrong(Icols=['I'],m={'I':1e5},out=1)
    # ibs = b0.get_refl(idx)

#multislice
if 'M' in opts:
    ax,bz = b0.crys.params[:2]
    dz=bz/1
    Nx,Nz=5,int(thick/dz)
    x,z,f = b0.crys.get_potential_grid_p1(ndeg=2**6)

    tilt = 0#np.linspace(0,0.08,40)[0]
    ms0 = ms.Multi2D([x,z,f],ax,bz,keV=keV,Nx=Nx,nz=Nz,
            dz=dz,eps=eps,sg=-1,
            iZs=1,iZv=np.inf,tilt=tilt)
    # ms0.V_show()
    # ms0.Qz_show([-1],opts='nO',xylims=['x',-6,6],lw=2)
    ms0.Bz_show(np.arange(-1,2)*Nx,lw=2,cmap=cmap,
        ax=axB,opt='',xylims=[0,thick,0,1])

fig.show()
