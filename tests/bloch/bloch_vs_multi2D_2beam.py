from utils import*
from blochwave import bloch2D as bl          ;imp.reload(bl)
from multislice import mupy_utils as mut     ;imp.reload(mut)
import multislice.multi_2D as ms             ;imp.reload(ms)
plt.close('all')
opts = 'BM' #(Bloch) (Multislice)

K = 5.024937810560445   # exact Bragg condition for hk=(1,0) beam
keV = cst.lam2keV(1/K)

fig,axB=dsp.create_fig()
cmap='jet'
thick=1000  #thickness
eps = 0.01  #strength of potential
n = 10      #u=[1,n]
e0=1
# e0=0.98,0.965
#### Bloch
if 'B' in opts:
    file = 'dat/p1.py'
    # p1 = mut.import_wallpp(file)
    # p1.plot_unit_cells(opts='uAa',nh=3,nk=3)
    b0  = bl.Bloch2D(file,keV=keV,u=[1,n],thick=thick,eps=e0*eps,
        Nmax=10,Smax=0.1,solve=1)
    # b0.show_beams(opts='B',fz=np.abs)
    # b0.show_ewald()
    # b0.show_beams(opts='SVB')
    b0.show_beams_vs_thickness(thicks=(0,thick,1000),strong=['I'],m={'I':100},
        linespec='--',marker='',cm=cmap,ax=axB,opt='')
    # b0.G('Sw')

    idx = b0.get_Istrong(Icols=['I'],out=1)
    ibs = b0.get_refl(idx)
    print('pendulosung thickness')
    print(b0.get_Xig())




#multislice
if 'M' in opts:
    file = 'dat/2_beam.npy'
    ax,bz = b0.crys.params[:2]
    pattern = b0.crys.pattern
    ax1 = np.sqrt(1**2+n**2)
    bz1 = ax1

    if 'V' in opts:
        import wallpp.plane_group as pg     ;imp.reload(pg)
        x1,y1,Za = pattern[0]
        x0,z0 = np.meshgrid(np.arange(-2,n+1),np.arange(-1,n+2))
        x0,z0 = np.stack([x0*ax+x1,z0*bz+y1])
        rot = lambda t:np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])
        Xa = rot(np.arctan(1/n)).dot(np.stack([x0.flatten(),z0.flatten()]))

        ndeg=2**6
        x,z = np.meshgrid(np.arange(n*ndeg)*ax1/(n*ndeg),np.arange(n*ndeg)*bz1/(n*ndeg))
        f = np.sum(np.array([pg.fv(np.stack([x.flatten(),z.flatten()]).T,X0,int(Za)) for X0 in Xa.T]),axis=0)
        f = np.reshape(f,(n*ndeg,n*ndeg))
        np.save(file,[x,z,f])
        print(colors.yellow+file+colors.green+" saved"+colors.black)
        if 'v' in opts:scat=[Xa[0],Xa[1],15,'k']

    x,z,f = np.load(file,allow_pickle=True)
    dz = bz1/(n)
    Nz = np.int(thick/dz) #1000*n
    tilt = np.linspace(0,0.08,40)[28]

    ms0 = ms.Multi2D([x[0],z.T[0],f],ax1,bz1,keV=keV,Nx=1,nz=Nz,
            dz=dz,eps=eps,sg=-1,
            iZs=1,iZv=np.inf,tilt=tilt)
    # ms0.V_show()
    # ms0.Qz_show([-1],opts='nO',xylims=['x',-6,6],lw=2)
    ms0.Bz_show(np.arange(0,2,1)*n,lw=2,ax=axB,opt='',cmap=cmap)

fig.show()
