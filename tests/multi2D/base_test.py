import importlib as imp
from utils import*
import multislice.multi_2D as MS2D  ;imp.reload(MS2D)
import wallpp.plane_group as pg     ;imp.reload(pg)
plt.close('all')
path='../multislice/docs_fig/multi2D/'

keV = 100#200
ndeg = 2**8
eps = 0.01
#### Potential of a structure from wallpaper library
pptype,a,b,angle = 'p1',10,4,90
# pattern = np.array([[2,2,4]])
pattern = np.array([[2,2,1]])
p1      = pg.Wallpaper(pptype,a,b,angle,pattern,ndeg=ndeg)
pattern = p1.get_potential_grid_p1()


Nx=1
iZv = np.inf
def base():
    mp0 = MS2D.Multi2D(pattern,a,b,keV=keV,
            Nx=Nx,dz=b,nz=100,ppopt='',#XQZTP
            TDS=False,
            iZs=1,iZv=10,opts='q',eps=eps)

    # mp0.Bz_show(iBs=np.arange(6)*Nx,lw=2,xylims=[0,400,0,1])
    # mp1.propagate(300,iZs=2)
    # mp0.Xxz_show(pOpt='pt',axPos=[0.14,0.12 ,0.85,0.85],
    #     opt='p',name='../docs_fig/multislice/multi_2D.png')
    # mp1.propagate(10,iZs=2)
    # mp1.save('test.pkl')
    # mp2 = MS2D.load('test.pkl')
    # mp2.Bz_show()
    mp0.Qz_show([0,-1],opts='S')
    mp0.Qxz_show(iZs=5)
    # mp0.Qz_show(slice(0,10,2))
    return mp0

def ms_sg():
    '''Difference between choosing propagator with + and with -'''
    mp0 = MS2D.Multi2D(pattern,a,b,keV,
            Nx=Nx,dz=b,nz=0,ppopt='',#XQZTP
            iZs=1,iZv=10,eps=eps)
    mp1 = MS2D.Multi2D(pattern,a,b,keV,
            Nx=Nx,dz=b,nz=0,ppopt='',#XQZTP
            iZs=1,iZv=10,eps=eps)
    mp0._set_propagator(sg=1)
    mp1._set_propagator(sg=-1)
    mp0.propagate(1000)
    mp1.propagate(1000)
    # plts=[[mp0.q,mp1.getI()-mp0.getI(),'b','']]#,[mp1.q,mp1.getI(),'r','-']]
    plts=[[mp0.z,mp0.getB(1),'b','+'],[mp1.z,mp1.getB(1),'r','-']]
    dsp.stddisp(plts)

def small_thick():
    mp1 = MS2D.Multi2D(pattern,a,b,keV,
            Nx=Nx,dz=0.1*b,nz=1000,ppopt='',#XQZTP
            iZs=1,iZv=iZv,eps=eps)
    # mp1.Bz_show(iBs=np.arange(6)*Nx,lw=2,xylims=[0,400,0,1])
    return mp1

def tilt(t=0.1,eps=eps):
    mpt = MS2D.Multi2D(pattern,a,b,keV,tilt=t,
            Nx=Nx,dz=0.1*b,nz=1000,ppopt='',#XQZTP
            iZs=1,iZv=1,eps=eps,v=0)
    # mp1.Bz_show(iBs=np.arange(6)*Nx,lw=2,xylims=[0,400,0,1])
    return mpt

def tilts_test(tilts=np.linspace(0,0.9,100),
    eps=0.01,nz=3000):
    mp2 = np.empty(tilts.size,dtype=object)
    for i,t in enumerate(tilts):
        print('theta=%.2f' %t)
        mp2[i] = MS2D.Multi2D(pattern,a,b,keV,tilt=-t,
                Nx=1,dz=0.1*b,nz=nz,ppopt='',#XQZTP
                iZs=10,iZv=100,eps=eps,sg=-1,v=1)
    return tilts,mp2

mp0=base()
# ms_sg()
# mp1=small_thick()
# mpt=tilt(t=0.1)
# ts,msts=tilts_test(tilts=np.linspace(0,0.2,51),eps=0.1,nz=10010)
# ts,msts=tilts_test(tilts=np.linspace(0,0.3,21)[6:8],eps=0.02,nz=50000)
# tilts_show(ts,msts,iBs=1,iZs=-1,name=path+'rocking_eps1.svg',opt='p')
# tilts_show(msts,iBs=np.arange(-3,3),iZs=-1,name=path+'rocking_eps1.svg',opt='p')
# MS2D.tilts_show(ts,msts,iBs=1,iZs=slice(0,10000,1000),name=path+'rocking_eps1.svg',opt='p')
# m0 = msts[25];m0.Bz_show(np.arange(-5,5))
# msts[10].Ewald_show(lw=2,xylims=[-0.01,0.4,-0.001,0.005])
# msts=tilts_test(tilts=np.linspace(0,0.2,51),eps=0.10,name=path+'rocking_eps2.svg',opt='ps')
