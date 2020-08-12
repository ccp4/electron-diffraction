from utils import*
import multislice.multi_2D as MS2D
import wallpp.plane_group as pg
import importlib as imp
imp.reload(MS2D)
imp.reload(pg)
plt.close('all')

ndeg = 2**8
#get potential of a structure from wallpaper library
pptype,a,b,angle = 'p1',10,4,90
pattern = np.array([[2,2,4]])
p1      = pg.Wallpaper(pptype,a,b,angle,pattern,ndeg=ndeg)
pattern = p1.get_potential_grid()


Nx=1
def base():
    mp0 = MS2D.Multi2D(pattern,a,b,keV=200,
            Nx=Nx,dz=b,nz=100,ppopt='',#XQZTP
            iZs=1,iZv=1,eps=0.1)

    mp0.Bz_show(iBs=np.arange(6)*Nx,lw=2,xylims=[0,400,0,1])
    # mp1.propagate(300,iZs=2)
    # mp0.Xxz_show(pOpt='pt',axPos=[0.14,0.12 ,0.85,0.85],
    #     opt='p',name='../docs_fig/multislice/multi_2D.png')
    # mp1.propagate(10,iZs=2)
    # mp1.save('test.pkl')
    # mp2 = MS2D.load('test.pkl')
    # mp2.Bz_show()
    # mp0.Qz_show(slice(0,10,2))
    return mp0

def small_thick():
    mp1 = MS2D.Multi2D(pattern,a,b,keV=200,
            Nx=Nx,dz=0.1*b,nz=1000,ppopt='',#XQZTP
            iZs=1,iZv=1,eps=0.1)
    mp1.Bz_show(iBs=np.arange(6)*Nx,lw=2,xylims=[0,400,0,1])
    return mp1

mp0=base()
mp1=small_thick()
