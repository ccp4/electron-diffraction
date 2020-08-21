from utils import*
import multislice.multi_2D as MS2D
import wallpp.plane_group as pg
import importlib as imp
imp.reload(MS2D)
imp.reload(pg)
plt.close('all')

keV = 50#200
ndeg = 2**8
eps = 0.01
#### Potential of a structure from wallpaper library
pptype,a,b,angle = 'p1',10,4,90
# pattern = np.array([[2,2,4]])
pattern = np.array([[2,2,1]])
p1      = pg.Wallpaper(pptype,a,b,angle,pattern,ndeg=ndeg)
pattern = p1.get_potential_grid()


Nx=1
iZv = np.inf
def base():
    mp0 = MS2D.Multi2D(pattern,a,b,keV,
            Nx=Nx,dz=b,nz=100,ppopt='',#XQZTP
            iZs=1,iZv=iZv,eps=eps)

    # mp0.Bz_show(iBs=np.arange(6)*Nx,lw=2,xylims=[0,400,0,1])
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
    mp1 = MS2D.Multi2D(pattern,a,b,keV,
            Nx=Nx,dz=0.1*b,nz=1000,ppopt='',#XQZTP
            iZs=1,iZv=iZv,eps=eps)
    # mp1.Bz_show(iBs=np.arange(6)*Nx,lw=2,xylims=[0,400,0,1])
    return mp1

def tilt(t=0.1):
    mp1 = MS2D.Multi2D(pattern,a,b,keV,tilt=t,
            Nx=Nx,dz=0.1*b,nz=1000,ppopt='',#XQZTP
            iZs=1,iZv=iZv,eps=eps)
    # mp1.Bz_show(iBs=np.arange(6)*Nx,lw=2,xylims=[0,400,0,1])
    return mp1

# mp0=base()
# mp1=small_thick()
tilts=np.linspace(0,0.2,50)
iBs = np.arange(5)
mp2,It = np.empty(tilts.size,dtype=object),np.zeros((tilts.size,iBs.size))
for i,t in enumerate(tilts):
    mp2[i]  = tilt(t)
    It[i,iBs] = mp2[i].getI()[iBs+1]
# plts=[[mp0.q,mp0.getI(),'b',''],[mp2.q,mp2.getI(),'r',r'$\theta_i=%.1f^{\circ}$' %mp2.tilt]]
cs = dsp.getCs('Spectral',iBs.size)
plts=[ [tilts,It[:,iB],cs[iB],'%d' %(iB+1)] for iB in iBs]
dsp.stddisp(plts)
