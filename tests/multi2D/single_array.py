import importlib as imp
from utils import*                  ; imp.reload(dsp)
import multislice.multi_2D as MS2D  ;imp.reload(MS2D)
plt.close('all')

keV  = 100  #200
nx,nz = 2**14,1000 #
eps   = 1#0.5#0.1 #
opts  = 'v'

ax,bz = 100, 4
Ai = np.array([0.1,0.25,0.26,0.27,1.5])
# fv = lambda X,X0,Za : V0*np.exp(-(np.linalg.norm(X-X0,axis=1)/Ai[Za])**2)

x0,z0 = np.linspace(0,ax,nx),np.linspace(0,bz,nz)


xa,za,Za = ax/2,bz/2,3
x,z = np.meshgrid(x0,z0)
r = np.sqrt((x-xa)**2+(z-za)**2)
fv = np.zeros((nz,nx))
print('..computing fv...')
# fv = np.exp(-(r/Ai[Za])**2)
r0,V0 = 0.5, 0.25
fv[np.abs(r)<r0] = V0



if 'v' in opts:
    print('..plotting fv...')
    dr=bz/2
    idx = np.abs(x0-ax/2)<dr
    xylims = np.hstack([za+dr*np.array([-1,1]),xa+dr*np.array([-1,1])])
    dsp.stddisp(im=[z[:,idx],x[:,idx],fv[:,idx]],labs=['$z$','$x$'],imOpt='c',axPos='V',xylims=xylims,opt='p')


Nz = 2     #nb atoms
pattern = [x0,z0,fv]
mp0 = MS2D.Multi2D(pattern,ax,bz,keV,
        Nx=1,dz=bz,nz=Nz,ppopt='',
        iZs=1,iZv=10,eps=eps)

mp0.Tz_show(iSz=slice(0,None,1),Vopt='VT',lw=2)
# mp0.Xxz_show()#iZs=1,iXs=2)
# mp0.Bz_show(iBs=np.arange(1,10),tol=1e-2,cmap='jet',lw=2)
mp0.Qz_show(iZs=1,lw=2,opts='S')
