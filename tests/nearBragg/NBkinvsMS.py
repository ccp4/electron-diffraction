from utils import*
from scipy.signal import find_peaks
import importlib as imp
import wallpp.plane_group as pg   ;imp.reload(pg)
import multislice.nearBragg as nb ;imp.reload(nb)
import multislice.multi_2D as ms  ;imp.reload(ms)
plt.close('all')
path='../../multislice/docs_fig/multi2D/NBkinvsMS'
opts='I'

# Problem
keV   = 200        # wavelength 200keV
Nx    = 10         # Number of transverse unit cells
Nz    = 30          #
ndeg  = 2**8       # number of pixels
pptype,ax,bz,angle,Za = 'p1',5,10,90,  2
pattern = np.array([[2.5,5,Za]])
eps = 0.1
opt='p'

### Structure definition
p1 = pg.Wallpaper(pptype,ax,bz,angle,pattern,ndeg=ndeg)
potential = p1.get_potential_grid()


### Run multislice
ms0 = ms.Multi2D(potential,ax,bz,keV=keV,
    Nx=1,nz=Nz,dz=bz,eps=eps,
    iZs=1,iZv=1)
dq=ms0.q[1]-ms0.q[0]
# Ims = ms0.psi_qz[-1,:].copy()/(ms0.nx**2*ms0.dq**2)
# print('MS Inorm:%.2f' %(Ims.sum()*ms0.dq))


### Run near bragg
npts=2049 #4097 #61
q0s = np.linspace(-3,3,npts)
nbG = nb.NearBragg(pattern.T,ax,bz,keV=keV,Nx=Nx,Nz=Nz,eps=eps,
    method='Greens2',q0s=q0s,iZv=1)
# nbG.F_show(qopt=0,lw=2)
# nbG.stat_show(lw=2,opt=opt,name=path+'_proba.svg')

#postprocessing
# Ncoh,Nkin = np.sum(nbG.S[0,:]),np.sum(nbG.S[1,:])
# dq = q0s[1]-q0s[0]
# Nscat = nbG.I.sum()*dq
# Ikin = nbG.I.copy()/Nscat*Nkin
# I = Ikin.copy()
# I[int(q0s.size/2)]+=Ncoh/dq
# print('Ncoh=%.4f, Nkin=%.4f,Nscat=%.4f, norm(I)=%.2f' %(Ncoh,Nkin,Nscat,np.sum(I)*dq))

#remove origin and normalize with highest peak
Ims100 = ms0.psi_qz[-1,:].copy()/ms0.nx**2*ms0.dq**2;
Ims100[0]=0
Ims100/=Ims100.max();


Inb100 = nbG.I.copy()
Inb100/=Inb100.max()
idx,props = find_peaks(Inb100, height=0.7, threshold=None, distance=10, prominence=0.75)#, width=None, wlen=None, rel_height=0.5, plateau_size=None)
hh = props['peak_heights'];hh[hh==1]=0
idx0 = idx[hh.argmax()]
Inb100/=Inb100[idx0]


pltsI=[]
# pltsI = [[q0s,Ikin,'g--','$I_{kin}$']]
pltsI+= [[q0s  ,Inb100 ,'g-','$NB$']]
pltsI+= [[ms0.q,Ims100 ,'bo','$MS$']]
dsp.stddisp(pltsI,labs=['$q(\AA^{-1})$','$I(q)$'],lw=2,
    name=path+'_I.svg',opt='p')
