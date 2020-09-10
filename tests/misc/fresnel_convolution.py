from scipy.integrate import nquad,trapz,quad
import scipy.fftpack as fft
from scipy.signal import fftconvolve
from scipy.interpolate import interp1d
from utils import*
import utils.FourierUtils as fu
import wallpp.plane_group as pg
import multislice.multi_2D as MS2D
import importlib as imp
imp.reload(pg)
imp.reload(dsp)
plt.close('all')
path = '../multislice/docs_fig/'

#parameters :
## structure
pptype,a,b,angle = 'p1',20,50,90
pattern = np.array([[10,0,5]])      #pattern
ndeg    = 2**9                      #sampling per unit cell
p1      = pg.Wallpaper(pptype,a,b,angle,pattern,ndeg=ndeg)
pattern = p1.get_potential_grid()

## multislice
ppopt=''       #P(Potential) T(Transmission) B(Beams) C(Chanel) Q(Psi_q2) X(Psi2)
ropts=''       #M(multislice) F(Fresnel)
dz = b
# angle = 1#0.005#1 #5*np.pi/180

mp0 = MS2D.multi2D(pattern,a,b,keV=200,
        Nx=1,dz=dz,nz=1,ppopt='',#XQZTP
        iZs=1,iZv=1)
if 'T' in ppopt:mp0.Tz_show()

k0 = mp0.k0
angle=mp0.q.max()/k0
##### propagator
Dz = dz
dx = 0.3*angle*Dz # 0.1305*dz
x,T  = mp0.x,mp0.T[0]
b0 = x.max()

xx = np.linspace(-dx,dx,1000)
cc,ss = np.cos(np.pi*k0/(2*Dz)*xx**2),np.sin(np.pi*k0/(2*Dz)*xx**2)

plts0 = [[xx,cc,'b','$P_r$'],[xx,ss,'r','$P_i$']]
dsp.stddisp(plts0,labs=['$x$','$Fresnel Propagator$'],lw=2,name=path+'fresnelP.svg',opt='ps')

# transmission
x0,Tr,Ti = np.hstack([-dx,x,b0+dx]),np.hstack([1,T.real,1]),np.hstack([0,T.imag,0])
plts = [[x0,Tr,'b','$T_r$'],[x0,Ti,'r','$T_i$']]
dsp.stddisp(plts,labs=['$x$','$Transmission$'],lw=2,name=path+'fresnelT.svg',opt='ps')
Tr,Ti = interp1d(x0,Tr),interp1d(x0,Ti)

### integrands
TrGr = lambda x,xp:Tr(xp-x)*cc
TrGi = lambda x,xp:Tr(xp-x)*ss
TiGr = lambda x,xp:Ti(xp-x)*cc
TiGi = lambda x,xp:Ti(xp-x)*ss
# xx = np.linspace(-dx,dx,1000)
# TrGr = lambda j:T.real*np.cos(np.pi*k0/dz*(x-x[j])**2)
# TrGi = lambda j:T.real*np.sin(np.pi*k0/dz*(x-x[j])**2)
# TiGr = lambda j:T.imag*np.cos(np.pi*k0/dz*(x-x[j])**2)
# TiGi = lambda j:T.imag*np.sin(np.pi*k0/dz*(x-x[j])**2)
psi_fr = np.zeros(T.shape,dtype=complex)
for j in range(0,x.size,1): #loop over observation point
    # plts0 = [
    #     [xx,TrGi(xx,x[j]),'b','$T_rG_r$'],[xx,TrGi(xx,x[j]),'m','$T_rG_i$'],
    #     [xx,TiGi(xx,x[j]),'c','$T_iG_i$'],[xx,TiGr(xx,x[j]),'r','$T_iG_r$'],
    #     ]
    # dsp.stddisp(plts0,labs=['$x$','$Fresnel$'],lw=2,title='x=%.2f' %x[j])

    # rr = quad(Trfr,0,b0,args=(x[j],Dz))[0]
    # ri = quad(Trfi,0,b0,args=(x[j],Dz))[0]
    # ir = quad(Tifr,0,b0,args=(x[j],Dz))[0]
    # ii = quad(Tifi,0,b0,args=(x[j],Dz))[0]
    # rr,ri= 0,0
    # ir = quad(TiGr,0,b0,args=(x[j],Dz))[0]
    # ii = quad(TiGi,0,b0,args=(x[j],Dz))[0]
    # rr = trapz(TrGr(j),x)
    # ri = trapz(TrGi(j),x)
    # ir = trapz(TiGr(j),x)
    # ii = trapz(TiGi(j),x)
    rr = trapz(TrGr(xx,x[j]),xx)
    ri = trapz(TrGi(xx,x[j]),xx)
    ir = trapz(TiGr(xx,x[j]),xx)
    ii = trapz(TiGi(xx,x[j]),xx)
    psi_fr[j]=rr-ii+1J*(ri+ir)
####real space
psi_cv = mp0.Psi_x
frr = psi_fr.real/psi_fr.real[0]
fri = psi_fr.imag/psi_fr.real[0]
ftr = psi_cv.real/psi_cv.real[0]
fti = psi_cv.imag/psi_cv.real[0]
plts = [[x,frr,'b','frr'],[x,fri,'r','fri']]
plts+= [[x,ftr,'b--','Re(fft)'],[x,fti,'r--','Im(fft)']]
dsp.stddisp(plts,labs=['$x$','$\Psi(x)$'],lw=2,name=path+'fresnelX.svg',opt='ps')
##reciprocal space
# Fcv = Pq*Tq
# Ffr = fft.fft(psi_fr)
# plts = [[q,Ffr.real,'b','Frr'],[q,Ffr.imag,'r','Fri']]
# plts+= [[q,Fcv.real,'b--','fftr'],[q,Fcv.imag,'r--','ffti']]
# dsp.stddisp(plts,labs=['$x$','$\Psi(q_x)$'],lw=2)
