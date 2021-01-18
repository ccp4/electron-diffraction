from scipy.integrate import nquad,trapz,quad
import scipy.fftpack as fft
from scipy.signal import fftconvolve
from scipy.interpolate import interp1d
import importlib as imp
from utils import*
import utils.FourierUtils as fu     ;
import wallpp.plane_group as pg     ;imp.reload(pg)
import multislice.multi_2D as MS2D  ;imp.reload(MS2D)
imp.reload(dsp)
plt.close('all')
path = '../multislice/docs_fig/'


#parameters :
## structure
pptype,a,b,angle = 'p1',20,50,90
pattern = np.array([[10,0,4]])      #pattern
ndeg    = 2**8                      #sampling per unit cell
p1      = pg.Wallpaper(pptype,a,b,angle,pattern,ndeg=ndeg)
pattern = p1.get_potential_grid()

## multislice
ppopt=''       #P(Potential) T(Transmission) B(Beams) C(Chanel) Q(Psi_q2) X(Psi2)
ropts=''       #M(multislice) F(Fresnel)
dz = b
# angle = 1#0.005#1 #5*np.pi/180

mp0 = MS2D.Multi2D(pattern,a,b,keV=200,
        Nx=4,dz=dz,nz=1,ppopt='',#XQZTP
        iZs=1,iZv=1,sg=-1,eps=0.1,copt=1)
if 'T' in ppopt:mp0.Tz_show()

k0 = mp0.k0
angle=mp0.q.max()/k0
##### propagator
Dz = dz
dx = angle*Dz # 0.1305*dz
x,T  = mp0.x,mp0.T[0]
b0 = x.max()

#real space propagator
xx = np.linspace(-dx,dx,10000)
# cc,ss = np.cos(np.pi*k0/Dz*xx**2),np.sin(np.pi*k0/Dz*xx**2)
ss,cc = -k0/dz*np.cos(np.pi*k0/Dz*xx**2),k0/dz*np.sin(np.pi*k0/Dz*xx**2)
plts0 = [[xx,cc,'b','$P_r$'],[xx,ss,'r','$P_i$']]
dsp.stddisp(plts0,labs=['$x$','$Fresnel Propagator$'],lw=2,name=path+'fresnelP.svg',opt='ps')

k=2
P = fft.fft(np.exp(-1J*k*xx**2))
q = fft.fftfreq(xx.size,xx[1]-xx[0])
q=fft.fftshift(q)
P=fft.fftshift(P)
plts=[[q,P.real,'b'],[q,P.imag,'r']]
dsp.stddisp(plts,lw=2)
mp0.Pq_show(lw=2)
exit()

# transmission
x0,Tr,Ti = np.hstack([-dx,x,b0+dx]),np.hstack([1,T.real,1]),np.hstack([0,T.imag,0])
plts = [[x0,Tr,'b','$T_r$'],[x0,Ti,'r','$T_i$']]
# dsp.stddisp(plts,labs=['$x$','$Transmission$'],lw=2,name=path+'fresnelT.svg',opt='ps')
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
rr,ii = np.zeros(T.shape),np.zeros(T.shape)
ri,ir = np.zeros(T.shape),np.zeros(T.shape)
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
    rr[j] = trapz(TrGr(xx,x[j]),xx)
    ri[j] = trapz(TrGi(xx,x[j]),xx)
    ir[j] = trapz(TiGr(xx,x[j]),xx)
    ii[j] = trapz(TiGi(xx,x[j]),xx)
    psi_fr[j]=rr[j]-ii[j]+1J*(ri[j]+ir[j])
####real space
psi_cv = mp0.Psi_x
frr = psi_fr.real/psi_fr.real[0]
fri = psi_fr.imag/psi_fr.real[0]
ftr = psi_cv.real/psi_cv.real[0]
fti = psi_cv.imag/psi_cv.real[0]
plts0  = [[x,rr,'b'  ,'rr'],[x,ri,'r'  ,'ri']]
plts0 += [[x,ii,'b--','ii'],[x,ir,'r--','ir']]
plts = [[x,frr,'b','frr'],[x,fri,'r','fri']]
plts+= [[x,ftr,'b--','Re(fft)'],[x,fti,'r--','Im(fft)']]
dsp.stddisp(plts0,labs=['$x$','$\Psi(x)$'],lw=2,name=path+'fresnelTG.svg',opt='p')
dsp.stddisp(plts,labs=['$x$','$\Psi(x)$'],lw=2,name=path+'fresnelX.svg',opt='p')
##reciprocal space
# Fcv = Pq*Tq
# Ffr = fft.fft(psi_fr)
# plts = [[q,Ffr.real,'b','Frr'],[q,Ffr.imag,'r','Fri']]
# plts+= [[q,Fcv.real,'b--','fftr'],[q,Fcv.imag,'r--','ffti']]
# dsp.stddisp(plts,labs=['$x$','$\Psi(q_x)$'],lw=2)
