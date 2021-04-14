import importlib as imp
from scipy.integrate import nquad,trapz,quad
from scipy.interpolate import interp1d
from utils import*                  ;imp.reload(dsp)
import wallpp.plane_group as pg     ;imp.reload(pg)
import multislice.multi_2D as MS2D  ;imp.reload(MS2D)
plt.close('all')
path = '../multislice/docs_fig/'


## structure
pptype,a,b,angle = 'p1',50,5,90
pattern = np.array([[25,2,1]])   #pattern
ndeg    = 2**11                  #sampling per unit cell
p1      = pg.Wallpaper(pptype,a,b,angle,pattern,ndeg=ndeg)
pattern = p1.get_potential_grid()
eps=1      #transmission function exp(j*eps*sig*Vz)


def manual_convolution(x,T,k0,angle,Dz):
    dx,b0 = angle*Dz,x.max()

    #real space propagator
    xx = np.linspace(-dx,dx,10000)
    P  = -1J*np.exp(1J*np.pi/4)*k0/Dz*np.exp(1J*np.pi*k0/Dz*xx**2)
    cc,ss = P.real,P.imag
    plts0 = [[xx,cc,'b','$P_r$'],[xx,ss,'r','$P_i$']]
    dsp.stddisp(plts0,labs=['$x$','Fresnel Propagator'],lw=2,xylims=['x',-2,2],name=path+'fresnelP.svg',opt='ps')

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
    psi_fr = np.zeros(T.shape,dtype=complex)
    rr,ii = np.zeros(T.shape),np.zeros(T.shape)
    ri,ir = np.zeros(T.shape),np.zeros(T.shape)
    for j in range(0,x.size,1): #loop over observation points
        rr[j] = trapz(TrGr(xx,x[j]),xx)
        ri[j] = trapz(TrGi(xx,x[j]),xx)
        ir[j] = trapz(TiGr(xx,x[j]),xx)
        ii[j] = trapz(TiGi(xx,x[j]),xx)
        psi_fr[j]=rr[j]-ii[j]+1J*(ri[j]+ir[j])
    # plts0  = [[mp0.x,rr,'b'  ,'rr']      ,[mp0.x,ri,'r'  ,'ri']]
    # plts0 += [[mp0.x,ii,'b--','ii']      ,[mp0.x,ir,'r--','ir']]
    # dsp.stddisp(plts0,labs=['$x$','$\Psi(x)$'],lw=2,name=path+'fresnelTG.svg',opt='p')
    return psi_fr

#multislice 1 slice
print('...ms...')
mp0 = MS2D.Multi2D(pattern,a,b,keV=200,
        Nx=1,dz=b,nz=1,ppopt='',#XQZTP
        iZs=1,iZv=1,sg=-1,eps=eps,copt=1)
psi_cv = mp0.Psi_x

#manual convolution
print('...convolution...')
angle = 10*mp0.q.max()/mp0.k0/2
psi_fr = manual_convolution(x=mp0.x,T=mp0.T[0], k0=mp0.k0,angle=angle,Dz=b)


####real space
frr = psi_fr.real/psi_fr.real[0]
fri = psi_fr.imag/psi_fr.real[0]
ftr = psi_cv.real/psi_cv.real[0]
fti = psi_cv.imag/psi_cv.real[0]

plts   = [[mp0.x,frr,'b','frr']     ,[mp0.x,fri,'r','fri']]
plts  += [[mp0.x,ftr,'b--','Re(MS)'],[mp0.x,fti,'r--','Im(MS)']]
dsp.stddisp(plts,labs=['$x$','$\Psi(x)$'],lw=2,xylims=['x',20,30],name=path+'fresnelX.svg',opt='ps')
