from utils import*
import scipy.signal as sig
import importlib as imp
import wallpp.lattice as lat        ;imp.reload(lat)
import wallpp.plane_group as pg     ;imp.reload(pg)
import multislice.multi_2D as ms    ;imp.reload(ms)
plt.close('all')
path='../multislice/docs_fig/multi2D/'
opts='M'

K = 5.008343065817388   # exact Bragg condition for 2nd beam
keV   = cst.lam2keV(1/K)#200       # wavelength 200keV
Nx,Nz = 1,2000
ndeg  = 2**6       # number of pixels
npx   = Nx*ndeg    # number of pixels
pptype,ax,bz,angle = 'p1',2,2,90
Za = 2
pattern = np.array([[1,1,1]])
Nh = 20


#### 2-beam setup
rot = lambda t:np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])
a1,a2 = lat.get_lattice_vec(lat_type='rect',a=ax,b=bz)
b1,b2 = lat.reciprocal_lattice_2D(a1,a2)
Ai = np.array([0.1,0.25,0.26,0.27,1.5])
fj = lambda ti,j:np.sqrt(np.pi*Ai[j])*np.exp(-(np.pi*Ai[j]*K*np.sin(ti))**2)


n = 10
x0,z0 = np.meshgrid(np.arange(-1,n),np.arange(n+1))
x0,z0 = np.stack([x0*ax+pattern[0][0],z0*bz+pattern[0][1]])
Xa = rot(1/n).dot(np.stack([x0.flatten(),z0.flatten()]))

ax1 = np.sqrt(1**2+n**2)*ax
bz1 = ax1

x,z   = np.meshgrid(np.arange(n*ndeg)*ax1/(n*ndeg),np.arange(n*ndeg)*bz1/(n*ndeg))
f = np.sum(np.array([pg.fv(np.stack([x.flatten(),z.flatten()]).T,X0,Za) for X0 in Xa.T]),axis=0)
f = np.reshape(f,(n*ndeg,n*ndeg))

if 'P' in opts:dsp.stddisp(im=[x,z,f],scat=[Xa[0],Xa[1],15,'k'],)

epss = np.arange(1,4)*0.01#array([0.01,0.05,0.1,0.2])

ms0 = []
for eps in epss:
    print(colors.red,eps,colors.black)
    ms0 += [ms.Multi2D([x[0],z.T[0],f],ax1,bz1,keV=keV,Nx=Nx,nz=Nz,
            dz=bz1/(2*n),eps=eps,sg=-1,
            iZs=1,iZv=np.inf)]

plts0,plts1 = [],[]
xi_0,xi_1   = np.zeros(epss.shape),np.zeros(epss.shape)
cs = dsp.getCs('hsv',epss.size)
for i,eps in enumerate(epss):
    Ib = ms0[i].Bz_show([10,20],out=1)
    idx0 = sig.find_peaks(-Ib[0],width=5)[0][0]
    idx1 = sig.find_peaks(-Ib[1],width=5,prominence=0.1)[0][0]
    xi_0[i] = ms0[i].z[idx0]
    xi_1[i] = ms0[i].z[idx1]
    plts0+=[[ms0[i].z,Ib[0],cs[i],r'$\epsilon=%.2f$' %eps]]
    plts1+=[[ms0[i].z,Ib[1],cs[i],r'$\epsilon=%.2f$' %eps]]
    plts0+=[[xi_0[i],Ib[0][idx0],[cs[i],'s'],'']]
    plts1+=[[xi_1[i],Ib[1][idx1],[cs[i],'s'],'']]
plts2 = [[epss,xi_0,'b-o','$g_1$'],[epss,xi_1,'g-o','$g_2$']]
dsp.stddisp(plts0 ,labs=[r'$z(\AA)$',r'$I$'],lw=2,name=path+'2_beam_I1.svg' ,opt='ps', xylims=[0,200,-0.0001,0.005])
dsp.stddisp(plts1 ,labs=[r'$z(\AA)$',r'$I$'],lw=2,name=path+'2_beam_I2.svg' ,opt='ps')
dsp.stddisp(plts2 ,labs=[r'Potential strength $V_{\epsilon}$',r'$\zeta(\AA)$'],lw=2,name=path+'2_beam_xi.svg' ,opt='ps')
