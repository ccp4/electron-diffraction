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

epss = np.arange(1,6,2)*0.01#array([0.01,0.05,0.1,0.2])

ms0 = []
for eps in epss:
    print(colors.red,eps,colors.black)
    ms0 += [ms.Multi2D([x[0],z.T[0],f],ax1,bz1,keV=keV,Nx=Nx,nz=Nz,
            dz=bz1/(2*n),eps=eps,
            iZs=1,iZv=np.inf)]

plts,xi_g,xi_1 = [],np.zeros(epss.shape),np.zeros(epss.shape)
cs = dsp.getCs('hsv',epss.size)
for i,eps in enumerate(epss):
    Ib = ms0[i].Bz_show([-20,-10],out=1)
    idx0 = sig.find_peaks(-Ib[0],width=5,prominence=0.1)[0][0]
    idx1 = sig.find_peaks(-Ib[1],width=5)[0][0]
    xi_g[i] = ms0[i].z[idx0]
    xi_1[i] = ms0[i].z[idx1]
    plts+=[[ms0[i].z,-Ib[0],cs[i],r'$\epsilon=%.2f$' %eps]]
    plts+=[[ms0[i].z,-Ib[1],[cs[i],'--'],'']]
    plts+=[[xi_g[i],-Ib[0][idx0],[cs[i],'s'],'']]
    plts+=[[xi_1[i],-Ib[1][idx1],[cs[i],'s'],'']]
plts2 = [[epss,xi_g,'b-o','$g$'],[epss,xi_1,'g-o','$2$']]
dsp.stddisp(plts ,lw=2,name=path+'2_beam_I.svg' ,opt='ps')
dsp.stddisp(plts2,lw=2,name=path+'2_beam_xi.svg',opt='ps')
