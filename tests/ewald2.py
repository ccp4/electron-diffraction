from utils import*
import scipy.signal as sig
import importlib as imp
import wallpp.lattice as lat        ;imp.reload(lat)
import wallpp.plane_group as pg     ;imp.reload(pg)
import multislice.multi_2D as ms    ;imp.reload(ms)
plt.close('all')
path='../multislice/docs_fig/multi2D/'
opts='M'

keV   = 200
Nx,Nz = 1,50,#10,2  # unit cells
ndeg  = 2**5       # number of pixels
npx   = Nx*ndeg    # number of pixels
pptype,ax,bz,angle = 'p1',2,2,90
Za = 2
pattern = np.array([[1,1,1]])
Nh = 20

p1 = pg.Wallpaper(pptype,ax,bz,angle,pattern,ndeg=ndeg)
potential = p1.get_potential_grid()


# #### 2-beam setup
# rot = lambda t:np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])
# a1,a2 = lat.get_lattice_vec(lat_type='rect',a=ax,b=bz)
# b1,b2 = lat.reciprocal_lattice_2D(a1,a2)
# Ai = np.array([0.1,0.25,0.26,0.27,1.5])
# fj = lambda ti,j:np.sqrt(np.pi*Ai[j])*np.exp(-(np.pi*Ai[j]*K*np.sin(ti))**2)
#
# n = 10
# x0,z0 = np.meshgrid(np.arange(-1,n),np.arange(n+1))
# x0,z0 = np.stack([x0*ax+pattern[0][0],z0*bz+pattern[0][1]])
# Xa = rot(1/n).dot(np.stack([x0.flatten(),z0.flatten()]))
#
# ax1 = np.sqrt(1**2+n**2)*ax
# bz1 = ax1
#
# x,z   = np.meshgrid(np.arange(n*ndeg)*ax1/(n*ndeg),np.arange(n*ndeg)*bz1/(n*ndeg))
# f = np.sum(np.array([pg.fv(np.stack([x.flatten(),z.flatten()]).T,X0,Za) for X0 in Xa.T]),axis=0)
# f = np.reshape(f,(n*ndeg,n*ndeg))
#
# if 'P' in opts:dsp.stddisp(im=[x,z,f],scat=[Xa[0],Xa[1],15,'k'],)

epss = np.arange(1,6)*0.01#array([0.01,0.05,0.1,0.2])
Nx,Nz = 1,2000
ms0 = []
for eps in epss:
    print(colors.red,eps,colors.black)
    ms0 += [ms.Multi2D(potential,ax,bz,keV=keV,Nx=Nx,nz=Nz,
            dz=bz/2,eps=eps,
            iZs=1,iZv=np.inf)]

plts,xi_g = [],np.zeros(epss.shape)
cs = dsp.getCs('hsv',epss.size)
for i,eps in enumerate(epss):
    Ib = ms0[i].Bz_show([1],out=1)[0]
    idx = sig.find_peaks(-Ib,width=5)[0][0]
    xi_g[i] = ms0[i].z[idx]
    plts+=[[ms0[i].z,-Ib,cs[i],r'$\epsilon=%.2f$' %eps]]
    plts+=[[xi_g[i],-Ib[idx],[cs[i],'s'],'']]

dsp.stddisp(plts,lw=2)
dsp.stddisp([epss,xi_g,'b-o'],lw=2)
