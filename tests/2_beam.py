from utils import*
import importlib as imp
import wallpp.lattice as lat        ;imp.reload(lat)
import wallpp.plane_group as pg     ;imp.reload(pg)
import multislice.multi_2D as ms    ;imp.reload(ms)
plt.close('all')
path='../multislice/docs_fig/multi2D/'
opts='ME'

K = 5.008343065817388
keV   = cst.lam2keV(1/K)#200       # wavelength 200keV
Nx,Nz = 1,50,#10,2  # unit cells
ndeg  = 2**6       # number of pixels
npx   = Nx*ndeg    # number of pixels
pptype,ax,bz,angle = 'p1',2,2,90
Za = 2
pattern = np.array([[1,1,1]])
Nh = 20
rot = lambda t:np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])

a1,a2 = lat.get_lattice_vec(lat_type='rect',a=ax,b=bz)
b1,b2 = lat.reciprocal_lattice_2D(a1,a2)
Ai = np.array([0.1,0.25,0.26,0.27,1.5])
fj = lambda ti,j:np.sqrt(np.pi*Ai[j])*np.exp(-(np.pi*Ai[j]*K*np.sin(ti))**2)


if 'M' in opts:
    # p1 = pg.Wallpaper(pptype,ax,bz,angle,pattern,ndeg=ndeg)
    # potential = p1.get_potential_grid()
    n = 10
    x0,z0 = np.meshgrid(np.arange(-1,n),np.arange(n+1))
    x0,z0 = np.stack([x0*ax+pattern[0][0],z0*bz+pattern[0][1]])
    Xa = rot(1/n).dot(np.stack([x0.flatten(),z0.flatten()]))

    ax1 = np.sqrt(1**2+n**2)*ax
    bz1 = ax1

    x,z   = np.meshgrid(np.arange(n*ndeg)*ax1/(n*ndeg),np.arange(n*ndeg)*bz1/(n*ndeg))
    f = np.sum(np.array([pg.fv(np.stack([x.flatten(),z.flatten()]).T,X0,Za) for X0 in Xa.T]),axis=0)
    f = np.reshape(f,(n*ndeg,n*ndeg))

    dsp.stddisp(im=[x,z,f],
        scat=[Xa[0],Xa[1],15,'k'],)

    Nx,Nz = 1,500
    ms0 = ms.Multi2D([x[0],z.T[0],f],ax1,bz1,keV=keV,Nx=Nx,nz=Nz,
            dz=bz1/(2*n),eps=0.01,
            iZs=1,iZv=1)
    ms0.Bz_show(np.arange(-4,5)*n)
    ms0.Qz_show(slice(0,-1,10))

    if 'D' in opts:
        iM = -n*4
        fig,ax = ms0.Bz_show(iM)
        m = ms0.psi_qz[:,iM].max()/ms0.nx**2/ms0.dq

        fj = lambda ti,j:np.sqrt(np.pi*Ai[j])*np.exp(-(np.pi*Ai[j]*K*np.sin(ti))**2)
        pg.fv()

        kx = 1#4/ax
        ti   = np.arcsin(kx/ms0.k0)
        xi_g = np.pi/(ms0.sig*ms0.eps*fj(ti,Za))*ax*bz
        w_g  = ms0.k0*(1-np.cos(ti))
        F = np.sin(t/xi_g)**2
        plts = [[ms0.z,F,'b--']]
        dsp.stddisp(plts,ax=ax0,lw=2,labs=[r'$\AA(nm)$','$I$'])


if 'E' in opts:
    g,theta = 1,0.1
    # Ks = np.linspace(0,2,100)
    # F = lambda K:K*(1-np.sin(2*np.arcsin(g/(2*K))))-g*np.sin(theta)
    # dsp.stddisp([Ks,F(Ks),'b'],lw=2)
    K  = 1/(2*np.sin(theta))#1/cst.keV2lam(keV)
    # K = cst.keV2A(12)
    dt = 30*np.pi/180 #10deg each side
    t  = 3*np.pi/2+dt*np.linspace(-1,1,1000)

    nh,nk = 20,10
    h,k = np.meshgrid(np.arange(-nh,nh+1),np.arange(nk))
    X = h*b1[0]+k*b2[0]
    Z = h*b1[1]+k*b2[1]
    # ss = X.shape
    X,Z  = rot(-theta).dot(np.stack([X.flatten(),Z.flatten()]))
    # X,Z = X.reshape(ss)
    H = Nz*bz
    zeta = np.linspace(0,0.25/bz,100)
    Fz = lambda i : 1/(1.1*ax)*np.sinc(zeta*H)**2+i/ax

    plts=[]
    # plts +=[[Fz(i),zeta,'b--',''] for i in range(nh)]
    plts += [[K*np.cos(t),K*np.sin(t)+K,'r','']]
    scat  = [X,Z,15,'k']
    ax1 = np.sqrt(1**2+10**2)*ax
    bz1 = ax1

    dsp.stddisp(plts,scat=scat,labs=[r'$k_x(\AA^{-1})$',r'$k_z(\AA^{-1})$'],
        lw=2,xylims=[0,3,0,3])#,xyTicks=[1/ax1,1/bz1])


if 'K' in opts:
    z = ms0.z
    K = ms0.k0
    kx = np.arange(1,Nh)/ax
    ti = np.arcsin(kx/K)
    xi = K*(1-np.cos(ti))

    iM = 3
    m = ms0.psi_qz[:,iM+1].max()/ms0.nx**2/ms0.dq
    F = lambda i:fj(ti[i],Za)*np.sin(np.pi*z*xi[i])/(np.pi*xi[i])#*ms0.sig/(ax*bz)
    fig,ax0 = ms0.Bz_show(np.arange(1,Nh)*Nx,lw=2)
    cs = dsp.getCs('jet',Nh-1)
    plts = [[z,(F(i)/F(iM).max())**2*m,[cs[i],'--'],'$%d$' %i] for i in range(Nh-1)]
    dsp.stddisp(plts,ax=ax0,lw=2,labs=[r'$\AA(nm)$','$I$'])
    # plts  = [[t,fj(t,Za)**2,'g--','fe**2']]
