from utils import*
import importlib as imp
import wallpp.plane_group as pg     ;imp.reload(pg)
import multislice.multi_2D as ms    ;imp.reload(ms)
import multislice.nearBragg as nb   ;imp.reload(nb)
plt.close('all')
path='../multislice/docs_fig/multi2D/'
opts='C'

# Problem
keV   = 200        # wavelength 200keV
Nx,Nz = 10,10,#10,2  # unit cells
ndeg  = 2**7       # number of pixels
npx   = Nx*ndeg    # number of pixels
pptype,ax,bz,angle = 'p1',10,5,90
Za = 2
pattern = np.array([[5,2.5,Za]])
eps=1.0
plts=[]


if 'I' in opts :
    sig_e = 0.5 #A
    S = Nx*ax
    Nslices = 1
    N = Nslices*Nz
    z,Ncoh = np.zeros((N,1)),np.ones((N,1))
    for i in range(1,Nz):
        za = np.array([2.5])
        Na = np.array([1])
        for iS in range(Nslices):
            j = i*Nslices+iS    #; print(j)
            z[j] = ax*i+za[iS]
            Ncoh[j]=Ncoh[j-1]*(1-sig_e*Na[iS]/S)

    plts = [z,Ncoh,'b']
    dsp.stddisp(plts,labs=['$z$','$I_{coh}$'],lw=2)

if 'C' in opts:
    # epss=[0.1] #[0.05,0.1,0.25,0.5,1,2]
    # for i,eps in enumerate(epss):
    p1 = pg.Wallpaper(pptype,ax,bz,angle,pattern,ndeg=ndeg)
    potential = p1.get_potential_grid()
    ms0 = ms.Multi2D(potential,ax,bz,keV=keV,Nx=Nx,nz=Nz,
            dz=bz,eps=eps,
            iZs=1,iZv=1)
    # ms0.Tz_show(                                         name=path+'multi2d_Z%d_pot%d_P.svg' %(Za,i),opt='s');plt.close('all')
    # ms0.Bz_show(np.arange(1,5)*Nx                       ,name=path+'multi2d_Z%d_pot%d_Z.svg' %(Za,i),opt='s');plt.close('all')
    # ms0.Xxz_show(iZs=slice(0,-1,10),iXs=slice(0,-1,4)   ,name=path+'multi2d_Z%d_pot%d_Q.png' %(Za,i),opt='s',pOpt='pt');plt.close('all')
    Ims = ms0.psi_qz[-1,:]; Ims/=ms0.nx**2*ms0.dq**2;
    # Im0=Ims[0];Ims[0]=0
    # Ims/=ms0.psi_qz[-1,:].max();
    plts+=[[ms0.q,Ims,'ro','MS']]

    qmax = 3.5#ms0.q.max()
    npx = 4096
    # nbG = nb.NearBragg(pattern.T,ax,bz,keV=keV,Nx=Nx,Nz=Nz,eps=0.1,
    #         npx=npx,qmax=qmax,method='Greens')
    # Inb = nbG.I/nbG.I.max()
    # Inb = nbG.I
    # Inb[int(Nx/2*ndeg-1)]=0;Inb/=Inb.max()
    # plts+=[[nbG.getQ(),Inb,'g-+','Greens']]

    nbFr = nb.NearBragg(pattern.T,ax,bz,keV=keV,Nx=Nx,Nz=Nz,eps=0.1,
            npx=npx,qmax=qmax,method='Fraunhofer')
    Ifr = nbFr.I
    Ifr /= nbFr.I.sum()
    #idx = ifr.argmax()
    # Ifr[idx]=0;Ifr/=Ifr.max()
    # Ifr = nbFr.I/nbFr.I.max()
    plts+=[[nbFr.getQ(),Ifr,'b','Fraunhofer']]

    dsp.stddisp(plts,labs=['q','I'])
