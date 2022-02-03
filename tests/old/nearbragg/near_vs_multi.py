from utils import*
import importlib as imp
import wallpp.plane_group as pg   ;imp.reload(pg)
import multislice.multi_2D as ms  ;imp.reload(ms)
import multislice.nearBragg as nb ;imp.reload(nb)
plt.close('all')
# path='../multislice/docs_fig/multi2D/MSvsNBstrong_'
path='../multislice/docs_fig/multi2D/'
opts='I'

# Problem
keV   = 200        # wavelength 200keV
Nx    = 10         # Transverse unit cells
Nz    = 100        #
ndeg  = 2**7       # number of pixels
pptype,ax,bz,angle,Za = 'p1',10,5,90,  2
pattern = np.array([[5,2.5,Za]])
eps = 0.1
opt = 'p'


if 'I' in opts :
    Nz = 1000
    sig_e = 0.05 #A
    Nslices = 1
    N = Nslices*Nz
    z,Ncoh = np.zeros((N,1)),np.ones((N,1))
    for i in range(1,Nz):
        za = np.array([2.5])
        Na = np.array([1])      #nb atoms
        for iS in range(Nslices):
            j = i*Nslices+iS    #; print(j)
            z[j] = bz*i+za[iS]
            Ncoh[j]=Ncoh[j-1]*(1-sig_e*Na[iS]/ax)

    le = ax*bz/sig_e; print(le) #Angstrom
    Npoi = np.exp(-z/le)
    z/=10
    plts  = [[z,Ncoh,'b','$P_{coh}$']]
    plts += [[z,Npoi,'b--','$Poisson$']]
    dsp.stddisp(plts,labs=[r'$z(nm)$','$P_{coh}$'],lw=2,
        name=path+'Pcohz.svg',opt='p')

if 'C' in opts:
    p1 = pg.Wallpaper(pptype,ax,bz,angle,pattern,ndeg=ndeg)
    potential = p1.get_potential_grid()
    ms0 = ms.Multi2D(potential,ax,bz,keV=keV,
        Nx=1,nz=Nz,dz=bz,eps=eps,
        iZs=1,iZv=1)
    # ms0.Xxz_show(iZs=slice(0,-1,10),iXs=slice(0,-1,4)   ,name=path+'multi2d_Z%d_pot%d_Q.png' %(Za,i),opt='s',pOpt='pt');plt.close('all')
    # ms0.Tz_show(Vopt='V',lw=2,
    #     name=path+'V.svg',opt=opt)
    # ms0.Qz_show(slice(0,11,2),'S',lw=2,xylims=['x',-2.5,2.5],
    #     name=path+'Iz.svg',opt=opt,axPos='T',setPos=1)
    # ms0.Bz_show(np.arange(1,6),lw=2,
    #     name=path+'Bz.svg',opt=opt)

    #### Near Bragg Nz=10
    qmax,npx  = 1.0,2500
    idx = (ms0.q<=qmax) &(ms0.q>=0)
    q0s = ms0.q[idx]
    # nbG = nb.NearBragg(pattern.T,ax,bz,keV=keV,Nx=Nx,Nz=Nz,eps=eps,
    #     method='Fresnel',q0s=q0s)
    #     # npx=npx,qmax=qmax,)
    nbG = nb.NearBragg(pattern.T,ax,bz,keV=keV,Nx=Nx,Nz=Nz,eps=eps,
        method='Dual',q0s=q0s,)

    # Ims10 = ms0.psi_qz[9,:].copy()/ms0.nx**2*ms0.dq**2;#Im0=Ims[0];
    # Ims10[0]=0
    # Ims10/=Ims10.max();
    # Inb10 = nbG.I.copy()
    # idx   = Inb10.argmax();Inb10[idx]=0;#remove forward propagation
    # Inb10/=Inb10.max() #normalize
    #
    # # plts+=[[nbFr.getQ(),Ifr,'b-','$NB_{kin}$']]
    # plts10 =[[nbG.getQ(),Inb10,'g-+','$NB_{kin}$']]
    # plts10+=[[ms0.q     ,Ims10,'ro','MS']]
    # dsp.stddisp(plts10,labs=['q','I'],lw=2,xylims=['x',0,2.01],
    #     # inset={'axpos':[0.55,0.3,0.4,0.4],'xylims':[1,2,-0.01,0.01],'ms':5,'lw':0.5},
    #     opt='p',name=path+'nz10.svg')
    #
    # #### Near Bragg Nz = 100
    # nbG.replicate(ax,bz,Nx,Nz=100)
    # nbG._Fresnel()

    Ims100 = ms0.psi_qz[-1,:].copy()/ms0.nx**2*ms0.dq**2;
    Ims100[0]=0
    Ims100/=Ims100.max();

    Inb100 = nbG.I.copy()
    idx = Inb100.argmax();Inb100[idx]=0;#remove forward propagation
    Inb100/=Inb100.max() #normalize

    plts100 =[[nbG.getQ(),Inb100,'g-+','$NB_{kin}$']]
    plts100+=[[ms0.q     ,Ims100,'ro','MS']]
    dsp.stddisp(plts100,labs=['$q(\AA^{-1})$','I'],lw=2,xylims=['x',0,1.01],
        # inset={'axpos':[0.55,0.3,0.4,0.4],'xylims':[1,2,-0.01,0.01],'ms':5,'lw':0.5},
        opt=opt,name=path+'I.svg')
