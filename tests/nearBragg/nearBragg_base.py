from utils import*
import importlib as imp
import multislice.nearBragg as nb   ;imp.reload(nb)
import wallpp.plane_group as pg     ;imp.reload(pg)
plt.close('all')
path='../../nearBragg/figures/'
dat='dat/'


Natoms,dim = 3,2
ax,bz = 5,3        # lattice constants A
lam   = 0.05       # wavelength A
npx   = 4096       # number of pixels
Nx,Nz = 20,10,#10,2  # unit cells
qmax = 4
opts='H' #HF

if 1:
    # x,z,Za = [0],[0],[1]
    x,z,Za = [0.25,0.75,0.5], [0.25,0.25,0.5], [1,1,2]
    # X   = np.random.rand(Natoms,dim)
    # Za  = np.random.randint(len(Zs),size=Natoms)
    np.save(dat+'atoms.npy',[x,z,Za])

pattern = np.load(dat+'atoms.npy')

nbG = nb.NearBragg(pattern,ax,bz,lam=lam,npx=npx,Nx=Nx,Nz=Nz,qmax=qmax,method='Greens')

if 'H' in opts:
    nbH = nb.NearBragg(pattern,ax,bz,lam=lam,npx=npx,Nx=Nx,Nz=Nz,qmax=qmax,method='Holton',path='dat/')
    objs=[(nbH,'k--+','Holton'),(nbH,'k-','Greens')]
    nb.display(objs,('q','I'),name=path+'comparisonHolton.svg',opt='sp')
    nbH.Pattern_show(name=path+'comparisonHoltonPattern.svg',opt='sp')

if 'F' in opts:
    nbFf = nb.NearBragg(pattern,ax,bz,lam=lam,npx=npx,Nx=Nx,Nz=Nz,qmax=qmax,method='Fraunhofer')
    nbFr = nb.NearBragg(pattern,ax,bz,lam=lam,npx=npx,Nx=Nx,Nz=Nz,qmax=qmax,method='Fresnel')
    # objs=[(nbG,'k-','Greens'),(nbFr,'k--o','Fresnel'),(nbFf,'k--x','Fraunhofer')]
    # nb.display(objs,('q','I'))
    m=nbG.I.max()
    plts =[[nbFf.getQ(),nbG.I/m,'g','Greens']]
    # plts+=[[nbFf.getQ(),nbFf.I/m,'b','Fresnel']]
    plts+=[[nbFf.getQ(),nbFr.I/m,'r','Fraunhofer']]
    dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$','$I$'],xylims=[-2,2,0,1],
        name=path+'path_length.svg',opt='sp')

    plts=[[nbFf.getQ(),np.log10(abs(nbG.I-nbFr.I)/nbG.I),'r','']]
    labs = [r'$q(\AA^{-1})$','$\log_{10} |I_{fr}-I_{Gr}/I_Gr|$']
    dsp.stddisp(plts,labs=labs,xylims=[-2,2,-10,5],xyTicks=[0.5,2],
        name=path+'path_length_diff.svg',opt='ps')
