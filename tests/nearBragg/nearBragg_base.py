from utils import*
from scipy.signal import find_peaks
import importlib as imp
import multislice.nearBragg as nb   ;imp.reload(nb)
import wallpp.plane_group as pg     ;imp.reload(pg)
plt.close('all')
path=dsp.get_figpath(__file__,'/../../nearBragg/figures/')
dat='dat/'

keV=200
Natoms,dim = 3,2
ax,bz = 10,5        # lattice constants A
npx   = 4096       # number of pixels
Nx,Nz = 10,100,#10,2  # unit cells
qmax = 4
opts='F' #HFld
eps=0.1

if 1:
    x,z,Za = [5],[2.5],[2]
    # pattern = np.array([[5,2.5,Za]])
    # x,z,Za = [0.25,0.75,0.5], [0.25,0.25,0.5], [1,1,2]
    # X   = np.random.rand(Natoms,dim)
    # Za  = np.random.randint(len(Zs),size=Natoms)
    np.save(dat+'atoms.npy',[x,z,Za])

pattern = np.load(dat+'atoms.npy')

nbG = nb.NearBragg(pattern,ax,bz,keV=keV,npx=npx,Nx=Nx,Nz=Nz,qmax=qmax,method='Greens')
plts = [[nbG.getQ(),nbG.I,'g-+','Greens']]

if 'H' in opts:
    nbH = nb.NearBragg(pattern,ax,bz,keV=keV,npx=npx,Nx=Nx,Nz=Nz,qmax=qmax,method='Holton',path='dat/')
    objs=[(nbH,'k--+','Holton'),(nbH,'k-','Greens')]
    nb.display(objs,('q','I'),name=path+'comparisonHolton.svg',opt='sp')
    nbH.Pattern_show(name=path+'comparisonHoltonPattern.svg',opt='sp')

if 'F' in opts:
    nbFf = nb.NearBragg(pattern,ax,bz,keV=keV,eps=eps,npx=npx,Nx=Nx,Nz=Nz,qmax=qmax,method='Fraunhofer')
    # nbFr = nb.NearBragg(pattern,ax,bz,keV=keV,eps=eps,npx=npx,Nx=Nx,Nz=Nz,qmax=qmax,method='Fresnel')
    # objs=[(nbG,'k-','Greens'),(nbFr,'k--o','Fresnel'),(nbFf,'k--x','Fraunhofer')]
    # nb.display(objs,('q','I'))
    nbG.I/=nbG.I.max()
    # nbFr.I/=nbFr.I.max()
    nbFf.I/=nbFf.I.max()

    if 'd' in opts :
        tol = 5e-5
        # idx = nbG.I>1e-3
        idx = find_peaks(nbG.I,height=tol,prominence=tol,distance=10)[0]
        plts += [[nbG.getQ()[idx],nbG.I[idx],'gs','']]
    # plts+=[[nbFr.getQ(),nbFr.I,'r','Fresnel']]
    plts+=[[nbFf.getQ(),nbFf.I,'b','Fraunhofer']]
    dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$','$I$'],xylims=[-2,2,0,1],
        name=path+'path_length.svg',opt='p')

    #plot differences
    if 'd' in opts:
        dI = 1e-15*np.ones(nbG.I.shape)
        dI[idx] = abs((nbG.I-nbFr.I)/nbG.I)[idx]
        lablog,xylims,xyTicks='',[-2,2,0,1],[0.5,0.1]
        if 'l' in opts:
            dI =np.log10(dI)
            lablog,xylims,xyTicks='\log_{10}',[-2,2,-10,5],[0.5,2]
        plts=[[nbFf.getQ(),dI,'r','']]
        labs = [r'$q(\AA^{-1})$','$%s |I_{fr}-I_{Gr}/I_Gr|$' %lablog]
        dsp.stddisp(plts,labs=labs,xylims=xylims,xyTicks=xyTicks,
            name=path+'path_length_diff.svg',opt='p')
