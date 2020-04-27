from utils import*
from mpl_toolkits.mplot3d import Axes3D
from scattering_factors import*
figpath=get_figpath(__file__,"/figures/")


def structure_factor3D(pattern,lat_vec,hkl=None,hklMax=10,sym=1,v='q'):
    '''Computes structure factor in 3D from :
    - `pattern` : Nx4 array - N atoms with each row : x,y,z,Z
    - `lat_vec` : 3x3 array - reciprocal lattice vectors
    - `hkl`     : list of 3 miller indices h,k,l each as 3Ndarray
        - `hklMax`  : int - max miller index in each direction from -hklMax to hklMax
    '''
    #unpack
    if not hkl : hkl = get_miller3D(hklMax,sym)
    hx,ky,lz = hkl
    ra,fa = pattern[:,:3],pattern[:,3]
    atoms = list(np.array(np.unique(fa),dtype=int));#print(atoms)
    # get scattering factor
    b1,b2,b3 = lat_vec
    k_x,k_y,k_z = hx*b1[0]+ky*b2[0]+lz*b3[0],hx*b1[1]+ky*b2[1]+lz*b3[1],hx*b1[2]+ky*b2[2]+lz*b3[2]
    q = np.sqrt(k_x**2+k_y**2+k_z**2)/(2*pi)
    q,fq = get_elec_atomic_factors(atoms,q)
    if 'q' in v:qmax=q.max();print('qmax=%.4f A^-1\nmax_res=%.4f A' %(qmax,1/qmax))
    #structure factor
    Fhkl,n_atoms = np.zeros(hx.shape,dtype=complex),len(atoms)
    for i,atom in zip(range(n_atoms),atoms):
        F_i = np.zeros(hx.shape,dtype=complex)
        idx = fa==atom
        for ri in ra[idx,:]:
            F_i += np.exp(2*pi*1J*(ri[0]*hx+ri[1]*ky+ri[2]*lz))
        Fhkl += F_i*fq[i]
    return hkl,Fhkl

def structure_factor2D(pattern,lat_vec,hk=None,hkMax=10,sym=1,v='q'):
    #unpack
    if not hk : hk = get_miller2D(hkMax,sym)
    hx,ky = hk
    ra,fa = pattern[:,:2],pattern[:,2]
    atoms = list(np.array(np.unique(fa),dtype=int));#print(atoms)
    #scattering factor
    b1,b2 = lat_vec
    k_x,k_y = hx*b1[0]+ky*b2[0],hx*b1[1]+ky*b2[1]
    q = np.sqrt(k_x**2+k_y**2)/(2*pi)
    if 'q' in v : qmax=q.max();print('qmax=%.4f A^-1\nmax_res=%.4f A' %(qmax,1/qmax))
    #compute structure factor
    Fhl,n_atoms = np.zeros(h.shape,dtype=complex),len(atoms)
    for i,atom in zip(range(n_atoms),atoms):
        F_i = np.zeros(h.shape,dtype=complex)
        idx = fa==atom
        for ri in ra[idx,:]:
            F_i += np.exp(2*pi*1J*(ri[0]*hx+ri[1]*ky))
        Fhl += F_i*fq[i]
    return Fhl

##########################################################################
###def : display
def plot_structure3D(hkl,Fhkl,**kwargs):
    '''scatter plot of the intensity'''
    hx,ky,lz = hkl
    hx,ky,lz = hx.flatten(),ky.flatten(),lz.flatten()
    S        = np.abs(Fhkl.flatten())**2
    c        = np.log10(S+1)
    N = hx.max()
    lims = [0,N,0,N,0,N]
    stddisp(scat=[hx,ky,lz,c],imOpt='c',rc='3d',legOpt=0,xylims=lims,**kwargs)
    #fig = plt.figure()
    #ax = plt.subplot(111,projection='3d')
    #c = ax.scatter(hx,ky,lz,c=np.log10(S+1),cmap='jet')
    #cb=fig.colorbar(c,ax=ax)

def plot2Dcutplane(Shkl,n='100',title='auto',**kwargs):
    '''Not working for arbitrary cut yet '''
    N = Shkl.shape[0]
    if title=='auto' : title='Silicon([%s]) $S_{hk}$' %n
    if   n=='100' : Scut = Shkl[0,:,:]
    elif n=='110' : Scut = Shkl[np.identity(N,dtype=bool),:]
    stddisp(im=Scut,imOpt='c',legOpt=0,title=title,**kwargs)

def plot1Dcutline(Shkl,u='001',lat_vec=None,dopt='fq',**kwargs):
    '''
    u : str cutline
    dopt:f(fe(q)),q(q units)
    '''
    N,labx = Shkl.shape[0],'$index$'
    ql = np.arange(N)
    if u == '001' :
        Scut = Shkl[0,0,:]
    else :
        Scut = Shkl[0,0,:]
    plts = [[ql,Scut,'bs-','$S_{%s}$' %u]]
    if 'q' in dopt or 'f' in dopt :
        b1,b2,b3 = lat_vec
        ql *= b3[2]/(2*pi)
        labx = '$q(A^{-1})$'
    if 'f' in dopt:
        #k_x,k_y,k_z = hx*b1[0]+ky*b2[0]+lz*b3[0],hx*b1[1]+ky*b2[1]+lz*b3[1],hx*b1[2]+ky*b2[2]+lz*b3[2]
        #q = np.sqrt(k_x**2+k_y**2+k_z**2)/(2*pi)
        qz = np.linspace(0,N,1000)*b3[2]/(2*pi)
        qz,fq = get_elec_atomic_factors([14],qz)
        plts += [[qz,64*fq[0]**2,'g--','$64f_e^2$'],[qz,32*fq[0]**2,'c--','$32f_e^2$']]
    stddisp(plts,labs=[labx,'$S_q(A^2)$'],**kwargs)

##########################################################################
##### def: misc
def get_miller3D(hklMax=10,sym=1):
    '''hklMax : Number of Beams such that Fourier components = 2^(N+1) - 1'''
    Nhkl = range(-hklMax*sym,hklMax+1)
    hkl = np.meshgrid(Nhkl,Nhkl,Nhkl)
    return hkl
def get_miller2D(hkMax=10,sym=1):
    Nhk = range(-hkMax*sym,hkMax+1)
    hk = np.meshgrid(Nhk,Nhk)
    return hk

##########################################################################
#### def : test
##########################################################################
def test_structure_factor2D(N=5):
    pptype,a,b,angle='p1',50,50,90      #A
    pattern = np.array([[0.5,0.5,1]])

    p1=pg.Wallpaper(pptype,a,b,angle,pattern=pattern,gen=False,pOpt='')
    b1,b2 = p1.get_reciprocal_lattice_2D()

    hl = np.arange(-2**N+1,2**N)
    h,l = np.meshgrid(hl,hl)
    Fhl = get_structure_factor2D(pattern,h,l,[b1,b2])

    qx,qy = (h*b1[0]+l*b2[0])/(2*pi),(h*b1[1]+l*b2[1])/(2*pi)
    fig,(axM,axP) = create_fig('12',pad=0.5,rc='12')
    axM.pcolor(np.abs(Fhl)**2);
    axP.pcolor(qx,qy,np.abs(Fhl)**2);
    #axP.pcolor(kx,ky,np.angle(Fhl));
    fig.show()
    return h,l,qx,qy,Fhl

def test_structure_factor3D(N=5,sym=False):
    from crystals import Crystal
    Si = Crystal.from_database('Si'); #vol     = Si.volume
    lat_vec = Si.reciprocal_vectors
    rj      = np.array([a.coords_fractional for a in Si.atoms])
    pattern = np.concatenate((rj,14*np.ones((rj.shape[0],1))),axis=1)
    hkl,Fhkl  = structure_factor3D(pattern,lat_vec,hklMax=N,sym=sym)
    Shkl = np.abs(Fhkl)**2
    #plot_structure3D(hkl,Fhkl,name=figpath+'Si_Shkl.png',opt='s',cmap='jet',ms=50,title='Si $c=log_{10}(S+1)$',view=[9,-72])
    Shkl = np.log10(Shkl+1)
    plot2Dcutplane(Shkl,n='100',name=figpath+'Si_S100.png',opt='s',cmap='jet')
    plot2Dcutplane(Shkl,n='110',name=figpath+'Si_S110.png',opt='s',cmap='jet')
    #plot1Dcutline(Shkl,u='001',dopt='')
    return hkl,Fhkl


if __name__ == "__main__" :
    plt.close('all')
    #hk,Fhk = test_structure_factor2D(N=8)
    hkl,Fhkl = test_structure_factor3D(N=4)
