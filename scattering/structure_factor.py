import io
from math import pi
from crystals import Crystal
from subprocess import check_output
import numpy as np
import pandas as pd
import gemmi
import utils.displayStandards as dsp
import utils.physicsConstants as cst
from . import scattering_factors as scatf

__all__=['structure_factor3D','plot_structure3D','get_miller3D',
        'get_pendulossung']

#2D for factors
ai = 1/np.array([0.1,0.25,0.26,0.27,1.5])**2
fj = lambda q,j,eps:eps*(np.pi/ai[j])*np.exp(-(np.pi*q)**2/ai[j])

def get_gemmi_structure_factor(cif_file,gemmi_prog,keV=200,dmin=0.5,aniso=True):
    '''Uses gemmi program to compute structure factors.Efficient but ouptuts only unique reflections
    Parameters : 
    -----------
        cif_file : cif file 
        gemmi_prog : location of gemmi program
        keV,dmin,aniso : energy,resolution, whether to model ADP
    Returns:
    --------
        DataFrame containing structure factors 
        
    '''    
    cmd="{gemmi_program} sfcalc --for=electron --wavelength={lam} --dmin={dmin} {aniso} {cif_file}".format(
        gemmi_program=gemmi_prog,    
        lam=cst.keV2lam(keV),
        cif_file=cif_file,
        dmin=dmin,
        aniso=['--noaniso',''][aniso]
        )
    out = check_output(cmd,shell=True).decode()    
    df_FhklG=pd.read_csv(io.StringIO(out),names=['hkl','amp','phi'],sep="\t",index_col=0)
    df_FhklG['F'] = df_FhklG.amp*np.exp(1J*np.deg2rad(df_FhklG.phi))
    df_FhklG['I'] = np.abs(df_FhklG.F)**2
    hkl = [eval(','.join(hkl.strip()[1:-1].split(" "))) for hkl in df_FhklG.index]
    df_FhklG['h']= [h[0] for h in hkl]
    df_FhklG['k']= [h[1] for h in hkl]
    df_FhklG['l']= [h[2] for h in hkl]
    df_FhklG.index=[str(h) for h in hkl]
    df_FhklG['hkl']=df_FhklG.index
    return df_FhklG
    
def get_small_structure_factor(cif_file,Nmax=2,aniso=True):
    '''Uses gemmi to compute structure factors. Not very efficient but accounts for ADPs
    Parameters : 
    -----------
        cif_file : cif file 
        
    Returns:
    --------
        DataFrame containing structure factors 
        
    '''
    small = gemmi.read_small_structure(cif_file)    
    small.change_occupancies_to_crystallographic()    
    calc_e = gemmi.StructureFactorCalculatorE(small.cell)
    if not aniso:
        for site in small.sites:
            site.aniso.u11 = site.aniso.u22 = site.aniso.u33 = 0
    # print(small.sites[0].aniso.nonzero())
    hkl = get_miller3D(Nmax,sym=1,flat=True)
    df  = pd.DataFrame(hkl,columns=['h','k','l'])
    df.index=[str(tuple(h)) for h in hkl]
    df['hkl']=df.index     
    df['F'] = [calc_e.calculate_sf_from_small_structure(small, (r.h,r.k,r.l)) for h,r in df.iterrows()]
    df['I'] = np.abs(df.F)**2
    return df


def get_structure_factor(cif_file,hkl=None,hklMax=10,sym=1,ed=True,aniso=False,v=''):
    '''Computes structure factor in 3D from :
    - `cif_file` : cif file
    - `hkl`     : list of miller indices h,k,l (shape Nx3 => [h;k;l])
    - `hklMax`  : int - max miller index in each direction from -hklMax to hklMax
    - `ed` : True for electron
    - `v` : verbosity
    '''
    ####init
    crys = Crystal.from_cif(cif_file)
    ra = np.array([a.coords_fractional for a in crys.atoms] )
    fa = np.array([a.atomic_number for a in crys.atoms] ).astype(int)
    lat_vec = np.array(crys.reciprocal_vectors)/(2*np.pi)
    df_Fhkl = pd.DataFrame()

    #### reflection list
    if not hkl :
        hkl = get_miller3D(hklMax,sym)
        hkl = np.array([h.flatten() for h in hkl]).T
    df_Fhkl[['h','k','l']] = hkl
    h,k,l = hkl.T

    #### wave vector
    #### [h,k,l].dot([b1;b2;b3] => hb1+kb2+lb3)
    q = np.linalg.norm(hkl.dot(lat_vec),axis=1)
    df_Fhkl['q']=q

    #### scattering factor
    atoms=np.unique(fa).tolist()
    if ed:
        q,fq = scatf.get_elec_atomic_factors(atoms,q)
        fq = np.array(fq).T
    else:
        elements = np.unique([a.element for a in crys.atoms])
        q,fq = scatf.get_fx(elements,q)
    if 'q' in v:qmax=q.max();print('qmax=%.4f A^-1\nmax_res=%.4f A' %(qmax,1/qmax))
 
    #structure factor
    Fhkl = np.zeros(q.shape,dtype=complex)
    for i,atom in enumerate(atoms):
        F_i = np.zeros(q.shape,dtype=complex)
        idx = fa==atom
        #loop over atoms with atomic number "atom"
        for ri in ra[idx,:]:
            x,y,z=ri
            F_i += np.exp(-2*pi*1J*(x*h+y*k+z*l))
        Fhkl += F_i*fq[:,i]

    
    df_Fhkl['F'] = Fhkl

    df_Fhkl.index  = [str((h,k,l)) for h,k,l in hkl]
    df_Fhkl['hkl'] = df_Fhkl.index 
    df_Fhkl['I']   = np.abs(df_Fhkl.F)**2
    
    return df_Fhkl


def structure_factor3D(pattern,lat_vec,hkl=None,hklMax=10,sym=1,ed=True,v=''):
    '''Computes structure factor in 3D from :
    - `pattern` : Nx4 array - N atoms with each row : fractional coordinates and Za
    - `lat_vec` : 3x3 array - reciprocal lattice vectors (2*pi/a convention)
    - `hkl`     : list of 3 miller indices h,k,l each as 3Ndarray
    - `hklMax`  : int - max miller index in each direction from -hklMax to hklMax
    '''
    #unpack
    if not hkl : hkl = get_miller3D(hklMax,sym)
    hx,ky,lz = hkl
    # print('rearranging h,k,l so h[(hi,k,l)]=hi')
    hx,ky,lz = np.transpose(hx,[1,0,2]),np.transpose(ky,[1,0,2]),np.transpose(lz,[1,0,2])
    hkl = [hx,ky,lz]
    ra,fa = pattern[:,:3],pattern[:,3]
    atoms = list(np.array(np.unique(fa),dtype=int));#print(atoms)
    # get scattering factor
    b1,b2,b3 = lat_vec/(2*pi)
    k_x,k_y,k_z = hx*b1[0]+ky*b2[0]+lz*b3[0],hx*b1[1]+ky*b2[1]+lz*b3[1],hx*b1[2]+ky*b2[2]+lz*b3[2]
    q = np.sqrt(k_x**2+k_y**2+k_z**2)
    if ed:
        q,fq = scatf.get_elec_atomic_factors(atoms,q)
    else:
        q,fq = scatf.get_fx(atoms,q)
    if 'q' in v:qmax=q.max();print('qmax=%.4f A^-1\nmax_res=%.4f A' %(qmax,1/qmax))
    #structure factor
    Fhkl,n_atoms = np.zeros(hx.shape,dtype=complex),len(atoms)
    for i,atom in zip(range(n_atoms),atoms):
        F_i = np.zeros(hx.shape,dtype=complex)
        idx = fa==atom
        #loop over atoms with of atomic number "atom"
        for ri in ra[idx,:]:
            # print(ri)
            # print(ri,np.exp(-2*pi*1J*(ri[0]*hx+ri[1]*ky+ri[2]*lz)))
            F_i += np.exp(-2*pi*1J*(ri[0]*hx+ri[1]*ky+ri[2]*lz))
        Fhkl += F_i*fq[i]

    # cs=dsp.getCs('Spectral',n_atoms)
    # plts=[[q.flatten(),fq[i].flatten(),[cs[i],'+'],atom] for i,atom in enumerate(atoms)]
    # dsp.stddisp(plts,lw=2)
    return hkl,Fhkl

def structure_factor2D(pattern,lat_vec,hk=None,hkMax=10,sym=1,v=0,eps=1):
    '''get structure factor
    - pattern : atomic positions in real space (fractional coordinates)
    - lat_vec : reciprocal lattice vectors (2pi/a convention)
    '''
    #unpack
    if not hk : hl = get_miller2D(hkMax,sym)
    hx,ky = hl
    hx,ky = hx.T,ky.T
    hl = [hx,ky]
    ra,fa = pattern[:,:2],pattern[:,2]
    atoms = list(np.array(np.unique(fa),dtype=int));#print(atoms)
    #scattering factor
    b1,b2 = lat_vec
    k_x,k_y = hx*b1[0]+ky*b2[0],hx*b1[1]+ky*b2[1]
    q = np.sqrt(k_x**2+k_y**2)/(2*pi)
    fq = [fj(q,Za,eps) for Za in atoms]
    if v : qmax=q.max();print('qmax=%.4f A^-1\nmax_res=%.4f A' %(qmax,1/qmax))
    #compute structure factor
    Fhl,n_atoms = np.zeros(hx.shape,dtype=complex),len(atoms)
    for i,atom in zip(range(n_atoms),atoms):
        F_i = np.zeros(hx.shape,dtype=complex)
        idx = fa==atom
        for ri in ra[idx,:]:
            F_i += np.exp(2*pi*1J*(ri[0]*hx+ri[1]*ky))
        Fhl += F_i*fq[i]
    return hl,Fhl

###
def get_pendulossung(name='Si',miller=[0,0,0],keV=200,opt='p'):
    ''' give the theorertical 2-beam approximation Pendullosung thickness
    - `name`    : compound
    - `miller`  : [h,k,l]
    - `keV`     : wavelength (keV)
    RETURN
    - xi : Pendullosung thickness (A)
    '''
    # compute structure factor
    from crystals import Crystal
    crys = Crystal.from_database(name)
    lat_vec  = crys.reciprocal_vectors
    pattern  = np.array([list(a.coords_fractional)+[a.atomic_number] for a in crys.atoms])
    hkl,Fhkl = structure_factor3D(pattern,lat_vec,hklMax=max(miller),sym=0,v='')
    ax,by,cz = crys.lattice_parameters[:3]
    Vcell    = crys.volume
    # compute Pendullosung
    h,k,l = miller
    Ug,K  = np.abs(Fhkl[h,k,l])/Vcell,1/scatf.wavelength(keV)
    xi  = K/Ug
    if 'p' in opt:
        print(green+"\tPendullosung thickness "+name+'[%d%d%d]' %(h,k,l)+black)
        print('%-5s= %.2f\n%-5s= %.2f\n%-5s= %.2f'%('K',K,'Ug',Ug,'xi',xi))
    return xi

##########################################################################
###def : display
##########################################################################
def plot_structure3D(hkl,Fhkl,log=0,**kwargs):
    '''scatter plot of the intensity'''
    hx,ky,lz = hkl
    hx,ky,lz = hx.flatten(),ky.flatten(),lz.flatten()
    S        = np.abs(Fhkl.flatten())**2
    if log : S = np.log10(S+1)
    N = hx.max()
    lims = [0,N,0,N,0,N]
    dsp.stddisp(scat=[hx,ky,lz,S],rc='3d',xylims=lims,**kwargs)

def plot2Dcutplane(Shkl,n='100',title='auto',**kwargs):
    '''Not working for arbitrary cut yet '''
    N = Shkl.shape[0]
    if title=='auto' : title='Silicon([%s]) $S_{hk}$' %n
    if   n=='100' : Scut = Shkl[0,:,:]
    elif n=='110' : Scut = Shkl[np.identity(N,dtype=bool),:]
    dsp.stddisp(im=Scut,imOpt='c',legOpt=0,title=title,**kwargs)

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
        qz,fq = scatf.get_elec_atomic_factors([14],qz)
        plts += [[qz,64*fq[0]**2,'g--','$64f_e^2$'],[qz,32*fq[0]**2,'c--','$32f_e^2$']]
    dsp.stddisp(plts,labs=[labx,'$S_q(A^2)$'],**kwargs)

##########################################################################
##### def: misc
def get_miller3D(hklMax=10,sym=1,flat=False):
    '''hklMax : Number of Beams such that Fourier components = 2^(N+1) - 1'''
    Nhkl = range(-hklMax*sym,hklMax+1)
    hkl = np.meshgrid(Nhkl,Nhkl,Nhkl)
    if flat :
        return np.array([h.flatten() for h in hkl]).T
    return hkl
def get_miller2D(hkMax=10,sym=1):
    Nhk = range(-hkMax*sym,hkMax+1)
    hk = np.meshgrid(Nhk,Nhk)
    return hk

##########################################################################
#### def : test
##########################################################################
def _test_structure_factor2D(N=5):
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

def _test_structure_factor3D(N=5,sym=False,**kwargs):
    from crystals import Crystal
    Si = Crystal.from_database('Si'); #vol     = Si.volume
    lat_vec = Si.reciprocal_vectors
    rj      = np.array([a.coords_fractional for a in Si.atoms])
    pattern = np.concatenate((rj,14*np.ones((rj.shape[0],1))),axis=1)
    hkl,Fhkl  = structure_factor3D(pattern,lat_vec,hklMax=N,sym=sym)
    Shkl = np.abs(Fhkl)**2
    #plot_structure3D(hkl,Fhkl,name=figpath+'Si_Shkl.png',opt='s',cmap='jet',ms=50,title='Si $c=log_{10}(S+1)$',view=[9,-72])
    Shkl = np.log10(Shkl+1)
    plot2Dcutplane(Shkl,n='100',name=figpath+'Si_S100.png',**kwargs)
    plot2Dcutplane(Shkl,n='110',name=figpath+'Si_S110.png',**kwargs)
    #plot1Dcutline(Shkl,u='001',dopt='')
    return hkl,Fhkl


if __name__ == "__main__" :
    figpath=get_figpath(__file__,"/figures/")
    plt.close('all')
    #hk,Fhk = _test_structure_factor2D(N=8)
    #hkl,Fhkl = _test_structure_factor3D(N=4,opt='p',cmap='jet')
