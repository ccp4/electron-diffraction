import importlib as imp
import json,tifffile,os,glob,pickle5,subprocess,crystals
import numpy as np,pandas as pd
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from crystals import Crystal
from utils import displayStandards as dsp   #;imp.reload(dsp)
from utils import glob_colors as colors     #;imp.reload(colors)
# from utils import handler3D as h3D          #;imp.reload(h3D)
# from . import rotate_exp as exp ;imp.reload(exp)
import gemmi

def sweep_var(Simu:object,
        params:Sequence[str],vals:Sequence[Sequence],
        tag:str='',path:str='',
        **kwargs,
    ):
    """runs a set of similar simulations Simu with a varying parameter

    Parameters
    ----------
    Simu
        Simulator object
    params
        list of parameters to be changed (can be a single string)
    vals
        values to assign to the parameters (nsimus x nparams)
    tag
        Name prefix of simulations
    path
        path to the simulations folder
    kwargs
        arguments passed to Simu

    Returns
    -------
        pd.DataFrame
            info about the set of simulations with params pkl
    """
    if not os.path.exists(path):
        print('creating directory:',path)
        print(subprocess.check_output('mkdir -p %s' %path,shell=True).decode())

    if isinstance(params,str):
        params=[params]
        vals=[[v] for v in vals]

    if not isinstance(vals,list):
        vals=list(vals)
    if not isinstance(vals[0],list):
        vals = [list(v) for v in vals]

    if not len(params)==len(vals[0]):
        raise Exception('len(params)=%d but it must equal len(vals[0])=%d' %(len(params),len(vals[0])))

    cols = params+['pkl']
    df = pd.DataFrame(columns=cols)#+pp.info_cols)

    nsimus = len(vals)
    pad = int(np.ceil(np.log10(nsimus)))
    for i,val in enumerate(vals):
        kwargs.update(dict(zip(params,val)))
        name = '%s_%s_%s' %('-'.join(params),tag,str(i).zfill(pad))
        sim_obj = Simu(path=path,name=name,**kwargs)
        df.loc[name,cols] = ''
        # df.loc[name,cols] = list(val)+[sim_obj._get_pkl()]
        df.loc[name,params] = list(val)
        df.loc[name,'pkl']  = sim_obj._get_pkl()
    df_file = os.path.join(path,'df_%s.pkl' %tag)
    df.to_pickle(df_file)
    print(colors.green+'Dataframe saved : '+colors.yellow+df_file+colors.black)
    return df

def get_pkl(file=None,path='',name='unkown.pkl'):
    if not file:
        file=os.path.join(path,name+'.pkl')
    return file

def save(obj,file=None,path='',name='unkown.pkl'):
    """save an object"""
    file = get_pkl(file,path,name)
    with open(file,'wb') as out :
        pickle5.dump(obj, out, pickle5.HIGHEST_PROTOCOL)
    print(colors.green+"object saved\n"+colors.yellow+file+colors.black)
def save_pkl(obj,file):
    save(obj,file=file)

def load_pkl(file):
    """load an object"""
    with open(file,'rb') as f : obj = pickle5.load(f)
    return obj

def rot(a,axis='x',deg=True):
    if deg:a = np.deg2rad(a)
    c,s = np.cos(a),np.sin(a)
    if   axis=='x':R = np.array([[1,0,0],[0,c,s],[0,-s,c]])
    elif axis=='y':R = np.array([[c,0,s],[0,1,0],[-s,0,c]])
    elif axis=='z':R = np.array([[c,s,0],[-s,c,0],[0,0,1]])
    return R


def rotation_matrix(k,a,deg=True):
    '''Rotation matrix from vector k and angle a'''
    #https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    if deg:a = np.deg2rad(a)
    kx,ky,kz = k
    K = np.array([[0,-kz,ky],[kz,0,-kx],[-ky,kx,0]])
    R = np.eye(3) + np.sin(a)*K + (1-np.cos(a))*K.dot(K)

    return R

def axis_and_angle_from_R(R,deg=True):
    ''' vector and angle from rotation matrix R'''
    g,C = np.linalg.eig(R)
    idx = np.where(np.abs(np.real(g)-1)<1e-5)[0][0]
    k = np.real(C[:,idx]).flatten()

    v0 = np.array([0,0,1])
    if np.linalg.norm(k-v0)<1e-2:
        v0=np.array([0,1,0])
    vr = R.dot(v0)
    v0_p = v0-v0.dot(k)*k
    vr_p = vr-vr.dot(k)*k
    alpha = np.arccos(np.dot(v0_p,vr_p)/(np.linalg.norm(v0_p)*np.linalg.norm(vr_p)))

    if np.sign(np.cross(k,v0_p).dot(vr_p))==-1:
        alpha=2*np.pi-alpha
    if deg:
        alpha=np.rad2deg(alpha)
    return k,alpha


#############################################################################
#### orientation vectors
#############################################################################
def get_uvw(u,npts=20,osc=0.2):
    """get a list of orientation vectors in a continuously rotated fashion
    about u += osc with npts point.
    The rotation axis is u x ez
    """
    ez = [0,0,1]
    rot_axis = np.cross(ez,u)

    u1  = rotation_matrix(rot_axis,-osc,deg=True).dot(u)
    u2  = rotation_matrix(rot_axis, osc,deg=True).dot(u)
    uvw = get_uvw_cont(u1,u2,nframes=npts)
    return uvw

def get_uvw_cont(u0:Sequence,u1:Sequence,nframes:int=2,show:bool=False,**kwargs):
    """get a list of orientation vectors in a continuously rotated fashion between u0 and u1.

    Parameters
    ----------
    u0
        initial vector
    u1
        final vector
    nframes
        number of vectors inbetween
    show
        show the vectors
    kwargs
        :func:~utilities.show_uvw
    """
    u0,u1=np.array(u0,dtype=float),np.array(u1,dtype=float)
    u0/=np.linalg.norm(u0)
    u1/=np.linalg.norm(u1)
    e1 = u1 - u0.dot(u1)*u0
    e1/=np.linalg.norm(e1)

    omega_max = np.arccos(u1.dot(u0))
    omega_r = np.linspace(0,omega_max,nframes)

    ct,st = np.cos(omega_r),np.sin(omega_r)
    u0e1  = np.hstack([np.array(u0)[:,None],np.array(e1)[:,None]])
    xy    = np.vstack([ct,st])
    uvw   = u0e1.dot(xy).T

    if show:
        if not 'eij' in kwargs.keys():kwargs['eij']={}
        kwargs['eij'].update(dict(zip(['u_0','u_1'],[u0,u1])))
        show_uvw(uvw,**kwargs)
    return uvw

def get_uvw_rock(e0,e1=[0,1],deg:float=1,dt:Optional[float]=None,npts:int=20,**kwargs):
    """create orientation vectors around u0 to simulate a rocking curve

    Parameters
    ----------
    e0
        main beam orientation
    e1
        rocking curve direction components in spherical basis ( u1[x,y,z] = u1[0]*e_theta+u1[1]*e_phi)
    deg
        angular range in degrees
    dt
        angular step in degrees
    npts
        number of points (angular spread npts*deg)
    kwargs
        passed to :func:~utilities.get_uvw_cont
    """
    e0=np.array(e0,dtype=float)/np.linalg.norm(e0)

    theta,phi = np.deg2rad(theta_phi_from_u(e0))    #euler angles
    ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
    e_theta = np.array([ct*cp,ct*st,-st])
    e_phi = np.array([-sp,cp,0])
    e1  = e1[0]*e_theta+e1[1]*e_phi
    e1 /= np.linalg.norm(e1)

    if dt:
        deg = dt*npts
        print('angular spread=%.1f degrees' %deg)
    deg_r = np.deg2rad(deg)
    u0 = e0-e1*deg_r
    u1 = e0+e1*deg_r

    if not 'eij' in kwargs.keys():kwargs['eij']={}
    kwargs['eij'].update(dict(zip(['e_0'],[e0])))
    uvw = get_uvw_cont(u0,u1,nframes=npts,**kwargs)
    return uvw

def get_uvw_CBED(u0:Sequence,deg:float,npts:int,
    cart:bool=True,show:bool=False,**kwargs
    ):
    """get orientation vectors in a CBED pattern

    Parameters
    ----------
    u0
        main orientation
    deg
        angular range in degrees
    npts
        number of points around u0 in each direction
    cart
        compute using cartesian transform
    show
        show uvw
    kwargs
        passed to :func:~utilities.show_uvw
    """
    u0 = np.array(u0,dtype=float)
    e0 = u0/np.linalg.norm(u0)
    e1 = np.cross(e0,[0,1,0])
    e2 = np.cross(e0,e1)
    eij = np.vstack([e0,e1,e2])
    dt = np.deg2rad(deg)

    if cart:
        dx = np.tan(dt)
        x0 = np.linspace(-dx,dx,npts)
        x,y = np.meshgrid(x0,x0)
        x,y = x.flatten(),y.flatten()
        uvw = np.vstack([x,y]).T.dot(np.vstack([e1,e2]))+e0
        # uvw += e0
        uvw = (uvw.T/np.linalg.norm(uvw,axis=1)).T
    else:
        if isinstance(npts,int):npts=[npts]*2
        ds = np.linspace(0,dt,npts[0])
        theta,phi = np.meshgrid(ds,np.linspace(0,2*np.pi,npts[1]))
        theta,phi = theta.flatten(),phi.flatten()
        ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
        uvw0 = np.vstack([ct,st*cp,st*sp]).T
        uvw = uvw0.dot(eij)
        # uvw /=np.linalg.norm(uvw,axis=1)

    if show:show_uvw(uvw,eij=dict(zip(['e_0','e_1','e_2'],eij)),**kwargs)
    return uvw

def show_uvw(uvw:Iterable[Sequence],ax=None,
        eij:dict={'e_z':[0,0,1]},h3d:bool=False,**kwargs
    ):
    """show orientation vectors

    Parameters
    ----------
    uvw
        Array of orientations
    eij
        additional vectors to plot
    h3d
        use vector handlers
    """
    ez =  np.cross(uvw[0],uvw[-1])
    cs = dsp.getCs('viridis',len(list(eij.keys())) )

    plts = [[[0,u[0]],[0,u[1]],[0,u[2]],'b'] for u in uvw ]
    plts+= [[[0,ei[0]],[0,ei[1]],[0,ei[2]] , c ,'$%s$' %k] for c,(k,ei) in zip(cs,eij.items())]
    fig,ax = dsp.stddisp(plots=plts,labs=['x','y','z'],rc='3d',lw=4,xylims=1,**kwargs)
    if h3d:h3d = h3D.handler_3d(fig,persp=False,xm0=1)
    return fig,ax

def u_from_theta_phi(theta,phi):
    theta,phi = np.deg2rad([theta,phi])
    ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
    u = [st*cp,st*sp,ct]
    return u

def theta_phi_from_u(u):
    u = u/np.linalg.norm(u)
    x,y,z  = u
    theta = np.rad2deg(np.arccos(z))
    phi   = np.rad2deg(np.arctan2(y,x))
    return theta,phi

def get_uvw_from_theta_phi(theta,phi,e1=[0,0,1],**kwargs):
    u0 = np.array(u_from_theta_phi(theta,phi))
    u1 = np.cross(e1,u0);u1/=np.linalg.norm(u1)
    uvw = get_uvw_cont(u0,u1,**kwargs)
    return uvw

def uvw_add_points(uvw0,npts=1,plot=0,h3d=1,**kwargs):
    ntot = (uvw0.shape[0]-1)*(npts+1) + 1
    uvw = np.zeros((ntot,3))
    for i in range(uvw0.shape[0]-1):
        e0,e1 = uvw0[i:i+2,:]
        # e0,u1 = uvw0[i:i+2,:]
        # e1 = u1 - e0.dot(u1)*e0
        # e1/=np.linalg.norm(e1)
        # alpha = np.rad2deg(np.arccos(u1.dot(e0)))#;print(alpha)
        uvw_i = get_uvw_cont(e0,e1,nframes=npts+2)#omega=np.linspace(0,alpha,npts+2))
        idx = i*(npts+1)
        uvw[idx:idx+npts+1,:] = uvw_i[:-1,:]
    uvw[-1,:] = uvw_i[-1,:]

    if plot:
        plts  = [[[0,u[0]],[0,u[1]],[0,u[2]],'r'] for u in uvw ]
        fig,ax = dsp.stddisp(plots=plts,labs=['x','y','z'],rc='3d',lw=4,opt='')
        plts = [[[0,u[0]],[0,u[1]],[0,u[2]],'b'] for u in uvw0 ]
        dsp.stddisp(ax=ax,plots=plts,labs=['x','y','z'],rc='3d',lw=2,**kwargs)
        # if h3d:h3d = h3D.handler_3d(fig,persp=False,xm0=1)
        return fig,ax
    return uvw



######################################################
#### misc
######################################################
def _find_files(path,f_type):
    files = glob.glob(path+'/*.%s' %f_type)
    if len(files)==1:
        file=files[0]
        return file
    else:
        msg ='''
    only 1 %s file should be found in path folder :
     %s
    but %d were found :
     %s''' %(f_type,os.path.realpath(path),len(files),str(files))
        raise Exception(msg)

def find_cif_file(path,cif_file=None):
    """find cif file in path"""
    if not cif_file:
        cif_file = _find_files(path,'cif')
    return cif_file

def pdb2npy(pdb,npy_file=''):
    crys=crystals.Crystal.from_pdb(pdb)
    Zxyz=np.array([
        np.array([a.atomic_number,a.coords_fractional%1],dtype=object)
            for a in crys],dtype=object)
    if not npy_file: npy_file=pdb+'.npy'
    np.save(npy_file,np.array([Zxyz ,crys.lattice_vectors],dtype=object))
    print(colors.green+'file saved : '+colors.yellow+npy_file+colors.black)
    return npy_file

def import_npy(npy_file):
    (Zxyz,lat) = np.load(npy_file,allow_pickle=True)
    crys=crystals.Crystal(
        crystals.AtomicStructure([
            crystals.Atom(element=Z,coords=xyz) for Z,xyz in Zxyz]),
        lat)
    return crys

def crys2felix(crys,opt='wr',out=None):
    # crys = Crystal.from_cif(cif_file)
    cif = """_chemical_formula_sum '%s'
loop_
_cell_length_a %.5f
_cell_length_b %.5f
_cell_length_c %.5f
_cell_angle_alpha %.3f
_cell_angle_beta %.3f
_cell_angle_gamma %.3f
_symmetry_Int_Tables_number %d
_symmetry_space_group_name_Hall '%s'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_B_iso_or_equiv
_atom_site_occupancy
""" %(
        crys.chemical_formula,*crys.lattice_parameters,
        crys.international_number,crys.hall_symbol)
    cif += '\n'.join([
        "%s %s %.5f %.5f %.5f 0.1 1" %("%s%d" %(a.element,i), a.element,*a.coords_fractional)
            for i,a in enumerate(crys.atoms) if all(a.coords_fractional<0.99)
        ])

    if not out:out='felix.cif'
    if 'w' in opt:
        with open(out,'w') as f : f.write(cif)
        print(colors.yellow+out+colors.green+" successfully created."+colors.black)
    if 'r' in opt:
        if 'w' in opt:
            print(colors.blue+"Content of "+ colors.yellow +out+colors.black)
            with open(out,'r') as f : print(''.join(f.readlines()))
        else:
            print(cif)

def gemmi_sf(pdb_file:str='',dmin=2):
    st = gemmi.read_structure(pdb_file)
    dc = gemmi.DensityCalculatorX()
    dc.d_min = dmin
    dc.addends.subtract_z()
    dc.set_grid_cell_and_spacegroup(st)
    dc.set_refmac_compatible_blur(st[0])
    dc.put_model_density_on_grid(st[0])
    grid = gemmi.transform_map_to_f_phi(dc.grid); print('sf : ',grid.shape)
    Nmax = (np.array(grid.shape)//2-1).min()
    Fhkl = np.zeros((2*Nmax+1,)*3,dtype=complex)
    hklF = np.meshgrid(*[np.arange(-Nmax,Nmax+1)]*3)
    hkls = np.array([i.flatten() for i in hklF]).T
    idxs = np.array([i.flatten() for i in np.meshgrid(*[np.arange(2*Nmax+1)]*3)]).T
    #TODO : multiple size grid
    # nu2,nv2,nw2 = np.array(grid.shape)//2-1
    # Fhkl = np.zeros((2*nu2+1,2*nv2+1,2*nw2+1),dtype=complex)
    # hkls=np.array([i.flatten() for i in np.meshgrid(range(-nu2,nu2+1),range(-nv2,nv2+1),range(-nw2,nw2+1))]).T
    # idxs=np.array([i.flatten() for i in np.meshgrid(range(2*nu2+1),range(2*nv2+1),range(2*nw2+1))]).T
    # hkls = pd.DataFrame()
    print(colors.blue+'...filling...'+colors.black)
    for idx,hkl in zip(idxs,hkls):
        Fhkl[tuple(idx)] = dc.mott_bethe_factor(hkl) * grid.get_value(*hkl)
    return hklF,Fhkl



def import_crys(file:str=''):
    """import a Crystal

    Parameters
    ----------
    file
        can be a cif file or a builtin
    """
    builtins=list(Crystal.builtins)

    if file.split('.')[-1]=='cif':
        crys = Crystal.from_cif(file)
    elif file.split('.')[-1]=='npy':
        crys = import_npy(file)
    elif sum(np.array(list(Crystal.builtins))==file):
        crys = Crystal.from_database(file)
    else:
        if file:
            msg = 'cannot import %s\n' %file
            msg+= 'Provide a .cif file or use a builtin structure \n%s' %builtins.__str__()
            raise Exception(msg)
        else:
            print('available builtins:\n')
            print(builtins.__str__())
            return
    print(colors.green+'imported file : '+colors.yellow+file+colors.black)
    return crys

def convert2tiff(tiff_file,im0,n0=0,rot=0,n=256,Imax=5e4):
    '''
    used to rotate a multislice image to a standard pets tiff file image
    '''
    alpha = np.deg2rad(rot)
    ct,st = np.cos(alpha),np.sin(alpha)
    u_n   = np.arange(-n,n+1)
    hi,ki = np.meshgrid(u_n,u_n)
    Hi,Ki = ct*hi+st*ki,st*hi-ct*ki
    H,K = np.array(np.round(Hi)+n0,dtype=int),np.array(np.round(Ki)+n0,dtype=int)
    # dsp.stddisp(im=[I[H,K]],caxis=[0,50],pOpt='im',cmap='gray')
    I = np.array(im0*Imax,dtype='uint16')
    tifffile.imwrite(tiff_file,I[H,K]) #,np.flipud(I))
    print(colors.yellow+tiff_file+colors.green+' saved'+colors.black)




################################################################################
#### Reciprocal space related
################################################################################
def get_reciprocal(abc):
    a1,a2,a3 = abc
    b1 = np.cross(a2,a3)/(a1.dot(np.cross(a2,a3)))
    b2 = np.cross(a3,a1)/(a2.dot(np.cross(a3,a1)))
    b3 = np.cross(a1,a2)/(a3.dot(np.cross(a1,a2)))
    abc_star = np.vstack([b1,b2,b3])#*(2*np.pi)
    return abc_star

def get_lattice(lat_vec,Nmax=5):
    '''Get a lattice up to Nmax points :
    - lat_vec : reciprocal lattice vectors
    - Nmax : max order '''
    a1,b1,c1 = lat_vec
    N = np.arange(-Nmax,Nmax+1)
    h,k,l = np.meshgrid(N,N,N)
    h,k,l = h.flatten(),k.flatten(),l.flatten()
    qx = h*a1[0]+k*b1[0]+l*c1[0]
    qy = h*a1[1]+k*b1[1]+l*c1[1]
    qz = h*a1[2]+k*b1[2]+l*c1[2]
    return (h,k,l),(qx,qy,qz)
def get_lattice2D(lat_vec,Nmax=5):
    '''Get a lattice up to Nmax points :
    - lat_vec : reciprocal lattice vectors
    - Nmax : max order '''
    a1,b1 = lat_vec
    N = np.arange(-Nmax,Nmax+1)
    h,k = np.meshgrid(N,N)
    h,k = h.flatten(),k.flatten()
    qx = h*a1[0]+k*b1[0]
    qz = h*a1[1]+k*b1[1]
    return (h,k),(qx,qz)

def get_ewald(K,nts=100,nps=200):
    ''' Get ewald sphere coordinates
    - K : reciprocal space beam vector
    '''
    Kx,Ky,Kz = K

    theta,phi = np.linspace(0,np.pi,nts),np.linspace(0,2*np.pi,nps)
    x = Kx+self.K0*np.outer(np.sin(theta),np.cos(phi))
    y = Ky+self.K0*np.outer(np.sin(theta),np.sin(phi))
    z = Kz+self.K0*np.outer(np.cos(theta),np.ones(phi.shape))
    return x,y,z

def get_excitation_errors(K,lat_vec,Nmax=5,Smax=0.02):
    ''' get excitation errors for lattice lat_vec and beam K
    - K : reciprocal space beam vector
    - lat_vec : reciprocal lattice vectors
    - Nmax : max order of reflections(resolution)
    - Smax : maximum excitation error to be included
    '''
    K0 = np.linalg.norm(K)
    Kx,Ky,Kz = K

    (h,k,l),(qx,qy,qz) = get_lattice(lat_vec,Nmax)

    # Sw = np.abs(np.sqrt((Kx-qx)**2+(Ky-qy)**2+(Kz-qz)**2) - K0)
    Sw = np.abs(np.sqrt((Kx+qx)**2+(Ky+qy)**2+(Kz+qz)**2) - K0)
    if Smax:
        idx = Sw<Smax
        h,k,l = np.array([h[idx],k[idx],l[idx]],dtype=int)
        qx,qy,qz,Sw = qx[idx],qy[idx],qz[idx],Sw[idx]
    d = dict(zip(['h','k','l','qx','qy','qz','Sw'],[h,k,l,qx,qy,qz,Sw]))
    df = pd.DataFrame.from_dict(d)
    df.index=[str(tuple(h)) for h in df[['h','k','l']].values]
    return df

def get_structure_factor(cif_file,dfout=0,**sf_args):
    '''computes structure factor 3D
    - sf_args : see (structure_factor3D)
    returns :
    - (qx,qy,qz),Fhkl
    '''
    crys = import_crys(cif_file)
    pattern = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in crys.atoms] )
    lat_vec = np.array(crys.reciprocal_vectors)
    (h,k,l),Fhkl = structure_factor3D(pattern, lat_vec, **sf_args)
    qx = h/crys.lattice_parameters[0]
    qy = k/crys.lattice_parameters[1]
    qz = l/crys.lattice_parameters[2]
    if dfout:
        data = [h.flatten(),k.flatten(),l.flatten(),qx.flatten(),qy.flatten(),qz.flatten(),Fhkl.flatten()]
        d = dict(zip(['h','k','l','qx','qy','qz','Fhkl'],data))
        return pd.DataFrame.from_dict(d)
    else:
        return (qx,qy,qz),Fhkl

def get_excited_beams(cif_file,K,lat_vec,Nmax,Smax):
    df_Sw   = get_excitation_errors(K,lat_vec,Nmax,Smax)
    df_Fhkl = get_structure_factor(cif_file,dfout=1,hklMax=Nmax+1)

    hkl0 = df_Sw[['h','k','l']].values
    hklF = df_Fhkl[['h','k','l']].values
    ridx = []
    for r in hkl0:
        idx = np.where(np.linalg.norm(r-hklF,axis=1)==0)[0]
        if idx.size:
            ridx+=[idx[0]]

    Fhkl_0 = df_Fhkl.iloc[ridx]['Fhkl'].values #;print(df_Fhkl.iloc[ridx]);#print(ridx);print(df_Fhkl.shape)
    df_Sw['Fhkl'] = Fhkl_0
    return df_Sw

def get_kinematic_intensities(cif_file,K,thick,Nmax=5,Smax=None):
    crys    = import_crys(cif_file)
    lat_vec = np.array(crys.reciprocal_vectors)/(2*np.pi)
    df_Sw   = get_excited_beams(cif_file,K,lat_vec,Nmax,Smax)

    #kinematic intensities
    K0  = np.linalg.norm(K)
    sig = cst.keV2sigma(cst.lam2keV(1/K0))      #;print(sig)
    Sw  = df_Sw.Sw.values
    gs = (sig*np.pi*thick*Sw)*np.sinc(thick*Sw)

    Fhkl_0 = df_Sw['Fhkl'].values
    I  = (np.abs(Fhkl_0)*gs)**2
    df_Sw['I'] = I
    return df_Sw

def project_beams(K,qxyz,e0=[1,0,0],v=0):
    '''pass from reciprocal cartesian coordinates to panel coordinates
    then uses the `aper` parameter to get values in pixels
    e0:rotation axis in the panel frame
    Note:
        In the images frame, the beam in perpendicular to the images(z axis)
        and the rotation axis is usually along the x axis of the images
    '''
    #build the frames basis
    e3 = K/np.linalg.norm(K)
    if not np.linalg.norm(e3-e0):
        if not np.linalg.norm(e3-np.array([1,0,0])):
            e0=[0,1,0]
        else:
            e0=[1,0,0]

    e0 = np.array(e0)/np.linalg.norm(e0)
    e2 = np.cross(e3,e0)
    e1 = np.cross(e2,e3)

    px,py = qxyz.dot(e1),qxyz.dot(e2)
    if v:
        return px,py,e0.dot(e1)
    else:
        return px,py

def project_beams2D(K,qxy):
    u = K/np.linalg.norm(K)
    v = np.array([u[1],-u[0]])
    px = qxy.dot(v)
    return px


def remove_friedel_pairs(reflF):
    hkls = np.array([np.array(h[1:-1].split(','),dtype=int)  for h in reflF])
    refl=[]
    for h,hkl in zip(reflF,hkls):
        if not str(tuple(-hkl)) in refl:
            refl.append(h)
    print('removing Friedel pairs')
    return refl

def show_tiff(file,cutoff=100):
    import tifffile
    im=tifffile.imread(file)
    dsp.stddisp(im=[im],title=file,
        cmap='gray',caxis=[0,cutoff],pOpt='tX')


def import_fcf(file_path = 'shelx_thick_10A.fcf'):
    # Initialize an empty list to store your data
    data = []
    # Open the file and read it line by line
    with open(file_path, 'r') as file:
        # A flag to mark when the relevant data block starts
        data_block_start = False
        for line in file:
            # Check for the start of the data block
            if line.strip() == '_refln_index_h':
                data_block_start = True
                continue

            # If in the data block and the line starts with an underscore, it's a column header, so skip it
            if data_block_start and line.strip().startswith('_'):
                continue

            # If in the data block and the line doesn't start with an underscore, it's a data row
            if data_block_start and not line.strip().startswith('_'):
                # Check for the end of the data block (empty line or new section starting)
                if line.strip() == '' or line.strip().startswith('_'):
                    data_block_start = False
                    continue
                # Split the line by whitespace and add the resulting list to the data
                data.append(line.split())

    # Define the column names
    columns = ['h', 'k', 'l', 'Fc-squared', 'Fo-squared', 'sigma(Fo-squared)', 'status flag']

    # # Create a DataFrame
    df = pd.DataFrame(data, columns=columns)

    # Convert numerical columns to appropriate types
    df[['h', 'k', 'l']] = df[['h', 'k', 'l']].astype(int)
    df[['Fc-squared', 'Fo-squared', 'sigma(Fo-squared)']] = df[['Fc-squared', 'Fo-squared', 'sigma(Fo-squared)']].astype(float)

    df.index=[str((h,k,l)) for h,k,l in df[['h','k','l']].values]
    df['hkl']=df.index
    return df

def to_shelx(df_hkl,file):
    '''converts to a.hkl file ready to use by shelx
    hkl: dataframe containing h,k,l,I,sig
    '''
    hkl=df_hkl.copy()
    hkl=hkl.drop(str((0,0,0)))
    #### max the intensity at 9999.99 for shelx (lol)
    maxI=hkl.I.max()
    # if maxI>=1e4:
    hkl.I=hkl.I*9999.99/maxI
    formats = {
        'h':'{:>4}', 'k':'{:>3}', 'l':'{:>3}',
        'I':'{:>7.2f}', 'sig':'{:>7.2f}',
        }
    formatters = {k: v.format for k, v in formats.items()}
    content=hkl.to_string(formatters=formatters, index=False, header=False)
    with open(file, 'w') as f:
        f.write(content)

    print(colors.green+'file saved : \n'+colors.yellow,file,colors.black)

def read_space_group(struct_file):
    dataset={}
    if struct_file.split('.')[-1]=='cif':
        with open(struct_file,'r') as f:
            lines=f.readlines()
        spg_info = {
            'lattice_system'      :['_space_group_crystal_system','_symmetry_cell_setting'],
            'international_number':['_space_group_IT_number'     ,],
            'international_symbol':['_space_group_name_H-M_alt'  ,'_symmetry_space_group_name_H-M'],
        }

        # (key,vals),l=next(spg_info),next(lines)
        for l in lines:
            if l.strip() in ['_space_group_symop_operation_xyz','_symmetry_equiv_pos_as_xyz']:
                break
            for key,vals in spg_info.items():
                for v in vals:
                    if v==l.split(' ')[0]:
                        val = l.replace(v,'').strip(" '\n")
                        # print('match found:',l.split(' ')[0],v,val)
                        dataset[key]=val
                        break

    elif struct_file.split('.')[-1]=='pdb':
        with open(struct_file,'r') as f:
            lines=f.readlines()
        for l in lines:
            if 'SYMMETRY OPERATORS FOR SPACE GROUP' in l:
                dataset['international_symbol'] = l.split(':')[1].strip()
                break

    if dataset:
        print(colors.green+'Reading space group info manually in cif file'+colors.black)
        if 'international_symbol' in dataset.keys() and 'international_number' not in dataset.keys():
            with open('static/spg/spg_groups_hm.json','r') as f:
                spg_hm=json.load(f)
            hm = dataset['international_symbol']
            if hm in spg_hm.keys():
                dataset['international_number'] = spg_hm[hm]
                print(colors.blue+'retrieved space group number from dictionary'+colors.black)

    return dataset


def compute_B(params):
    a,b,c,alpha,beta,gamma = params
    angles=np.deg2rad([alpha,beta,gamma])
    ca,cb,cc = np.cos(angles)
    sa,sb,sc = np.sin(angles)
    V=a*b*c*np.sqrt(1-ca**2-cb**2-cc**2 + 2*ca*cb*cc)
    B = np.array([
        [1/a,0,0],
        [-cc/(a*sc),1/(b*sc),0],
        [b*c/V*(cc*(ca-cb*cc)/sc - cb*sc), a*c/(V*sc)*(ca-cb*cc),a*b*sc/V ],
        ])
    return B
