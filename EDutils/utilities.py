import importlib as imp
import tifffile,os,glob,pickle5,subprocess
import numpy as np,pandas as pd
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from crystals import Crystal
from utils import displayStandards as dsp   #;imp.reload(dsp)
from utils import glob_colors as colors     #;imp.reload(colors)
# from utils import handler3D as h3D          #;imp.reload(h3D)
# from . import rotate_exp as exp ;imp.reload(exp)

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
        name = '%s_%s%s' %('-'.join(params),tag,str(i).zfill(pad))
        sim_obj = Simu(path=path,name=name,**kwargs)
        df.loc[name,cols] = list(val)+[sim_obj._get_pkl()]
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

def load_pkl(file):
    """load an object"""
    with open(file,'rb') as f : obj = pickle5.load(f)
    return obj


#############################################################################
#### orientation vectors
#############################################################################
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
    - Nmax : max order of reflections(resolution)
    - Smax : maximum excitation error to be included
    '''
    K0 = np.linalg.norm(K)
    Kx,Ky,Kz = K

    (h,k,l),(qx,qy,qz) = get_lattice(lat_vec,Nmax)

    Sw = np.abs(np.sqrt((Kx-qx)**2+(Ky-qy)**2+(Kz-qz)**2) - K0)
    if Smax:
        idx = Sw<Smax
        h,k,l = np.array([h[idx],k[idx],l[idx]],dtype=int)
        qx,qy,qz,Sw = qx[idx],qy[idx],qz[idx],Sw[idx]
    d = dict(zip(['h','k','l','qx','qy','qz','Sw'],[h,k,l,qx,qy,qz,Sw]))
    return pd.DataFrame.from_dict(d)

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
    e3 = K/np.linalg.norm(K)
    e0 = e0/np.linalg.norm(e0)
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
