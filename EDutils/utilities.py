import importlib as imp
import tifffile,os,glob,pickle5
import numpy as np,pandas as pd
from utils import displayStandards as dsp   #;imp.reload(dsp)
from utils import glob_colors as colors     #;imp.reload(colors)
from utils import handler3D as h3D          #;imp.reload(h3D)

from . import rotate_exp as exp ;imp.reload(exp)

def sweep_var(Simu,param,vals,tag='',path='',**kwargs):
    '''runs a set of similar simulations Simu with a varying parameter
    args :
        - Simu      : Simulator constructor
        - param,vals: the parameter and values to sweep
        - tag       : Name prefix of simulations
        - path      : path to the simulations folder
    returns :
        - df : pd.DataFrame of the simulation
    - kwargs : Simu constructor arguments
    '''
    cols = [param,'pkl']
    df = pd.DataFrame(columns=cols)#+pp.info_cols)
    pad = int(np.ceil(np.log10(vals.size)))
    for i,val in enumerate(vals):
        kwargs[param] = val
        name = '%s_%s%s' %(tag,param,str(i).zfill(pad))
        sim_obj = Simu(path=path,name=name,**kwargs)
        df.loc[name,cols] = [val,sim_obj.get_pkl()]
    df_file = os.path.join(path,'df_%s.pkl' %tag)
    df.to_pickle(df_file)
    print(colors.green+'Dataframe saved : '+colors.yellow+df_file+colors.black)
    return df

def load_pkl(file):
    with open(file,'rb') as f : obj = pickle5.load(f)
    return obj


class Rocking(exp.Rocking):
    def test():print('ok')

#############################################################################
#### misc
#############################################################################
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
    uvw = get_uvw(u0,u1,**kwargs)
    return uvw

def get_uvw(u0,u1,omega,nframes=2,plot=False,h3d=False):
    ''' get a list of continuously rotated vector between u0 and u1.
    u0,u1 : initial and final vector
    omega : array or tuple or int : rotation range to cover
    nframes : number of vectors (omega is int or tuple)
    '''
    if type(omega) in [float,int]:omega=(0,omega)
    if isinstance(omega,tuple):omega=np.linspace(omega[0],omega[1],nframes)
    u1/=np.linalg.norm(u1)
    u1 = u1 - u0.dot(u1)*u0
    u1/=np.linalg.norm(u1)

    omega_r = np.deg2rad(omega)
    ct,st = np.cos(omega_r),np.sin(omega_r)
    u0u1  = np.hstack([np.array(u0)[:,None],np.array(u1)[:,None]])
    xy    = np.vstack([ct,st])
    uvw   = u0u1.dot(xy).T
    if plot:show_uvw(uvw,h3d)
    return uvw

def uvw_add_points(uvw0,npts=1,plot=0,h3d=1):
    ntot = (uvw0.shape[0]-1)*(npts+1) + 1
    uvw = np.zeros((ntot,3))
    for i in range(uvw0.shape[0]-1):
        e0,u1 = uvw0[i:i+2,:]
        e1 = u1 - e0.dot(u1)*e0
        e1/=np.linalg.norm(e1)
        alpha = np.rad2deg(np.arccos(u1.dot(e0)))#;print(alpha)
        uvw_i = get_uvw(e0,e1,omega=np.linspace(0,alpha,npts+2))
        idx = i*(npts+1)
        uvw[idx:idx+npts+1,:] = uvw_i[:-1,:]
    uvw[-1,:] = uvw_i[-1,:]

    if plot:
        plts  = [[[0,u[0]],[0,u[1]],[0,u[2]],'r'] for u in uvw ]
        plts += [[[0,u[0]],[0,u[1]],[0,u[2]],'b'] for u in uvw0 ]
        fig,ax = dsp.stddisp(plots=plts,labs=['x','y','z'],rc='3d',lw=4)
        if h3d:h3d = h3D.handler_3d(fig,persp=False,xm0=1)
    return uvw

def show_uvw(uvw,h3d=False,**kwargs):
    ez =  np.cross(uvw[0],uvw[-1])

    plts = [[[0,u[0]],[0,u[1]],[0,u[2]],'b'] for u in uvw ]
    plts +=[[ [0,ez[0]],[0,ez[1]],[0,ez[2]] , 'r' ]]
    fig,ax = dsp.stddisp(plots=plts,labs=['x','y','z'],rc='3d',lw=4,**kwargs)
    if h3d:h3d = h3D.handler_3d(fig,persp=False,xm0=1)
    return fig,ax

def get_uvw_rock(u0,u1=[0,1],omega=np.linspace(0,0.5,32)):
    '''create orientation vectors around u0 to simulate rocking curve
    - u0 : main beam orientation
    - u1 : so that u1|x,y,z = u1[0]*e_theta+u1[1]*e_phi
    - omega : angular spread
    '''
    theta,phi = np.deg2rad(theta_phi_from_u(u0))
    ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
    e_theta = np.array([ct*cp,ct*st,-st])
    e_phi = np.array([-sp,cp,0])
    u1  = u1[0]*e_theta+u1[1]*e_phi
    u1/=np.linalg.norm(u1)
    uvw = get_uvw(u0,u1,omega,plot=0)
    return uvw

def convert2tiff(tiff_file,im0,n0,rot,n=256,Imax=5e4):
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

def make_pets(pts_file,aperpixel,deg=0.0652,ref_cell=None,tag=''):
    center,pixel_size='AUTO',0.01,
    phi,omega= deg/2 ,0
    ax,by,cz,alpha,beta,gamma = ref_cell
    # ax,by,cz,alpha,beta,gamma = 8.1218, 9.977, 17.725, 90.0, 90.0, 90.0
    pts = '''lambda 0.025080
geometry continuous
omega  %d
phi  %.5f
virtualframes   7 5 1

aperpixel    %.6f
noiseparameters      2.5000      1.0000
saturationlimit   20000

center    256.0 256.0
centermode friedelpairs 0

beamstop no

dstarmax  3.0
dstarmaxps  3.0
i/sigma    7.00    5.00
reflectionsize  7

referencecell     %.5f    %.5f     %.5f    %.5f   %.5f    %.5f 1

#List of images
#columns: file name,alpha,beta,delta omega,frame scale,calibration correction(px/rec.angstrom),ellipt.distortion correction(amplitude),ellipt.distortion correction(phase), use for calculations(0/1)
imagelist
'''%(omega,phi,aperpixel,  ax,by,cz,alpha,beta,gamma)
    out = os.path.dirname(pts_file)
    tif_files = np.sort(glob.glob(out+'/tiff/%s*.tiff' %tag))
    alphas = np.arange(tif_files.size)*deg
    for i,tif_file in enumerate(tif_files):
        tif_file = os.path.basename(tif_file)
        pts+='%s %.4f 0.0 0.0 1.0 0 0 0  1\n' %('tiff\\'+tif_file,alphas[i])
    pts+='endimagelist\n'
    # pts_file = out+name+'.pts'
    with open(pts_file,'w') as f:
        f.write(pts)
        print(colors.green+'file saved : '+colors.yellow+pts_file+colors.black)
