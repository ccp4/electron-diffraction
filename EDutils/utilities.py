import importlib as imp
import tifffile,os,glob,pickle5
import numpy as np,pandas as pd
from utils import displayStandards as dsp   #;imp.reload(dsp)
from utils import glob_colors as colors     #;imp.reload(colors)
from utils import handler3D as h3D          #;imp.reload(h3D)

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

class Rocking:
    def __init__(self,Simu,param,vals,ts,tag,path,**kwargs):
        ''' simulate rocking curve
        - ts : rocking parameter list (theta,Sw,tilt)
        - kwargs : Simu constructor arguments
        '''
        self.path = path
        self.tag  = tag
        self.vals = vals
        self.df = sweep_var(Simu,param,vals,tag=tag,path=path,**kwargs)
        self.ts       = ts
        self.df['ts'] = ts
        self.save(v=1)

    def save(self,v=0):
        '''save this object'''
        file = os.path.join(self.path,'rock_%s.pkl' %(self.tag))
        with open(file,'wb') as out :
            pickle5.dump(self, out, pickle5.HIGHEST_PROTOCOL)
        if v:print(colors.green+"object saved\n"+colors.yellow+file+colors.black)

    def load(self,i=0,ts=None):
        if type(ts) in [float,int]:i=np.argmin(abs(self.ts-ts))
        file = self.df.iloc[i].pkl
        sim_obj = load_pkl(file)
        return sim_obj

    def do(self,f,**args):
        for i in range(self.df.shape[0]):
            obj = self.load(i)
            obj.__getattribute__(f)(**args)
            obj.save()

    def plot_rocking(self,iZs=-1,zs=None,bargs={},cmap='viridis',**kwargs):
        '''plot rocking curve for set of beams at thickness zs
        - iZs : int or list - slice indices (last slice by default )
        - zs  : float or list - selected thickness (in A) to show(takes preference over iZs if set)
        - bargs : see beam_vs_thickness
        '''
        nbs,refl,nzs,iZs = self._init_rocking(iZs,zs,**bargs)      #;print(refl)
        I = np.zeros((self.ts.size,nbs,nzs))
        for i,name in enumerate(self.df.index):
            sim_obj = self.load(i)
            refl,z,re,im,Ib = sim_obj.beam_vs_thickness(refl=refl) #;print(refl)
            I[i] = np.array(Ib[:,iZs])

        plts = []
        if nbs>=nzs:
            cs,ms = dsp.getCs(cmap,nbs), dsp.markers
            legElt = { '%s' %refl0:[cs[i],'-'] for i,refl0 in enumerate(refl)}
            for iz,iZ in enumerate(z[iZs]):
                legElt.update({'$z=%d A$' %(iZ):['k',ms[iz]+'-']})
                plts += [[self.ts,I[:,i,iz],[cs[i],ms[iz]+'-'],''] for i in range(nbs)]
        else:
            cs,ms = dsp.getCs(cmap,nzs),  dsp.markers
            legElt = { '%s' %refl0:['k','-'+ms[i]] for i,refl0 in enumerate(refl)}
            for iz,iZ in enumerate(z[iZs]):
                legElt.update({'$z=%d A$' %(iZ):[cs[iz],'-']})
                plts += [[self.ts,I[:,i,iz],[cs[iz],ms[i]+'-'],''] for i in range(nbs)]

        dsp.stddisp(plts,labs=[r'$\theta$(deg)','$I$'],legElt=legElt,**kwargs)

    def _init_rocking(self,iZs,zs,**bargs):
        sim_obj = self.load(0)
        refl,z,re,im,Ib = sim_obj.beam_vs_thickness(**bargs)
        if isinstance(zs,float) or isinstance(zs,float):zs = [zs]
        if isinstance(zs,list) or isinstance(zs,np.ndarray):
            iZs = [np.argmin(abs(z-z0)) for z0 in zs]
        if isinstance(iZs,int):iZs=[iZs]
        nbs = np.array(refl).size
        nzs = z[iZs].size
        return nbs,refl,nzs,iZs


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

def get_uvw(u0,u1,omega,nframes=2,plot=False,h3d=False):
    ''' get a list of continuously rotated vector between u0 and u1.
    u0,u1 : initial and final vector
    omega : array or tuple or int : rotation range to cover
    nframes : number of vectors (omega is int or tuple)
    '''
    if type(omega) in [float,int]:omega=(0,omega)
    if isinstance(omega,tuple):omega=np.linspace(omega[0],omega[1],nframes)
    omega_r = np.deg2rad(omega)
    ct,st = np.cos(omega_r),np.sin(omega_r)
    u0u1  = np.hstack([np.array(u0)[:,None],np.array(u1)[:,None]])
    xy    = np.vstack([ct,st])
    uvw   = u0u1.dot(xy).T
    if plot:show_uvw(uvw,h3d)
    return uvw

def show_uvw(uvw,h3d=False):
    ez =  np.cross(uvw[0],uvw[-1])

    plts = [[[0,u[0]],[0,u[1]],[0,u[2]],'b'] for u in uvw ]
    plts +=[[ [0,ez[0]],[0,ez[1]],[0,ez[2]] , 'r' ]]
    fig,ax = dsp.stddisp(plots=plts,labs=['x','y','z'],rc='3d',lw=4)
    if h3d:h3d = h3D.handler_3d(fig,persp=False,xm0=1)


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
