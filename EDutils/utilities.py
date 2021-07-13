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

    def set_tag(self,tag):
        nts   = self.ts.size
        pad   = int(np.ceil(np.log10(nts)))
        names = [ '%s_%s%s' %(tag,'u',str(i).zfill(pad)) for i in range(nts)]
        self.df.index = names
        cmd ="cd %s; rename 's/%s/%s/' %s*.pkl df_%s.pkl rock_%s.pkl" %(self.path,self.tag,tag,self.tag,self.tag,self.tag)
        Popen(cmd,shell=True)
        self.tag = tag
        for i,name in enumerate(names):
            sim_obj = load_pkl(os.path.join(self.path,name+'.pkl'))
            sim_obj.set_name(name=name,path=self.path)
            sim_obj.save()
            self.df.loc[name,'pkl'] = sim_obj.get_pkl()
        df_file = os.path.join(self.path,'df_%s.pkl' %tag)
        self.df.to_pickle(df_file)
        print(colors.green+'Dataframe saved : '+colors.yellow+df_file+colors.black)
        self.save()

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

    def integrate_rocking(self,cond='',refl=[]):
        refl,nbs = self._get_refl(cond=cond,refl=refl)
        z,nzs = self.get_z()

        nts = self.ts.size
        self.Iz     = dict(zip([str(h) for h in refl], np.zeros((nbs,nzs)) ))
        self.Iz_kin = dict(zip([str(h) for h in refl], np.zeros((nbs,nzs)) ))
        for i in range(nts):
            sim_obj = self.load(i)
            idx  = sim_obj.get_beam(refl=refl,cond='')
            if idx:
                hkl0 = [str(tuple(h)) for h in sim_obj.get_hkl()[idx]]
                for idB,hkl_0 in zip(idx,hkl0):
                    self.Iz[hkl_0] += sim_obj.Iz[idB,:]
                    self.Iz_kin[hkl_0] += sim_obj.Iz_kin[idB,:]
        self.save()

    ###############
    #### Display
    ###############
    def plot_integrated(self,cond='',refl=[],cm='Spectral',**kwargs):
        if not 'Iz' in self.__dict__:self.integrate_rocking(cond=cond,refl=refl)
        refl,nbs = self._get_refl(cond=cond,refl=refl)
        refl = [str(tuple(h)) for h in refl]
        z = self.load(0).z
        cs = dsp.getCs(cm,nbs)
        plts = [[z,self.Iz[h],cs[i],'%s' %h] for i,h in enumerate(refl)]
        dsp.stddisp(plts,labs=[r'$z(\AA)$','$I_{int}$'],**kwargs)

    def plot_rocking(self,iZs=-1,zs=None,bargs={},cmap='viridis',**kwargs):
        '''plot rocking curve for set of selected beams at thickness zs
        - iZs : int or list - slice indices (last slice by default )
        - zs  : float or list - selected thickness (in A) to show(takes preference over iZs if set)
        - bargs : see beam_vs_thickness
        '''
        iZs,nzs  = self._get_iZs(iZs,zs)    #;print(iZs)
        refl,nbs = self._get_refl(**bargs)#;print(refl)
        nts = self.ts.size
        # I = np.zeros((self.ts.size,nbs,nzs))
        I = dict(zip([str(h) for h in refl], [np.nan*np.ones((nts,nzs))]*nbs))
        for i in range(nts):
            sim_obj = self.load(i)
            idx  = sim_obj.get_beam(refl=refl,cond='')
            if idx:
                hkl0 = [str(tuple(h)) for h in sim_obj.get_hkl()[idx]]
                for idB,hkl_0 in zip(idx,hkl0):
                    I[hkl_0][i,:] = np.array(sim_obj.Iz[idB,iZs])

        refl = list(I.keys())
        z = self.load(0).z
        plts = []
        if nbs>=nzs:
            cs,ms = dsp.getCs(cmap,nbs), dsp.markers
            legElt = { '%s' %refl0:[cs[i],'-'] for i,refl0 in enumerate(refl)}
            for iz,iZ in enumerate(z[iZs]):
                legElt.update({'$z=%d A$' %(iZ):['k',ms[iz]+'-']})
                plts += [[self.ts,I[refl0][:,iz],[cs[i],ms[iz]+'-'],''] for i,refl0 in enumerate(refl)]
        else:
            cs,ms = dsp.getCs(cmap,nzs),  dsp.markers
            legElt = { '%s' %refl0:['k','-'+ms[i]] for i,refl0 in enumerate(refl)}
            for iz,iZ in enumerate(z[iZs]):
                legElt.update({'$z=%d A$' %(iZ):[cs[iz],'-']})
                plts += [[self.ts,I[refl0][:,iz],[cs[iz],ms[i]+'-'],''] for i,refl0 in enumerate(refl)]

        dsp.stddisp(plts,labs=[r'$\theta$(deg)','$I$'],legElt=legElt,**kwargs)

    def QQplot(self,zs=None,iZs=10,refl=[],cmap='Spectral',**kwargs):
        if not 'Iz' in self.__dict__:self.integrate_rocking(cond='I>1e-4')#,refl=refl)
        iZs,nzs  = self._get_iZs(iZs,zs)    #;print(iZs)
        z  = self.load(0).z.copy()[iZs]

        if not any(refl):
            refl = list(self.Iz.keys());refl.remove(str((0,)*3))
        else:
            refl = [str(tuple(h)) for h in refl]
        Iz_dyn = np.array([self.Iz[    h].copy()[iZs] for h in refl])
        Iz_kin = np.array([self.Iz_kin[h].copy()[iZs] for h in refl])
        # print(Iz_dyn.shape)
        iB = np.argsort(np.sum(Iz_dyn,axis=1))[-1]
        Iz_dyn/= Iz_dyn[iB,:]
        Iz_kin/= Iz_kin[iB,:]
        cs = dsp.getCs(cmap,nzs) #; print(len(cs),Iz_dyn.shape,Iz_kin.shape)

        plts=[[I_kin,I_dyn,[cs[i],'o'],r'$z=%d \AA$' %z0] for i,(z0,I_dyn,I_kin) in enumerate(zip(z,Iz_dyn.T,Iz_kin.T))]
        plts+=[ [[0,1],[0,1],[(0.5,)*3,'--'],''] ]
        dsp.stddisp(plts,labs=['$I_{kin}$','$I_{dyn}$'],sargs={'alpha':0.5},**kwargs)


    def Sw_vs_theta(self,refl=[[0,0,0]],cond='',thick=None,fz=abs,opts='',
        iTs=slice(0,None),ts=None,
        cm='Spectral',figname='',**kwargs):
        '''Displays Sw and I for a range of angle simulations at given thickness
        - refl,cond : selected reflections
        - thick : thickness
        - fz : functor for Sw
        - Iopt : plot I
        '''
        Iopt = 'I' in opts
        iTs,nts = self._get_iTs(iTs,ts)
        xlab,ts = r'$\theta$',self.ts.copy()[iTs]
        if 'f' in opts:
            xlab,ts = 'frame',np.arange(1,self.ts.size+1)[iTs]       #;print(iTs,ts)

        if thick and Iopt:
            if not self.load(0).thick==thick:self.do('set_thickness',thick=thick)
        refl,nbs = self._get_refl(cond,refl)            #;print(nbs)

        Sw = pd.DataFrame(np.ones((nts,nbs)),columns=[str(h) for h in refl])
        if Iopt:I  = pd.DataFrame(np.zeros((nts,nbs)),columns=[str(h) for h in refl])
        for i,name in enumerate(self.df.index[iTs]):
            b = self.load(i) #;print(i)
            idx = b.get_beam(refl=refl,cond=cond)
            hkl0 = [str(tuple(h)) for h in b.get_hkl()[idx]]
            Sw.loc[i,hkl0] = b.df_G.loc[idx,'Sw'].values
            if Iopt:I.loc[i,hkl0] = b.df_G.loc[idx,'I'].values

        #locate minimum excitation errors
        iSmin = np.argmin(Sw.values.T,axis=1) #locate minimums
        Sw[Sw==1]=np.nan
        SwE = fz(Sw.values.T)

        if 'i' in opts:
            dsp.stddisp(im=[SwE],pOpt='im',labs=[xlab,'$beam$'],caxis=[fz(1e-2),fz(1e-6)],
                cmap='Reds',title='Excitation error',name=figname+'_Sw.svg',**kwargs)
            # print(refl)

        else:
            cs,txts = dsp.getCs(cm,nbs),[]
            plts = [[ts,Sw0,[cs[i],'-o'],''] for i,Sw0 in enumerate(SwE)]
            if 't' in opts:txts = [[ts[idx],SwE[i,idx],'%s' %str(h),cs[i]] for i,(h,idx) in enumerate(zip(refl,iSmin))]
            dsp.stddisp(plts,texts=txts,labs=[xlab,'$S_w$'],name=figname+'_Sw.svg',**kwargs)

        if Iopt:
            IE = I.values.T
            if 'i' in opts:
                dsp.stddisp(im=[IE],pOpt='im',labs=[xlab,'$beam$'],caxis=[0,0.2],
                    cmap='YlGnBu',title='Intensity',name=figname+'_Sw.svg',**kwargs)
            else:
                plts = [[ts,I0,[cs[i],'-o'],''] for i,I0 in enumerate(IE)]
                txts = [[ts[idx],IE[i,idx],'%s' %str(h),cs[i]] for i,(h,idx) in enumerate(zip(refl,iSmin))]
                dsp.stddisp(plts,texts=txts,labs=[xlab,'$I$'],
                    title=r'thickness=$%d A$' %thick,name=figname+'_I.svg',**kwargs)
            # return IE
        return Sw

    ##############################
    #### gets
    ##############################
    def _get_refl(self,cond,refl):
        if cond:
            refl = []
            for i,name in enumerate(self.df.index):
                b = self.load(i)
                idx = b.get_beam(cond=cond)
                hkl = b.get_hkl()[idx]
                refl += [tuple(h) for h in hkl]
            refl = np.unique(refl,axis=0)           #;print(refl)
        if not isinstance(refl[0],tuple):refl=[tuple(h) for h in refl]
        nbs = len(refl)#.size;print(nbs,refl)
        return refl,nbs

    def _get_iTs(self,iTs,ts):
        t = self.ts
        if isinstance(ts,float) or isinstance(ts,int):ts = [ts]
        if isinstance(ts,list) or isinstance(ts,np.ndarray):
            iTs = [np.argmin(abs(t-t0)) for t0 in ts]
        if isinstance(iTs,int):iTs=[iTs]
        if not type(iTs) in [list,np.ndarray,slice]:iTs = slice(0,None,1)
        nts = t[iTs].size
        return iTs,nts

    def _get_iZs(self,iZs,zs):
        z = self.load(0).z
        if isinstance(zs,float) or isinstance(zs,int):zs = [zs]
        if isinstance(zs,list) or isinstance(zs,np.ndarray):
            iZs = [np.argmin(abs(z-z0)) for z0 in zs]
        if isinstance(iZs,int):iZs=[iZs]
        nzs = z[iZs].size
        return iZs,nzs

    def get_z(self):
        z = self.load(0).z
        nzs = z.size
        return z,nzs


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
    u0 = u_from_theta_phi(theta,phi)
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
