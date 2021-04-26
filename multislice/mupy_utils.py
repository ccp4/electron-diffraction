import importlib as imp
import pandas as pd, numpy as np
import os,matplotlib,cbf,re,glob
from subprocess import check_output
from utils import displayStandards as dsp   ; imp.reload(dsp)
from utils import glob_colors as colors
from crystals import Crystal
from matplotlib import rc
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from . import multislice as mupy            ; imp.reload(mupy)
from . import rotating_crystal as rcc       #; imp.reload(rcc)
from . import postprocess as pp             #; imp.reload(pp)
from . import multi_3D as MS3D              ;imp.reload(MS3D)
from . import pymultislice                  ;imp.reload(pymultislice)

cs = {1:colors.unicolor(0.75),3:(1,0,0),
      6:colors.unicolor(0.25),7:(0,0,1),8:(1,0,0),16:(1,1,0),14:(1,0,1),17:(0,1,1)}

def multi3D(name='./unknwon',filename='',load_opt=0,**kwargs):
    if filename :
        name=filename.replace('_3D.pkl','')
    else :
        filename=name+'_3D.pkl'

    if load_opt :
        if os.path.exists(filename):
            print(colors.green+'loading '+colors.yellow+filename+colors.black)
            return pymultislice.load(filename)
        else:
            print(colors.red+'file not found : '+colors.yellow+filename+colors.black)
    return MS3D.Multi3D(name=name,**kwargs)



def get_reciprocal(abc):
    a1,a2,a3 = abc
    b1 = np.cross(a2,a3)/(a1.dot(np.cross(a2,a3)))
    b2 = np.cross(a3,a1)/(a2.dot(np.cross(a3,a1)))
    b3 = np.cross(a1,a2)/(a3.dot(np.cross(a1,a2)))
    abc_star = np.vstack([b1,b2,b3])#/(2*np.pi)
    return abc_star

def ewald_sphere(lat_params,lam=0.025,tmax=7,T=0.2,nx=20,ny=10,**kwargs):
    '''lam(wavelength Angstrum), tmax(angle degrees), T(Thickness mum)'''
    a1,a2,a3 = lat_params
    b1,b2,b3 = 1/a1,1/a2,1/a3
    tmax, K=np.deg2rad(tmax), 1/lam
    dqT = 10/(T*1e4)*np.array([-1,1])
    #reciprocal lattice
    h,k = np.meshgrid(b1*np.arange(-nx,nx+1),b2*np.arange(-1,ny))
    #sphere
    t = 3*np.pi/2+np.linspace(-tmax,tmax,100)
    plts = [[K*np.cos(t),K*np.sin(t)+K,'r',r'$\lambda=%.2f\AA$' %lam]]
    #rel rods
    plts += [[ [h,h],k+dqT ,'b-',''] for h,k in zip(h.flatten(),k.flatten())]

    #display
    scat=[h,k,15,'b']
    fig,ax = dsp.stddisp(plts,scat=scat,labs=['$q_x$','$q_y$'],
        lw=3,#,xylims=[-nx*b1,nx*b1,-b2,ny*b2],xyTickLabs=[[],[]],
        **kwargs)


def sweep_var(datpath,param,vals,df=None,ssh='',tail='',pargs=True,do_prev=0,
    **kwargs):
    '''
    runs a set of similar simulations with one varying parameter
    - datpath       : path to the simulation folder
    - param,vals    : the parameters and values to sweep
    - df :
        - int - create and save the new dataframe if 1
        - pd.Dataframe - to update(since parsed as a reference)
    - pargs         : bool - pass param and values to Multislice  if True
    - do_prev       : Used for iterative fourier transform
    - kwargs : see help(Multislice)
    '''
    do_df,save = isinstance(df,pd.core.frame.DataFrame),0
    if isinstance(df,int):
        if df : df,do_df,save = pd.DataFrame(columns=[param,'host','state']+pp.info_cols),1,1
    nvals,prev = len(vals),None
    for i,val in zip(range(nvals),vals):
        if pargs:kwargs[param]=val
        if do_prev and i: prev = multi.outf['image']
        multi=mupy.Multislice(datpath,prev=prev,
            ssh=ssh,tail=tail+param+str(i).zfill(int(np.ceil(nvals/10))),
            **kwargs)
        if do_df:
            df.loc[multi.outf['obj']] = [np.nan]*len(df.columns)
            df.loc[multi.outf['obj']][[param,'host','state']] = [val,ssh,'start']
    if save :
        df.to_pickle(datpath+'df.pkl')
        print(colors.green+'Dataframe saved : '+colors.yellow+datpath+'df.pkl'+colors.black)
        return df



################################################################################
#### coordinate file generation from cif files
################################################################################
def import_crys(file):
    if file.split('.')[-1]=='cif':
        crys = Crystal.from_cif(file)
    elif sum(np.array(list(Crystal.builtins))==file):
        crys = Crystal.from_database(file)
    else:
        raise Exception('cannot import %s' %file)
    return crys

def gen_xyz2(file,xyz,lat_params,n=[0,0,1],theta=0,pad=0,fmt='%.4f',opts=''):
    if 'v' in opts:print('...import file...')
    crys = import_crys(file)
    n_u = n
    pattern = np.array([[a.atomic_number]+list(a.coords_cartesian)+[a.occupancy,1.0] for a in crys.atoms])
    lat_vec = np.array(crys.lattice_vectors)
    coords = pattern[:,1:4]
    Za,occ,bfact = pattern[:,[0,4,5]].T
    ax,by,cz = lat_params

    #### replicate
    if 'v' in opts:print('...finding unit cells index...')
    l,m,n = find_xyz(lat_vec,lat_params,n_u,theta,plot='p' in opts)
    a1,a2,a3 = np.array(lat_vec)
    coords = np.vstack([coords + i*a1+j*a2+k*a3 for i,j,k in zip(l,m,n)])
    Za     = np.tile(Za   ,[n.size])
    occ    = np.tile(occ  ,[n.size])
    bfact  = np.tile(bfact,[n.size])

    ###orient
    coords = rcc.orient_crystal(coords,n_u=n_u,theta=theta)

    #### apply padding
    if isinstance(pad,int) or isinstance(pad,float):pad=[pad]*2
    if sum(pad)>0:
        coords[:,0] += ax*pad[0]
        coords[:,1] += by*pad[1]
        ax *= 1+2*pad[0]
        by *= 1+2*pad[1]

    pattern = np.hstack([Za[:,None],coords,occ[:,None],bfact[:,None]])
    #### save
    if 'v' in opts:print('...saving to file ...')
    dir=''.join(np.array(n,dtype=str))
    header = 'one unit cell of %s at %s, %.1f degree\n' %(dsp.basename(file),str(n),theta)
    header+= ' '.join([fmt]*3) %(ax,by,cz)
    np.savetxt(xyz,pattern,footer='-1',header=header,fmt='%d '+' '.join([fmt]*5),comments='')
    print(colors.green+"coords file saved : \n"+colors.yellow+xyz+colors.black)
    npy_file = xyz.replace('.xyz','.npy')
    np.save(npy_file,[(ax,by,cz),pattern])
    print(colors.green+'binary coords file saved :'+colors.yellow+npy_file+colors.black)
    # if 'v' in opts:
    print('number of coordinates = %d' %pattern.shape[0])
    print('number of unit cells  = %d' %l.shape[0])
    # if 'p' in dopt :
    #     with open(name,'r') as f : print(''.join(f.readlines()))

def find_xyz(lat_vec,lat_params,n_u,theta,plot=0):
    ax,by,cz = lat_params
    rot_vec = rcc.orient_crystal(lat_vec,n_u=n_u,T=True,theta=theta)
    # ra1,ra2,ra3 = rot_vec

    #### brute force unit cells generation
    N = np.ceil(1.5*max(lat_params)/min(np.linalg.norm(lat_vec,axis=0)))
    u1 = np.arange(-N,N+1)
    l,m,n = np.meshgrid(u1,u1,u1)
    #### l*a1+m*a2+n*a3
    lmn = np.vstack([l.flatten(),m.flatten(),n.flatten()]).T
    xyz  = lmn.dot(rot_vec)
    x,y,z = xyz.T
    idx = (x>0) & (y>0) & (z>0) & (x<ax) & (y<by) & (z<cz)
    lmn = lmn[idx]

    if plot:
        l,m,n = lmn.T
        x,y,z = xyz[idx].T
        # print(l.min(),l.max(),m.min(),m.max(),n.min(),n.max())

        #### generate super cell corners
        l0,m0,n0 = np.meshgrid([0,1],[0,1],[0,1])
        ax,by,cz = np.diag(lat_params)
        sc = np.vstack([ax*i+by*j+cz*k for i,j,k in zip(l0.flatten(),m0.flatten(),n0.flatten())])
        x0,y0,z0 = sc.T
        # nlm1 = np.array([np.round(rot_vec.dot(v)/lat_p2) for v in sc])
        #### get actual corners
        # x1,y1,z1 = nlm1.dot(rot_vec).T

        scat = ()
        scat+=([x0,y0,z0,50,'b','o'],)
        scat+=([x,y,z,50,'r','s'],)
        dsp.stddisp(rc='3d',scat=scat)
    return lmn.T


def gen_xyz(file,n=[0,0,1],rep=[1,1,1],pad=0,xyz='',**kwargs):
    ''' convert cif file into autoslic .xyz input file
    - file : cif_file
    - rep : super cell repeat
    - n : reorient of the z axis into n
    - pad : amount of padding on each side (in unit of super cell size)
    '''
    if isinstance(file,str):
        tail = ''.join(np.array(n,dtype=str)) + '.xyz'
        if file.split('.')[-1]=='cif':
            crys = Crystal.from_cif(file)
            if not xyz : xyz = file.replace('.cif',tail)
        elif sum(np.array(list(Crystal.builtins))==file):
            crys = Crystal.from_database(file)
            if not xyz : xyz = file+tail
        lat_vec    = np.array(crys.lattice_vectors)
        lat_params = crys.lattice_parameters[:3]
        pattern    = np.array([[a.atomic_number]+list(a.coords_cartesian)+[a.occupancy,1.0] for a in crys.atoms])
    else:
        if not xyz: raise Exception('xyz filename required')
        pattern = file
    pattern,lat = make_xyz(xyz,pattern,lat_vec,lat_params,n=n,pad=pad,rep=rep,**kwargs)
    npy_file = xyz.replace('.xyz','.npy')
    np.save(npy_file,[lat,pattern])
    print(colors.green+'binary coords file saved :'+colors.yellow+npy_file+colors.black)

def import_cif(file,xyz='',n=[0,0,1],rep=[1,1,1],pad=0,dopt='s',lfact=1.0,tail=''):
    ''' convert cif file into autoslic .xyz input file
    - file : cif_file
    - rep : super cell repeat
    - n : reorient of the z axis into n
    - pad : amount of padding on each side (in unit of super cell size)
    '''
    crys       = Crystal.from_cif(file)
    lat_vec    = np.array(crys.lattice_vectors)
    lat_params = crys.lattice_parameters[:3]
    pattern    = np.array([[a.atomic_number]+list(lfact*a.coords_cartesian)+[a.occupancy,1.0] for a in crys.atoms])

    if xyz:make_xyz(xyz,pattern,lat_vec,lat_params,n=n,pad=pad,rep=rep,fmt='%.4f',dopt=dopt)
    return pattern #,lat_params #pattern,crys # file

def make_xyz(name,pattern,lat_vec,lat_params,n=[0,0,1],theta=0,rep=[1,1,1],pad=0,fmt='%.4f',dopt='s'):
    '''Creates the.xyz file from a given compound and orientation
    - name       : Full path to the file to save
    - pattern    : Nx6 ndarray - Z,x,y,z,occ,wobble format
    - lat_vec    : 3x3 ndarray - lattice vectors [a1,a2,a3]
    - lat_params : ax,by,cz
    - pad  : amount of padding on each side (in unit of super cell size)
    - n    : beam direction axis
    - dopt : p(print file),s(save)
    '''
    compound = dsp.basename(name)
    Za,occ,bfact = pattern[:,[0,4,5]].T
    coords = pattern[:,1:4]
    Nx,Ny,Nz = rep
    ax0,by0,cz0 = lat_params
    ax,by,cz = lat_params
    #replicate
    if sum(rep)>3 :
        ni,nj,nk = np.meshgrid(range(Nx),range(Ny),range(Nz))
        ni,nj,nk = ni.flatten(),nj.flatten(),nk.flatten()
        a1,a2,a3 = np.array(lat_vec)
        coords = np.vstack([coords + i*a1+j*a2+k*a3 for i,j,k in zip(ni,nj,nk)])
        Za     = np.tile(Za   ,[ni.size])
        occ    = np.tile(occ  ,[ni.size])
        bfact  = np.tile(bfact,[ni.size])
        ax = Nx*ax0
        by = Ny*by0
        cz = Nz*cz0

    #orient
    coords = rcc.orient_crystal(coords,n_u=n)
    if theta :
        t = np.deg2rad(theta)
        st,ct = np.sin(t),np.cos(t)
        coords = np.array([[ct,st,0],[-st,ct,0],[0,0,1]]).dot(coords.T).T
    #apply padding
    if isinstance(pad,int) or isinstance(pad,float):pad=[pad]*3
    if sum(pad)>0:
        coords[:,0] += Nx*ax0*pad[0]
        coords[:,1] += Ny*by0*pad[1]
        coords[:,2] += Nz*cz0*pad[2]
        ax *= 1+2*pad[0]
        by *= 1+2*pad[1]
        cz *= 1+2*pad[2]
    pattern = np.hstack([Za[:,None],coords,occ[:,None],bfact[:,None]])
    # ax,by,cz = lat_params
    lat_params = (ax,by,cz)
    #write to file
    if 's' in dopt :
        dir=''.join(np.array(n,dtype=str))
        header = 'one unit cell of %s\n' %(compound)
        header+= ' '.join([fmt]*3) %(ax,by,cz)
        np.savetxt(name,pattern,footer='-1',header=header,fmt='%d '+' '.join([fmt]*5),comments='')
        print(colors.green+"coords file saved : \n"+colors.yellow+name+colors.black)
        if 'p' in dopt :
            with open(name,'r') as f : print(''.join(f.readlines()))
    return pattern,[ax,by,cz]

def make_mulslice_datfile(dat_file,cif_file):
    ''' create a data file used by atompot(mulslice) from a cif file
    - dat_file : name of file to save
    - cif_file : the file to import
    '''
    crys = Crystal.from_cif(cif_file)
    pattern = np.array([np.hstack([a.coords_cartesian,a.atomic_number]) for a in crys.atoms] )
    deck = ' '.join(np.array(crys.lattice_parameters[:3],dtype=str))+'\n'
    deck+='0\n'
    Za = np.array(np.unique(pattern[:,-1]),dtype=int)
    for Z in Za:
        deck+=str(Z)+'\n'#; print()
        xyz = np.array(pattern[pattern[:,-1]==1][:,:3],dtype=str)
        xyz = '\n'.join([' '.join(l) for l in xyz]) #; print(xyz)
        deck+=xyz+'\n'
    deck+='\n\nComment here'
    with open(dat_file,'w') as f : f.write(deck)
    # with open(datpath,'r') as f :print(f.readlines())


################################################################################
#### display coordinate file
################################################################################
def show_grid(file,opts='',popts='pv',figs='21',**kwargs):
    '''
    file : path to .xyz file
    opt : str format 'x1x2' - 'xy','xz','yz','zx',...
    '''
    if 'v' in popts:print('...loading file...')
    npy_file = file.replace('.xyz','.npy')
    if os.path.exists(npy_file):
        lat_params,pattern = np.load(npy_file,allow_pickle=True)
    else:
        with open(file,'r') as f:
            l=list(map(lambda s:s.strip().split(' '),f.readlines()))
            lat_params = np.array(l[1],dtype=float)
            pattern   = np.array(l[2:-1],dtype=float)
            np.save(npy_file,[lat_params,pattern])
            print(colors.green+'binary coords file saved :'+colors.yellow+npy_file+colors.black)

    if isinstance(opts,list):
        fig,axs = dsp.create_fig(figsize=figs,rc=[1,len(opts)])
        for axi,opts_i in zip(axs,opts):
            plot_grid(pattern,lat_params,opts_i,popts=popts,ax=axi,setPos=0,opt='',**kwargs)
        fig.show()
    else:
        plot_grid(pattern,lat_params,opts,**kwargs)

def plot_grid(pattern,lat_params,opts,popts='pv',**kwargs):
    #a1,a2,a3,alpha,beta,gamma=crys.lattice_parameters
    xij = {'x':1,'y':2,'z':3}
    vopt = 'v' in popts
    if opts:
        x1,x2 = opts
        i,j = xij[x1],xij[x2]
        Z = np.array(pattern[:,0],dtype=int)
        # print('...assigning colors...')
        C = np.array([cs[E] for E in Z])
        pps = [dsp.Rectangle((0,0),lat_params[i-1],lat_params[j-1],linewidth=2,edgecolor='b',alpha=0.1)]
        scat=[]
        if 'p' in popts:
            scat = [pattern[:,i],pattern[:,j],C]
        if 'h' in popts :
            if vopt:print('...finding convex hull...')
            hull = ConvexHull(pattern[:,[i,j]])#,incremental=True)
            idx = hull.vertices #;print(idx)
            points = np.array([pattern[idx,i],pattern[idx,j]]).T
            pps+=[dsp.matplotlib.patches.Polygon(points,linewidth=2,facecolor='r',edgecolor='r',alpha=0.2)]
        if vopt:print('...plotting...')
        return dsp.stddisp(labs=['$%s$'%x1,'$%s$'%x2],patches=pps,scat=scat,
            xylims=[0,lat_params[i-1],0,lat_params[j-1]],**kwargs)

    # return lat_params,pattern

def import_xyz(xyz_file):
    '''import .xyz file
    - xyz_file : file to import
    returns :
    - pattern    : [Za,coords,occ,bfactor]
    - lat_params : [ax,by,cz] for the supercell
    '''
    with open(xyz_file,'r') as f:l=list(map(lambda s:s.strip().split(' '),f.readlines()))
    lat_params = np.array(l[1],dtype=float)
    pattern = np.array(l[2:-1],dtype=float)
    return pattern,lat_params

def show_grid3(xyz_file,ms=1,**kwargs):
    '''display an .xyz file content in a 2D plane
    xyz_file : .xyz file
    '''
    pattern,lat_params = import_xyz(xyz_file)
    Z = np.array(pattern[:,0],dtype=int)
    C = [cs[E] for E in Z]

    scat = [pattern[:,1],pattern[:,2],pattern[:,3],ms,C]
    return dsp.stddisp(scat=scat,labs=['$x$','$y$','$z$'],rc='3d',**kwargs)


################################################################################
################################################################################
class Image_viewer:
    '''a copy of dials.image_viewer'''
    def __init__(self,file,sym=1,pad=None):
        basename      = file.split('_')
        self.im       = int(basename[-1].replace('.cbf',''))
        self.basename = '_'.join(basename[:-1])
        file_ids      = [int(f.split('_')[-1].replace('.cbf','')) for f in glob.glob(self.basename+'*.cbf')]
        if pad :
            self.pad=pad
        else:
            self.pad=int(np.log10(max(file_ids)))+1       #;print(self.pad)

        ls = np.array(check_output("head -n25 %s" %file,shell=True).decode().split('\n'))
        px_str      = ls[['Pixel_size' in l for l in ls]][0]
        D_str       = ls[['Detector_distance' in l for l in ls]][0] #;print(D_str)
        lam_str     = ls[["Wavelength" in l for l in ls]][0]

        self.pxy    = np.array(re.findall("\d+e-\d+", px_str),dtype=float)
        self.D      = 1 #float(re.findall("\d+.\d+",D_str)[0])
        self.lam    = float(re.findall("\d+.\d+",lam_str)[0])

        shape = np.array(cbf.read(file).data.shape)
        self.image = cbf.read(file).data

        if sym:
            Nx,Ny = np.array(shape/2,dtype=int)
            self.Nx,self.Ny = Nx,Ny
            self.pX,self.pY = np.meshgrid(np.arange(-Nx,Nx+1),np.arange(-Ny,Ny+1))
        else:
            Nx,Ny = shape
            # self.pX,self.pY = np.meshgrid(np.arange(Nx),np.arange(Ny))
            self.pX,self.pY = np.meshgrid(np.arange(Nx),np.flipud(np.arange(Ny)))
        self.qx,self.qy = self.px2s(self.pX,self.pY)

    def s2px(self,xh):
        '''xh : recorded coordinates on image'''
        px,py = self.pxy
        pX = self.lam*self.D/px*xh[:,0] + self.Nx
        pY = -self.lam*self.D/py*xh[:,1] + self.Ny
        return pX,pY

    def px2s(self,pX,pY):
        px,py = self.pxy
        qx,qy = px/self.D/self.lam*pX, py/self.D/self.lam*pY
        return qx,qy

    def show_image(self,im=None,lab='q',stack=1,**kwargs):
        rc('text', usetex=False)
        if not im:im = self.im
        file = "%s_%s.cbf" %(self.basename,str(im).zfill(self.pad))   #;print(filename)
        # with open(file,'wb') as f:print(file)
        content = cbf.read(file)
        image = content.data
        if stack>1:
            for i in range(1,stack):
                filename = "%s_%s.cbf" %(self.basename,str(im+i).zfill(self.pad))   #;print(filename)
                image+=cbf.read(filename).data
            if 'caxis' in list(kwargs.keys()):
                kwargs['caxis'] = list(np.array(kwargs['caxis'])*stack) #; print(kwargs['caxis'])
        if 'q' in lab:
            labs = ['$q_x(A^{-1})$','$q_y(A^{-1})$']
            X,Y = self.qx,-self.qy
        elif 'p' in lab:
            labs = ['','']#['$p_x$','$p_y$']
            X,Y = self.pX,self.pY
        elif 'x' in lab:
            labs = ['$x(mm)$','$y(mm)$']
            px,py = self.pxy
            X,Y = self.pX*px*1e3,self.pY*py*1e3
        return dsp.stddisp(im=[X,Y,image],labs=labs,
            pOpt='ptX',title="%s" %os.path.basename(file),**kwargs)



#######
class Viewer_cbf:
    '''similar to adxv. Works with raw cbf'''
    def __init__(self,exp_path,figpath,i=0):
        ''' View cbf files
        - exp_path : path to images
        - figpath : place to save the figures
        - i : starting image
        '''
        self.exp_path = exp_path
        self.figpath = figpath
        self.figs  = np.sort(os.listdir(exp_path));print(self.figs)
        self.nfigs = self.figs.size
        self.fig,self.ax = dsp.stddisp()
        cid = self.fig.canvas.mpl_connect('key_press_event', self)
        self.i=i     #starting image
        self.inc=1   #increment(use 'p' or 'm' to change)
        self.mode=1
        rc('text', usetex=False)
        self.import_exp()

    def import_exp(self):
        print(colors.yellow+self.figs[self.i]+colors.black)
        try:
            content = cbf.read(self.exp_path+self.figs[self.i])
        except:#UnicodeDecodeError
            self.i=self.i+self.mode
            print(colors.red+'error reading file'+colors.black)
            self.import_exp()
            return
        numpy_array_with_data = content.data
        header_metadata = content.metadata
        print(colors.blue,header_metadata,colors.black)
        self.ax.cla()


        dsp.stddisp(fig=self.fig,ax=self.ax,im=[numpy_array_with_data],
            cmap='gray',caxis=[0,50],pOpt='t',title="image %d:%s" %(self.i,self.figs[self.i]),opt='')
        self.fig.canvas.draw()

    def __call__(self, event):
        # print(event.key)
        if event.key in ['up','right']:
            self.i=min(self.i+self.inc,self.nfigs-1)
            self.mode=1
            self.import_exp()
        elif event.key in ['left','down']:
            self.i=max(0,self.i-self.inc)
            self.mode=-1
            self.import_exp()

        if event.key=='s':
            dsp.saveFig(self.figpath+'exp_%s.png' %str(self.i).zfill(3),ax=self.ax)

        if event.key=='p':
            self.inc=min(self.inc+1,100);print(self.inc)
        if event.key=='m':
            self.inc=max(1,self.inc-1);print(self.inc)
