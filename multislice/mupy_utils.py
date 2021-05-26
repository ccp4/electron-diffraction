import importlib as imp
import pandas as pd, numpy as np
import os,matplotlib,cbf,tifffile,re,glob,easygui
from subprocess import check_output
from matplotlib import rc
from crystals import Crystal
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from utils import displayStandards as dsp   ; imp.reload(dsp)
from utils import glob_colors as colors,handler3D as h3d
from utils import physicsConstants as cst
from scattering.structure_factor import structure_factor3D
from . import rotating_crystal as rcc       #; imp.reload(rcc)
from . import postprocess as pp             #; imp.reload(pp)
from . import multi_3D as MS3D              #;imp.reload(MS3D)
from . import pymultislice                  #;imp.reload(pymultislice)

cs = {1:colors.unicolor(0.75),3:(1,0,0),
      6:colors.unicolor(0.25),7:(0,0,1),8:(1,0,0),16:(1,1,0),14:(1,0,1),17:(0,1,1)}

################################################################################
#### Multislice related
################################################################################
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


def sweep_var(datpath,param,vals,df=None,ssh='',tail='',pargs=True,do_prev=0,
    **kwargs):
    from . import multislice as mupy            ; imp.reload(mupy)
    '''
    runs a set of similar simulations with one varying parameter
    - datpath    : path to the simulation folder
    - param,vals : the parameters and values to sweep
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

def gen_xyz2(file,xyz,lat_params,n=[0,0,1],theta=0, ff=0,fmt='%.4f',opts=''):
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
    l,m,n = find_xyz(lat_vec,lat_params,n_u,theta,plot='p' in opts,v='v' in opts)
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
    header = 'one unit cell of %s at %s, %.1f degree\n' %(dsp.basename(file),str(n_u),theta)
    header+= ' '.join([fmt]*3) %(ax,by,cz)
    np.savetxt(xyz,pattern,footer='-1',header=header,fmt='%d '+' '.join([fmt]*5),comments='')
    print(colors.green+"coords file saved : \n"+colors.yellow+xyz+colors.black)
    npy_file = xyz.replace('.xyz','.npy')
    np.save(npy_file,[(ax,by,cz),pattern])
    print(colors.green+'binary coords file saved :'+colors.yellow+npy_file+colors.black)
    # if 'v' in opts:
    print('number of coordinates = %d' %pattern.shape[0])
    print('number of unit cells  = %d' %l.shape[0])

def find_xyz(lat_vec,lat_params,n_u,theta,plot=0,v=0):
    ax,by,cz = lat_params
    rot_vec = rcc.orient_crystal(lat_vec,n_u=n_u,T=True,theta=theta)
    # ra1,ra2,ra3 = rot_vec
    #### brute force unit cells generation
    N = np.ceil(1.5*max(lat_params)/min(np.linalg.norm(lat_vec,axis=0)))
    u1 = np.arange(-N,N+1)
    l,m,n = np.meshgrid(u1,u1,u1)
    #### l*a1+m*a2+n*a3
    lmn = np.vstack([l.flatten(),m.flatten(),n.flatten()]).T
    if v:print('...\torienting full mesh N=%d...' %N)
    xyz  = lmn.dot(rot_vec)
    x,y,z = xyz.T
    if v:print('...\tkeeping valid unit cells...')
    idx = (x>0) & (y>0) & (z>0) & (x<ax) & (y<by) & (z<cz)
    lmn = lmn[idx]

    if plot:
        print('N',N)
        print('rot_vec',rot_vec)
        l,m,n = lmn.T
        x,y,z = xyz[idx].T
        print('lminmax',l.min(),l.max(),m.min(),m.max(),n.min(),n.max())

        #### generate super cell corners
        l0,m0,n0 = np.meshgrid([0,1],[0,1],[0,1])
        ax,by,cz = np.diag(lat_params)
        sc = np.vstack([ax*i+by*j+cz*k for i,j,k in zip(l0.flatten(),m0.flatten(),n0.flatten())])
        x0,y0,z0 = sc.T
        #### get actual corners
        lat_p2 = np.linalg.norm(rot_vec,axis=0)**2
        nlm1 = np.array([np.round(rot_vec.dot(v)/lat_p2) for v in sc])
        print('corners' ,nlm1)
        x1,y1,z1 = nlm1.dot(rot_vec).T

        scat = ()
        scat+=([x0,y0,z0,50,'b','o'],)
        scat+=([x0,y0,z0,50,'g','d'],)
        scat+=([x,y,z,20,'r','s'],)
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
    pattern[:,1:4] = rcc.orient_crystal(pattern[:,1:4],n_u=n) #,lat_params #pattern,crys # file
    return pattern

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

def get_unit_cell(cell_mesh,n):
    alpha,lw,c = 0.1,2,'c'
    x,y,z = cell_mesh
    planes_coords,ncells = [], x.shape[0]
    for i in range(ncells):
        planes_coords+=[[x[i,:,:],y[i,:,:],z[i,:,:]],
                        [x[:,i,:],y[:,i,:],z[:,i,:]],
                        [x[:,:,i],y[:,:,i],z[:,:,i]]]
    surfs = []
    for plane_coords in planes_coords :
        x,y,z  = plane_coords
        p_shape = x.shape
        coords = np.array([x.flatten(),y.flatten(),z.flatten()]);#print(coords.shape)
        x,y,z  = rcc.orient_crystal(coords,n_u=n,T=False);
        x,y,z  = np.reshape(x,p_shape),np.reshape(y,p_shape),np.reshape(z,p_shape)
        surfs += [[x,y,z,c,alpha,lw,c]]
    return surfs

def get_arrow_3d(n,x0,rc=0.1,h=0.2):
    '''for trihedron'''
    nu,nv = 15,2; #print(u)
    shape = (nu,nv)
    u = np.linspace(0,2*np.pi,nu)
    v = np.linspace(0,np.pi/2,nv)
    x = np.outer(np.cos(u),np.sin(v))*h/2
    y = np.outer(np.sin(u),np.sin(v))*h/2
    z = np.outer(np.ones(u.shape),np.cos(v))*h
    coords = np.array([x.flatten(),y.flatten(),z.flatten()]);#print(coords.shape)
    #uvec = np.array(n);print(n)
    x,y,z = rcc.orient_crystal(coords,n,[0,0,1],T=False)
    x,y,z = np.reshape(x,shape),np.reshape(y,shape),np.reshape(z,shape)

    O = np.array([0,0,0]);#print(x0)#,x0 + np.array([O,n]))
    xl,yl,zl = (x0 + np.array([O,n])).T ;#print(xl,yl,zl)
    return [x+x0[0]+n[0],y+x0[1]+n[1],z+x0[2]+n[2]],[xl,yl,zl]

def get_vec(n,crys,bopt) :
    if isinstance(n,int) : n = crys.lattice_vectors[n]
    if bopt : n = np.array(n).dot(np.array(crys.lattice_vectors))
    return n


################################################################################
#### display coordinate file
################################################################################
def show_cell(file,n=[0,0,1],bopt=1,x0=None,rep=[1,1,1],h3D=1,
    xylims=None,axPos=[],**kwargs):
    '''Show unit cell and coords from cif file '''
    crys = Crystal.from_cif(file)
    n    = get_vec(n,crys,bopt)

    #unit cells,lattice vectors,atom coordinates
    cell_mesh = crys.mesh(range(rep[0]+1),range(rep[1]+1),range(rep[2]+1))
    surfs     = get_unit_cell(cell_mesh,n)
    uvw       = rcc.orient_crystal(np.array(crys.lattice_vectors),n_u=n,T=True)
    pattern   = import_cif(file,n=n,rep=rep)
    E,X,Y,Z   = pattern[:,:4].T
    scat      = [X,Y,Z, [cs[int(e)] for e in E]]

    #display options
    if x0==None : x0 = -1
    if isinstance(x0,float) or isinstance(x0,int) : x0=np.array([x0]*3)
    if not xylims:
        c1,c2,c3 = (X.min()+X.max())/2,(Y.min()+Y.max())/2,(Z.min()+Z.max())/2
        w = 0.75*max(X.max()-X.min(),Y.max()-Y.min(),Z.max()-Z.min())
        xylims=[c1-w/2,c1+w/2,c2-w/2,c2+w/2,c3-w/2,c3+w/2]
        # print(xylims)
    if not axPos:axPos=[0,0,1,0.95]
    fig,ax=dsp.stddisp(scat=scat,ms=100,surfs=surfs,rc='3d',std=0)
    show_trihedron(ax,uvw=uvw,x0=[0,0,0],cs=['r','g','b'],labs=['$a$','$b$','$c$'],lw=2,rc=0.1,h=0.2)
    show_trihedron(ax,x0=x0,labs=['$x$','$y$','$z$'],lw=2,rc=0.1,h=0.2)
    dsp.stddisp(ax=ax,xylims=xylims,axPos=axPos,
        pOpt='Xpt',**kwargs)

    if h3D:hdl = h3d.handler_3d(fig,persp=False)

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

def plot_grid(pattern,lat_params,opts,popts='pv',xylims=[],**kwargs):
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
        if not xylims : xylims = [0,lat_params[i-1],0,lat_params[j-1]]
        return dsp.stddisp(labs=['$%s$'%x1,'$%s$'%x2],patches=pps,scat=scat,
            xylims=xylims,**kwargs)

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

def show_trihedron(ax,uvw=None,x0=[0,0,0],cs=None,clabs=None,lw=2,rc=0.1,h=0.2,**kwargs):
    '''
    x0 : position of trihedron
    uvw : Nx3 ndarray
    ll,rc,h : length,radius and height af arrow/cones
    cs,labs : colors and labels of cones (default black and None)
    '''
    if not isinstance(uvw,np.ndarray) : uvw = np.identity(3)
    txts,x0=[],np.array(x0)
    if not cs : cs=['k']*uvw.shape[0]
    if clabs :
        xtxt = x0 + 1.1*uvw
        for i in range(len(clabs)):
            xt0,yt0,zt0 = xtxt[i,:]
            txts += [[xt0,yt0,zt0,clabs[i],cs[i][0]]]
    plots,surfs = [],[]
    for u,cu in zip(uvw,cs) :
        (x,y,z),(xl,yl,zl) = get_arrow_3d(u,x0,rc=0.1,h=0.2)
        plots += [[xl,yl,zl,cu]]
        surfs += [[x,y,z,cu[0],None,lw,cu[0]]]
    dsp.stddisp(ax=ax,texts=txts,plots=plots,surfs=surfs,lw=lw,**kwargs)

def show_ewald_sphere(lat_params,lam=0.025,tmax=7,T=0.2,nx=20,ny=10,**kwargs):
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

################################################################################
#### Class viewers
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

class Base_Viewer:
    def __init__(self,figpath='',frame=None,thick=5,cutoff=50,i=0,v=1,pargs={}):#,fig=None,ax=None,):
        self.fig,self.ax = dsp.stddisp()
        cid = self.fig.canvas.mpl_connect('key_press_event', self)
        self.figpath = figpath
        self.nfigs   = self._get_nfigs()
        if frame:
            i = min(max(frame,0),self.nfigs)-1
        else:
            frame = i+1
        self.frame  = frame       #starting image
        self.i      = frame       #starting image
        self.inc    = 1           #increment(use 'p' or 'm' to change)
        self.cutoff = cutoff
        self.thick  = thick
        self.fieldNames  = ['thick','frame','cutoff','inc']
        self.mode = 1
        self.pargs  = pargs
        rc('text', usetex=False)
        self.show()
        if v:self.show_help()

    def __call__(self, event):
        # print(event.key)
        if event.key in ['up','right']:
            self.i=min(self.i+self.inc,self.nfigs-1)
            self.mode=1
        elif event.key in ['left','down']:
            self.i=max(0,self.i-self.inc)
            self.mode=-1
        self.frame=self.i+1

        if event.key=='s':
            dsp.saveFig(self.figpath+'_%s.png' %str(self.i).zfill(3),ax=self.ax)

        #increment rate
        if event.key=='p':
            self.inc=min(self.inc+1,100)        ;print('increment rate : %d' %self.inc)
        elif event.key=='m':
            self.inc=max(1,self.inc-1)          ;print('increment rate : %d' %self.inc)
        elif event.key=='ctrl+r':
            self.inc=1                          ;print('increment rate : %d' %self.inc)

        #brightness
        if event.key=='pageup':
            self.cutoff=min(self.cutoff+5,500)  ;print('cutoff : %d' %self.cutoff)
        elif event.key=='pagedown':
            self.cutoff=max(1,self.cutoff-5)    ;print('cutoff : %d' %self.cutoff)
        elif event.key=='r':
            self.cutoff=50                      ;print('cutoff : %d' %self.cutoff)

        #thickness
        if event.key=='ctrl+t':
            self.thick+=5                       ;print('thickness : %d' %self.thick)
        elif event.key=='ctrl+T':
            self.thick=max(self.thick-5,5)      ;print('thickness : %d' %self.thick)

        if event.key=='h':self.show_help()
        elif event.key=='enter':self.settings()
        keys = self.call(event)

        update_keys = keys+['enter','ctrl+t','ctrl+T','pageup','pagedown','r','left','right','down','up']
        if event.key in update_keys:self.show()

    def settings(self):
        fieldValues = ['%d' %self.__dict__[f] for f in self.fieldNames]
        fieldNames = self.fieldNames.copy()
        self.get_fieldValues(fieldNames,fieldValues)
        dict_fv = multenterbox("Change settings","settings", fieldValues,fieldNames)
        if dict_fv:
            for f in self.fieldNames:
                self.__dict__[f] = int(dict_fv[f])
            self.set_fieldValues(dict_fv)
        self.i = self.frame-1

    def show(self):
        tle = "frame %d, thickness=%d $\AA$" %(self.i+1,self.thick)
        self.ax.cla()
        self.get_im(fig=self.fig,ax=self.ax,title=tle,opt='',**self.pargs)
        self.fig.canvas.draw()

    def show_help(self):
        msg = '''
    'up or right'  : show frame+1
    'down or left' : show frame-1
    ##
    'p' : increase increment rate
    'm' : decrease increment rate
    ##
    'pageup'   : increae cutoff brightness
    'pagedown' : decrease cutoff brightness
    'r'        : reset cutoff brightness
    ##
    'ctrl+t' : increase thickness
    'ctrl+T' : decrease thickness
    ##
    'enter' : change settings
    's' : save image
    'h' : show help
        '''
        print(colors.green+'shortcuts : '+colors.blue+msg+colors.black)
    ###################################
    ##### virtual functions
    ###################################
    def _get_nfigs(self):
        print('base nfigs')
        return 0
    def get_im(self,**kwargs):
        print('base get_im')
    def get_fieldValues(self,fieldNames,fieldValues):return None
    def set_fieldValues(self,dict_fv):return None
    def call(self,event):return []


class Frames_Viewer(Base_Viewer):
    '''Viewer for pets'''
    def __init__(self,pets,thick,Imag,kargs,**sargs):
        self.pets   = pets
        self.kargs  = kargs
        super().__init__(thick=thick,cutoff=Imag,**sargs)

    def get_fieldValues(self,fieldNames,fieldValues):
        fieldNames+=list(self.kargs.keys())
        fieldValues+=[str(f) for f in self.kargs.values()]
    def set_fieldValues(self,dict_fv):
        float_keys = np.setdiff1d(list(self.kargs.keys()),'opts')
        for f in float_keys:
            self.kargs[f] = eval(dict_fv[f])
        self.kargs['opts']=dict_fv['opts']
    def _get_nfigs(self):
        return self.pets.nFrames
    def get_im(self,**kwargs):
        self.pets.show_frame(frame=self.i+1,thick=self.thick,Imag=self.cutoff,
            **self.kargs,**kwargs)
    def call(self,event):
        chars = 'KkPhr'
        keys = ['ctrl+K','ctrl+k', 'ctrl+H','ctrl+h','ctrl+g']

        vals = [c in self.kargs['opts'] for c in chars]
        for i,(c,k) in enumerate(zip(chars,keys)):
            if event.key==k:
                vals[i] = not vals[i]

        self.kargs['opts']='q'+''.join([c for c,val in zip(chars,vals) if val])

        return keys

class Viewer(Base_Viewer):
    '''similar to adxv. Works with raw cbf/tiff'''
    def __init__(self,exp_path,v=1,**sargs):
        ''' View cbf files
        - exp_path : path to images
        - figpath : place to save the figures
        - i : starting image
        '''
        d_fmt = {'cbf':self.load_cbf,'tiff':self.load_tif,'tif':self.load_tif}
        self.supported_fmts = d_fmt.keys()

        self.exp_path = os.path.realpath(exp_path)
        self.fmt      = self.find_format(v)
        self.figs     = np.sort(glob.glob(self.exp_path+'/*.%s' %self.fmt))#;print(self.figs)
        self.load     = d_fmt[self.fmt]

        super().__init__(v=v,**sargs)

    def get_im(self,**kwargs):
        fig = self.figs[self.i]
        figname = os.path.basename(fig)
        print(colors.yellow+fig+colors.black)
        im = self.load(fig)
        dsp.stddisp(im=[im],cmap='gray',caxis=[0,self.cutoff],pOpt='t',**kwargs)

    def _get_nfigs(self):
        return self.figs.size

    ###############################################################
    #### misc
    ###############################################################
    def find_format(self,v=1):
        fmts = np.unique([f.split('.')[-1] for f in os.listdir(self.exp_path)])
        fmts = [fmt for fmt in fmts if fmt in self.supported_fmts]
        if not len(fmts):
            raise Exception('no supported format found in %s. Supported formats :' %(self.exp_path),self.supported_fmts)
        fmt = fmts[0]
        if len(fmts)>1:
            print('warning multiple formats found',fmts)
            print('using %s' %fmt)
        if v:print('%s format detected' %fmt)
        return fmt

    def load_cbf(self,fig):
        try:
            content = cbf.read(fig)
        except:#UnicodeDecodeError
            self.i=self.i+self.mode
            print(colors.red+'error reading file'+colors.black)
            self.import_exp()
            return
        numpy_array_with_data = content.data
        header_metadata = content.metadata
        # print(colors.blue,header_metadata,colors.black)
        return numpy_array_with_data

    def load_tif(self,fig):
        return tifffile.imread(fig)


def multenterbox(msg,title,fieldValues,fieldNames):
    fieldValues = easygui.multenterbox(msg, title, fieldNames,fieldValues)
    while True:
        if fieldValues is None:
            break
        errs = list()
        for n, v in zip(fieldNames, fieldValues):
            if v.strip() == "":errs.append('"{}" is a required field.'.format(n))
        if not len(errs):
            break
        fieldValues = easygui.multenterbox("\n".join(errs), title, fieldNames, fieldValues)
    return dict(zip(fieldNames,fieldValues))
