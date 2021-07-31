import os,pickle5
import numpy as np, pandas as pd
from matplotlib import patches
import matplotlib.pyplot as plt
from utils import displayStandards as dsp
from utils import glob_colors as colors
from scipy.interpolate import interp2d,RegularGridInterpolator
from scipy.spatial import cKDTree
# from matplotlib import collections as coll
# from scipy.spatial import Delaunay
# import matplotlib.tri as mtri
from . import lattice
from . import config as cg
# from . import nodesTri as nt

class Wallpaper :
    '''Class for creating wallpapers'''
    def __init__(self,pp_type,a=1,b=1,alpha=90,pattern=[[0,0,0]],
        gen=True,pOpt='',nh=3,nk=3,ndeg=100,
        name='',path='') :
        ''' initialize wallpaper
        pp_type   : plane group type (str)
        a,b,alpha : lattice parameters (double)
        pattern   : Nx3 array : x y f coordinates in unit cells
        gen       : generate wallpaper applying symmetry operations (see generate_wallpp for ndeg,nh,nk)
        pOpt      : plot options                                    (see plot_unit_cells)
                  : 'w' to see the unit cell if warning is issued
            Note : gen,ndeg,nh,nk will work only if pattern is defined
        ndeg   : number of discretization points along one the direction for building the potential map
        '''
        self.pp_type  = pp_type
        self.name     = name if name else pp_type #cg.df_wallpp.loc[pp_type]['long_name']
        self.path     = path
        self.df = cg.df_wallpp.loc[pp_type]
        self._init_lattice(a,b,alpha)
        self._init_unit_cell()
        # self.ndeg = ndeg
        self.pattern0 = np.array(pattern).copy()
        self.pattern  = cg.generate(self.pp_type,self.lat_vec,pattern)
        self._interp_pattern()
        # if gen:
        if pOpt:
            self.plot_wallpp(pOpt,nh,nk)#,xylims=[0,1,0,0.5])

    def _init_lattice(self,a,b,alpha) :
        self.a,self.b,self.alpha = a,b,alpha
        self.lat_type = cg.df_wallpp.loc[self.pp_type]['lat_type']
        self.lat      = lattice.Lattice2D(self.lat_type,a,b,alpha)
        self.params   = self.lat.params
        self.lat_vec    = self.lat.lat_vec
        a1,a2           = self.lat_vec
        self.unit_cell  = np.array([0*a1,a1,a1+a2,a2])
        self.rec_vec    = self.lat.rec_vec

    def _init_unit_cell(self) :
        verts    = self.df.asym_cell
        a,b,angle = self.params
        self.asym_cell = verts.dot(self.lat_vec)

    def _interp_pattern(self):
        pattern = self.pattern0.copy()
        lat1 = np.linalg.inv(self.lat_vec)
        pattern[:,:2] = pattern[:,:2].dot(lat1)                         #fractional coordinates
        pattern1 = cg.generate_symmetries(self.pp_type,pattern=pattern) #generate
        x1,y1 = pattern1[:,:2].T
        z1 = pattern1[:,2:]
        # dsp.stddisp(scat=[x1,y1,z1[:,1]],pOpt='Xet');plt.show()

        print(colors.blue+'...nearest neighbor ...'+colors.black)
        ndeg = 100
        x0 = np.linspace(0,1,ndeg)
        x,y = np.meshgrid(x0,x0)
        xy  = np.array([x.flatten(),y.flatten()]).T
        idx = np.array([np.linalg.norm(X-pattern1[:,:2],axis=1).argmin() for X in xy])
        idx = np.reshape(idx,(ndeg,ndeg))
        # print(idx[0],idx[-1],idx.shape,idx.max())

        print(colors.blue+'...2D interpolation ...'+colors.black)
        # print(z1.shape,idx.max())
        # data = np.reshape(z1[idx,0],x.shape)
        # print(data.shape,x.shape)
        self.f = [RegularGridInterpolator((x0,x0), np.reshape(z[idx],x.shape),method='nearest') for z in z1.T]

        # self.f = [interp2d(x0,x0,z[idx]) for z in z1.T]
        # self.f = [interp2d(x1,y1,z) for z in z1.T]


    def show_ref(self,ndeg=100):
        print(colors.blue+'...plotting ref ...'+colors.black)
        x = np.linspace(0,1,ndeg)
        im = np.stack([f(x,x) for f in self.f],axis=2)
        im = np.array(im*255,dtype=int)#;print(im.max())
        # print(im.shape,x.shape)
        plt.imshow(im);plt.show()
        # dsp.stddisp(im=[im])

    def show_xygrid(self,x,y,**kwargs):
        # x,y=x-np.floor(x),y-np.floor(y)
        x0,y0 = x.copy(),y.copy()
        x,y = np.meshgrid(x0,y0)
        xy = x.flatten(),y.flatten()
        xy = np.array([x-np.floor(x),y-np.floor(y)]) .T
        im = np.stack([np.reshape(f(xy),x.shape) for f in self.f],axis=2)

        # im = np.reshape(im,x.shape)
        im = np.array(im*255,dtype=int)#;print(im.max())
        plt.imshow(im);plt.show()

        # dsp.stddisp(im=[x0,y0,im],**kwargs)#;plt.show()


    def plot(self,opts='aug',nh=3,nk=3,pOpt='Xe',**kwargs):
        '''plot wallpaper
        - opts : a(asym_cell),u(unit_cell),g(lattice grid)
        - nh,nk : multiplicity
        '''
        fig,ax = dsp.stddisp(opt='')
        hk,R     = self.lat.get_lattice(nh,nk)
        R        = np.hstack([R,np.zeros((R.shape[0],self.pattern.shape[1]-2))])
        pattern1 = np.vstack([self.pattern +r for r in R])
        x,y = pattern1[:,:2].T
        z   = pattern1[:,2:]
        dsp.stddisp(ax=ax,scat=[x,y,z],opt='')

        pp = []
        if 'u' in opts:
            pp += [patches.Polygon(self.unit_cell,alpha=0.25,color='b')]
        if 'a' in opts:
            pp += [patches.Polygon(self.asym_cell,alpha=0.25,color='r')]
        if 'g' in opts:
            self.lat.plot_lattice(nh,nk,ax=ax,gridcolor='k',opt='')

        dsp.stddisp(ax=ax,patches=pp,pOpt=pOpt,**kwargs)

    ######################################################################
    ##### misc
    ######################################################################
    def save(self,file=None):
        file = os.path.join(self.path,self.name+'_wallpp.pkl')
        with open(file,'wb') as f:pickle5.dump(self,f,pickle5.HIGHEST_PROTOCOL)
        print(colors.green+'file saved : ' + colors.yellow+ file+colors.black)

######################################################################
##### misc
######################################################################
def load(path='',name='',file=None):
    if name:file = os.path.join(path,name+'_wallpp.pkl')
    with open(file,'rb') as f : obj = pickle5.load(f)
    return obj
