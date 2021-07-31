from utils import displayStandards as dsp
from utils import glob_colors as colors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
import easygui,os
from . import config as cg
from . import lattice as lat

if plt.rcParams['keymap.save']:
    plt.rcParams['keymap.save']=[]

class Wallpaper_gui():
    def __init__(self,image,path='',pg='p1',npyfile=None,config=None,opt=0,
            params={'a':1},asym_verts=[],P1=[0,0]):
        self.path  = path
        self.name  = ''.join(os.path.basename(image).split('.')[:-1])#;print(self.name)
        self.image = image
        self.im0   = plt.imread(self.image)
        self._set_keys()
        self.pg  = pg
        if not self.load_config(config):
            self.params = params
            self.P1     = np.array(P1)
            self.rot    = 0
        self._set_plane_group()
        self.set_npy_file(npyfile)

        self.xylims=[]
        self.fig,self.ax = dsp.stddisp(title=self.pg)#opt='')
        nx,ny = self.im0.shape[:2]
        self.xylims=np.array([0,nx,0,ny])
        self.fig.canvas.mpl_connect('key_press_event',self)
        self.show()
        if opt:self.show_pattern(opt)

    ################################################################
    ### callbacks
    ################################################################
    def __call__(self,event):
        key = event.key         #;print(key)
        kopts=self.opts_keys
        if key in self.d_keys.keys():
            self.__getattribute__(self.d_keys[key])()
            self.show()

    def _set_keys(self):
        opts = ['vec','asym']
        self.opts=dict(zip(opts,[True]*len(opts)))
        self.opts_keys=dict(zip(['2','3'],opts))
        self.d_keys={
            'h':'help',
            'ctrl+enter'        :'set_vectors',
            'shift+enter'       :'set_origin',
            'shift+ctrl+enter'  :'set_plane_group',
            'ctrl+2'            :'show_pattern',
            'ctrl+3'            :'show_unit_cell',
            'ctrl+4'            :'show_wallpp',
            'ctrl+='            :'increase_lims',
            'ctrl+-'            :'decrease_lims',
            'ctrl+o'            :'open_image',
            'ctrl+O'            :'load_config',
            'ctrl+S'            :'save_config',
            'ctrl+s'            :'save',
            'S'                 :'save_pattern',
            's'                 :'save_fig'
        }
        self.df_keys = pd.DataFrame.from_dict({v:k for k,v in self.d_keys.items()},orient='index',columns=['key'])

    def increase_lims(self):self.xylims *= 2
    def decrease_lims(self):self.xylims = self.xylims/2

    def set_vectors(self):
        print('select points defining lattice parameters')
        ns = {'oblique':3,'rect':3,'hex':2,'square':2}
        lat_type=self.lat_type
        n = ns[lat_type]
        pts = np.array(plt.ginput(n=n,timeout=1000))
        if len(pts)==n:
            if lat_type=='square' or lat_type=='hex':
                a1 = pts[1]-pts[0]
                a = np.linalg.norm(a1)
                params = {'a':int(a)}
            else:
                P1,P2,P3 = pts
                a1,a2 = P2-P1,P3-P1
                a,b = np.linalg.norm(a1),np.linalg.norm(a2)
                params = {'a':int(a),'b':int(b)}
                if lat_type=='oblique':
                    alpha = np.arccos(np.dot(a1,a2)/(a*b))
                    params['alpha'] = np.rad2deg(alpha)
            rot = np.rad2deg(np.arctan2(a1[1],a1[0]))#*np.sign(a1[1])
            eps = rot-90*np.round(rot/90)           #;print(rot,eps)
            if abs(eps)<8:rot=90*np.round(rot/90)   #;print(rot)
            self.rot = -np.deg2rad(rot)
            self.params = params
            self._set_plane_group()

    def set_origin(self):
        print('select an origin ')
        pts = np.array(plt.ginput(n=1,timeout=10000),dtype=int)
        self.P1=pts[0]

    def set_plane_group(self):
        pg = easygui.choicebox('choose plane group','plane',cg.pp_types)
        if pg :
            self.pg=pg
            self._set_plane_group()
            self.set_npy_file()

    def _set_plane_group(self):
        self.lat_type = cg.df_wallpp.loc[self.pg,'lat_type']
        asym_cell_frac = cg.df_wallpp.loc[self.pg,'asym_cell']
        self.lat_vec = lat.set_lattice_vec(self.lat_type,**self.params)
        self.asym_verts = asym_cell_frac.dot(self.lat_vec)
        if self.rot:
            ct,st = np.cos(self.rot),np.sin(self.rot)
            rot2=np.array([[ct,st],[-st,ct]])
            self.lat_vec = rot2.dot(self.lat_vec.T).T
            self.asym_verts = rot2.dot(self.asym_verts.T).T

    def set_pattern(self):
        nx = min(self.im0.shape[:2])
        patch = plt.Polygon(self.asym_verts+self.P1)
        x,y = np.meshgrid(np.arange(nx),np.arange(nx))
        x0,y0 = x.flatten(),y.flatten()
        # print(patch.contains_point([320,90]))
        # im0 = np.sum(self.im0[:nx,:nx,:],axis=2)
        im0 =  self.im0[:nx,:nx,:]#/255#-self.im0[:nx,:nx,0]
        idx = patch.contains_points(list(np.array([x0,y0]).T),radius=0)#;print(idx)
        x0,y0,z0,z1,z2 = x0[idx],y0[idx],im0[:,:,0].flatten()[idx],im0[:,:,1].flatten()[idx],im0[:,:,2].flatten()[idx]
        # z = z1
        # u = np.linalg.norm(np.vstack([x0,y0]).T-np.array([320,90]),axis=1)
        # print(x0[u.argmin()],y0[u.argmin()])

        self.vec = self.lat_vec
        x0-=self.P1[0]
        y0-=self.P1[1]
        if self.rot:
            ct,st = np.cos(-self.rot),np.sin(-self.rot)
            rot2  = np.array([[ct,st],[-st,ct]])
            x0,y0 = rot2.dot(np.vstack([x0,y0]))
            self.vec = rot2.dot(self.vec.T).T
        self.pattern = np.array([x0,y0,z0,z1,z2]).T

    ################################################################
    ### show
    ################################################################
    def show(self):
        self.ax.cla()
        self.ax.imshow(self.im0,origin='lower')
        pp = []
        if self.opts['vec']:
            a1,a2 = self.lat_vec
            verts = self.P1+np.array([0*a1,a1,a1+a2,a2])
            pp += [patches.Polygon(verts,alpha=0.2,color='b')]
        if self.opts['asym']:
            verts = self.P1+self.asym_verts
            pp += [patches.Polygon(verts,alpha=0.2,color='r')]

        dsp.stddisp(ax=self.ax,patches=pp,pOpt='tX',xylims=self.xylims,opt='',title=self.pg)
        self.fig.canvas.draw()

    def show_asym(self)     :self.show_pattern(0)
    def show_unit_cell(self):self.show_pattern(1)
    def show_wallpp(self)   :self.show_pattern(2)
    def show_pattern(self,opt=0):
        npyfile=os.path.join(self.path,'tmp.npy')
        self.save(npyfile)
        pattern,vec = np.load(npyfile,allow_pickle=True)
        a1,a2 = vec
        verts = [0*a1,a1,a1+a2,a2]
        pp = [patches.Polygon(verts,alpha=0.5,color='b')]
        # pattern = np.array(pattern).T
        if opt:pattern = cg.generate(self.pg,vec,pattern)
        if opt>1:
            pattern=cg.replicate(2,2,self.pg,self.params,pattern)
        x,y = pattern[:,:2].T
        F = pattern[:,2:] #;print(F.shape)
        dsp.stddisp(scat=[x,y,F],sargs={'s':40,'marker':'s'},patches=pp,
            title=self.pg,equal=1)

    ################################################################
    ### help,save,load,open
    ################################################################
    def help(self):
        print(colors.green+'help:')
        print(colors.blue,self.df_keys)
        print(colors.green+'good luck'+colors.black)

    def save_config(self,config=None):
        if not config:config=os.path.join(self.path,self.name+'%s_config.npy' %self.pg)
        np.save(config,[self.pg,self.params,self.rot,self.P1])
        print(colors.green+'config saved : '+colors.yellow+config+colors.black)

    def load_config(self,config=None):
        if config:
            if not isinstance(config,str):config=os.path.join(self.path,self.name+'%s_config.npy' %self.pg)
            if os.path.exists(config):
                self.pg,self.params,self.rot,self.P1 = np.load(config,allow_pickle=True)
                self._set_plane_group()
                print(colors.green+'config loaded : '+colors.yellow+config+colors.black)
                return 1

    def open_image(self):
        image = easygui.fileopenbox('choose image','image',self.image)
        if image:
            self.image=image
            self.im0 = plt.imread(self.image)
            self.set_npy_file()

    def set_npy_file(self,npyfile=None):
        if not npyfile:
            self.npyfile = os.path.join(self.path,self.name+'_%s.npy' %self.pg)
        else:
            self.npyfile=npyfile

    def save_pattern(self):
        npyfile = easygui.filesavebox('save pattern','save',self.npyfile)
        self.set_npy_file(npyfile)
        if npyfile:self.save()

    def save(self,npyfile=''):
        if not npyfile:npyfile = self.npyfile
        self.set_pattern()
        np.save(npyfile,[self.pattern,self.vec])
        print(colors.green+'file saved : '+colors.yellow+npyfile+colors.black)
    def save_fig(self):
        figname = os.path.join(self.path,self.name+'.png')
        dsp.saveFig(figname,ax=self.ax)
