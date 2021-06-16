import importlib as imp
import easygui,os,glob,copy,pickle
from subprocess import Popen,PIPE
import numpy as np,matplotlib.pyplot as plt,pandas as pd
from utils import displayStandards as dsp
from utils import glob_colors as colors
# from multislice.pets import Pets
# from multislice.multislice import Multislice
# print('...bloch,multi...')
from blochwave.bloch import load_Bloch #Bloch
from multislice import mupy_utils as mut    #;imp.reload(mut)
from multislice import postprocess as pp    #;imp.reload(pp)
from multislice import multislice as ms     #;imp.reload(ms)
from blochwave import bloch as bl           ;imp.reload(bl)
from multislice import pets as pt           #;imp.reload(pt)
from . import __version__
from . import display as EDdisp             #;imp.reload(EDdisp)

pd.set_option('precision',3)
pd.set_option('max_rows',100)
if plt.rcParams['keymap.save']:
    plt.rcParams['keymap.save']=[]
    plt.rcParams['keymap.grid']=[]
    plt.rcParams['keymap.xscale']=[]
    plt.rcParams['keymap.yscale']=[]
    # plt.rc('text', usetex=False)



class Viewer:
    def __init__(self,config=None,help=False,**kwargs):
        '''Dynamic beam viewer :
        - config : str - saved configuration
        - help : bool - show help at startup
        - kwargs : see init_args'''
        self.__version__=__version__
        # print('...setting up keys...')
        self.set_keys()
        if help:self.show_help()

        self.path,self.cif_file,self.bloch,self.pets = None,None,None,None
        self.multi_path,self.tifpath,self.multi = None,None,None
        self.rpl,self.im,self.Icols = None,None,[]
        self.nfigs,self.kargs,self.u = 0,{},[0,0,1]
        self.pattern_args = {'Iopt':'csngt','Imax':5e4,'gs':0.025,'rmax':25,'Nmax':512}
        self.multi_args = {'data':'','mulslice':False,'NxNy':512,
            'slice_thick':1,'i_slice':10,'ssh':'','opt':'srf'}
        self.xyz_params = {'theta':0,'lat_params':[20,20,100],'pad':1.0,'rep':None}

        # print('...load parameters...')
        if config:self.load_config(config)
        else:self.init_args(**kwargs)
        for k,v in self.dsp_chars.items():self.dsp_d[k] = v in self.pets_opts
        for k,v in self.opt_chars.items():self.opt_d[k] = v in self.pets_opts
        self.cutoff_r,self.mag_r,self.dthick_r=self.cutoff,self.mag,self.dthick

        self.load_objects()
        self.set_mode(self.mode)
        if not (self.tifpath or self.multi):self.set_mode('rotate')

        # print('...Complete initialization...')
        self.multi_args['tag']  = self.tag+self.frame_str()
        self.multi_args['name'] = self.multi_path

        #print('... init figure ...')
        self.fig,self.ax = dsp.stddisp()
        cid = self.fig.canvas.mpl_connect('key_press_event', self)
        self.call = self.call_update

        # print('... update and show ...')
        self.xyz_params['file'] = self.cif_file
        self.xyz_params['n']    = self.u
        self.xyz_params['xyz']  = self.get_xyz()
        if not self.thick : self.thick=100
        self.update_thickness()
        self.update(fsolve=1)
        self.show()
        # plt.show()


    def init_args(self,path=None,cif_file=None,bloch=None,pets=None,multi=None,tag='',
        orient=[0,0,5,5],u=None,Smax=0.01,Nmax=5,thick=None,dthick=5,thicks=(0,1000,1000),
        F='Sw',xylims=1.5,init_opts='R',cutoff=100,dcutoff=10,mag=500,cmap='Greens',
        frame=1,incrate=1,drot=1,drotS=1,pets_opts='Pr',rot=0,rotS=0,kargs={},**kwargs):
        ''' configuration Initializer
    #### path setup (at least one of these need to be defined)
        - path      : str - path to the main folder
        - cif_file  : str - full path to .cif file (will be automatically fetched in path if not defined)
        - bloch     : Bloch object or str (will be automatically created if not defined)
        - multi     : str - path to multislice simulations (default multi)
        - pets      : Pets object or str to .pts file -
    #### Simulation setup
        - orient : [theta,phi,dtheta,dphi] in degrees
        - u      : list3 - set beam direction(takes preference over orient)
        - Smax,Nmax,thick,F : see blochwave.Bloch
        - init_opts: str -
            - h(show help at start) R(rotate mode) B(save Bloch)
            - k(show_hkl) i(show_i) z(show_z) u(show_u) v(show u in zone axis reciprocal space)
        - cutoff : brightness cutoff
        - mag    : circle magnification
        - kwargs : see Bloch.show_beams
        - kargs  : see Pets.show_frame
        '''
        #generic
        self.cutoff  = cutoff
        self.dcutoff = dcutoff
        self.mag     = mag
        self.xylims  = xylims
        self.Smax    = Smax
        self.Nmax    = Nmax
        self.thick   = thick
        self.dthick  = dthick
        #bloch related
        self.theta,self.phi,self.dtheta,self.dphi = orient
        self.dtheta_dphi = [self.dtheta,self.dphi]
        self.thicks      = thicks
        self.beams_args  = kwargs
        self.set_theta_phi_from_u(u)
        self.F,self.fopts = F,self.df_F.loc[self.df_F.F==F].iloc[0].fopt

        #frames related
        self.frame   = frame
        self.incrate = incrate
        self.rot     = rot
        self.rotS    = rotS
        self.drot,self.drotS = drot,drotS
        #pets related
        self.pets_opts = pets_opts

        #multislice related
        self.tag  = tag
        self.cmap = cmap

        #init_opts
        self.mode       = ['frames','rotate']['R' in init_opts]
        self.show_hkl   = 'k' in init_opts
        self.show_i     = 'i' in init_opts
        self.show_Ig    = 'I' in init_opts
        self.show_Vg    = 'V' in init_opts
        self.show_z     = 'z' in init_opts
        self.show_u     = 'u' in init_opts
        self.show_v     = 'v' in init_opts
        self.save_bloch = 'B' in init_opts
        self.hold_mode  = ' ' in init_opts

        # print('...init path...')
        self.init_path(path,cif_file,bloch,pets,multi)

    def load_config(self,config):
        with open(config,'rb') as f:d_config = pickle.load(f)
        for k,v in d_config.items(): self.__dict__[k] = v

    def init_path(self,path,cif_file,bloch,pets,multi):
        if all([type(o)==type(None) for o in [path,cif_file,bloch,pets] ]):
            args = ['path','cif_file','bloch','pets']
            raise Exception('at least one of these must be defined : ',args)

        if type(path)==str:self.path = path
        if type(cif_file)==str:self.cif_file = cif_file
        if pets:self.init_pets(pets)
        if bloch:self.init_bloch(bloch)
        if multi:self.init_multi(multi)

        if not type(self.path)==str:
            self.path = os.path.dirname(self.cif_file)
        if not type(self.cif_file)==str:
            self.cif_file = mut._get_cif_file(self.path,cif_file)

        if bloch and pets :
            bloch_cif,pets_cif = os.path.basename(self.bloch.cif_file),os.path.basename(self.pets.cif_file)
            if not bloch_cif==pets_cif:
                raise Exception('bloch and pets seems to have different cif_files : bloch=%s, pets=%s'%(bloch_cif,pets_cif))

    def load_objects(self):
        if not type(self.bloch)==bl.Bloch:
            bloch_path = os.path.join(self.path,'bloch')
            if not os.path.exists(bloch_path):
                p = Popen("mkdir %s" %bloch_path,shell=True,stderr=PIPE,stdout=PIPE);p.wait()
            # print('creating Bloch')
            self.bloch = bl.Bloch(self.cif_file,path=bloch_path)#,solve=1)

        if not type(self.pets)==pt.Pets:
            pets_path = os.path.join(self.path,'pets')
            if not os.path.exists(pets_path):pets_path=self.path
            pts_files = glob.glob(os.path.join(pets_path,'*.pts'))
            if len(pts_files)>=1:
                print(colors.green+'found .pts file. Loading %s ' %pts_files[0]+colors.black)
                self.init_pets(pts_files[0])
            else:
                print(colors.red+'no pts files found'+colors.black)
        if not type(self.multi_path)==str:
            multi_path = os.path.join(self.path,'multi')
            if os.path.exists(multi_path):
                self.multi_path = multi_path
                msg ='''multislice folder found '''
                print(colors.green+msg+colors.black)
            else:
                print(colors.red+'no multislice folder'+colors.black)

        if self.multi_path:
            self.update_nfigs()
            self.multi = pp.load(self.multi_path,tag=self.tag+str(self.frame).zfill(3))
            if not self.thick and self.multi:
                self.thick = self.multi.thickness
                # self.update_frame()
        # else:
        #     self.update(fsolve=1)

        self.figpath = os.path.join(self.path,'figures')
        if not os.path.exists(self.figpath):
            p = Popen("mkdir %s" %self.figpath,shell=True,stderr=PIPE,stdout=PIPE);p.wait()

    def init_bloch(self,bloch):
        if type(bloch)   == bl.Bloch:self.bloch = bloch
        elif type(bloch) == str:self.bloch = load_Bloch(file=bloch)
        else:raise Exception('bloch args need be NoneType, str or Bloch object ')
        if not type(self.path)     == str : self.path = self.bloch.path
        if not type(self.cif_file) == str : self.cif_file = self.bloch.cif_file

    def init_multi(self,multi_path):
        if type(multi_path) == str:self.multi_path = multi_path

    def init_pets(self,pets):
        if type(pets)==pt.Pets : self.pets = pets
        elif type(pets)==str : self.pets = pt.Pets(pets,gen=1)
        else:raise Exception('pets args need be NoneType, str or Pets object ')
        if not type(self.path)     == str : self.path = self.pets.path
        if not type(self.cif_file) == str : self.cif_file = self.pets.cif_file
        self.nfigs    = self.pets.nFrames
        self.tifpath  = os.path.join(self.pets.path,'tiff')
        # self.fmt,self.load = mut.find_format(self.tifpath)
        # self.figs = np.sort(glob.glob(self.tifpath+'/*.%s' %self.fmt))#;print(self.figs)

    def open_config(self,config_file=''):
        ext = config_file.split('.')[-1]
        if not ext=='pkl':config_file=''
        elif not os.path.exists(config_file):
            config_file = os.path.join(self.path,config_file)

        if not os.path.exists(config_file):
            config_file = os.path.join(self.path,'')
            config_file = easygui.fileopenbox('open config','config',config_file,['*.pkl'])

        ext = config_file.split('.')[-1]
        if os.path.exists(config_file) and ext=='pkl' :
            print(colors.green+'Loading '+colors.yellow+config_file+colors.black)
            Viewer(config=config_file)

    def save_config(self,config_file=''):
        keys = ['path','cif_file','multi_path','tifpath','mode','tag','u',      #paths
        'cutoff','dcutoff','mag','xylims','Smax','Nmax','thick','dthick',       #generic
        'theta','phi','dtheta','dphi','dtheta_dphi','thicks',                   #bloch
        'frame','incrate','rot','rotS','drot','drotS','pets_opts','cmap',       #frames
        'show_hkl','show_i','show_Ig','show_Vg','show_z','show_u','show_v',     #shows
        'save_bloch','hold_mode','beams_args','F','fopts',
        'multi_args','xyz_params',
        ]
        if not config_file:
            config_file = os.path.join(self.path,'config.pkl')
            config_file = easygui.filesavebox('save config','save',config_file)

        if config_file:
            d_config = {k:self.__dict__[k] for k in keys}
            with open(config_file,'wb') as f:pickle.dump(d_config,f,pickle.HIGHEST_PROTOCOL)
            print(colors.yellow+config_file+colors.green+' saved'+colors.black)

    def update_nfigs(self):
        if not self.tifpath and self.multi_path:
            self.nfigs = len(glob.glob(os.path.join(self.multi_path,'*.pkl')))
    ######################################################################
    ##### Update and show functions
    ######################################################################
    def __call__(self,event):
        ''' event callback function
        normal mode : call = call_update
        find mode   : call = shortcut_call
        '''
        self.call(event)

    def set_mode(self,mode):
        self.mode = mode
        self.show_im = {'frames':self.show_frames,'rotate':self.show_bloch}[self.mode]
        if self.mode=='frames':
            if self.multi_path and 'M' in self.pets_opts:
                self.update_frame()
                self.dthick_r,self.dthick=self.dthick,self.multi.dzs
                self.cutoff,self.mag = self.cutoff_r,self.mag_r
                print('frames mode. Changing dthick=%d' %self.dthick)
        else:
            self.dthick=self.dthick_r
            self.cutoff_r,self.mag_r = self.cutoff,self.mag
            self.cutoff,self.mag = 5,500
            self.set_theta_phi_from_u(self.u)
            self.update()

    def call_update(self, event):
        # print(event.key)
        key = event.key
        alias_match = self.df_keys.loc[(self.df_keys.alias==key) & (self.df_keys['mode']==self.mode)]
        if any(alias_match.index):
            # print(alias_match)
            key = alias_match.iloc[0].key
        # print(key)
        self.get_keys(key)
        self.do_update(key)

    def set_theta_phi_from_u(self,u=None):
        if type(u) in [list,np.ndarray]:
            self.u = u/np.linalg.norm(u)
            x,y,z  = self.u
            self.theta = np.rad2deg(np.arccos(z))
            self.phi   = np.rad2deg(np.arctan2(y,x))
            print('u set to : ', str(self.u))
            print('theta=%.2f, phi=%.2f' %(self.theta,self.phi))

    def update_u(self):
        if self.mode == 'rotate':
            theta,phi = np.deg2rad([self.theta,self.phi])
            ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
            self.u = [st*cp,st*sp,ct]
        else:
            if self.tifpath:
                self.u = self.pets.get_beam(self.frame)
            elif self.multi_path:
                if 'u' in self.multi.__dict__.keys():
                    self.u = self.multi.u
                else:
                    print(colors.red+'warning : multislice does not contain u switching to rotate mode'+colors.black)
                    self.set_mode('rotate')
                    self.u=[0,0,1]
        self.xyz_params['n']=self.u
        self.xyz_params['xyz']=self.get_xyz()

    def update_frame(self):
        if 'M' in self.pets_opts and self.multi_path:
            # if self.tifpath:
            tag = self.tag+str(self.frame).zfill(3)
            multi = pp.load(self.multi_path,tag=tag)
            # else:
                # pkls = glob.glob(os.path.join(self.multi_path,'*.pkl'))
                # multi = pp.load_multi_obj(pkls[self.frame])
            if multi:
                state = self.multi.check_simu_state()
                if state=='done':
                    self.multi=multi
                    self.multi.set_thicks()
                    self.multi.datpath=os.path.join(self.multi_path,'')
                    self.dthick = self.multi.dzs
                    self.update_thickness()
                else:
                    print(colors.blue+'multislice not done yet : '+colors.green+state+colors.black)
                    self.set_mode('rotate')
                    if not self.tifpath:return
            else:
                print(colors.red+'warning : file multislice not found switching to rotate mode'+colors.black)
                self.set_mode('rotate')
                if not self.tifpath:return

        if 'P' in self.pets_opts and self.tifpath:
            self.rpl = self.pets.rpl.loc[self.pets.rpl.F==self.frame]
        # if 'B' in self.pets_opts :
        self.update(fsolve='B' in self.pets_opts)

    def update_thickness(self):
        if 'M' in self.pets_opts and self.multi_path and self.mode=='frames':
            if self.thick<self.multi.thickness :
                self.iz,self.thick = self.multi.get_iz(self.thick,v=2)
                self.iz = int(self.iz)
            else:
                self.iz=None
                self.thick=self.multi.thickness
            qx,qy,I = self.multi.pattern(out=1,iz=self.iz,rot=self.rotS,v=1,**self.pattern_args)
            self.im=[qx,qy,I]
        if self.F in ['S','I','Ig']:self.bloch.set_thickness(self.thick)

    def update(self,fsolve=0):
        self.update_u()
        self.bloch.update_Nmax(self.Nmax)
        self.bloch.set_beam(u=self.u)                       #;print(self.u)
        self.bloch._set_excitation_errors(Smax=self.Smax)   #;print(self.Smax)
        self.bloch._set_Vg()
        self.bloch._set_kinematic()
        if fsolve:
            self.bloch._solve_Bloch(opts='0v')
            self.bloch.set_thickness(thick=self.thick)
            if self.save_bloch:
                if self.mode=='frames':self.bloch.save(self.path+'/%s_bloch.pkl' %str(self.frame).zfill(4))
                else:self.bloch.save()
                self.save_bloch=0

    def show_vals(self):
        if self.show_hkl : print(self.bloch.get_kin())
        if self.show_u   : print('[uvw]:');print(self.bloch.u)
        if self.show_v   : print('[uvw]_rec:');print(self.bloch.Kuvw0)
        if self.show_z   : print('zones:');print(self.bloch.get_zones())

        self.Icols = []
        if self.show_i  : self.Icols += ['I']
        if self.show_Ig : self.Icols += ['Ig']
        if self.show_Vg : self.Icols += ['Vg']
        if any(self.Icols) : self.bloch.get_Istrong(Icols=self.Icols)

    def show_frames(self):
        title = self.tag + ', frame %d' %self.frame
        if any([c in self.pets_opts for c in 'MK']):
            title += ', thickness=%.1f$\AA$' %self.thick

        hkl_idx = self.bloch.get_Istrong(Icols=self.Icols,out=1)
        EDdisp.show_frame(opts=self.pets_opts,mag=self.mag,rot=self.rot,
            df_pets=self.rpl,im_multi=self.im,df_bloch=self.bloch.df_G,hkl_idx=hkl_idx,
            ax=self.ax,title=title,xylims=self.xylims,
            cmap='Reds',caxis=[0,self.cutoff])

    def show_bloch(self):
        self.bloch.show_beams(ax=self.ax,fig=self.fig,F=self.F,fopts=self.fopts,
            mag=self.mag,cutoff=self.cutoff,opts='x',
            cmap=self.cmap,xylims=self.xylims,opt='',gridOn=0)#,**self.beams_args)

    def show(self):
        self.show_vals()
        self.ax.cla()
        self.show_im()
        self.fig.canvas.draw()

    ######################################################################
    ##### Multislice stuff
    ######################################################################
    def solve_multi(self):
        btns = ['Cancel','Change params','Gen xyz','Show xyz','Run']
        btn,tag = 'Start',self.tag
        # while (btn in btns[1:-2] or btn=='Start'):
        self.xyz_params['xyz']  = self.get_xyz()
        self.multi_args['data'] = self.get_xyz()
        self.multi_args['tag']  = self.tag+self.frame_str()
        msg = self.get_multi_params()
        btn = easygui.buttonbox(msg,'Multislice',btns)
        if   btn=='Change params' : self.set_multi()
        elif btn=='Gen xyz'       : self._gen_xyz(force=1);self.solve_multi()
        elif btn=='Show xyz'      : self.show_xyz()
        elif btn=='Run'           : self._run_multi()

        if not self.multi.tail=='_'+self.tag+self.frame_str():
            # print(self.multi.tail,self.tag+self.frame_str())
            self.update_frame()
            self.show()

    def set_multi(self):
        self.multi_path = os.path.join(self.path,'multi')
        if not os.path.exists(self.multi_path):
            p = Popen("mkdir %s" %self.multi_path,shell=True,stderr=PIPE,stdout=PIPE);p.wait()

        if not self.multi_args['ssh']:self.multi_args['ssh']=None
        if not self.tag:self.tag='None'
        mult_keys = ['NxNy','mulslice','slice_thick','i_slice','ssh','opt']
        xyz_keys  = ['rep','lat_params','pad']
        fieldNames  = ['tag'] + mult_keys+xyz_keys
        fieldValues = [self.tag]
        fieldValues += [str(self.multi_args[k]) for k in mult_keys]
        fieldValues += [str(self.xyz_params[k]) for k in xyz_keys]
        vals = EDdisp.multenterbox('solve multislice','multi',fieldValues,fieldNames)
        if vals:
            self.xyz_params['lat_params'] = eval(vals['lat_params'])
            self.xyz_params['pad']        = eval(vals['pad'])
            self.xyz_params['rep']        = eval(vals['rep'])
            self.tag = vals['tag']
            if self.tag=='None':self.tag=''
            # self.frame = eval(vals['frame'])
            for f in mult_keys:
                if isinstance(self.multi_args[f],str):
                    self.multi_args[f] = vals[f]
                else:
                    self.multi_args[f] = eval(vals[f])
        self.solve_multi()

    def show_multi(self):
        self.show_xyz()
        msg = self.get_multi_params()
        cli = '''Multislice. Use :
        ctrl+M to solve
        ctrl+m to change params
        '''
        print(colors.green+cli+colors.black)
        print(msg)

    def show_xyz(self):
        self._gen_xyz();
        mut.show_grid(self.get_xyz(),opts=['xy','xz'],popts='hv',figs='f',title=os.path.basename(self.get_xyz()))
    def get_multi_params(self):
        df_multi = pd.DataFrame.from_dict(self.multi_args,orient='index',columns=['values'])
        df_xyz   = pd.DataFrame.from_dict(self.xyz_params,orient='index',columns=['values'])
        msg ='''
    xyz params :
        %s

    Multislice params :
        %s
        '''%(str(df_xyz),str(df_multi))
        return msg

    def _gen_xyz(self,force=0):
        xyz = self.get_xyz()
        if not os.path.exists(xyz) or force:
            print('generating xyz')
            # print(self.cif_file,xyz,self.u,
            #     self.lat_params,self.pad,self.rotS,'v')
            xyz_params = self.xyz_params.copy()
            if self.xyz_params['rep']:
                xyz_params.pop('lat_params');
                mut.gen_xyz(**xyz_params)
            else:
                xyz_params.pop('rep');
                mut.gen_xyz2(opts='v',**xyz_params)

    def _run_multi(self):
        self._gen_xyz()
        # print('Solving multi : ', self.multi_args)
        self.multi = ms.Multislice(u=self.u,**self.multi_args)
        self.update_nfigs()

    def show_tags(self):
        if self.multi_path:
            pkls = glob.glob(os.path.join(self.multi_path,'*.pkl'))
            tails = ['_'.join(pkl.split('_')[1:-1]) for pkl in pkls]
            tails = np.unique([pkl[:-3] for pkl in pkls])
            print(colors.green+'available tails:'+colors.yellow,tails,colors.black)

    def get_xyz(self):
        return os.path.join(self.multi_path,self.tag+self.frame_str()+'.xyz')
    def frame_str(self,pad=3):return '%s' %str(self.frame).zfill(pad)

    ######################################################################
    ##### Misc
    ######################################################################
    def get_figname(self):
        if self.mode=='frame':
            png_file = '_%s.png' %str(self.frame).zfill(3)
        else:
            png_file = '%s.png' %str(self.bloch.name)
        return os.path.join(self.figpath,png_file)

    ######################################################################
    ##### Keys
    ######################################################################
    def do_update(self,key):
        if   key == self.dict_k['save_bloch']:self.update(fsolve=1)
        elif key in self.modes_keys : print('mode %s' %self.mode)
        elif key in self.frame_keys   and self.mode=='frames': print('frame %d' %self.frame)
        elif key in self.dframe_keys  and self.mode=='frames': print('increment rate : %d' %self.incrate)
        elif key in self.orient_keys  and self.mode=='rotate': print('theta=%.4f, phi=%.4f' %(self.theta,self.phi))
        elif key in self.dorient_keys and self.mode=='rotate': print('dtheta=%.1f, dphi=%.1f' %(self.dtheta,self.dphi))
        elif key in self.rot_keys     and self.mode=='frames': print('rot=%.1f, rotS=%.1f' %(self.rot,self.rotS))
        elif key in self.thick_keys      : print('thickness : %d' %self.thick);self.update_thickness()
        elif key in self.dthick_keys     : print('dthick=%d' %self.dthick)
        elif key in self.brightness_keys : print('brightness cutoff : %.1f' %self.cutoff)
        elif key in self.dbrightness_keys: print('dbrightness cutoff : %.1f' %self.dcutoff)
        elif key in self.rbrightness_keys: print('cutoff : %.1f, dcutoff : %.1f' %(self.cutoff,self.dcutoff))
        elif key in self.vals_keys       : self.show_vals()
        #
        if (key in self.orient_keys or key in ['enter']) and self.mode=='rotate' :
            self.update(fsolve=self.F in ['S','I'])
        elif key in self.F_keys or key in self.Fnum_keys and self.mode=='rotate':
            print('displaying F="%s"' %self.F)
            if not self.bloch.solved:
                if self.F in ['S','I'] : self.update(fsolve=1)
        elif (key in self.frame_keys or key in ['enter'] ) and self.mode=='frames':
            self.update_frame()
        elif key in self.pets_keys and self.mode=='frames':
            print('frame options %s' %self.pets_opts)

        show_keys = key in self.show_keys
        show_keys = show_keys or (key in self.frame_keys+self.dsp_keys+self.opt_keys) and self.mode=='frames'
        show_keys = show_keys or (key in self.orient_keys+self.F_keys+self.Fnum_keys) and self.mode=='rotate'
        if show_keys : self.show()


    def get_keys(self,key):
        # print(key)
        if   key == self.dict_k['settings']     : self.settings()
        elif key == self.dict_k['find_shortcut']: self.find_shortcut()
        elif key == self.dict_k['new_fig'] :
            config_file = os.path.join(self.path,'new_fig.pkl')
            self.save_config(config_file)
            Viewer(config=config_file)
            # v1.fig,v1.ax = dsp.stddisp()
            # cid = v1.fig.canvas.mpl_connect('key_press_event', v1)
        elif key == self.dict_k['help_cli']      : self.show_help_cli()
        elif key == self.dict_k['help_gui']      : self.show_help_gui()
        elif key == self.dict_k['help_simple']   : print(self.df_keys)
        elif key == self.dict_k['save_fig']      : dsp.saveFig(self.get_figname(),ax=self.ax)
        elif key == self.dict_k['save_config']   : self.save_config(os.path.join(self.path,'config.pkl'))
        elif key == self.dict_k['save_config_as']: self.save_config()
        elif key == self.dict_k['reload_config'] : self.open_config('config.pkl')
        elif key == self.dict_k['open_config']   : self.open_config()
        elif key == self.dict_k['save_bloch']    : self.save_bloch=1 #not self.save_bloch;print('%s saving bloch' %['not ',''][self.save_bloch])
        elif key == self.dict_k['solve_bloch']   : self.update(fsolve=1)
        elif key == self.dict_k['set_multi']     : self.set_multi()
        elif key == self.dict_k['solve_multi']   : self.solve_multi()
        elif key == self.dict_k['show_multi']    : self.show_multi()
        elif key == self.dict_k['show_tags']     : self.show_tags()
        elif key == self.dict_k['hold']          : self.hold_mode=not self.hold_mode;print('hold mode ' + ['off','on'][self.hold_mode])
        #modes
        elif key==self.dict_k['frames_mode'] : self.set_mode('frames')
        elif key==self.dict_k['rotate_mode'] : self.set_mode('rotate')

        if self.mode=='rotate':
            # theta,phi
            if   key == self.dict_k['elev_up']   : self.theta=(self.theta+self.dtheta)%180
            elif key == self.dict_k['elev_down'] : self.theta=(self.theta-self.dtheta)%180
            elif key == self.dict_k['azim_up']   : self.phi=(self.phi+self.dphi)%360
            elif key == self.dict_k['azim_down'] : self.phi=(self.phi-self.dphi)%360
            # dtheta,dphi,dthick
            elif key == self.dict_k['delev_up']   : self.dtheta=min(self.dtheta+1,90)
            elif key == self.dict_k['delev_down'] : self.dtheta=max(self.dtheta-1,0)
            elif key == self.dict_k['dazim_up']   : self.dphi=min(self.dphi+1,90)
            elif key == self.dict_k['dazim_down'] : self.dphi=max(self.dphi-1,0)
        elif self.mode=='frames':
            # frame
            if   key==self.dict_k['frame_up']   : self.frame=min(self.frame+self.incrate,self.nfigs)
            elif key==self.dict_k['frame_down'] : self.frame=max(1,self.frame-self.incrate)
            # increment rate
            elif key==self.dict_k['inc_frame_up']  :self.incrate=min(self.incrate+1,100)
            elif key==self.dict_k['inc_frame_down']:self.incrate=max(1,self.incrate-1)
            #in plane rotate
            if   key == self.dict_k['rot_up']   : self.rot=(self.rot+self.drot)%360
            elif key == self.dict_k['rot_down'] : self.rot=(self.rot-self.drot)%360
            elif key == self.dict_k['rotS_up']  : self.rotS=(self.rotS+self.drotS)%360
            elif key == self.dict_k['rotS_down']: self.rotS=(self.rotS-self.drotS)%360

        #brightness
        if   key==self.dict_k['brightness up']     : self.cutoff=min(self.cutoff+self.dcutoff,1e5)
        elif key==self.dict_k['brightness down']   : self.cutoff=max(1e-5,self.cutoff-self.dcutoff)
        if   key==self.dict_k['dbrightness up']    : self.dcutoff*=2
        elif key==self.dict_k['dbrightness down']  : self.dcutoff/=2
        elif key==self.dict_k['brightness reset']  : self.cutoff  = 50
        elif key==self.dict_k['dbrightness reset'] : self.dcutoff = 10
        # thickness
        elif key == self.dict_k['thick_up']     : self.thick+=self.dthick
        elif key == self.dict_k['thick_down']   : self.thick=max(self.thick-self.dthick,self.dthick)
        elif key == self.dict_k['dthick_up']    : self.dthick=min(self.dthick+1,1000)
        elif key == self.dict_k['dthick_down']  : self.dthick=max(self.dthick-1,1)
        # vals
        elif key == self.dict_k['show_u']   : self.show_u   = not self.show_u  ;print('%sshowing [uvw]'     %['not ',''][self.show_u])
        elif key == self.dict_k['show_v']   : self.show_v   = not self.show_v  ;print('%sshowing [uvw]_rec' %['not ',''][self.show_v])
        elif key == self.dict_k['show_hkl'] : self.show_hkl = not self.show_hkl;print('%sshowing hkl'       %['not ',''][self.show_hkl])
        elif key == self.dict_k['show_I']   : self.show_i   = not self.show_i  ;print('%sshowing Istrong'   %['not ',''][self.show_i])
        elif key == self.dict_k['show_Ig']  : self.show_Ig  = not self.show_Ig ;print('%sshowing Ig_strong' %['not ',''][self.show_Ig])
        elif key == self.dict_k['show_Vg']  : self.show_Vg  = not self.show_Vg ;print('%sshowing Vg_strong' %['not ',''][self.show_Vg])
        elif key == self.dict_k['show_z']   : self.show_z   = not self.show_z  ;print('%sshowing zones'     %['not ',''][self.show_z])
        elif key == self.dict_k['show_beams_vs_thickness']:
            if not self.bloch.solved : self.update(fsolve=1)
            self.bloch.show_beams_vs_thickness(self.thicks,strong=self.Icols)
        #Rotate mode
        if self.mode=='rotate':
            if (key in self.F_keys or key in self.Fnum_keys):
                loc = self.df_F.loc[(self.df_F.key==key) | (self.df_F.alias==key)].iloc[0]#;print(loc)
                self.F,self.fopts = loc.F,loc.fopt
        # Frames mode
        elif self.mode=='frames':
            do_up=False
            if key in self.opt_keys :
                self.opt_d[key] = not self.opt_d[key]
            elif key in self.dsp_keys:
                if self.hold_mode:
                    self.dsp_d[key] = not self.dsp_d[key]
                else:
                    self.dsp_d = self.dsp_d.fromkeys(self.dsp_d,False)
                    self.dsp_d[key] = True
                do_up = do_up or (self.dsp_d['M'] and not 'M' in self.pets_opts and self.multi_path and not self.im)
                do_up = do_up or (self.dsp_d['P'] and not 'P' in self.pets_opts and self.tifpath )
                do_up = do_up or (self.dsp_d['B'] and not 'B' in self.pets_opts and not self.bloch.solved)
            elif key==self.dict_k['show_Exp'] and self.tifpath:
                self.pets.show_exp(frame=self.frame)
            self.pets_opts =''.join([self.opt_chars[c] for c,val in self.opt_d.items() if val])
            self.pets_opts+=''.join([self.dsp_chars[c] for c,val in self.dsp_d.items() if val])
            if do_up:self.update_frame()


    def set_keys(self):
        self.dict_k,self.dict_a = {},{}

        generic       = ['new_fig'   ,'settings' ,'help_cli','help_gui','help_simple' ]
        self.gen_keys = ['ctrl+enter','enter'    ,'ctrl+2'  ,'ctrl+1'  ,'ctrl+ '      ]
        generic      += ['save_config','save_config_as','open_config','reload_config',]
        self.gen_keys+= ['ctrl+s'     ,'ctrl+S'        ,'ctrl+o'     ,'backspace'    ,]
        generic      += ['save_fig','save_bloch','solve_bloch']
        self.gen_keys+= ['s'       ,'ctrl+alt+S','ctrl+B'     ]
        generic      += ['set_multi','solve_multi','show_multi','show_tags']
        self.gen_keys+= ['ctrl+m'   ,'ctrl+M'     ,'ctrl+alt+m','ctrl+O'    ]
        generic      += ['hold','find_shortcut']
        self.gen_keys+= [' '   ,'ctrl+P'       ]
        self.dict_k.update(dict(zip(generic,self.gen_keys)))

        #modes
        modes = ['frames_mode','rotate_mode']
        self.modes_keys = ['ctrl+alt+'+c for c in '12']
        self.dict_k.update(dict(zip(modes,self.modes_keys)))

        #orientation
        rots  = ['elev_up','elev_down','azim_up','azim_down']
        self.orient_keys  = ['up','down','right','left']
        self.dict_k.update(dict(zip(rots,self.orient_keys)))
        drots = ['d'+s for s in rots]
        self.dorient_keys = ['ctrl+'+s for s in self.orient_keys]
        self.dict_k.update(dict(zip(drots,self.dorient_keys)))

        #rotate in frame
        rots  = ['rot_up','rot_down','rotS_up','rotS_down']
        self.rot_keys  = ['r','R','ctrl+r','ctrl+R']
        self.dict_k.update(dict(zip(rots,self.rot_keys)))
        # drots = ['d'+s for s in rots]
        # self.dorient_keys = ['ctrl+'+s for s in self.orient_keys]
        # self.dict_k.update(dict(zip(drots,self.dorient_keys)))

        #frames
        frames = ['frame_'+s for s in ['up','down']] #'right','down','left']]
        dframes = ['inc_'+s for s in frames]
        self.frame_keys    = ['up','down']
        self.dframe_keys   = ['ctrl+'+s for s in self.frame_keys]
        self.frame_keys_a  = ['right','left']
        # self.dframe_keys_a = ['ctrl+'+s for s in self.frame_keys_a]
        self.dict_k.update(dict(zip(frames,self.frame_keys)))
        self.dict_k.update(dict(zip(dframes,self.dframe_keys)))
        self.dict_a.update(dict(zip(frames,self.frame_keys_a)))
        self.frame_keys = self.frame_keys.copy()+self.frame_keys_a

        #thickness
        thicks  = ['thick_up','thick_down']
        dthicks = ['d'+s for s in thicks]
        self.thick_keys  = ['shift+ctrl+up','shift+ctrl+down']
        self.dthick_keys = ['ctrl+t','ctrl+T']
        self.dict_k.update(dict(zip(thicks,self.thick_keys)))
        self.dict_k.update(dict(zip(dthicks,self.dthick_keys)))

        #brightness
        bright  = ['brightness '+s for s in ['down','up']]
        dbright = ['d'+s for s in bright]
        rbright = [s+'brightness reset' for s in ['','d']]
        self.brightness_keys  = ['pagedown','pageup']
        self.dbrightness_keys = ['ctrl+'+k for k in self.brightness_keys]
        self.rbrightness_keys = ['shift+ctrl+'+k for k in self.brightness_keys]
        self.dict_k.update(dict(zip(bright,self.brightness_keys)))
        self.dict_k.update(dict(zip(dbright,self.dbrightness_keys)))
        self.dict_k.update(dict(zip(rbright,self.rbrightness_keys)))

        #vals
        vals = ['hkl','I','z','u','v','beams_vs_thickness','Ig','Vg']
        vals = ['show_'+s for s in vals]
        self.vals_keys = ['ctrl+'+c for c in 'hizuvbIV']
        self.dict_k.update(dict(zip(vals,self.vals_keys)))

        #display Bloch
        Fs    = ['L','Sw','Vg','S','I','Ig']#,'E']
        fopt  = ['m','L' ,'l' ,'m','m','m' ]#,'m']
        bloch = ['Lattice','Excitation error','Potential','Scattering','Intensities','Kinematic']
        bloch = ['show '+s for s in bloch]
        self.F_keys    = [''+c for c in 'LGVSIK']
        self.Fnum_keys = [''+c for c in '123456']
        self.df_F = pd.DataFrame.from_dict(dict(zip(['names','key','alias','F','fopt'],[bloch,self.F_keys,self.Fnum_keys,Fs,fopt])))
        self.df_F = self.df_F.set_index('names')

        #display Frames
        self.dict_k['show_Exp']='0'
        petsD = ['Processed','Multislice','Bloch','Kinematic','Potential','Excitation error','Lattice']
        petsD = ['show_'+s for s in petsD]
        self.dsp_keys   = ['P','M','B','K','V','S','L']
        self.dsp_keys_a = [''+c for c in '1234567']
        self.dsp_chars = dict(zip(self.dsp_keys,'PMBKVSL'))
        self.dsp_d = dict(zip(self.dsp_keys,np.zeros((len(self.dsp_keys)),dtype=bool) ))
        self.dict_k.update(dict(zip(petsD,self.dsp_keys)))
        self.dict_a.update(dict(zip(petsD,self.dsp_keys_a)))
        #options
        petsO = ['hkl_exp','hkl_kin','rings','legend']
        petsO = ['show_'+s for s in petsO]
        self.opt_keys = ['j','k','g','l']
        self.opt_chars = dict(zip(self.opt_keys,'hkrl'))
        self.opt_d = dict(zip(self.opt_keys,np.zeros((len(self.opt_keys)),dtype=bool) ))
        self.dict_k.update(dict(zip(petsO,self.opt_keys)))
        self.pets_keys = self.dsp_keys+self.opt_keys

        #### Keys that trigger the show
        self.show_keys = self.modes_keys.copy()
        self.show_keys += self.rot_keys+self.thick_keys+self.brightness_keys+self.rbrightness_keys
        self.show_keys +=[self.dict_k[k] for k in ['settings','solve_bloch']]

        #### Data frame for the help command
        self.df_keys = pd.DataFrame.from_dict(self.dict_k,orient='index',columns=['key'])
        self.df_keys['alias'] = ['']*self.df_keys.shape[0]
        for k,v in self.dict_a.items():self.df_keys.loc[k,'alias'] = v
        self.df_keys = self.df_keys.append(self.df_F[['key','alias']])
        self.df_keys['mode']  = ['any']*self.df_keys.shape[0]
        self.df_keys.loc[rots+drots+bloch,'mode'] = 'rotate'
        self.df_keys.loc[frames+dframes+petsD+petsO,'mode'] = 'frames'

        self.df_generic = self.df_keys.loc[generic]
        self.df_modes   = self.df_keys.loc[modes]
        self.df_orient  = self.df_keys.loc[rots+drots]
        self.df_frames  = self.df_keys.loc[frames+dframes]
        self.df_thicks  = self.df_keys.loc[thicks+dthicks]
        self.df_bright  = self.df_keys.loc[bright]
        self.df_vals    = self.df_keys.loc[vals]
        self.df_bloch   = self.df_keys.loc[bloch]
        self.df_pets    = self.df_keys.loc[petsD+petsO]
        # self.df_pets['alias'].iloc[:4] = self.Pnum_keys

    def settings(self):
        fieldNames  =['Smax','Nmax','thick','dthick','xylims','mag','cutoff']
        if   self.mode=='rotate':
            self.dtheta_dphi=[self.dtheta,self.dphi]
            fieldNames+=['theta','phi','dtheta_dphi','thicks','F','cmap']
        elif self.mode=='frames':
            fieldNames+=['frame','incrate','rot','rotS','pets_opts']
        fieldValues = [str(self.__dict__[f]) for f in fieldNames]
        dict_fv = EDdisp.multenterbox("Change settings","settings", fieldValues,fieldNames)
        if dict_fv:
            for f in fieldNames:
                if isinstance(self.__dict__[f],str):
                    self.__dict__[f] = dict_fv[f]
                else:
                    self.__dict__[f] = eval(dict_fv[f])

        #update dtheta,dphi
        if self.mode=='rotate':
            if type(self.dtheta_dphi) in [float,int]:
                self.dtheta_dphi=[self.dtheta_dphi]*2
            self.dtheta,self.dphi = self.dtheta_dphi

    ###################################
    ##### Shotcut call
    ###################################
    def is_key(self,key):return self.df_keys.loc[(self.df_keys.key==key) | (self.df_keys.alias==key)]
    def is_cmd(self,cmd):return self.df_keys.iloc[[cmd in k for k in  self.df_keys.index]]
    def find_shortcut(self):
        print(colors.green+'find shortcut mode :'+colors.black)
        msg = '''
        Type a shortcut to check its effect.
        Use escape to switch back to normal mode'''
        print(msg)
        self.call = self.shortcut_call

    def shortcut_call(self,event):
        key = event.key
        combs = ['control','shift','alt',
         'ctrl+shift','ctrl+alt','ctrl+alt+shift',
         'shift+control','shift+super','shift+ctrl+super','shift+alt+control',
         'alt+control','alt+shift']
        if key=='escape' :
            self.call=self.call_update
            print(colors.green+'Back to normal mode'+colors.black)
        elif not key in combs :
            dfk = self.is_key(key)
            if any(dfk.index):print(colors.green+'found key'+colors.black);print(dfk)
            else:print(colors.red+'no key binding for %s' %key+colors.black)

    ###################################
    ##### Help
    ###################################
    def show_help_cli(self):
        print(colors.green+'Shortcuts :  '+colors.black)
        print(colors.blue+'Generic :     '+colors.black);print(self.df_generic)
        print(colors.blue+'Modes :       '+colors.black);print(self.df_modes)
        print(colors.blue+'Frames (Frames mode) :      '+colors.black);print(self.df_frames)
        print(colors.blue+'Orientation (Rotate mode) : '+colors.black);print(self.df_orient)
        print(colors.blue+'Thickness :   '+colors.black);print(self.df_thicks)
        print(colors.blue+'Brightness :  '+colors.black);print(self.df_bright)
        print(colors.blue+'Show values : '+colors.black);print(self.df_vals)
        print(colors.blue+'Display (Rotate mode): '+colors.black);print(self.df_bloch)
        print(colors.blue+'Display (Frames mode): '+colors.black);print(self.df_pets)

    def show_help_gui(self):
        # msg = str(self.df_keys)
        msg =''
        msg+= '\t\t\t SHORTCUTS : \n'
        msg+= 'Generic :     '+'\n'
        msg+=str(self.df_generic)+'\n'
        msg+= 'Modes :       '+'\n'
        msg+=str(self.df_modes)+'\n'
        msg+= 'Frames (Frames mode) :      '+'\n'
        msg+=str(self.df_frames)+'\n'
        msg+= 'Orientation (Rotate mode) : '+'\n'
        msg+=str(self.df_orient)+'\n'
        msg+= 'Thickness :   '+'\n'
        msg+=str(self.df_thicks)+'\n'
        msg+= 'Brightness :  '+'\n'
        msg+=str(self.df_bright)+'\n'
        msg+= 'Show values : '+'\n'
        msg+=str(self.df_vals)+'\n'
        msg+= 'Display (Rotate mode): '+'\n'
        msg+=str(self.df_bloch)+'\n'
        msg+= 'Display (Frames mode): '+'\n'
        msg+=str(self.df_pets)+'\n'
        easygui.msgbox(msg, title='help')
