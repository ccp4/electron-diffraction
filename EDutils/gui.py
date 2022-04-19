import importlib as imp
import easygui,os,glob,copy,pickle
from subprocess import Popen,PIPE
import numpy as np,matplotlib.pyplot as plt,pandas as pd
from utils import displayStandards as dsp
from utils import glob_colors as colors
# from multislice.pets import Pets
# from multislice.multislice import Multislice
# print('...bloch,multi...')
from blochwave.util import load_bloch #Bloch
from multislice import mupy_utils as mut    #;imp.reload(mut)
from multislice import postprocess as pp    #;imp.reload(pp)
from multislice import multislice as ms     #;imp.reload(ms)
from blochwave import bloch as bl           ;imp.reload(bl)
from blochwave import util as bl_ut         ;imp.reload(bl_ut)
from . import __version__
from . import pets as pt                    ;imp.reload(pt)
from . import display as EDdisp             ;imp.reload(EDdisp)
from . import gui_config as cfg             ;imp.reload(cfg)
from . import utilities as ut               ;imp.reload(ut)

# currently bugs with pandas 1.4
#pd.set_option('precision',3)
#pd.set_option('max_rows',100)
if plt.rcParams['keymap.save']:
    plt.rcParams['keymap.save']=[]
    plt.rcParams['keymap.grid']=[]
    plt.rcParams['keymap.xscale']=[]
    plt.rcParams['keymap.yscale']=[]
    # plt.rc('text', usetex=False)



class Gui:
    def __init__(self,config=None,help=True,pargs={},**kwargs):
        '''Dynamic beam viewer :
            - config : str - saved configuration
            - help : bool - show help at startup
            - kwargs : see init_args
        '''
        self.__version__=__version__
        # print('...setting up keys...')
        # self.set_keys()
        self.keymap = cfg.KeyBinder()
        if help:self.keymap.show_help_cli()

        self.path,self.cif_file,self.bloch,self.pets = None,None,None,None
        self.multi_path,self.tifpath,self.multi = None,None,None
        self.rpl,self.im,self.Icols,self.rock,self.dyn_opt = None,None,[],None,1
        self.cond = '(Sw<1e-3) & (Vga>1e-5)'
        self.nfigs,self.kargs,self.u = 0,{},[0,0,1]
        self.pattern_args = {'Iopt':'csngt','Imax':5e4,'gs':0.025,'rmax':25,'Nmax':512}
        self.multi_args = {'data':'','mulslice':False,'NxNy':512,'Nhk':0,'hk_pad':None,
            'slice_thick':1,'i_slice':10,'ssh':'','opt':'srf'}
        self.xyz_params = {'theta':0,'lat_params':[20,20,100],'pad':1.0,'rep':None}
        self.dsp_d = dict(zip(self.keymap.dsp_keys,np.zeros((len(self.keymap.dsp_keys)),dtype=bool) ))
        self.opt_d = dict(zip(self.keymap.opt_keys,np.zeros((len(self.keymap.opt_keys)),dtype=bool) ))
        self.show_Sw,self.Swtol=False,2
        self.m={'I':1e3,'Ig':1e2,'Vg':10,'Swl':self.Swtol}

        # print('...load parameters...')
        if config:self.load_config(config)
        else:self.init_args(**kwargs)
        self.set_dsp_d()

        self.cutoff_r,self.mag_r,self.dthick_r=self.cutoff,self.mag,self.dthick
        if 'tag' in self.multi_args.keys():
            self.multi_args['tail']=self.multi_args['tag']
            self.multi_args.pop('tag')

        self.load_objects()
        self.multi_args['tail'] = self.tag+self.frame_str()
        self.multi_args['name'] = self.multi_path

        self.set_mode(self.mode)
        if not (self.tifpath or self.multi or self.rock):self.set_mode('rotate')

        # print('...Complete initialization...')

        #print('... init figure ...')
        self.pargs=pargs
        self.fig,self.ax = dsp.stddisp(**pargs )
        cid = self.fig.canvas.mpl_connect('key_press_event', self)
        self.call = self.call_update

        # print('... update and show ...')
        if self.multi_path:
            self.xyz_params['xyz']  = self.get_xyz()
        self.xyz_params['file'] = self.cif_file
        self.xyz_params['n']    = self.u
        if not self.thick : self.thick=100
        if not self.rock:
            self.update(fsolve=1)
            self.update_thickness()
        self.show()
        # plt.show()


    def init_args(self,path=None,cif_file=None,bloch=None,rock=None,pets=None,multi=None,tag='',
        orient=[0,0,5,5],u=None,Smax=0.01,Nmax=5,thick=None,dthick=5,thicks=(0,1000,1000),
        xylims=1.5,cutoff=50,dcutoff=10,mag=500,cmap=None,cmap_beams='jet',
        frame=1,incrate=1,drot=1,drotS=1,pets_opts='Pr',rot=0,rotS=0,dyn_opt=1,
        init_opts='R'):
        # F='Sw',kargs={},**kwargs):
        ''' configuration Initializer
    #### path setup (at least one of these need to be defined)
        - path      : str - path to the main folder
        - cif_file  : str - full path to .cif file (will be automatically fetched in path if not defined)
        - rock      : Rocking object or str for a rocking object
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
        self.cmap    = cmap
        self.rot ,self.drot  = rot, drot
        self.rotS,self.drotS = rotS, drotS

        #rotate related
        self.theta,self.phi,self.dtheta,self.dphi = orient
        self.dtheta_dphi = [self.dtheta,self.dphi]
        #frames related
        self.frame   = frame
        self.incrate = incrate
        #bloch related
        self.cmap_beams = cmap_beams
        self.thicks  = thicks #; self.F,self.fopts = F,self.df_F.loc[self.df_F.F==F].iloc[0].fopt
        #pets related
        self.pets_opts = pets_opts
        self.dyn_opt = dyn_opt
        #multislice related
        self.tag  = tag

        #init_opts
        self.mode       = ['frames','rotate']['R' in init_opts]
        self.show_hkl   = 'k' in init_opts
        self.show_i     = 'i' in init_opts
        self.show_Ig    = 'I' in init_opts
        self.show_Vg    = 'V' in init_opts
        self.show_Sw    = 'W' in init_opts
        self.show_z     = 'z' in init_opts
        self.show_u     = 'u' in init_opts
        self.show_v     = 'v' in init_opts
        self.save_bloch = 'B' in init_opts
        self.hold_mode  = ' ' in init_opts

        # print('...init path...')
        self.set_theta_phi_from_u(u)
        self.init_path(path,cif_file,bloch,rock,pets,multi)

    def load_config(self,config):
        with open(config,'rb') as f:d_config = pickle.load(f)
        for k,v in d_config.items(): self.__dict__[k] = v

    def init_path(self,path,cif_file,bloch,rock,pets,multi):
        if all([type(o)==type(None) for o in [path,cif_file,bloch,pets] ]):
            args = ['path','cif_file','bloch','pets']
            raise Exception('at least one of these must be defined : ',args)

        if type(path)==str:self.path = path
        if type(cif_file)==str:self.cif_file = cif_file
        if pets:self.init_pets(pets)
        if rock:
            self.rock = ut.load_pkl(rock)
            self.bloch=self.rock.load(self.frame)
            self.pets_opts,self.mode,self.nfigs = 'B','frames',self.rock.ts.size
            bloch = None
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

        if not type(self.pets)==pt.Pets :
            if type(self.pets)==int:
                print(colors.green+'ignoring pets'+colors.black)
                self.pets = None
            else:
                pets_path = os.path.join(self.path,'pets')
                if not os.path.exists(pets_path):pets_path=self.path
                pts_files = glob.glob(os.path.join(pets_path,'*.pts'))
                if len(pts_files)>=1:
                    print(colors.green+'found .pts file. Loading %s ' %pts_files[0]+colors.black)
                    if not self.pets==-1:
                        self.init_pets(pts_files[0])
                else:
                    print(colors.red+'no pts files found'+colors.black)

        if not type(self.multi_path)==str :
            if type(self.multi_path)==int:
                print(colors.green+'ignoring multi'+colors.black)
                self.multi_path = None
            else:
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
            if self.multi and not self.thick :
                self.thick = self.multi.thickness


        self.figpath = os.path.join(self.path,'figures')
        if not os.path.exists(self.figpath):
            p = Popen("mkdir %s" %self.figpath,shell=True,stderr=PIPE,stdout=PIPE);p.wait()

    def init_bloch(self,bloch):
        if type(bloch)   == bl.Bloch:self.bloch = bloch
        elif type(bloch) == str:self.bloch = load_bloch(file=bloch)
        else:raise Exception('bloch args need be NoneType, str or Bloch object ')
        if not type(self.path)     == str : self.path = self.bloch.path
        if not type(self.cif_file) == str : self.cif_file = self.bloch.cif_file
        self.u=self.bloch.u

    def init_multi(self,multi_path):
        if type(multi_path) in [str,int]:self.multi_path = multi_path

    def init_pets(self,pets):
        if type(pets)==pt.Pets : self.pets = pets
        elif type(pets)==str : self.pets = pt.Pets(pets,gen=1,cif_file=self.cif_file,dyn=self.dyn_opt)
        elif type(pets)==int : self.pets = pets;return
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
        'show_hkl','show_i','show_Ig','show_Vg','show_Sw','Swtol','show_z','show_u','show_v',     #shows
        'save_bloch','hold_mode','cmap_beams',#'F','fopts',
        'multi_args','xyz_params','cond',
        ]
        if not config_file:
            config_file = os.path.join(self.path,'config.pkl')
            config_file = easygui.filesavebox('save config','save',config_file)

        if config_file:
            d_config = {k:self.__dict__[k] for k in keys}
            with open(config_file,'wb') as f:pickle.dump(d_config,f,pickle.HIGHEST_PROTOCOL)
            print(colors.yellow+config_file+colors.green+' saved'+colors.black)

    ######################################################################
    ##### call
    ######################################################################
    def __call__(self,event):
        ''' event callback function
        normal mode : call = call_update
        find mode   : call = shortcut_call
        '''
        self.call(event)

    def call_update(self, event):
        # print(event.key)
        key = event.key
        dfk = self.keymap.df_keys
        alias_match = dfk.loc[(dfk.alias==key) & ( (dfk['mode']==self.mode) | (dfk['mode']=='any'))]
        if any(alias_match.index):
            # print(alias_match)
            key = alias_match.iloc[0].key
        # print(key)
        self.get_keys(key)
        self.do_update(key)

    def find_shortcut(self):
        print(colors.green+'find shortcut mode :'+colors.black)
        msg = '''
        Type a shortcut to check its effect.
        Use escape to switch back to normal mode'''
        print(msg)
        self.call = self.shortcut_call

    def shortcut_call(self,event):
        key = event.key
        if key=='escape' :
            self.call=self.call_update
            print(colors.green+'Back to normal mode'+colors.black)
        else:
            self.keymap.shortcut_call(key)

    ######################################################################
    ##### update
    ######################################################################
    def update_u(self):
        if self.mode == 'rotate':
            theta,phi = np.deg2rad([self.theta,self.phi])
            ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
            self.u = [st*cp,st*sp,ct]
        else:
            if self.tifpath:
                self.u = self.pets.get_beam(self.frame)
            elif self.multi:
                if 'u' in self.multi.__dict__.keys():
                    self.u = self.multi.u
                else:
                    print(colors.red+'warning : multislice does not contain u switching to rotate mode'+colors.black)
                    self.set_mode('rotate');return
                    self.u=[0,0,1]

        if self.multi_path:
            self.xyz_params['n']   = self.u
            self.xyz_params['xyz'] = self.get_xyz()


    def update_frame(self):
        if self.dsp_d['M']:
            loaded_multi = self.load_multi()
            if loaded_multi:
                state = self.multi.check_simu_state()       #;print(state)
                if state=='done':
                    self.multi.set_thicks()
                    self.multi.datpath=os.path.join(self.multi_path,'')
                    self.dthick = self.multi.dzs
                    self.update_thickness()
                else:
                    print(colors.blue+'multislice not done yet : '+colors.green+state+colors.black)
                    self.set_mode('rotate')
                    if not self.tifpath:return
            else:
                self.dsp_d['M'] = 0
                print(colors.red+'warning : no multislice simu found. Back to rotate mode'+colors.black)
                self.set_mode('rotate')
                if not self.tifpath:return

        if self.dsp_d['P'] and self.tifpath:
            self.rpl = self.pets.rpl.loc[self.pets.rpl.F==self.frame]
        if self.rock:
            self.bloch=self.rock.load(i=self.frame-1)
            self.u = self.bloch.u
        else:
            self.update(fsolve=1)#'B' in self.pets_opts)
        return 1

    def update_thickness(self):
        if self.dsp_d['M'] and self.mode=='frames':
            if self.thick<self.multi.thickness :
                self.iz,self.thick = self.multi.get_iz(self.thick,v=2)
                self.iz = int(self.iz)
            else:
                self.iz=None
                self.thick=self.multi.thickness
            qx,qy,I = self.multi.pattern(out=1,iz=self.iz,rot=self.rotS,v=1,**self.pattern_args)
            self.im=[qx,qy,I]
        if self.dsp_d['B'] or self.dsp_d['K']:
            self.bloch.set_thickness(self.thick)
            # print('update thickness',self.thick,self.bloch.thick)

    def update(self,fsolve=0):
        self.update_u()
        self.bloch.update_Nmax(self.Nmax)
        self.bloch.set_beam(u=self.u)                       #;print(self.u)
        self.bloch._set_excitation_errors(Smax=self.Smax)   #;print(self.Smax)
        self.bloch._set_Vg()
        self.bloch._set_kinematic()
        self.bloch.df_G['I']=0
        if fsolve:
            self.bloch._solve_Bloch()#opts='0v')
            self.bloch.set_thickness(thick=self.thick)
            if self.save_bloch:
                if self.mode=='frames':
                    bloch_file = os.path.join(self.path,'bloch','%s%s.pkl' %(self.tag,str(self.frame).zfill(4)))
                    self.bloch.save(bloch_file)
                else:self.bloch.save()
                self.save_bloch=0

    def set_mode(self,mode):
        self.mode = mode
        # self.show_im = {'frames':self.show_frames,'rotate':self.show_bloch}[self.mode]
        if self.mode=='frames':
            if self.update_frame():
                if self.dsp_d['M']:
                    self.dthick_r,self.dthick=self.dthick,self.multi.dzs
                    self.cutoff,self.mag = self.cutoff_r,self.mag_r
                    print('frames mode. Changing dthick=%d' %self.dthick)
        else:
            self.dthick=self.dthick_r
            self.cutoff_r,self.mag_r = self.cutoff,self.mag
            self.cutoff,self.mag = 50,500
            self.set_theta_phi_from_u(self.u)
            self.update(fsolve=1)

    ######################################################################
    ##### show
    ######################################################################
    def show_vals(self):
        if self.show_hkl : print(self.bloch.get_kin())
        if self.show_u   : print('[uvw]:');print(self.bloch.u)
        if self.show_v   : print('[uvw]_rec:');print(self.bloch.Kuvw0)
        if self.show_z   : print('zones:');print(self.bloch.get_zones())

        self.Icols=[] #,cond=[],[]
        if self.show_i  : self.Icols+=['I'] #;cond += ['(I>1e-4)']
        if self.show_Sw : self.Icols+=['Sw']#;cond += ['(Sw<1e-2)']
        if self.show_Ig : self.Icols+=['Ig']#;cond += ['(Ig>1e-2)']
        if self.show_Vg : self.Icols+=['Vg']#;cond += ['(Vga>1e-8)']
        # if any(cond):cond = self.cond+' & '.join(cond);print(cond)
        if self.cond and self.Icols:
            idx = self.bloch.get_beam(cond=self.cond)#;print(idx)
            cols = ['h','k','l']+self.Icols
            print(self.bloch.df_G.loc[idx,cols])

    def show_frames(self):
        title = ''
        sargs = {}
        cs = 'I'
        if not self.hold_mode:
            i,k = [(i,k) for i,(k,v) in enumerate(self.dsp_d.items()) if v][0]
            title += '%s, ' %list(EDdisp.legElts.keys())[i]
            cs = ['S','I'][k=='M']
        else:
            sargs={'alpha':0.75}
        if self.mode=='frames':
            title += 'frame %s, ' %self.frame_str()
        if self.multi_path and self.tag:
            title +='tag=%s, ' %(self.tag.replace('_',' '))
        if any([c in self.pets_opts for c in 'MKB']):
            title += 'thickness=%.1f$\AA$, ' %self.thick

        # sefl.refl = self.bloch.get_beam(cond=self.cond)#;print(hkl_idx)
        self.n,self.tol=5,1e-3
        self.refl = self.bloch.get_beam(cond=lambda dfG:bl_ut.strong_beams(dfG,tol=self.tol,n=self.n),index=False)
        # hkl_idx = self.refl
        print(self.refl)
        df = self.bloch.df_G
        df = df.drop(str((0,0,0)))
        EDdisp.show_frame(opts=self.pets_opts,mag=self.mag,rot=self.rot,
            df_pets=self.rpl,im_multi=self.im,df_bloch=df,hkl_idx=self.refl,
            ax=self.ax,title=title,xylims=self.xylims,single_mode=not self.hold_mode,
            cmap=self.cmap,cutoff=self.cutoff,sargs=sargs,cs=cs,**self.pargs)

    def show(self):
        self.show_vals()
        self.ax.cla()
        self.show_frames()
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
        self.multi_args['tail']  = self.tag+self.frame_str()
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
        mult_keys = ['NxNy','mulslice','slice_thick','i_slice','ssh','opt','Nhk','hk_pad']
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

    def open_multi(self):
        tails = self.show_tags(out=1)
        msg = '''Load multislice
        Available tags :
        %s
        ''' %tails
        tag = easygui.enterbox(msg,'load multislice',self.tag)
        if tag:
            if not tag==self.tag:
                self.tag=tag
                self.load_multi()
                self.update_frame()
    def load_multi(self):
            tag = self.tag+str(self.frame).zfill(3)
            multi = pp.load(self.multi_path,tag=tag)
            loaded_multi=False
            if multi:
                self.multi=multi
                self.multi.datpath = os.path.join(self.multi_path,'')
                self.multi.set_thicks() #; print("dz:",self.multi.dzs)
                self.multi.save()
                if 'xyz_params' in self.multi.__dict__.keys():
                    self.xyz_params = self.multi.xyz_params
                else:
                    print('warning : xyz_params not in multi')
                multi_args = list(self.multi_args.keys())
                for k in ['data','name','ssh','mulslice','opt']:multi_args.remove(k)
                self.multi_args['data'] = self.multi.data[0]
                for k in multi_args :
                    if k in self.multi.__dict__.keys():
                        self.multi_args[k]=self.multi.__dict__[k]
                    else:
                        print('warning : %s not in multi' %k)
                loaded_multi=True
            return loaded_multi

    def _run_multi(self):
        self._gen_xyz()
        # print('Solving multi : ', self.multi_args)
        self.multi = ms.Multislice(u=self.u,xyz_params=self.xyz_params,
            **self.multi_args)
        self.update_nfigs()

    def show_tags(self,out=0):
        if self.multi_path:
            pkls = glob.glob(os.path.join(self.multi_path,'*.pkl'))
            tails = ['_'.join(os.path.basename(pkl).split('_')[1:-1]) for pkl in pkls]
            tails = np.unique([tail[:-3] for tail in tails]) #;print(tails)
        if out:
            return tails
        else:
            print(colors.green+'available tails:'+colors.yellow,tails,colors.black)
    def get_xyz(self):
        return os.path.join(self.multi_path,self.tag+self.frame_str()+'.xyz')


    ######################################################################
    ##### Keys
    ######################################################################
    def do_update(self,key):
        if   key == self.keymap.dict_k['save_bloch']:self.update(fsolve=1)
        elif key in self.keymap.modes_keys : print('mode %s' %self.mode)
        elif key in self.keymap.frame_keys   and self.mode=='frames': print('frame %d' %self.frame)
        elif key in self.keymap.dframe_keys  and self.mode=='frames': print('increment rate : %d' %self.incrate)
        elif key in self.keymap.orient_keys  and self.mode=='rotate': print('theta=%.4f, phi=%.4f' %(self.theta,self.phi))
        elif key in self.keymap.dorient_keys and self.mode=='rotate': print('dtheta=%.1f, dphi=%.1f' %(self.dtheta,self.dphi))
        elif key in self.keymap.rot_keys     and self.mode=='frames': print('rot=%.1f, rotS=%.1f' %(self.rot,self.rotS))
        elif key in self.keymap.thick_keys      : print('thickness : %d' %self.thick);self.update_thickness()
        elif key in self.keymap.dthick_keys     : print('dthick=%d' %self.dthick)
        elif key in self.keymap.brightness_keys : print('brightness cutoff : %.1f' %self.cutoff)
        elif key in self.keymap.dbrightness_keys: print('dbrightness cutoff : %.1f' %self.dcutoff)
        elif key in self.keymap.rbrightness_keys: print('cutoff : %.1f, dcutoff : %.1f' %(self.cutoff,self.dcutoff))
        elif key in self.keymap.vals_keys       : self.show_vals()
        elif key in self.keymap.pets_keys       : print('frame options %s' %self.pets_opts)
        elif key == self.keymap.dict_k ['settings']:
            if self.mode=='rotate'   : self.update(fsolve=1)
            elif self.mode=='frames' : self.update_frame()
        #
        if self.mode=='rotate' :
            if key in self.keymap.orient_keys :
                self.update(fsolve=self.dsp_d['B'] and not self.bloch.solved )
        elif self.mode=='frames':
            if (key in self.keymap.frame_keys ) :
                self.update_frame()

        show_keys = key in self.keymap.show_keys
        show_keys = show_keys or (key in self.keymap.frame_keys  and self.mode=='frames')
        show_keys = show_keys or (key in self.keymap.orient_keys and self.mode=='rotate')
        if show_keys : self.show()#;print('showing')


    def get_keys(self,key):
        # print(key)
        dict_k = self.keymap.dict_k
        if   key == dict_k['settings']     : self.settings()
        elif key == dict_k['find_shortcut']: self.find_shortcut()
        elif key == dict_k['new_fig'] :
            config_file = os.path.join(self.path,'new_fig.pkl')
            self.save_config(config_file)
            Viewer(config=config_file)
            # v1.fig,v1.ax = dsp.stddisp()
            # cid = v1.fig.canvas.mpl_connect('key_press_event', v1)
        elif key == dict_k['help_cli']      : self.keymap.show_help_cli()
        elif key == dict_k['help_gui']      : self.keymap.show_help_gui()
        elif key == dict_k['help_simple']   : self.keymap.show_help_simple()
        elif key == dict_k['save_fig']      : dsp.saveFig(self.get_figname(),ax=self.ax)
        #config
        elif key == dict_k['save_config']   : self.save_config(os.path.join(self.path,'config.pkl'))
        elif key == dict_k['save_config_as']: self.save_config()
        elif key == dict_k['reload_config'] : self.open_config('config.pkl')
        elif key == dict_k['open_config']   : self.open_config()
        #bloch
        elif key == dict_k['save_bloch']    : self.save_bloch=1 #not self.save_bloch;print('%s saving bloch' %['not ',''][self.save_bloch])
        elif key == dict_k['solve_bloch']   : self.update(fsolve=1)
        #multi
        elif key == dict_k['open_multi']    : self.open_multi()
        elif key == dict_k['solve_multi']   : self.solve_multi()
        elif key == dict_k['set_multi']     : self.set_multi()
        elif key == dict_k['run_multi']     : self._run_multi()
        elif key == dict_k['show_multi']    : self.show_multi()
        elif key == dict_k['show_tags']     : self.show_tags()
        elif key == dict_k['log_info']      :
            if self.multi:print(self.multi.check_simu_state())
        #modes
        elif key == dict_k['frames_mode'] : self.set_mode('frames')
        elif key == dict_k['rotate_mode'] : self.set_mode('rotate')
        #in plane rotate
        if   key == dict_k['rot_up']   : self.rot=(self.rot+self.drot)%360
        elif key == dict_k['rot_down'] : self.rot=(self.rot-self.drot)%360
        elif key == dict_k['rotS_up']  : self.rotS=(self.rotS+self.drotS)%360
        elif key == dict_k['rotS_down']: self.rotS=(self.rotS-self.drotS)%360

        if self.mode=='rotate':
            # theta,phi
            if   key == dict_k['elev_up']   : self.theta=(self.theta+self.dtheta)%180
            elif key == dict_k['elev_down'] : self.theta=(self.theta-self.dtheta)%180
            elif key == dict_k['azim_up']   : self.phi=(self.phi+self.dphi)%360
            elif key == dict_k['azim_down'] : self.phi=(self.phi-self.dphi)%360
            # dtheta,dphi,dthick
            elif key == dict_k['delev_up']   : self.dtheta=min(self.dtheta+1,90)
            elif key == dict_k['delev_down'] : self.dtheta=max(self.dtheta-1,0)
            elif key == dict_k['dazim_up']   : self.dphi=min(self.dphi+1,90)
            elif key == dict_k['dazim_down'] : self.dphi=max(self.dphi-1,0)
        elif self.mode=='frames':
            # frame
            if   key == dict_k['frame_up']   : self.frame=min(self.frame+self.incrate,self.nfigs)
            elif key == dict_k['frame_down'] : self.frame=max(1,self.frame-self.incrate)
            # increment rate
            elif key == dict_k['inc_frame_up']  :self.incrate=min(self.incrate+1,100)
            elif key == dict_k['inc_frame_down']:self.incrate=max(1,self.incrate-1)

        #brightness
        if   key==dict_k['brightness up']     : self.cutoff=min(self.cutoff+self.dcutoff,1e5)
        elif key==dict_k['brightness down']   : self.cutoff=max(1e-5,self.cutoff-self.dcutoff)
        if   key==dict_k['dbrightness up']    : self.dcutoff*=2
        elif key==dict_k['dbrightness down']  : self.dcutoff/=2
        elif key==dict_k['brightness reset']  : self.cutoff  = 50
        elif key==dict_k['dbrightness reset'] : self.dcutoff = 10
        # thickness
        elif key == dict_k['thick_up']     : self.thick+=self.dthick
        elif key == dict_k['thick_down']   : self.thick=max(self.thick-self.dthick,self.dthick)
        elif key == dict_k['dthick_up']    : self.dthick=min(self.dthick+1,1000)
        elif key == dict_k['dthick_down']  : self.dthick=max(self.dthick-1,1)
        # vals
        elif key == dict_k['show_u']   : self.show_u   = not self.show_u  ;print('%sshowing [uvw]'     %['not ',''][self.show_u])
        elif key == dict_k['show_v']   : self.show_v   = not self.show_v  ;print('%sshowing [uvw]_rec' %['not ',''][self.show_v])
        elif key == dict_k['show_hkl'] : self.show_hkl = not self.show_hkl;print('%sshowing hkl'       %['not ',''][self.show_hkl])
        elif key == dict_k['show_I']   : self.show_i   = not self.show_i  ;print('%sshowing Istrong'   %['not ',''][self.show_i])
        elif key == dict_k['show_Ig']  : self.show_Ig  = not self.show_Ig ;print('%sshowing Ig_strong' %['not ',''][self.show_Ig])
        elif key == dict_k['show_Vg']  : self.show_Vg  = not self.show_Vg ;print('%sshowing Vg_strong' %['not ',''][self.show_Vg])
        elif key == dict_k['show_Sw']  : self.show_Sw  = not self.show_Sw ;print('%sshowing Sw_strong' %['not ',''][self.show_Sw])
        elif key == dict_k['show_z']   : self.show_z   = not self.show_z  ;print('%sshowing zones'     %['not ',''][self.show_z])
        elif key == dict_k['show_beams_vs_thickness']:
            if not self.bloch.solved : self.update(fsolve=1)
            # refl = self.bloch.get_beam(cond=lambda dfG:bl_rock.strong_beams(dfG,tol=0.01,n=5),index=False)
            # self.refl = self.bloch.get_beam(cond=lambda dfG:bl_ut.strong_beams(dfG,tol=0.01,n=10),index=False)

            self.bloch.show_beams_vs_thickness(self.thicks,cm=self.cmap_beams,refl=self.refl)
        elif key== dict_k['show_Exp'] and self.tifpath:
            self.pets.show_exp(frame=self.frame)
        # Display
        elif key == dict_k['hold']:
            self.hold_mode=not self.hold_mode
            if not self.hold_mode :
                k='K'
                if any(self.dsp_d.values()):k=self.pets_opts[0]
                self.dsp_d = self.dsp_d.fromkeys(self.dsp_d,False)
                self.dsp_d[k]=True
            print('hold mode ' + ['off','on'][self.hold_mode])

        do_up=False
        if key in self.keymap.opt_keys :
            self.opt_d[key] = not self.opt_d[key]
        elif key in self.keymap.dsp_keys:
            if self.hold_mode:
                self.dsp_d[key] = not self.dsp_d[key]
            else:
                self.dsp_d = self.dsp_d.fromkeys(self.dsp_d,False)
                self.dsp_d[key] = True
            do_up = do_up or (self.dsp_d['M'] and not 'M' in self.pets_opts and self.multi_path and not self.im)
            do_up = do_up or (self.dsp_d['P'] and not 'P' in self.pets_opts and self.tifpath )
            do_up = do_up or (self.dsp_d['B'] and not 'B' in self.pets_opts and not self.bloch.solved)
        self.set_pets_opts()
        if do_up:
            if self.mode=='frames' :self.update_frame()
            else:self.update(fsolve=1)

    ######################################################################
    ##### Misc
    ######################################################################
    def set_pets_opts(self):
        if not self.tifpath and self.dsp_d['P'] : self.dsp_d['P']=False;print(colors.red+'Processed data can not be displayed'+colors.black)
        if not self.multi   and self.dsp_d['M'] : self.dsp_d['M']=False;print(colors.red+'multislice can not be displayed'    +colors.black)
        if not self.hold_mode:
            if not any(self.dsp_d.values()):self.dsp_d['K']=True
        self.pets_opts =''.join([self.keymap.dsp_chars[c] for c,val in self.dsp_d.items() if val])
        self.pets_opts+=''.join([self.keymap.opt_chars[c] for c,val in self.opt_d.items() if val])
    def set_dsp_d(self):
        for k,v in self.keymap.dsp_chars.items():self.dsp_d[k] = v in self.pets_opts
        for k,v in self.keymap.opt_chars.items():self.opt_d[k] = v in self.pets_opts

    def set_theta_phi_from_u(self,u=None):
        if type(u) in [list,np.ndarray]:
            self.u = u/np.linalg.norm(u)
            x,y,z  = self.u
            self.theta = np.rad2deg(np.arccos(z))
            self.phi   = np.rad2deg(np.arctan2(y,x))
            print('u set to : ', str(self.u))
            print('theta=%.2f, phi=%.2f' %(self.theta,self.phi))

    def frame_str(self,pad=3):return '%s' %str(self.frame).zfill(pad)
    def get_figname(self):
        if self.mode=='frame':
            png_file = '_%s.png' %str(self.frame).zfill(3)
        else:
            png_file = '%s.png' %str(self.bloch.name)
        return os.path.join(self.figpath,png_file)
    def update_nfigs(self):
        if not self.tifpath and self.multi_path:
            self.nfigs = len(glob.glob(os.path.join(self.multi_path,'*.pkl')))

    def settings(self):
        fieldNames  =['Smax','Nmax','thick','dthick','cond',
            'xylims','mag','cutoff','Swtol','rot','rotS',
            'pets_opts','cmap','cmap_beams']
        if   self.mode=='rotate':
            self.dtheta_dphi=[self.dtheta,self.dphi]
            fieldNames+=['theta','phi','dtheta_dphi','thicks']
        elif self.mode=='frames':
            fieldNames+=['frame','incrate']
        if not self.cmap:self.cmap='None'
        fieldValues = [str(self.__dict__[f]) for f in fieldNames]
        dict_fv = EDdisp.multenterbox("Change settings","settings", fieldValues,fieldNames)
        if dict_fv:
            for f in fieldNames:
                if isinstance(self.__dict__[f],str):
                    self.__dict__[f] = dict_fv[f]
                else:
                    self.__dict__[f] = eval(dict_fv[f])

        if self.cmap == 'None':self.cmap=None
        self.set_dsp_d()
        #update dtheta,dphi
        if self.mode=='rotate':
            if type(self.dtheta_dphi) in [float,int]:
                self.dtheta_dphi=[self.dtheta_dphi]*2
            self.dtheta,self.dphi = self.dtheta_dphi
        self.m={'I':1e3,'Ig':1e2,'Vg':10,'Swl':self.Swtol}
