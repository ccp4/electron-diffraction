import importlib as imp
import easygui,os,glob
from subprocess import Popen,PIPE
import numpy as np,matplotlib.pyplot as plt,pandas as pd
from utils import displayStandards as dsp
from utils import glob_colors as colors
from blochwave.bloch import Bloch,load_Bloch#;imp.reload(mut)
from multislice.pets import Pets            #;imp.reload()
from multislice import mupy_utils as mut             ;imp.reload(mut)

pd.set_option('precision',3)
if plt.rcParams['keymap.save']:
    plt.rcParams['keymap.save'].remove('s')
    plt.rcParams['keymap.save'].remove('ctrl+s')
    plt.rcParams['keymap.grid'].remove('g')
    plt.rc('text', usetex=False)

class Viewer:
    def __init__(self,path=None,cif_file=None,bloch=None,pets=None,
        orient=[0,0,5,5],u=None,Smax=0.01,Nmax=3,thick=100,
        F='Sw',xylims=1.5,init_opts='R',cutoff=50,
        frame=1,pets_opts='Pr',rot=0,kargs={},**kwargs):
        '''Dynamic beam viewer :
    #### path setup (at least one of these need to be defined)
        - path      : str - path to the main folder
        - cif_file  : str - full path to .cif file (will be automatically fetched in path if not defined)
        - bloch     : Bloch object or str (will be automatically created if not defined)
        - pets      : Pets object or str to .pts file -
    #### Simulation setup
        - orient : [theta,phi,dtheta,dphi] in degrees
        - u : set beam direction(takes preference over orient)
        - Smax,Nmax,thick,F : see blochwave.Bloch
        - init_opts: str -
            - h(show help at start) R(rotate mode) B(save Bloch)
            - k(show_hkl) i(show_i) z(show_z) u(show_u) v(show u in zone axis reciprocal space)
        - kwargs : see Bloch.show_beams
        - kargs  : see Pets.show_frame
        '''
        self.init_path(path,cif_file,bloch,pets)
        self.cif_file = mut._get_cif_file(self.path,cif_file)
        self.cutoff  = cutoff
        self.xylims  = xylims
        self.Smax  = Smax
        self.Nmax  = Nmax
        self.thick = thick
        self.dthick = 5
        #rotate related
        self.theta,self.phi,self.dtheta,self.dphi = orient
        self.dtheta_dphi = [self.dtheta,self.dphi]
        self.thicks = (0,1000,1000)
        self.beams_args = kwargs
        self.set_theta_phi_from_u(u)
        #frames related
        self.frame   = frame
        self.incrate = 1
        self.rot     = rot
        self.pets_opts = pets_opts+'q'
        self.kargs = kargs

        self.set_keys()
        self.F,self.fopts = F,self.df_F.loc[self.df_F.F==F].iloc[0].fopt
        self.mode,self.show_im = [['frames',self.show_exp],['rotate',self.show_sim]]['R' in init_opts]
        self.show_hkl   = 'k' in init_opts
        self.show_i     = 'i' in init_opts
        self.show_z     = 'z' in init_opts
        self.show_u     = 'u' in init_opts
        self.show_v     = 'v' in init_opts
        self.save_bloch = 'B' in init_opts
        if 'h' in init_opts:self.show_help()

        self.fig,self.ax = dsp.stddisp()
        cid = self.fig.canvas.mpl_connect('key_press_event', self)

        self.update(fsolve=1)
        self.show()

    def init_path(self,path,cif_file,bloch,pets):
        if all([type(o)==type(None) for o in [path,cif_file,bloch,pets] ]):
            args = ['path','cif_file','bloch','pets']
            raise Exception('at least one of these must be defined : ',args)
        self.path,self.cif_file,self.bloch,self.pets,self.tifpath = None,None,None,None,None

        if type(path)==str:self.path = path
        if type(cif_file)==str:self.cif_file = cif_file
        if pets:self.init_pets(pets)
        if bloch:self.init_bloch(bloch)

        if not type(self.path)==str:
            self.path = os.path.dirname(self.cif_file)
        if not type(self.cif_file)==str:
            self.cif_file = mut._get_cif_file(self.path,cif_file)
        if not type(self.bloch)==Bloch:
            self.bloch = Bloch(self.cif_file,path=self.path)
        if not type(self.pets)==Pets:
            pts_files = glob.glob(os.path.join(self.path,'*.pts'))
            if len(pts_files)>=1:
                print(colors.green+'loading found .pts file %s ' %pts_files[0]+colors.black)
                self.init_pets(pts_files[0])
            if not len(pts_files)==1:
                msg ='''%d pts files  found :%s''' %(len(pts_files),str(pts_files))

        if bloch and pets :
            if not self.bloch.cif_file == self.pets.cif_file:
                raise Exception('bloch and pets seems to have different cif_files')

        self.figpath = os.path.join(self.path,'figures')
        if not os.path.exists(self.figpath):
            p = Popen("mkdir %s" %self.figpath,shell=True,stderr=PIPE,stdout=PIPE);p.wait()

    def init_bloch(self,bloch):
        if type(bloch) == Bloch:self.bloch = bloch
        elif type(bloch) == str:self.bloch = load_Bloch(file=bloch)
        else:raise Exception('bloch args need be NoneType, str or Bloch object ')
        if not type(self.path)      == str : self.path = self.bloch.path
        if not type(self.cif_file)  == str : self.cif_file = self.bloch.cif_file

    def init_pets(self,pets):
        if type(pets)== Pets : self.pets = pets
        elif type(pets)==str : self.pets = Pets(pets,gen=1)
        else:raise Exception('pets args need be NoneType, str or Pets object ')
        if type(self.path)      == str : self.path = self.pets.path
        if type(self.cif_file)  == str : self.cif_file = self.pets.cif_file
        self.nfigs    = self.pets.nFrames
        self.tifpath  = os.path.join(self.pets.path,'tiff')
        self.fmt,self.load = mut.find_format(self.tifpath)
        self.figs = np.sort(glob.glob(self.tifpath+'/*.%s' %self.fmt))#;print(self.figs)

    def __call__(self, event):
        # print(event.key)
        self.get_keys(event.key)
        self.do_update(event.key)

    def get_figname(self):
        if self.mode=='frame':
            png_file = '_%s.png' %str(self.frame).zfill(3)
        else:
            png_file = '%s.png' %str(self.bloch.name)
        return os.path.join(self.figpath,png_file)

    ######################################################################
    ##### Update and show functions
    ######################################################################
    def update_u(self):
        if self.mode == 'rotate':
            theta,phi = np.deg2rad([self.theta,self.phi])
            ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
            self.u = [st*cp,st*sp,ct]
        else:
            self.u = self.pets.uvw[self.frame-1]
    def set_theta_phi_from_u(self,u=None):
        if type(u) in [list,np.ndarray]:
            self.u = u/np.linalg.norm(u)
            x,y,z  = self.u
            self.theta = np.rad2deg(np.arccos(z))
            self.phi   = np.rad2deg(np.arctan2(y,x))
            print('u set to : ', str(self.u))
            print('theta=%.2f, phi=%.2f' %(self.theta,self.phi))

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

    def show_vals(self):
        if self.show_hkl : print(self.bloch.get_kin())
        if self.show_i   : self.bloch.get_Istrong()
        if self.show_u   : print('[uvw]:');print(self.bloch.u)
        if self.show_v   : print('[uvw]_rec:');print(self.bloch.Kuvw0)
        if self.show_z   : print('zones:');print(self.bloch.get_zones())

    def show_exp(self):
        fig = self.figs[self.frame-1]
        figname = os.path.basename(fig)
        print(colors.yellow+fig+colors.black)
        im = self.load(fig)

        self.pets.show_frame(frame=self.frame,thick=self.thick,Imag=self.cutoff,
            opts=self.pets_opts,rot=self.rot,ax=self.ax,title=figname,
            xylims=self.xylims,**self.kargs)

    def show_sim(self):
        self.bloch.show_beams(ax=self.ax,fig=self.fig,F=self.F,fopts=self.fopts,
            opt='',xylims=self.xylims,gridOn=0,**self.beams_args)

    def show(self):
        self.show_vals()
        self.ax.cla()
        self.show_im()
        self.fig.canvas.draw()


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
        elif key in self.dthick_keys     : print('dthick=%d' %self.dthick)
        elif key in self.brightness_keys : print('brightness cutoff : %d' %self.cutoff)
        elif key in self.vals_keys       : self.show_vals()
        elif key in self.thick_keys :
            print('thickness : %d' %self.thick)
            if self.F in ['S','I','Ig']:self.bloch.set_thickness(self.thick)
        elif key in self.F_keys or key in self.Fnum_keys:
            print('displaying F="%s"' %self.F)
            if not self.bloch.solved:
                if self.F in ['S','I'] : self.update(fsolve=1)

        if (key in self.orient_keys or key in ['enter']) and self.mode=='rotate' :
            self.update(fsolve=self.F in ['S','I'])
        elif (key in self.frame_keys or key in ['enter'] ) and self.mode=='frames':
            self.update()

        show_keys = key in self.show_keys
        show_keys = show_keys or (key in self.frame_keys+self.pets_keys+self.Pnum_keys) and self.mode=='frames'
        show_keys = show_keys or (key in self.orient_keys+self.F_keys+self.Fnum_keys) and self.mode=='rotate'
        if show_keys : self.show()


    def get_keys(self,key):
        if   key == self.dict_k['settings'] : self.settings()
        elif key == self.dict_k['help']     : self.show_help()
        elif key == self.dict_k['save_fig'] : dsp.saveFig(self.get_figname(),ax=self.ax)
        elif key == self.dict_k['save_bloch'] : self.save_bloch=1;print('%s saving bloch' %['not ',''][self.save_bloch])
        #modes
        elif key==self.dict_k['frames_mode'] and self.tifpath :
            self.mode,self.show_im = 'frames',self.show_exp
            self.update()
        elif key==self.dict_k['rotate_mode'] :
            self.mode,self.show_im = 'rotate',self.show_sim
            self.set_theta_phi_from_u(self.u)
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
            if   key in [self.dict_k['frame_up']  ,self.dict_k['frame_right']] : self.frame=min(self.frame+self.incrate,self.nfigs)
            elif key in [self.dict_k['frame_down'],self.dict_k['frame_left']]  : self.frame =max(1,self.frame-self.incrate)
            # increment rate
            elif key==self.dict_k['inc_frame_up']  :self.incrate=min(self.incrate+1,100)
            elif key==self.dict_k['inc_frame_down']:self.incrate=max(1,self.incrate-1)

        #brightness
        if   key==self.dict_k['brightness up']    : self.cutoff=min(self.cutoff+5,500)
        elif key==self.dict_k['brightness down']  : self.cutoff=max(1,self.cutoff-5)
        elif key==self.dict_k['brightness reset'] : self.cutoff=50
        # thickness
        elif key == self.dict_k['thick_up']     : self.thick+=self.dthick
        elif key == self.dict_k['thick_down']   : self.thick=max(self.thick-self.dthick,self.dthick)
        elif key == self.dict_k['dthick_up']    : self.dthick=min(self.dthick+1,1000)
        elif key == self.dict_k['dthick_down']  : self.dthick=max(self.dthick-1,1)
        # vals
        elif key == self.dict_k['show_u'] : self.show_u = not self.show_u  ;print('%sshowing [uvw]' %['not ',''][self.show_u])
        elif key == self.dict_k['show_v'] : self.show_v = not self.show_v  ;print('%sshowing [uvw]_rec' %['not ',''][self.show_v])
        # Bloch related
        if self.mode=='rotate':
            if   key == self.dict_k['show_hkl'] : self.show_hkl = not self.show_hkl;print('%sshowing hkl' %['not ',''][self.show_hkl])
            elif key == self.dict_k['show_I']   : self.show_i   = not self.show_i  ;print('%sshowing Istrong' %['not ',''][self.show_i])
            elif key == self.dict_k['show_z']   : self.show_z   = not self.show_z  ;print('%sshowing zones' %['not ',''][self.show_z])
            elif key == self.dict_k['show_beams_vs_thickness']:
                if not self.bloch.solved :self.update(fsolve=1)
                self.bloch.show_beams_vs_thickness(self.thicks,strong=True)
            elif (key in self.F_keys or key in self.Fnum_keys):
                loc = self.df_F.loc[(self.df_F.key==key) | (self.df_F.alias==key)].iloc[0]#;print(loc)
                self.F,self.fopts = loc.F,loc.fopt
        # Pets display
        elif self.mode=='frames':
            vals = [c in self.pets_opts for c in self.pets_chars]
            if key in self.pets_keys :
                for i,(c,k) in enumerate(zip(self.pets_chars,self.pets_keys)):
                    if key==k:vals[i] = not vals[i]
            elif key in self.Pnum_keys:
                i = int(key[-1])
                vals[i]=not vals[i]
            self.pets_opts='q'+''.join([c for c,val in zip(self.pets_chars,vals) if val])

    def set_keys(self):
        self.dict_k = {}

        generic = ['settings','help','save_fig','save_bloch']
        self.dict_k.update(dict(zip(generic,['enter','h','s','S'])))

        #modes
        modes = ['frames_mode','rotate_mode']
        self.modes_keys = ['ctrl+alt+'+c for c in '12']
        self.dict_k.update(dict(zip(modes,self.modes_keys)))

        #orientation
        rots  = ['elev_up','elev_down','azim_up','azim_down']
        drots = ['d'+s for s in rots]
        self.orient_keys  = ['up','down','right','left']
        self.dorient_keys = ['ctrl+'+s for s in self.orient_keys]
        self.dict_k.update(dict(zip(rots,self.orient_keys)))
        self.dict_k.update(dict(zip(drots,self.dorient_keys)))

        #frames
        frames = ['frame_'+s for s in ['up','right','down','left']]
        dframes = ['inc_'+s for s in frames]
        self.frame_keys  = ['up','right','down','left']
        self.dframe_keys = ['ctrl+'+s for s in self.frame_keys]
        self.dict_k.update(dict(zip(frames,self.frame_keys)))
        self.dict_k.update(dict(zip(dframes,self.dframe_keys)))

        #thickness
        thicks  = ['thick_up','thick_down']
        dthicks = ['d'+s for s in thicks]
        self.thick_keys   = ['ctrl+t','ctrl+T']
        self.dthick_keys  = ['shift+ctrl+up','shift+ctrl+down']
        self.dict_k.update(dict(zip(thicks,self.thick_keys)))
        self.dict_k.update(dict(zip(dthicks,self.dthick_keys)))

        #brightness
        self.brightness_keys = ['pagedown','pageup','ctrl+r']
        bright = ['brightness '+s for s in ['down','up','reset']]
        self.dict_k.update(dict(zip(bright,self.brightness_keys)))

        #vals
        vals = ['hkl','I','z','u','v','beams_vs_thickness']
        vals = ['show_'+s for s in vals]
        self.vals_keys = ['ctrl+'+c for c in 'hizuvb']
        self.dict_k.update(dict(zip(vals,self.vals_keys)))

        #display Bloch
        Fs    = ['L','Sw','Vg','S','I','Ig']#,'E']
        fopt  = ['m','L' ,'l' ,'m','m','l' ]#,'m']
        names = ['Lattice','Excitation error','Potential','Scattering','Intensities','Kinematic']
        names = ['show '+s for s in names]
        self.F_keys    = ['ctrl+'+c for c in 'LGVSIK']
        self.Fnum_keys = ['ctrl+'+c for c in '123456']
        self.df_F = pd.DataFrame.from_dict(dict(zip(['names','key','alias','F','fopt'],[names,self.F_keys,self.Fnum_keys,Fs,fopt])))
        self.df_F = self.df_F.set_index('names')

        #display Pets
        self.pets_chars = 'EPSKhkr'
        pets = ['Exp','Processed','Simulated','Kinematic','hkl_exp','hkl_kin','rings']
        pets = ['show_'+s for s in pets]
        self.pets_keys = ['E','P','S','K','j','k','g']
        self.Pnum_keys = ['ctrl+'+c for c in '1234']
        self.dict_k.update(dict(zip(pets,self.pets_keys)))


        #### Keys that trigger the show
        self.show_keys = self.modes_keys.copy()
        self.show_keys += self.thick_keys+self.brightness_keys +['enter','S']

        #### Data frame for the help command
        self.df_keys = pd.DataFrame.from_dict(self.dict_k,orient='index',columns=['key'])
        self.df_keys['alias'] = ['']*self.df_keys.shape[0]

        self.df_generic = self.df_keys.loc[generic]
        self.df_modes   = self.df_keys.loc[modes]
        self.df_orient  = self.df_keys.loc[rots+drots]
        self.df_frames  = self.df_keys.loc[frames+dframes]
        self.df_thicks  = self.df_keys.loc[thicks+dthicks]
        self.df_bright  = self.df_keys.loc[bright]
        self.df_vals    = self.df_keys.loc[vals]
        self.df_bloch   = self.df_F[['key','alias']]
        self.df_pets    = self.df_keys.loc[pets]
        self.df_pets['alias'].iloc[:4] = self.Pnum_keys

    ###################################
    ##### Help
    ###################################
    def show_help(self):
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

    def settings(self):
        fieldNames  =['Smax','Nmax','thick','dthick','xylims']
        if   self.mode=='rotate':
            self.dtheta_dphi=[self.dtheta,self.dphi]
            fieldNames+=['theta','phi','dtheta_dphi','thicks','F']
        elif self.mode=='frames':
            fieldNames+=['frame','incrate','rot','pets_opts']
        fieldValues = [str(self.__dict__[f]) for f in fieldNames]
        dict_fv = multenterbox("Change settings","settings", fieldValues,fieldNames)
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

######################################################################
##### Misc
######################################################################
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
