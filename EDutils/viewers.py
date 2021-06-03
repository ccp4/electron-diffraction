import numpy as np,matplotlib.pyplot as plt,pandas as pd
from utils import displayStandards as dsp
from utils import glob_colors as colors
import easygui
from blochwave import bloch

pd.set_option('precision',3)

class Rotate_Viewer:
    def __init__(self,bloch,orient=[0,0,5,5],u=None,Smax=0.01,Nmax=3,thick=100,
            F='L',xylims=None,init='h',**kwargs):
        '''Dynamic beam viewer :
        - bloch : Bloch object
        - orient : [theta,phi,dtheta,dphi] in degrees
        - u : set beam direction(takes preference over orient)
        - Smax,Nmax,thick,F : see blochwave.Bloch
        - init:'h'(show help at start) k(show_hkl) i(show_i) z(show_z) u(show_u) v(show in u in zone axis reciprocal space)
        '''
        self.set_keys()
        self.bloch = bloch
        self.theta,self.phi,self.dtheta,self.dphi = orient
        self.dtheta_dphi = [self.dtheta,self.dphi]
        self.dthick = 5
        self.beams_args = kwargs
        self.Smax  = Smax
        self.Nmax  = Nmax
        self.thick = thick
        self.thicks = (0,1000,1000)
        self.xylims = xylims
        self.F,self.fopts = F,self.df_F.loc[self.df_F.F==F].iloc[0].fopt
        self.show_hkl = 'k' in init
        self.show_i   = 'i' in init
        self.show_z   = 'z' in init
        self.show_u   = 'u' in init
        self.show_v   = 'v' in init

        self.fig,self.ax = dsp.stddisp()
        cid = self.fig.canvas.mpl_connect('key_press_event', self)

        self.set_theta_phi_from_u(u)
        self.update(fsolve=1)
        self.show()
        if 'h' in init:self.show_help()


    def __call__(self, event):
        # print(event.key)
        self.get_keys(event.key)
        self.do_update(event.key)

    ###################################
    ##### Update and show functions
    ###################################
    def update_u(self):
        theta,phi = np.deg2rad([self.theta,self.phi])
        ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
        self.u = [st*cp,st*sp,ct]

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

    def show_vals(self):
        if self.show_hkl : print(self.bloch.get_kin())
        if self.show_i   : self.bloch.get_Istrong()
        if self.show_u   : print('[uvw]:');print(self.bloch.u)
        if self.show_v   : print('[uvw]_rec:');print(self.bloch.Kuvw0)
        if self.show_z   : print('zones:');print(self.bloch.get_zones())

    def show(self):
        self.show_vals()
        self.ax.cla()
        self.bloch.show_beams(ax=self.ax,fig=self.fig,F=self.F,fopts=self.fopts,
            opt='',xylims=self.xylims,**self.beams_args)
        self.fig.canvas.draw()

    def settings(self):
        self.dtheta_dphi=[self.dtheta,self.dphi]
        fieldNames  = ['theta','phi','dtheta_dphi','Smax','Nmax','thick','dthick','thicks','F','xylims']
        fieldValues = [str(self.__dict__[f]) for f in fieldNames]
        dict_fv = multenterbox("Change settings","settings", fieldValues,fieldNames)
        if dict_fv:
            for f in fieldNames:
                if isinstance(self.__dict__[f],str):
                    self.__dict__[f] = dict_fv[f]
                else:
                    self.__dict__[f] = eval(dict_fv[f])

        #update dtheta,dphi
        if type(self.dtheta_dphi) in [float,int]:
            self.dtheta_dphi=[self.dtheta_dphi]*2
        self.dtheta,self.dphi = self.dtheta_dphi


    ###################################
    ##### Keys
    ###################################
    def do_update(self,key):
        if key in self.dorient_keys :print('dtheta=%.1f, dphi=%.1f' %(self.dtheta,self.dphi))
        elif key in self.dthick_keys:print('dthick=%d' %self.dthick)
        elif key in self.vals_keys  :self.show_vals()
        elif key in self.thick_keys :
            print('thickness : %d' %self.thick)
            if self.F in ['S','I','Ig']:self.bloch.set_thickness(self.thick)
        elif key in self.F_keys or key in self.Fnum_keys:
            print('displaying F="%s"' %self.F)
            if not self.bloch.solved:
                if self.F in ['S','I'] : self.update(fsolve=1)
        elif key in self.orient_keys + ['enter']:
            print('theta=%.4f, phi=%.4f' %(self.theta,self.phi))
            self.update(fsolve=self.F in ['S','I'])

        # if key in self.update_keys:self.update()
        if key in self.show_keys:self.show()

    def get_keys(self,key):
        # dict_k = self.dict_k
        if   key == self.dict_k['settings'] : self.settings()
        elif key == self.dict_k['help']     : self.show_help()
        #modes
        elif key==self.dict_k['frames_mode']:self.mode='frames';print('passing into mode %s' %self.mode)
        elif key==self.dict_k['rotate_mode']:self.mode='rotate';print('passing into mode %s' %self.mode)
        # theta,phi
        elif key == self.dict_k['elev_up']    : self.theta=(self.theta+self.dtheta)%180
        elif key == self.dict_k['elev_down']  : self.theta=(self.theta-self.dtheta)%180
        elif key == self.dict_k['azim_up']    : self.phi=(self.phi+self.dphi)%360
        elif key == self.dict_k['azim_down']  : self.phi=(self.phi-self.dphi)%360
        # dtheta,dphi,dthick
        elif key == self.dict_k['dthick_up'  ] : self.dthick=min(self.dthick+1,1000)
        elif key == self.dict_k['dthick_down'] : self.dthick=max(self.dthick-1,1)
        elif key == self.dict_k['delev_up'   ] : self.dtheta=min(self.dtheta+1,90)
        elif key == self.dict_k['delev_down' ] : self.dtheta=max(self.dtheta-1,0)
        elif key == self.dict_k['dazim_up'   ] : self.dphi=min(self.dphi+1,90)
        elif key == self.dict_k['dazim_down' ] : self.dphi=max(self.dphi-1,0)
        # thickness
        elif key == self.dict_k['thick_up'  ]:self.thick+=self.dthick
        elif key == self.dict_k['thick_down']:self.thick=max(self.thick-self.dthick,self.dthick)
        # vals
        elif key == self.dict_k['show_hkl']:self.show_hkl = not self.show_hkl;print('%sshowing hkl' %['not ',''][self.show_hkl])
        elif key == self.dict_k['show_I'  ]:self.show_i   = not self.show_i  ;print('%sshowing Istrong' %['not ',''][self.show_i])
        elif key == self.dict_k['show_u'  ]:self.show_u   = not self.show_u  ;print('%sshowing [uvw]' %['not ',''][self.show_u])
        elif key == self.dict_k['show_v'  ]:self.show_v   = not self.show_v  ;print('%sshowing [uvw]_rec' %['not ',''][self.show_v])
        elif key == self.dict_k['show_z'  ]:self.show_z   = not self.show_z  ;print('%sshowing zones' %['not ',''][self.show_z])
        elif key == self.dict_k['show_beams_vs_thickness']:
            if not self.bloch.solved :self.update(fsolve=1)
            self.bloch.show_beams_vs_thickness(self.thicks,strong=True)
        # solution display
        elif key in self.F_keys or key in self.Fnum_keys:
            loc = self.df_F.loc[(self.df_F.key==key) | (self.df_F.alias==key)].iloc[0]
            self.F,self.fopts = loc.F,loc.fopt


    def set_keys(self):
        self.dict_k = {}

        generic = ['settings','help']
        self.dict_k.update(dict(zip(generic,['enter','h'])))

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

        #thickness
        thicks  = ['thick_up','thick_down']
        dthicks = ['d'+s for s in thicks]
        self.thick_keys   = ['ctrl+t','ctrl+T']
        self.dthick_keys  = ['shift+ctrl+up','shift+ctrl+down']
        self.dict_k.update(dict(zip(thicks,self.thick_keys)))
        self.dict_k.update(dict(zip(dthicks,self.dthick_keys)))

        #vals
        vals = ['hkl','I','z','u','v','beams_vs_thickness']
        vals = ['show_'+s for s in vals]
        self.vals_keys = ['ctrl+'+c for c in 'hizuvb']
        self.dict_k.update(dict(zip(vals,self.vals_keys)))

        #display
        Fs  = ['L','Sw','Vg','S','I','Ig']
        fopt  = ['m','L','l']+['m']*2+['l']
        names = ['Lattice','Excitation error','Potential','Scattering','Intensities','Kinematic']
        names = ['show '+s for s in names]
        self.F_keys    = ['ctrl+'+c for c in 'LGVSIK']
        self.Fnum_keys = ['ctrl+'+c for c in '123456']
        self.df_F = pd.DataFrame.from_dict(dict(zip(['names','key','alias','F','fopt'],[names,self.F_keys,self.Fnum_keys,Fs,fopt])))
        self.df_F = self.df_F.set_index('names')

        self.update_keys=self.orient_keys+['enter']
        self.show_keys=self.update_keys+self.F_keys+self.Fnum_keys+self.thick_keys

        self.df_keys = pd.DataFrame.from_dict(self.dict_k,orient='index',columns=['key binding'])
        self.df_generic = self.df_keys.loc[generic]
        self.df_modes   = self.df_keys.loc[modes]
        self.df_orient  = self.df_keys.loc[rots+drots]
        self.df_thicks  = self.df_keys.loc[thicks+dthicks]
        self.df_vals    = self.df_keys.loc[vals]

    ###################################
    ##### Help
    ###################################
    def show_help(self):
        print(colors.green+'Shortcuts :  '+colors.black)
        print(colors.blue+'Generic :     '+colors.black);print(self.df_generic)
        print(colors.blue+'Modes :       '+colors.black);print(self.df_modes)
        print(colors.blue+'Orientation : '+colors.black);print(self.df_orient)
        print(colors.blue+'Thickness :   '+colors.black);print(self.df_thicks)
        print(colors.blue+'Show values : '+colors.black);print(self.df_vals)
        print(colors.blue+'Display :     '+colors.black);print(self.df_F)


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
