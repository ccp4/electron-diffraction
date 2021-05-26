import numpy as np,matplotlib.pyplot as plt,pandas as pd
from utils import displayStandards as dsp
import easygui
from blochwave import bloch

pd.set_option('precision',3)

class Rotate_Viewer:
    def __init__(self,bloch,u_params=[0,0,5,5],u=None,Smax=0.01,Nmax=3,thick=100,
            F='L',xylims=None,**kwargs):
        '''Dynamic beam viewer :
        - u_params : [theta,phi,dtheta,dphi] in degrees
        '''
        self.bloch = bloch
        self.theta,self.phi,self.dtheta,self.dphi = u_params
        self.dtheta_dphi = [self.dtheta,self.dphi]
        self.dthick = 5
        self.beams_args = kwargs
        self.Smax  = Smax
        self.Nmax  = Nmax
        self.thick = thick
        self.xylims = xylims
        self.Fs = {'L':'L','G':'Sw','V':'Vg','S':'S','I':'I'}
        self.F  = F
        self.show_hkl = 0
        self.show_i   = 0
        self.show_u   = 0

        self.fig,self.ax = dsp.stddisp()
        cid = self.fig.canvas.mpl_connect('key_press_event', self)

        self.set_theta_phi_from_u(u)
        self.solved=0
        self.update()
        self.show()

    def __call__(self, event):
        # print(event.key)
        orient_keys  = ['left','up','down','right']
        dorient_keys = ['ctrl'+s for s in orient_keys]
        F_keys       = ['ctrl+'+c for c in 'LGVSI']
        thick_keys   = ['ctrl+t','ctrl+T']
        dthick_keys  = ['shift+ctrl+up','shift+ctrl+down']
        vals_keys    = ['ctrl+'+c for c in 'hiu']

        if event.key=='enter':self.settings()
        # theta,phi
        elif event.key =='up'   : self.theta=(self.theta+self.dtheta)%180
        elif event.key =='down' : self.theta=(self.theta-self.dtheta)%180
        elif event.key =='left' : self.phi=(self.phi+self.dphi)%360
        elif event.key =='right': self.phi=(self.phi-self.dphi)%360
        # dtheta,dphi,dthick
        elif event.key =='shift+ctrl+up'   : self.dthick=min(self.dthick+1,1000)
        elif event.key =='shift+ctrl+down' : self.dthick=max(self.dthick-1,1)
        elif event.key =='ctrl+up'   : self.dtheta=min(self.dtheta+1,90)
        elif event.key =='ctrl+down' : self.dtheta=max(self.dtheta-1,0)
        elif event.key =='ctrl+right': self.dphi=min(self.dphi+1,90)
        elif event.key =='ctrl+left' : self.dphi=max(self.dphi-1,0)
        # thickness
        elif event.key=='ctrl+t':
            self.thick+=self.dthick
        elif event.key=='ctrl+T':
            self.thick=max(self.thick-self.dthick,self.dthick)
        # misc
        elif event.key=='ctrl+h':self.show_hkl= not self.show_hkl;print('%sshowing hkl' %['not ',''][self.show_hkl])
        elif event.key=='ctrl+i':self.show_i  = not self.show_i  ;print('%sshowing Istrong' %['not ',''][self.show_i])
        elif event.key=='ctrl+u':self.show_u  = not self.show_u  ;print('%sshowing u' %['not ',''][self.show_u])
        elif event.key in F_keys:self.F = self.Fs[event.key[-1]]  #;print('displaying F="%s"' %self.F)

        #### do the updates
        # solve=self.F in 'SI' and event.key in F_keys
        if event.key in dorient_keys:
            print('dtheta=%.1f, dphi=%.1f' %(self.dtheta,self.dphi))
        elif event.key in dthick_keys:
            print('dthick=%d' %self.dthick)
        elif event.key in vals_keys:
            self.show_vals()
        elif event.key in thick_keys:
            print('thickness : %d' %self.thick)
            if self.F in 'SI':self.bloch.set_thickness(self.thick)
        elif event.key in F_keys:
            print('displaying F="%s"' %self.F)
            if not self.solved:self.update()
        elif event.key in orient_keys + ['enter']:
            print('theta=%.4f, phi=%.4f' %(self.theta,self.phi))
            self.solved=not self.F in 'SI'
            self.update()

        update_keys=orient_keys+F_keys+thick_keys+['enter']
        if event.key in update_keys:
            self.show()

    def update_u(self):
        theta,phi = np.deg2rad([self.theta,self.phi])
        ct,st,cp,sp = np.cos(theta),np.sin(theta),np.cos(phi),np.sin(phi)
        self.u = [st*cp,st*sp,ct]

    def set_theta_phi_from_u(self,u=None):
        if type(u) in [list,np.ndarray]:
            self.u = u/np.linalg.norm(u)
            x,y,z  = self.u
            self.theta = np.arccos(z)
            self.phi   = np.arctan2(y,x)

    def update(self):
        self.update_u()
        self.bloch.update_Nmax(self.Nmax)
        self.bloch.set_beam(u=self.u)                       #;print(self.u)
        self.bloch._set_excitation_errors(Smax=self.Smax)   #;print(self.Smax)
        self.bloch._set_Vg()
        if not self.solved:
            self.bloch._solve_Bloch(opts='0v')
            self.bloch.set_thickness(thick=self.thick)
            self.solved=1
        else:
            self.solved=0

    def show_vals(self):
        if self.show_hkl:print(self.bloch.get_kin())
        if self.show_i  :self.bloch.get_Istrong()
        if self.show_u  :print('u:');print(self.bloch.u)

    def show(self):
        self.show_vals()

        self.fopts = {'Vg':'m','L':'m','Sw':'L','S':'m','I':'m'}[self.F]
        self.ax.cla()
        self.bloch.show_beams(ax=self.ax,fig=self.fig,F=self.F,fopts=self.fopts,
            opt='',xylims=self.xylims,**self.beams_args)
        self.fig.canvas.draw()

    def settings(self):
        self.dtheta_dphi=[self.dtheta,self.dphi]
        fieldNames  = ['theta','phi','dtheta_dphi','Smax','Nmax','thick','dthick','F','xylims']
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
