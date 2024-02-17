import importlib as imp
import numpy as np,pandas as pd
from wallpp import wallpaper as wallpp  ;imp.reload(wallpp)
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
# from multislice import mupy_utils as mut;imp.reload(mut)
from scattering import structure_factor as sf       #;imp.reload(sf)
from utils import displayStandards as dsp           #;imp.reload(dsp)
from utils import physicsConstants as cst           #;imp.reload(cst)
from EDutils import utilities as ut                 #;imp.reload(ut)
from utils import glob_colors as colors
from multislice import postprocess as pp            #;imp.reload(pp)
from . import util as bloch_util                    ;imp.reload(bloch_util)

class Bloch2D():
    def __init__(self,wallpp_args,solve_args):
        self._set_structure(wallpp_args)
        self.solve(**solve_args)

    def solve(self,Smax=0.02,Nmax=10,keV=200,u=[0,1],
        thicks=(1,100,100),
        eps=1,v=True):
        '''Solve '''
        self.eps=eps
        self.update_Nmax(Nmax)
        self.set_beam(keV,u)

        if Smax :
            self._set_excitation_errors(Smax)
            self._set_Vg()
        self._solve_Bloch(v=v)
        self._set_beams_vs_thickness(thicks)#;print('thicks solved')

    def set_beam(self,keV:float=200,u:Sequence[float]=[0,1]):
        self.keV = keV
        self.lam = cst.keV2lam(self.keV)

        self.k0  = 1/self.lam
        self.sig = cst.keV2sigma(self.keV)
        self.u   = np.array(u)/np.linalg.norm(u)
        self.K   = self.k0*self.u

    def _set_structure(self,wallpp_args):
        '''wallpaper file (see import_wallpp)'''
        # self.file = file
        # if not crys:crys=file
        # self.crys = mut.import_wallpp(crys)
        self.crys = wallpp.Wallpaper(**wallpp_args)
        self.lat_vec0 = self.crys.lattice_vectors
        self.lat_vec  = self.crys.reciprocal_vectors#/(2*np.pi)
        self.pattern  = self.crys.pattern_fract


    def update_Nmax(self,Nmax:int):
        """
        Update resolution/max order

        Parameters
        ----------
        Nmax
            maximum h,k order
        """
        self.Nmax=Nmax
        (h,k),(qx,qz) = ut.get_lattice2D(self.lat_vec,self.Nmax)
        self.lattice = [(h,k),(qx,qz)]
        self.hkF,self.Fhk = sf.structure_factor2D(self.pattern, 2*np.pi*self.lat_vec,
            hkMax=2*self.Nmax,eps=self.eps)


    def _get_central_beam():
        return self.get_beam(refl=[str(0,0)],index=True)[0]
    def _set_zones(self):return

    def _set_excitation_errors(self,Smax=0.02):
        ''' get excitation errors for Sg<Smax
        - Smax : maximum excitation error to be included
        '''
        K,k0 = self.K,self.k0
        (h,k),(qx,qy) = self.lattice

        Kx,Ky = K
        # Sw = np.abs(np.sqrt((Kx-qx)**2+(Ky-qy)**2) - K0)
        Sw = (k0**2 - ((Kx+qx)**2+(Ky+qy)**2)**2 )/(2*k0)
        if Smax:
            idx = Sw<Smax
            h,k = np.array([h[idx],k[idx]],dtype=int)
            qx,qy,Sw = qx[idx],qy[idx],Sw[idx]
        d = dict(zip(['h','k','qx','qy','Sw'],[h,k,qx,qy,Sw]))

        self.Smax   = Smax
        self.nbeams = Sw.size
        self.df_G   = pd.DataFrame.from_dict(d)
        self.solved = False
        # print(self.df_G)

    def _set_Vg(self):
        #in 2D Fhk is just the Fourier transform of fv
        hk  = self.df_G[['h','k']].values.astype(np.int)
        Fhk = self.Fhk.copy()
        V0_idx = np.array([2*self.Nmax]*2)
        Fhk[tuple(V0_idx)] = 0
        Fhk = np.array([ Fhk[tuple(hk_G+V0_idx)] for hk_G in hk])

        #It needs to be converted into a scattering amplitude [A]
        #Fg[no dim] = vg[kVA^2], sig=[1/kVA], k0[A-1]

        Fg = Fhk*self.sig*self.k0/np.pi
        # Vg[A^-2]
        Vg_G =  Fg/self.crys.area

        # self.set_beam_positions()
        self.df_G['Vg'] = Vg_G
        self.df_G['Vga'] = abs(Vg_G)
        self.df_G['Swl'] = bloch_util.logM(self.df_G['Sw'])
        self.df_G['L']  = np.ones(Vg_G.shape)
        self._set_zones()
        self.df_G.index = [str(tuple(h)) for h in self.df_G[['h','k']].values.astype(int)]
        self.df_G.loc[str((0,0)),'Vg'] = 0
        # print(self.df_G)



    def _solve_Bloch(self,show_H=False,v=False):
        ''' Diagonalize the Hamiltonian'''
        # Ug is a (2*Nmax+1)^3 tensor :
        # Ug[l] = [U(-N,-N,l) .. U(-N,0,l) .. U(-N,N,l)
        #          U( 0,-N,l) .. U( 0,0,l) .. U( 0,N,l)
        #          U(-N,-N,l) .. U( N,0,l) .. U( N,N,l)]

        hk = self.df_G[['h','k']].values.astype(np.int)
        Sg  = self.df_G.Sw.values
        pre = 1#/np.sqrt(1-cst.keV2v(self.keV)**2)
        Ug = pre*self.Fhk/(self.crys.area*np.pi)*self.eps #/3

        #####################
        # setting average potential to 0
        # Ug[U0_idx] = U0
        # and Ug(iG,jG) are obtained from Ug[h,k,l] where h,k,l = hlk_iG-hkl_jG
        #####################
        U0_idx = [2*self.Nmax]*2
        Ug[tuple(U0_idx)] = 0
        # if Vopt0 :Ug[tuple(U0_idx)] = 0   #setting average potential to 0

        if v:
            print(colors.blue+'...assembling %dx%d matrix...' %((Sg.shape[0],)*2)+colors.black)
        H = np.diag(Sg+0J)
        for iG,hk_G in enumerate(hk) :
            # hk_J=(hk_G-hk)[0]
            # print(hk,hk_G-hk)#,tuple(hk_J+U0_idx))
            U_iG = np.array([Ug[tuple(hk_J+U0_idx)] for hk_J in hk_G-hk]) #;print(V_iG.shape)
            H[iG,:] += U_iG/(2*self.k0)  #off diagonal terms as potential

        H *= 2*np.pi #to get same as felix
        self.H=H
        if show_H:self.show_H()

        if v:print(colors.blue+'...diagonalization...'+colors.black)
        self.gammaj,self.CjG = np.linalg.eigh(self.H) #;print(red+'Ek',lk,black);print(wk)
        self.invCjG = np.linalg.inv(self.CjG)
        self.solved = True


    def get_beam(self,refl):
        if isinstance(refl,tuple):
            return np.where(self.df_G.index==str(refl))[0][0]
        else:
            return [np.where(self.df_G.index==str(h))[0][0] for h in refl]

    def _set_beams_vs_thickness(self,thicks=None):
        """ get Scattering matrix as function of thickness for all beams
        - thicks : tuple(ti,tf,step) or list or np.ndarray thicknesses
        """
        print(colors.blue+'... beam vs thickness ...'+colors.black)
        if not type(thicks)==type(None):
            self._set_thicks(thicks)

        # id0 = self.get_beam(refl=[(0,0)],index=True)[0]#;print(id0)
        id0 = self.get_beam(refl=(0,0))
        gammaj,CjG,invCjG = self.gammaj,self.CjG,self.invCjG

        # import time
        # t0  = time.time()
        #### fast efficient
        # M = np.exp(2J*np.pi*np.outer(gammaj,self.z))*(invCjG[:,id0][:,None])
        M = np.exp(1J*np.outer(gammaj,self.z))*(invCjG[:,id0][:,None])
        St = CjG.dot(M)

        #### slower
        # St = np.zeros((self.df_G.shape[0],self.z.size),dtype=complex)
        # for iT,thick in enumerate(self.z):
        #     S = CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*thick))).dot(invCjG)
        #     # S=CjG.T.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG)
        #     # St[:,iT]=CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG.T)[0,:]
        #     St[:,iT] = S[:,id0]
        # print('done %.2f' %(time.time()-t0))

        self.Sz = St
        self.Iz = np.abs(self.Sz)**2
        Sw,Ug = self.df_G[['Sw','Vg']].values.T
        # Sz_kin = np.array([np.pi/self.k0*Ug*t*np.sinc(Sw*t) for t in self.z]).T
        # self.Iz_kin = np.abs(Sz_kin)**2

#########################################################
#### misc
#########################################################
    def convert2tiff(self):
        printf('not relevant in 2D')
    def show_Fhk(self,opts='m',**kwargs):
        """Displays structure factor over grid

        Parameters
        -----------
        opts
            see bloch_util:get_fz
        """
        fz,fz_str = bloch_util.get_fz(opts)
        Fhk   = self.Fhk
        nx,ny  = np.array((np.array(Fhk.shape)-1)/2,dtype=int)
        i,j    = np.meshgrid(np.arange(-nx,nx+1),np.arange(-ny,ny+1))
        return dsp.stddisp(scat=[i,j,fz(Fhk)],title=tle,**kwargs)


    def _set_thicks(self,thicks):
        if isinstance(thicks,tuple):thicks = np.linspace(*thicks)
        elif isinstance(thicks,float) or isinstance(thicks,int) : thicks=[thicks]
        self.z = np.array(thicks)

    def set_beam_positions(self):
        self.df_G['px'] = ut.project_beams2D(K=self.K,qxy=self.get_G())


    def show_beams_vs_thickness(self,
        thicks:Optional[Sequence]=None,
        beam_args:dict={},
        refl:Iterable=None,
        iZs:Iterable=slice(0,None,1),
        **kwargs
    ):
        """ display beams vs thickness

        Parameters
        ----------
        thicks
            thickness range (see :meth:~Bloch.set_beams_vs_thickness)
        iZs
            thickness index
        refl
            selected beams
        beam_args
            To pass to get_beam
        kwargs
            to pass to :meth:~EDutil.postprocess.plot_beam_thickness
        """

        if type(refl)==type(None):
            idx = np.arange(self.nbeams)
        else:
            idx  = self.get_beam(refl)

        z   = self.z[iZs]
        # Iz  = self.get_beams_vs_thickness(idx=idx,iZs=iZs,dict_opt=False)
        Iz = self.Iz[idx,iZs].copy()
        hkl = self.df_G.index[idx]

        beams=[hkl,z,None,None,Iz]
        return pp.plot_beam_thickness(beams,**kwargs)
