"""Bloch wave solver"""
import importlib as imp
import numpy as np,pandas as pd,pickle5,os,glob,tifffile
from crystals import Crystal
from EDutils import display as EDdisp               #;imp.reload(EDdisp)
from multislice import mupy_utils as mut            #;imp.reload(mut)
from multislice import postprocess as pp            #;imp.reload(pp)
from scattering import structure_factor as sf       #;imp.reload(sf)
from scattering import scattering_factors as scatf  #;imp.reload(scatf)
from utils import displayStandards as dsp           #;imp.reload(dsp)
from utils import physicsConstants as cst           #;imp.reload(cst)
from utils import glob_colors as colors,handler3D as h3d
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
# from EDutils.viewers import Viewer            #;imp.reload(mut)
from EDutils import viewers            #;imp.reload(viewers)
from .util import*

class Bloch:
    """
    Bloch wave simulation class

    Parameters
    ---------
    cif_file
        structure (see import_crys)
    name
        optional name (used to save object see set_name)
    path
        A path for the simulation folder
    beam
        dictionary passed to :meth:`~Bloch.set_beam`
    keV
        electron wavelength (see :meth:`~Bloch.set_beam`)
    u
        beam direction [ux,uy,uz] in reciprocal space (see :meth:`~Bloch.set_beam`)
    Nmax
        max order of reflections/resolution (see :meth:`~Bloch.update_Nmax`)
    Smax
        maximum excitation error (see :meth:`~Bloch.solve`)
    solve
        solve the system by diagonalizing the Bloch matrix
    kwargs
        arguments to be passed to :meth:`~Bloch.solve`
    """
    def __init__(self,
        cif_file:str,name:str='',path:str='',
        beam:Optional[dict]={},keV:float=200,u:Sequence[float]=[0,0,1],
        Nmax:int=1,
        Smax:float=0.2,
        solve:bool=True,**kwargs,
    ):

        self.cif_file = cif_file
        self.crys     = mut.import_crys(cif_file)
        self.lat_vec0 = np.array(self.crys.lattice_vectors)
        self.lat_vec  = np.array(self.crys.reciprocal_vectors)/(2*np.pi)
        self.pattern  = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in self.crys.atoms] )
        self.solved = False
        self.Nmax   =  0
        self.thick  = 100
        self.thicks = self._set_thicks((0,1000,1000))

        self.update_Nmax(Nmax)
        beam_args={'keV':keV,'u':u}
        beam_args.update(beam)
        self.set_beam(**beam_args)
        self.set_name(name,path)
        self._set_excitation_errors(Smax)
        self._set_Vg()

        if solve:
            self.solve(Smax=Smax,**kwargs)
        # if show_thicks or 't' in opts:
        #     self.show_beams_vs_thickness(strong=['I'])
        # if 's' in opts:
        #     self.save()

    def set_name(self,name='',path=''):
        """
        Set the basename for saving a Bloch obj as path+name+'.pkl'

        Parameters
        ----------
        name
            The name of the simulation(can include the full path)
        path
            The path to the simulation folder


        .. NOTE::
            By default the simulation name is '<cif><zone axis>_<keV>keV_bloch'
        """
        basefile = self.cif_file.replace('.cif','')
        if not name:
            u_str = ''.join(['%d' %np.round(u) for u in self.Kabc0])
            name='%s%s_%dkeV_bloch' %(os.path.basename(basefile),u_str,np.round(self.keV))
        if not path:path=os.path.dirname(basefile)
        self.path = path                #; print(self.path)
        self.name = name                #;print('name:',self.name)

    def update_Nmax(self,Nmax:int):
        """
        Update resolution/max order (lattice and Fhkl)

        Parameters
        ----------
        Nmax
            maximum h,k,l order
        """
        if type(Nmax) in [int,np.int64] :
            if not Nmax==self.Nmax:
                self.Nmax=Nmax
                (h,k,l),(qx,qy,qz) = mut.get_lattice(self.lat_vec,self.Nmax)
                self.lattice = [(h,k,l),(qx,qy,qz)]
                self.hklF,self.Fhkl = sf.structure_factor3D(self.pattern, 2*np.pi*self.lat_vec,hklMax=2*self.Nmax)

    def set_beam(self,
        keV:float=200,
        # lam:float=None,
        u:Sequence[float]=[0,0,1],
        K:Optional[Sequence[float]]=None,
        u0:Optional[Sequence[int]]=None,
    ):
        """ Set the beam direction and wavelength

        Parameters
        ----------
        keV
            beam wavelength in keV
        u
            beam direction
        K
            update beam direction and wavelength (overwrites keV)
        u0
            beam direction zone axis notation


        .. note::

            The order of priority for setting up the beam is K,u0,u
        """
        if isinstance(K,list) or isinstance(K,np.ndarray):
            self.k0 = np.linalg.norm(K)
            self.u  =  K/self.k0
        else:
            if isinstance(u0,list) or isinstance(u0,np.ndarray):
                u = np.array(u0).dot(self.lat_vec)
            self.u  = np.array(u)/np.linalg.norm(u)
            self.k0 = 1/cst.keV2lam(keV)

        self.lam = 1/self.k0
        self.keV = cst.lam2keV(self.lam)
        self.sig = cst.keV2sigma(self.keV)
        # self.u   = self.K/self.k0
        self.K   = self.k0*self.u
        self.Kuvw = self.lat_vec.dot(self.K)  #projection in reciprocal crystal basis
        self.Kabc = self.lat_vec0.dot(self.K) #projection in crystal basis
        self.Kabc0 = self.Kabc/(abs(self.Kabc)[abs(self.Kabc)>0.01].min()) #;print(Kabc)
        self.Kuvw0 = self.Kuvw/(abs(self.Kuvw)[abs(self.Kuvw)>0.01].min()) #;print(Kabc)

    def set_thickness(self,thick:float=100):
        """set thickness and update beams

        Parameters
        ----------
        thick : thickness in A
        """
        if type(thick) in [int,float] :self.thick=thick
        self._set_kinematic()
        if self.solved:self._set_intensities()

    def solve(self,
        Smax:Optional[float]=None,hkl:Optional[Iterable[int]]=None,
        Nmax:Optional[int]=None,
        beam:Optional[dict]={},
        thick:float=None,thicks:Sequence[float]=None,
        opts:str='sv0',
    ):
        """ Diagonalize the Blochwave matrix

        Parameters
        ----------
        Smax
            maximum excitation error to be included
        hkl
            Beams to include
        Nmax
            max order of reflections/resolution (see :meth:`~Bloch.update_Nmax` )
        beam
            dictionary passed to :meth:`~Bloch.set_beam`

        thick
            thickness of crystal (can be modified without resolving)
        thicks
            range of thickness [z_min,z_max,z_step]
        opts
            s(save) t(show intensities) z(show beam vs thickness) 0(Vopt0) v(verbose) H(show H)


        .. note ::

            If hkl is specified, Smax is not taken into account
        """
        if Nmax : self.update_Nmax(Nmax)
        if beam : self.set_beam(**beam)
        if Smax or isinstance(hkl,np.ndarray):
            self._set_excitation_errors(Smax,hkl)
            self._set_Vg()
        self._solve_Bloch(show_H='H' in opts,Vopt0='0' in opts,v='v' in opts)
        ##### postprocess
        if 's' in opts:
            self.save()
        if thick or 't' in opts:
            self.set_thickness(thick)
        if not type(thicks)==type(None) or 'z' in opts:
            self._set_beams_vs_thickness(thicks)#;print('thicks solved')

    ################################################################################
    #### private
    ################################################################################
    def _solve_Bloch(self,show_H=False,Vopt0=True,v=False):
        ''' Diagonalize the Hamiltonian'''
        # Ug is a (2*Nmax+1)^3 tensor :
        # Ug[l] = [U(-N,-N,l) .. U(-N,0,l) .. U(-N,N,l)
        #          U( 0,-N,l) .. U( 0,0,l) .. U( 0,N,l)
        #          U(-N,-N,l) .. U( N,0,l) .. U( N,N,l)]

        hkl = self.df_G[['h','k','l']].values
        Sg  = self.df_G.Sw.values
        pre = 1/np.sqrt(1-cst.keV2v(self.keV)**2)

        Ug = pre*self.Fhkl/(self.crys.volume*np.pi) #/3

        #####################
        # setting average potential to 0
        # Ug[U0_idx] = U0
        # and Ug(iG,jG) are obtained from Ug[h,k,l] where h,k,l = hlk_iG-hkl_jG
        #####################
        U0_idx = [2*self.Nmax]*3
        Ug[tuple(U0_idx)] = 0
        # if Vopt0 :Ug[tuple(U0_idx)] = 0   #setting average potential to 0

        if v:print(colors.blue+'...assembling %dx%d matrix...' %((Sg.shape[0],)*2)+colors.black)
        H = np.diag(Sg+0J)
        for iG,hkl_G in enumerate(hkl) :
            U_iG = np.array([Ug[tuple(hkl_J+U0_idx)] for hkl_J in hkl_G-hkl]) #;print(V_iG.shape)
            # print('G : ',hkl_G)
            # print('U_G:',Ug[tuple(U0_idx+hkl_G)])
            # print('idx iG:',[tuple(hkl_J+U0_idx) for hkl_J in hkl_G-hkl])
            # print(U_iG)
            H[iG,:] += U_iG/(2*self.k0)  #off diagonal terms as potential

        self.H = H*2*np.pi #to get same as felix
        if show_H:self.show_H()

        if v:print(colors.blue+'...diagonalization...'+colors.black)
        self.gammaj,self.CjG = np.linalg.eigh(self.H) #;print(red+'Ek',lk,black);print(wk)
        self.invCjG = np.linalg.inv(self.CjG)
        self.solved = True

    def _set_excitation_errors(self,Smax=0.02,hkl=None):
        """ get excitation errors for Sg<Smax
        - Smax : maximum excitation error to be included
        - hkl : list of tuple or nbeams x 3 ndarray - beams to be included (for comparison with other programs)
        """
        K,K0 = self.K,self.k0
        if isinstance(hkl,list) or isinstance(hkl,np.ndarray):
            hkl = np.array(hkl)
            h,k,l = hkl.T
            qx,qy,qz = hkl.dot(self.lat_vec).T
        else:
            (h,k,l),(qx,qy,qz) = self.lattice

        Kx,Ky,Kz = K
        Sw = K0-np.sqrt((Kx-qx)**2+(Ky-qy)**2+(Kz-qz)**2)
        Swa = (K0**2-((Kx+qx)**2+(Ky+qy)**2+(Kz+qz)**2))/(2*K0)
        if Smax:
            idx = abs(Sw)<Smax
            h,k,l = np.array([h[idx],k[idx],l[idx]],dtype=int)
            qx,qy,qz,Sw,Swa = qx[idx],qy[idx],qz[idx],Sw[idx],Swa[idx]
        q = np.linalg.norm(np.array([qx,qy,qz]).T,axis=1)
        d = dict(zip(['h','k','l','qx','qy','qz','q','Sw','Swa'],[h,k,l,qx,qy,qz,q,Sw,Swa]))

        self.Smax = Smax
        self.nbeams = Sw.size
        self.df_G = pd.DataFrame.from_dict(d)
        self.df_G.index = [str(tuple(h)) for h in self.df_G[['h','k','l']].values]
        self.solved = False


    def _set_Vg(self,opt=0):
        if opt:
            q    = self.df_G['q']
            hkl  = self.df_G[['h','k','l']].values
            xi   = self.pattern[:,:3].T
            Za   = np.array(self.pattern[:,-1],dtype=int)
            fj   = scatf.get_fe(Za,q)
            # fj   = 1 #testing exp factor alone
            # Fhkl[nGs] = sum_{natoms} fj [nGs x natoms] * hkl[nGsx3].dot(xi[3,natoms])
            Fhkl = np.sum(fj*np.exp(-2J*np.pi*hkl.dot(xi)),axis=1)
            # Fhkl = np.sum(fj,axis=1)
        else:
            hkl  = self.df_G[['h','k','l']].values
            Fhkl = self.Fhkl.copy()
            V0_idx = np.array([2*self.Nmax]*3)
            Fhkl[tuple(V0_idx)] = 0
            Fhkl = np.array([ Fhkl[tuple(hkl_G+V0_idx)] for hkl_G in hkl])

        pre = 1/np.sqrt(1-cst.keV2v(self.keV)**2)
        Vg_G = Fhkl*pre/(self.crys.volume*np.pi)
        px,py,e0x = mut.project_beams(K=self.K,qxyz=self.get_G(),e0=[1,0,0],v=1)
        self.e0x = e0x
        self.df_G['px'] = px
        self.df_G['py'] = py
        self.df_G['Vg'] = Vg_G
        self.df_G['Vga'] = abs(Vg_G)
        self.df_G['Swl'] = logM(self.df_G['Sw'])
        self.df_G['L']  = np.ones(Vg_G.shape)
        self._set_zones()
        self.df_G.loc[str((0,0,0)),'Vg'] = 0

    def _set_zones(self):
        hkl     = self.df_G[['h','k','l']].values
        # Khkls   = np.round(hkl.dot(self.Kuvw/np.linalg.norm(self.Kuvw))*100)/100
        Khkls    = np.round(hkl.dot(self.Kuvw0)*100)/100 #integer only in zone axis orientation
        Khkl,ar = np.unique(Khkls,return_inverse=True) #; print(Khkls);print(zones);print(ar)
        zones = np.argsort(Khkl)
        self.df_G['zone'] = zones[ar]

    def _set_intensities(self):
        """get beam intensities at thickness"""
        id0 = self.is_hkl([0,0,0],v=0)
        gammaj,CjG = self.gammaj,self.CjG
        S = CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*self.thick))).dot(self.invCjG)
        # S = S[:,id0]
        S = S[id0,:]
        self.df_G['S'] = S
        self.df_G['I'] = np.abs(S)**2

    def _set_kinematic(self):
        Sw,Ug = self.df_G[['Sw','Vg']].values.T
        t,sig = self.thick, self.sig

        #[Ug]=[A-2], [k0]=[A^-1], [t]=[A], [Fhkl]=[fe]=[A]
        Sg = np.pi/self.k0*Ug*t*np.sinc(Sw*t)
        self.df_G['Sg'] = Sg
        self.df_G['Ig'] = np.abs(Sg)**2

    def _set_beams_vs_thickness(self,thicks=None):
        """ get Scattering matrix as function of thickness for all beams
        - thicks : tuple(ti,tf,step) or list or np.ndarray thicknesses
        """
        print(colors.blue+'... beam vs thickness ...'+colors.black)
        if not type(thicks)==type(None):
            self._set_thicks(thicks)

        id0 = self.get_beam(refl=[(0,0,0)],index=True)[0]#;print(id0)
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
        Sz_kin = np.array([np.pi/self.k0*Ug*t*np.sinc(Sw*t) for t in self.z]).T
        self.Iz_kin = np.abs(Sz_kin)**2

    def _set_thicks(self,thicks):
        if isinstance(thicks,tuple):thicks = np.linspace(*thicks)
        elif isinstance(thicks,float) or isinstance(thicks,int) : thicks=[thicks]
        self.z = np.array(thicks)

    def _get_pkl(self,file=None):
        if not file:file=os.path.join(self.path,self.name+'.pkl')
        return file

    def _get_slice(self,s):
        tle = ''
        if isinstance(s,tuple) or isinstance(s,int) or isinstance(s,str):
            if isinstance(s,int):s = 'l=%d' %s
            if isinstance(s,str):
                i,n = s.split('=')
                tle = 'plane %s, ' %s
                n   = 2*self.Nmax+int(n)
                if   i=='h':s = np.s_[n,:,:]
                elif i=='k':s = np.s_[:,n,:]
                elif i=='l':s = np.s_[:,:,n]
        return s,tle

    ################################################################################
    #### getter
    ################################################################################

    def get_beam(self,
        cond:str='',
        refl:Iterable[tuple or str]=[],
        index:bool=True,opt=1,
    ) -> Iterable[int]:
        """ select reflections according to a condition or indices

        Parameters
        ----------
        cond
            condition on the beams
        list-tuple/str refl
            list of beams either as tuple or str
        index
            if True returns indices as accessible through self.df_G.iloc[idx], otherwise through self.df_G.loc[idx]

        Returns
        -------
        Iterable[int] or Iterable[str]
            index of selected reflections


        .. note::
            reflections will first be seleced by applying cond
        """
        #### ensure self.df_G indices are str
        if not isinstance(self.df_G.index[0],str):
            self.df_G.index = [str(tuple(h)) for h in self.df_G[['h','k','l']].values]

        if cond:
            # print(cond)
            refl = list(self.df_G.loc[self.df_G.eval(cond)].index.values)
        if len(list(refl)):
            if not isinstance(refl[0],str):
                refl = [str(tuple(h)) for h in refl]

        #### keep valid reflections only
        refl = [h for h in refl if h in self.df_G.index]

        if index or opt:
            idx = [np.where(self.df_G.index==h)[0][0] for h in refl]
            return idx
        else:
            return refl

    def get_beams_vs_thickness(self,
        dict_opt:bool=False,
        idx:Iterable=slice(0,None,1),
        iZs:Iterable=slice(0,None,1),
    )->np.ndarray:
        """get beams as function of thickness

        Parameters
        ----------
        dict_opt
            returns I(z) as a dictionary
        idx
            beam indices (see meth:~Bloch.get_beam for beam selection)
        iZs
            selected thicknesses
        returns
        -------
        dict or np.ndarray
            beams as function of thickness
        """
        if 'Iz' not in self.__dict__:
            print('warning : thickness calc since Iz not an attribute')
            self._set_beams_vs_thickness()
        Iz = self.Iz[idx,iZs].copy()
        if dict_opt:
            return dict(zip(self.df_G.index[idx], Iz))
        else:
            return Iz


    def get_intensities(self):return self.df_G.I
    def get_hkl(self):return self.df_G[['h','k','l']].values
    def get_kin(self):return self.df_G[['h','k','l','Sw','Vg','Sg','Ig']]
    def get_zones(self):return self.df_G[['h','k','l','zone']].values
    def get_G(self):return self.df_G[['qx','qy','qz']].values
    def get_Xig(self,tol=1e4):
        self.df_G['Xi_g'] = self.k0/abs(self.df_G.Vg)
        return self.df_G.loc[self.df_G.Xi_g<tol,['h','k','l','Sw','Vg','Xi_g']]
    def get_Sw(self,Smax=1e-2):
        Smax = min(Smax,self.Smax)
        print(self.df_G[['h','k','l','Sw']].loc[self.df_G.Sw<=Smax])

    def get_Istrong(self,opts:str='p',
        sort:str='I',cols:Iterable[str]=[],
    ):
        """Get strong beams

        Parameters
        ----------
        cols
            columns to display in addition to I,Sw Vg
        sort
            sort along
        opts
            p(print values) o(return hkl index)


        .. note ::
            condition used : '(Swl<2) & (I>10) & (Vga>0.1)'
        """
        hkl = get_beam(cond='(Sw<1e-2) & (I>10) & (Vga>0.1)',index=False)

        if 'p' in opts:
            print(self.df_G[hkl].sort_values(sort,ascending=False)[['I','Sw','Vg']+cols])
        if 'o' in opts:
            return hkl

    def is_hkl(self,Ghkl,v=1):
        Gidx = np.where(np.linalg.norm(self.get_hkl() - Ghkl,axis=1)==0)[0]
        if v:print(self.df_G.iloc[Gidx])
        return Gidx[0]

    ################################################################################
    #### display
    ################################################################################
    def show_beams(self,opts='SVkr',cond='',refl=[],name='',**kwargs):
        """ display beams"""
        # cond = self.get_cond(cond)
        idx = self.get_beam(cond=cond,refl=refl)#;print(idx)
        if not name : name=os.path.join(self.path,self.name+'_setup.svg')
        EDdisp.show_frame(df_bloch=self.df_G,single_mode=False,
            opts=opts,hkl_idx=idx,name=name,
            **kwargs)

    def show_beams_vs_thickness(self,
        thicks:Optional[Sequence]=None,
        beam_args:dict=None,
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
        beam_args
            To pass to get_beam
        kwargs
            to pass to :meth:~EDutil.postprocess.plot_beam_thickness
        """

        if isinstance(thicks,tuple):
            self._set_beams_vs_thickness(thicks)

        ### beam selection
        if isinstance(beam_args,dict):
            beam_args.update({'index':True})
            idx = self.get_beam(**beam_args)
        else:
            idx = self.get_beam(cond='(Sw<1e-2) & (Vga>0)',index=True)

        z   = self.z[iZs]
        Iz  = self.get_beams_vs_thickness(idx=idx,iZs=iZs,dict_opt=False)
        hkl = self.df_G.index[idx]

        beams=[hkl,z,None,None,Iz]
        return pp.plot_beam_thickness(beams,**kwargs)


    def show_Idyn_vs_Ikin(self,thicks=None,iZs=10,cmap='Spectral',**kwargs):
        """Display correlation plot of dynamical vs kinematical intensities

        Parameters
        -----------
        thicks
            thicknesses
        hkl
            beams to include
        """
        if thicks:self._set_beams_vs_thickness(thicks)
        if isinstance(iZs,int):iZs=slice(0,None,iZs)
        # iB = self.get_Istrong(out=1)[0]#;print(iB,self.Iz[    iB,iZs],)
        iB = np.argsort(np.sum(self.Iz,axis=1))[-2]
        z      = self.z.copy()[       iZs]
        Iz_dyn = self.Iz.copy()[    :,iZs]/self.Iz[    iB,iZs]
        Iz_kin = self.Iz_kin.copy()[:,iZs]/self.Iz_kin[iB,iZs]#;print(Iz_dyn[Iz_dyn>1e-3])
        nzs = z.size
        cs = dsp.getCs(cmap,nzs) #; print(len(cs),Iz_dyn.shape,Iz_kin.shape)
        # scat=tuple( [[I_kin,I_dyn,cs[i]] for i,(I_dyn,I_kin) in enumerate(zip(self.Iz[:,::iZ],self.Iz_kin[:,::iZ]))])
        plts=[[I_kin,I_dyn,[cs[i],'o'],'$z=%d\AA$' %z] for i,(z,I_dyn,I_kin) in enumerate(zip(z,Iz_dyn.T,Iz_kin.T))]
        plts+=[ [[0,1],[0,1],[(0.5,)*3,'--'],''] ]
        dsp.stddisp(plts,labs=['$I_{kin}$','$I_{dyn}$'],sargs={'alpha':0.5},**kwargs)

    def show_Fhkl(self,s=None,opts='m',h3D=0,**kwargs):
        """Displays structure factor over grid

        Parameters
        -----------
        opts
            see get_fz
        s
            slice or str('k=0' => Fhkl(k=0)) or int('l==<int>')
        """
        fz,fz_str = get_fz(opts)
        s,s_str = self._get_slice(s)
        tle = 'Structure factor($\AA$), showing %s in %s'  %(fz_str,s_str)
        h,k,l = self.hklF
        Fhkl  = self.Fhkl #.copy()
        if isinstance(s,tuple):
            Fhkl   = Fhkl[s]
            nx,ny  = np.array((np.array(Fhkl.shape)-1)/2,dtype=int)
            i,j    = np.meshgrid(np.arange(-nx,nx+1),np.arange(-ny,ny+1))
            fig,ax = dsp.stddisp(scat=[i,j,fz(Fhkl)],title=tle,**kwargs)
        else:
            fig,ax = dsp.stddisp(scat=[h,k,l,fz(self.Fhkl)],title=tle,rc='3d',**kwargs)
            if h3D:h3d.handler_3d(fig,persp=False)

    def show_H(self,**kwargs):
        dsp.stddisp(im=[np.abs(self.H)],title='abs(H)', pOpt='im')


    ################################################################################
    #### misc
    ################################################################################

    def save(self,file=None,v=1):
        """save this object"""
        file = self._get_pkl(file)
        with open(file,'wb') as out :
            pickle5.dump(self, out, pickle5.HIGHEST_PROTOCOL)
        if v:print(colors.green+"object saved\n"+colors.yellow+file+colors.black)

    def convert2tiff(self,tiff_file:str=None,
        Imax:int=3e4,
        Nmax:int=512,aperpixel:Optional[float]=None,
        gs3:float=0.1,
        rmax:int=0,
        thick:float=None,
        show:bool=False,**kwargs,
    ):
        """Write intensities to a tiff file

        Parameters
        -----------
        tiff_file
            name of the file (default name if not provided)
        Imax
            max value for the intensities
        Nmax
            image resolution in pixel
        aperpixel
            reciprocal size for each pixel (aperture per pixel)
        gs3
            Gaussian broadening for each reflections
        rmax
            radius for the noise to be added
        show
            Show produced image
        kwargs
            args to be passed to the tiff viewer


        .. note ::
            The reciprocal space resolution is aperpixel*Nmax. As a result some of teh
            high resolution reflections might not appear if aperpixel or Nmax are too small.
            If aperpixel is not specified, it is automatically so it contains the image
            will contain all reflections.
        """
        Nmax=Nmax//2
        if thick:
            self.set_thickness(thick)
        px,py,I = self.df_G[['px','py','I']].values.T
        if not aperpixel:
            aperpixel = 1.1*max(px.max(),py.max())/Nmax
            print('aperpixel set to %.1E A^-1' %aperpixel)
        dqx,dqy = [aperpixel]*2

        #convert to pixel locations
        nx,ny = int(np.floor(gs3/dqx)),int(np.floor(gs3/dqy)) ;
        ix,iy = np.meshgrid(range(-nx,nx+1),range(-ny,ny+1))
        x,y = ix*dqx,iy*dqy
        i,j = np.array([np.round(px/dqx),np.round(py/dqy)],dtype=int)+Nmax

        #### Gaussian broadening
        Pb = np.exp(-(x**2+y**2)/(gs3/3)**2)#;dsp.stddisp(im=[x,y,Pb],pOpt='im')
        im0 = np.zeros((2*Nmax,)*2)
        for i0,j0,I0 in zip(i,j,I):
            i0x,j0y = i0+ix,j0+iy
            idx     = (i0x>=0) & (j0y>=0)  & (i0x<2*Nmax) & (j0y<2*Nmax)
            i0x,j0y = i0+ix[idx],j0+iy[idx]
            im0[i0x,j0y] += Pb[idx]/Pb[idx].sum()*I0

        if rmax:
            h,k = np.meshgrid(range(-Nmax,Nmax),range(-Nmax,Nmax))
            r = np.sqrt(h**2+k**2);r[r==0]=1
            im0 += rmax*np.random.rand(*im0.shape)#/(rmax*r)

        if not tiff_file:
            tiff_file = os.path.join(self.path,self.name+'.tiff')
        I = np.array(im0*Imax,dtype='uint16')
        tifffile.imwrite(tiff_file,np.flipud(I))
        print(colors.yellow+tiff_file+colors.green+' saved'+colors.black)
        if show:
            viewers.Viewer(self.path,frame=1,thick=self.thick,**kwargs)
            # dsp.stddisp(im=[I],cmap='')


    ###################################################################################
    #### convergence stuff
    ###################################################################################
    def convergence_test(self,
        Nmax:Iterable[int],
        Smax:Iterable[float],
        z=(0,100,100),
        hkl=None,
        opts:str='Ss'
    ):
        """creates a convergence test for a range of Nmax and Smax
        and displays it as funtion of thickness for beams hkl

        Parameters
        ----------
        opt:  S(solve) s(save)
        """
        solve,save = 'S' in opts,'s' in opts
        Nmax,Smax=[a.flatten() for a in np.meshgrid(Nmax,Smax)]
        z=np.array(z)
        df=pd.DataFrame(columns=['Nmax','Smax','nbeams','Iz'])
        # simulate
        self.update_Nmax(Nmax[0])
        self._set_excitation_errors(Smax=Smax[0])
        hkl0 = self.df_G.index
        if solve:
            for i,(Nm,Sm) in enumerate(zip(Nmax,Smax)):
                 self.solve(Nmax=Nm,Smax=Sm,thicks=z)
                 idx='(%d,%.3f)'%(Nm,Sm)
                 idx_b = self.get_beam_index(refl=hkl0)
                 df.loc[idx,['Nmax','Smax','nbeams']]=[Nm,Sm,self.nbeams]
                 df.loc[idx,'Iz']=self.Iz[idx_b,:].copy()
            self.df=df
            if save:self.save(self.path+self.name+'_cv.pkl')

            self.show_convergence(hkl)
        else:
            for i,(Nm,Sm) in enumerate(zip(Nmax,Smax)):
                self.update_Nmax(Nm)
                self._set_excitation_errors(Smax=Sm)
                idx='(%d,%.3f)'%(Nm,Sm)
                df.loc[idx,['Nmax','Smax','nbeams']]=[Nm,Sm,self.nbeams]
            print(df)

    def show_convergence(self,hkl):
        #beam selection
        valid_hkl=0
        if isinstance(hkl,list):valid_hkl=isinstance(hkl[0],str)
        self.update_Nmax(int(self.df.Nmax.min()))
        self._set_excitation_errors(Smax=self.df.Smax.min())
        if not valid_hkl:
            hkl = self.get_beam(cond='Sw<1e-2' ,opt=0)[0:2]
            hkl += [self.get_beam(cond='Sw>=1e-2',opt=0)[0]]
            print('default hkl : ',hkl)
        Icols=hkl
        idx_s = self.get_beam_index(refl=hkl)
        self.df[Icols]=None
        #plot
        cs=dsp.getCs('RdYlGn',self.df.shape[0])
        plts=[]
        ms = dsp.markers[:len(hkl)]
        for idx,m in zip(idx_s,ms):
            plts+=[[self.z,r.Iz[idx],[c,m+'-' ],
                 ] for c,(l,r) in zip(cs,self.df.iterrows())]
        legElt={'$%s$' %h:'k-'+ms[i] for i,h in enumerate(hkl)}
        legElt.update({'$%s$' %l:[c,'-'] for c,l in zip(cs,self.df.index)})
        dsp.stddisp(plts,labs=['$z$','$I$'],legElt=legElt,lw=2)
