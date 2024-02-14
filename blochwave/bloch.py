"""Bloch wave solver"""
import importlib as imp
import numpy as np,pandas as pd,pickle5,os,glob,tifffile,mrcfile
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from subprocess import check_output,Popen,PIPE
from crystals import Crystal
from utils import displayStandards as dsp           #;imp.reload(dsp)
from utils import physicsConstants as cst           #;imp.reload(cst)
from utils import glob_colors as colors,handler3D as h3d
from multislice import postprocess as pp            #;imp.reload(pp)
from scattering import structure_factor as sf       ;imp.reload(sf)
from scattering import scattering_factors as scatf  #;imp.reload(scatf)
from EDutils import viewers                         #;imp.reload(viewers)
from EDutils import utilities as ut                 #;imp.reload(ut)
from EDutils import pets as pt                      ;imp.reload(pt)
from EDutils import display as EDdisp               ;imp.reload(EDdisp)
from . import util as bloch_util                    ;imp.reload(bloch_util)
felix='%s/bin/felix' %os.path.dirname(__file__)

class Bloch:
    """
    Bloch wave simulation class

    Parameters
    ---------
    file
        structure file (3D:cif_file or wallpaper info in 2D)
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
    felix
        Use felix solver
    nbeams
        Number of beams to use (using felix only)
    v
        verbose
    kwargs
        arguments to be passed to :meth:`~Bloch.solve`
    """
    def __init__(self,
        cif_file:str,name:str='',path:str='',
        beam:Optional[dict]={},keV:float=200,u:Sequence[float]=[0,0,1],
        Nmax:int=1,dmin:int=None,
        Smax:float=0.2,
        solve:bool=True,init=True,
        felix:bool=False,nbeams:int=200,
        eps:float=1,f_sw=None,frame=None,
        v = True,
        **kwargs,
    ):
        self.solved = False
        self.Nmax   = 0
        self.dmin   = None
        self.thick  = 100
        self.thicks = self._set_thicks((0,1000,1000))
        self.eps=eps
        self.frame=frame

        self._set_structure(cif_file)
        beam_args={'keV':keV,'u':u}
        beam_args.update(beam)
        self.set_beam(**beam_args)
        self.set_name(name,path)
        self.nbeams=nbeams

        if solve :
            self.solve(Smax=Smax,Nmax=Nmax,dmin=dmin,f_sw=f_sw,**kwargs)
        else:
            if not felix and init:
                print(colors.blue+'...Nmax... '+colors.black)
                self.update_Nmax(Nmax,dmin)
                print(colors.blue+'...Excitation errors... '+colors.black)
                self._set_excitation_errors(Smax,f_sw=f_sw)
                print(colors.blue+'...Vg... '+colors.black)
                self._set_Vg()
        self.save(v=v)

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
        if not os.path.exists(self.path) and self.path:
            cmd='mkdir -p %s' %self.path
            p=check_output(cmd,shell=True).decode();print(p)

        if self.pdb_file:
            check_output('mv %s %s' %(self.cif_file,self.path),shell=True)

    def update_Nmax(self,Nmax:int,dmin:int=None):
        """
        Update resolution/max order (lattice and Fhkl)

        Parameters
        ----------
        Nmax
            maximum h,k,l order
        dmin
            minimum resolution
        """
        gemmi=not isinstance(dmin,type(None)) and self.pdb_file
        if Nmax:
            b1,b2,b3  = self.lat_vec
            gmax = Nmax*np.linalg.norm(b1+b2+b3) #A^-1
            dmin = max(0.1,min(3,gmax))
        self.dmin=dmin
        self.Nmax=Nmax

        if not os.path.exists(self.Fhkl_file()):
            if gemmi:
                print(colors.blue+'...gemmi structure factors...'+colors.black)
                hklF,Fhkl = ut.gemmi_sf(self.pdb_file,dmin)

                Nmax=np.array(Fhkl.shape[0])//4#.min()
                self.Nmax=Nmax
                print('Nmax gemmi:',Nmax)

            else:
                print(colors.blue+'...Structure factors...'+colors.black)
                idx = [i for i,x in enumerate(self.pattern[:,:3]) if all(x<0.99)] #should not be necessary
                self.pattern=self.pattern[idx,:]
                hklF,Fhkl = sf.structure_factor3D(self.pattern,
                    2*np.pi*self.lat_vec,hklMax=2*Nmax)
                df_Fhkl = sf.get_structure_factor(self.cif_file,hklMax=2*Nmax)
                df_Fhkl['Fga']  = np.real(np.abs(df_Fhkl.F))
                df_Fhkl['Uga']  = np.real(np.abs(df_Fhkl.F*cst.meff(self.keV)/(np.pi*self.crys.volume)))
                df_Fhkl['xi_g'] = self.k0/df_Fhkl.Uga
                df_Fhkl.to_pickle(self.get_Fhkl_pkl())
            #save
            np.save(self.Fhkl_file(),Fhkl)
            np.save(os.path.join(self.path,'hklF.npy'),hklF)
            print(colors.green+'structure factors updated.'+colors.black)

    def get_Fhkl_pkl(self):
        return self.Fhkl_file().replace('npy','pkl')

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

    def set_thickness(self,thick:float=100,v=True):
        """set thickness and update beams

        Parameters
        ----------
        thick : thickness in A
        """
        # print(type(thick))
        if type(thick) in [int,float] :self.thick=thick
        self._set_kinematic()
        if self.solved:
            self._set_intensities()
            if v: print('updated intensities')


    def solve(self,
        Smax:Optional[float]=None,hkl:Optional[Iterable[int]]=None,
        Nmax:Optional[int]=None,dmin:Optional[int]=None,
        beam:Optional[dict]={},dyngo_args:Optional[dict]={},
        thick:float=None,thicks:Sequence[float]=None,
        opts:str='sv0',
        felix:bool=False,
        nbeams:int=None,f_sw=None,
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
        dyngo_args
            dict contain dyngo type informations. If not empty, it solves  with Dyngo formulation which is slightly different.
        thick
            thickness of crystal (can be modified without resolving)
        thicks
            range of thickness [z_min,z_max,z_step]
        opts
            s(save) t(show intensities) z(show beam vs thickness) 0(Vopt0) v(verbose) H(show H)
        felix
            use felix solver
        nbeams
            Number of beams when using felix solver

        .. note ::

            If hkl is specified, Smax is not taken into account
        """
        if felix:
            self._prepare_Felix(nbeams=nbeams)
            self._run_felix() #wait='w' in opts)
            self._postprocess_felix(show_log='v' in opts)
            self._set_excitation_errors(hkl=self.df_G[['h','k','l']].values,felix=True)
            self._set_Vg(felix=True)
        else:
            if Nmax or dmin:self.update_Nmax(Nmax,dmin)
            if beam : self.set_beam(**beam)
            if Smax or isinstance(hkl,np.ndarray) :
                self._set_excitation_errors(Smax=Smax,hkl=hkl,f_sw=f_sw)
                self._set_Vg()
            self._solve_Bloch(show_H='H' in opts,Vopt0='0' in opts,v='v' in opts,
                dyngo_args=dyngo_args)
        ##### postprocess
        if thick or 't' in opts:
            self.set_thickness(thick)
        if not type(thicks)==type(None) and 'z' in opts:
            self._set_beams_vs_thickness(thicks)#;print('thicks solved')
        if 's' in opts:
            self.save(v='v' in opts)

    def update(self,keV=None,u=None,Smax=None,Nmax=None,dmin=None,
            gemmi=False,hkl=None,
            f_sw=None,
            **kwargs):
        if not gemmi:dmin=None
        self.set_beam(keV=keV,u=u)
        if Nmax or dmin :self.update_Nmax(Nmax,dmin)
        if Smax or isinstance(hkl,np.ndarray):
            self._set_excitation_errors(Smax,hkl,f_sw=f_sw)
            self._set_Vg()
        self.save()

    ################################################################################
    #### private
    ################################################################################
    def _prepare_Felix(self,npx=20,nbeams=None,thicks=(0,0,20),show_log=False):
        """"""
        if not nbeams:nbeams=self.nbeams
        print(colors.blue+"preparing simulation"+colors.black)
        cmd="""cd %s;
        if [ ! -d felix ];then mkdir felix;fi;
        """ %(self.path) #,os.path.realpath(felix_cif))
        # cp %s felix/felix.cif;
        # mv felix.inp felix;
        # p=Popen(cmd,shell=True)#;p.wait();o,e = p.communicate();if e:print(e)
        p=check_output(cmd,shell=True).decode();print(p)
        bloch_util.get_inp(npx,nbeams,self.u,self.keV,thicks,out=os.path.join(self.path,'felix','felix.inp'))
        ut.crys2felix(self.crys,opt='w',out=os.path.join(self.path,'felix','felix.cif'))

    def _run_felix(self,wait=True):
        if not os.path.exists(felix):
            print('Error : Running with Felix not available. You need to install felix in %s :' )
            print(felix)
            return

        print(colors.blue+"... running felix ..."+colors.black)
        cmd="""cd %s;
        cd felix;
        %s > felix.log 2>&1;
        """ %(self.path,felix)#os.path.dirname(__file__))
        p=Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE);
        if wait:
            p.wait();print(p.communicate())
        else:
            return p


    def _postprocess_felix(self,show_log=False):
        if show_log:
            print(colors.blue+"felix output"+colors.black)
            with open('%s/felix/felix.log' %self.path,'r') as f:print('\n'.join(f.readlines()))

        eigvals = os.path.join(self.path,'felix/eigenvals.txt')
        check_output("sed -i -E 's/ {1,}/,/g' %s;sed -i -E 's/^,//' %s" %(eigvals,eigvals) ,shell=True)
        df = pd.read_csv(eigvals,names=['h','k','l','gr','gi'])
        self.gammaj = (df.gr+1J*df.gi).values
        # g = np.loadtxt(os.path.join(self.path,'felix/eigenvals.txt'))
        # self.gammaj = g[:,3::2]+1J*g[:,4::2];g=g[:,0]
        C = np.loadtxt(os.path.join(self.path,'felix/eigenvec.txt'))
        self.CjG = C[:,3::2]+1J*C[:,4::2]
        self.invCjG=np.conj(self.CjG.T)

        self.df_G = df[['h','k','l']]
        self.df_G.index = [str(tuple(h)) for h in self.df_G.values]
        self.solved=True

    def _solve_Bloch(self,show_H=False,Vopt0=True,dyngo_args={},v=False):
        ''' Diagonalize the Hamiltonian'''
        # Ug is a (2*Nmax+1)^3 tensor :
        # Ug[l] = [U(-N,-N,l) .. U(-N,0,l) .. U(-N,N,l)
        #          U( 0,-N,l) .. U( 0,0,l) .. U( 0,N,l)
        #          U(-N,-N,l) .. U( N,0,l) .. U( N,N,l)]
        self.dyngo = any(dyngo_args)

        hkl = self.df_G[['h','k','l']].values
        Sg  = self.df_G.Sw.values
        if self.dyngo:
            Rmat=dyngo_args['Rmat']
            self.scale=dyngo_args['scale']
            surf_norm = Rmat.dot([0,0,1])
            gn = Rmat.dot(self.df_G[['qx','qy','qz']].T).T.dot(surf_norm)
            Knorm = self.k0*surf_norm[-1]
            self.Knorm=Knorm
            sqrtkg=np.sqrt(1+gn/Knorm)  #dyngo implementation
            Sg*=2*self.k0/sqrtkg

        pre = 1/np.sqrt(1-cst.keV2v(self.keV)**2)
        Ug = pre*self.get_Fhkl()/(self.crys.volume*np.pi)*self.eps #/3
        #####################
        # setting average potential to 0
        # Ug[U0_idx] = U0
        # and Ug(iG,jG) are obtained from Ug[h,k,l] where h,k,l = hlk_iG-hkl_jG
        #####################
        U0_idx = [2*self.Nmax]*3
        Ug[tuple(U0_idx)] = 0
        # if Vopt0 :Ug[tuple(U0_idx)] = 0   #setting average potential to 0
        if v:
            msg = '...assembling {N}x{N} matrix (structure factor shape : {Ug}) ...\
            '.format(N=Sg.shape[0],Ug=Ug.shape)
            print(colors.blue,msg,colors.black)
        H = np.diag(Sg+0J)
        for iG,hkl_G in enumerate(hkl) :
            U_iG = np.array([Ug[tuple(hkl_J+U0_idx)] for hkl_J in hkl_G-hkl]) #;print(V_iG.shape)
            # print('G : ',hkl_G)
            # print('U_G:',Ug[tuple(U0_idx+hkl_G)])
            # print('idx iG:',[tuple(hkl_J+U0_idx) for hkl_J in hkl_G-hkl])
            # print(U_iG)
            if self.dyngo:
                # qh = np.linalg.norm(self.lat_vec.T.dot((hkl_G-hkl).T),axis=0)
                qh = Rmat.dot((hkl_G-hkl).T).T.dot(surf_norm)
                sqrtkh = np.sqrt(1+qh/Knorm)
                H[iG,:] += U_iG/(sqrtkg[iG]*sqrtkh)  #so dyngo implementation
            else:
                H[iG,:] += U_iG/(2*self.k0)  #off diagonal terms as potential

        if self.dyngo:
            H/=(2*self.k0)

        H *= 2*np.pi #to get same as felix
        self.H=H
        if show_H:self.show_H()

        if v:print(colors.blue+'...diagonalization...'+colors.black)
        self.gammaj,self.CjG = np.linalg.eigh(self.H)
        # print(self.gammaj)
        self.invCjG = np.linalg.inv(self.CjG)
        self.solved = True

    def _set_structure(self,cif_file):
        pdb_file = ''
        if cif_file[-3:]=='pdb':
            pdb_file=cif_file
            cif_file = ut.pdb2npy(os.path.basename(cif_file[:-4]))

        self.cif_file = cif_file
        self.pdb_file = pdb_file
        self.crys     = ut.import_crys(self.cif_file)
        self.lat_vec0 = np.array(self.crys.lattice_vectors)
        self.lat_vec  = np.array(self.crys.reciprocal_vectors)/(2*np.pi)
        self.pattern  = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in self.crys.atoms] )

    def _set_excitation_errors(self,Smax=0.02,hkl=None,felix=False,f_sw=None):
        """ get excitation errors for Sg<Smax
        - Smax : maximum excitation error to be included
        - hkl : list of tuple or nbeams x 3 ndarray - beams to be included (for comparison with other programs)
            Note that the Smax filter will still be applied to these beams(use Smax=0 to force all beams)
        - felix : to compare with felix
        - f_sw : optional function to compute Sw f(hkl,frame) for
        """
        # print(colors.blue+'... setting excitation error ... '+colors.black)
        K,K0 = self.K,self.k0
        if isinstance(hkl,list) or isinstance(hkl,np.ndarray):
            hkl = np.array(hkl)
            h,k,l = hkl.T
            qx,qy,qz = hkl.dot(self.lat_vec).T
        else:
            (h,k,l),(qx,qy,qz) = self.get_lattice()
            # print(hkl.shape)
            # h,k,l = hkl
            hkl = np.array([h,k,l])
        if f_sw:
            args  = {'hkl':hkl.T,'frame':self.frame}
            Sw = f_sw(**args)
        else:
            Kx,Ky,Kz = K
            Sw = (K0**2-((Kx+qx)**2+(Ky+qy)**2+(Kz+qz)**2))/(2*K0)
        q = np.linalg.norm(np.array([qx,qy,qz]).T,axis=1)
        if felix:
            self.df_G[['qx','qy','qz','q','Sw','Swa']] = np.array([qx,qy,qz,q,Sw,abs(Sw)]).T
            self.df_G['I'] = 0
            self.df_G.loc[str((0,0,0)),'I'] = 1
            return

        # print('ok')
        if Smax:
            idx = abs(Sw)<Smax
            h,k,l = np.array([h[idx],k[idx],l[idx]],dtype=int)
            qx,qy,qz,q,Sw = qx[idx],qy[idx],qz[idx],q[idx],Sw[idx]#,Swa[idx]
            # Gz,GzG = Gz[idx],GzG[idx]
        d = dict(zip(['h','k','l','qx','qy','qz','q','Sw','Swa'],[h,k,l,qx,qy,qz,q,Sw,abs(Sw)]))


        self.Smax = Smax
        self.nbeams = Sw.size
        self.df_G = pd.DataFrame.from_dict(d)
        self.df_G.index = [str(tuple(h)) for h in self.df_G[['h','k','l']].values]
        self.solved = False
        self.df_G['I'] = 0
        self.df_G.loc[str((0,0,0)),'I'] = 1
        # print(self.df_G.iloc[60])


    def _set_Vg(self,felix=0):
        ##### opt was used for testing against FELIX
        if felix:
            q    = self.df_G['q']
            hkl  = self.df_G[['h','k','l']].values
            idx = [i for i,x in enumerate(self.pattern[:,:3]) if all(x<0.99)]
            xi   = self.pattern[idx,:3].T
            # print(xi)
            Za   = np.array(self.pattern[idx,-1],dtype=int)
            fj   = scatf.get_fe(Za,q) #;fj = 1 #testing exp factor alone
            #### Fhkl[nGs] = sum_{natoms} fj [nGs x natoms] * hkl[nGsx3].dot(xi[3,natoms])
            Fhkl = np.sum(fj*np.exp(-2J*np.pi*hkl.dot(xi)),axis=1)
            # print(self.crys,xi)
        else:
            hkl  = self.df_G[['h','k','l']].values
            Fhkl = self.get_Fhkl()
            V0_idx = np.array([2*self.Nmax]*3)
            Fhkl[tuple(V0_idx)] = 0
            # print(Fhkl.shape,hkl.max())
            Fhkl = np.array([ Fhkl[tuple(hkl_G+V0_idx)] for hkl_G in hkl])

        self.pre = 1/np.sqrt(1-cst.keV2v(self.keV)**2)
        # Vg_G = Fhkl/(self.crys.volume*np.pi)*self.pre*self.eps
        Vg_G = Fhkl
        px,py,e0x = ut.project_beams(K=self.K,qxyz=self.get_G(),e0=[1,0,0],v=1)
        self.e0x = e0x
        self.df_G['px'] = px
        self.df_G['py'] = py
        self.df_G['Fg' ]  = Fhkl
        self.df_G['Fg2']  = np.abs(self.df_G.Fg)**2
        self.df_G['Vg']   = np.abs(self.df_G.Fg)/self.crys.volume
        self.df_G['Ug']   = self.pre*self.df_G.Fg/(self.crys.volume*np.pi)*self.eps
        self.df_G['Uga']  = np.abs(self.df_G.Ug)
        self.df_G['xi_g'] = self.k0/self.df_G.Uga
        self.df_G['Swl']  = bloch_util.logM(self.df_G['Sw'])
        self.df_G['L']    = np.ones(Vg_G.shape)

        self._set_zones()
        self.df_G.loc[str((0,0,0)),'Vg'] = 0
        # print('set_Vg : Vg=',self.df_G['Vg'][80])

    def _set_zones(self):
        hkl     = self.df_G[['h','k','l']].values
        # Khkls   = np.round(hkl.dot(self.Kuvw/np.linalg.norm(self.Kuvw))*100)/100
        Khkls    = np.round(hkl.dot(self.Kuvw0)*100)/100 #integer only in zone axis orientation
        Khkl,ar = np.unique(Khkls,return_inverse=True) #; print(Khkls);print(zones);print(ar)
        zones = np.argsort(Khkl)
        self.df_G['zone'] = zones[ar]

    def _get_central_beam(self):
        return self.get_beam(refl=[str((0,0,0))],index=True)[0]

    def _set_intensities(self):
        """get beam intensities at thickness"""
        id0 = self._get_central_beam()
        gammaj,CjG = self.gammaj,self.CjG
        if self.dyngo:
            S = CjG.dot(np.diag(np.exp(1J*gammaj*self.thick*self.k0/self.Knorm))).dot(self.invCjG)
        else:
            S = CjG.dot(np.diag(np.exp(1J*gammaj*self.thick))).dot(self.invCjG)
        # S = S[:,id0]
        S = S[id0,:]
        self.df_G['S'] = S
        self.df_G['I'] = np.abs(S)**2
        # print(CjG[0]);print(self.df_G.sort_values('I')[['Sw','Ig','I']])
        if self.dyngo:
            self.df_G['I']*=self.scale**2
        # print(self.df_G['I'])

    def _set_I(self,iZ=-1):
        idx=self.get_beam(refl=self.df_G.index)
        self.thick=self.z[iZ]
        self.df_G.I=self.Iz[idx,iZ]

    def _set_kinematic(self):
        Sw,xi_g = self.df_G[['Sw','xi_g']].values.T
        t,sig = self.thick, self.sig

        #[Ug]=[A-2], [k0]=[A^-1], [t]=[A], [Fhkl]=3[fe]=[A]
        # print(Ug[0],t,Sw[0])
        Sg = np.pi*t/xi_g*np.sinc(Sw*t)
        self.df_G['Skin'] = Sg
        self.df_G['Ikin'] = np.abs(Sg)**2

    def _set_beams_vs_thickness(self,thicks=None,v=True):
        """ get Scattering matrix as function of thickness for all beams
        - thicks : tuple(ti,tf,step) or list or np.ndarray thicknesses
        """
        if v:print(colors.blue+'... beam vs thickness ...'+colors.black)
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
        Sw,xi_g = self.df_G[['Sw','xi_g']].values.T
        Sz_kin = np.array([np.pi*t/xi_g*np.sinc(Sw*t) for t in self.z]).T
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
        index:bool=True,
    ) -> Iterable[int]:
        """ select reflections according to a condition or indices.
        Returns strong beams if nothing is specified

        Parameters
        ----------
        cond
            condition on the beams (can be a function see meth:~`bloch_pp.strong_beams` for example)
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

        if not cond and not any(refl):cond={'tol':1e-3}
        if cond:
            if isinstance(cond,dict):
                cond_args=cond.copy()
                if not 'I' in self.df_G.columns:
                    if not 'Iz' in self.__dict__:
                        self._set_beams_vs_thickness(thicks=(0,1000,3))
                    self.df_G.I=self.Iz[:,-1]

                cond=lambda dfG:bloch_util.strong_beams(dfG,**cond_args)
            if isinstance(cond,str):
                refl = list(self.df_G.loc[self.df_G.eval(cond)].index.values)
            else:
                refl=cond(self.df_G)

        if len(list(refl)):
            if not isinstance(refl[0],str):
                refl = [str(tuple(h)) for h in refl]

        #### keep valid reflections only
        refl = [h for h in refl if h in self.df_G.index]

        if index:
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
            beam indices in dataFrame (see meth:~`Bloch.get_beam` for selected beams from their miller indices)
        iZs
            selected thicknesses
        returns
        -------
        dict or np.ndarray
            beams as function of thickness
        """
        if 'Iz' not in self.__dict__:
            self.thicks=(0,1000,1000)
            print('warning : computing Iz for thickness ', self.thicks)
            self._set_beams_vs_thickness()
        Iz = self.Iz[idx,iZs].copy()
        if dict_opt:
            return dict(zip(self.df_G.index[idx], Iz))
        else:
            return Iz


    def get_lattice(self):
        # U0_idx = [2*self.Nmax]*3
        # lat=np.load(os.path.join(self.path,'lattice.npy'))
        # print(lat)#self.Nmax,lat.shape)
        # if lat.shape[0]>4*self.Nmax:
        #     print(self.Nmax)
        #     snmax=slice(U0_idx-2*Nmax,U0_idx-2+Nmax)
        #     lat=lat[snmax,snmax,snmax]
        #     print(lat)
        lat = ut.get_lattice(self.lat_vec,self.Nmax)
        return lat
    def get_Fhkl(self):
        return np.load(self.Fhkl_file())
    def Fhkl_file(self):
        return os.path.join(self.path,'Fhkl_%d.npy' %self.Nmax)
        # return os.path.join(self.path,'Fhkl_%.2f.npy' %self.dmin)

    def get_intensities(self):return self.df_G.I
    def get_hkl(self):return self.df_G[['h','k','l']].values
    def get_kin(self):return self.df_G[['h','k','l','Sw','Vg','Sg','Ig']]
    def get_zones(self):return self.df_G[['h','k','l','zone']].values
    def get_G(self):return self.df_G[['qx','qy','qz']].values
    def get_Xig(self,tol=1e4):
        self.df_G['Xi_g'] = self.k0/abs(self.df_G.Vg)
        xig = self.df_G.loc[self.df_G.Xi_g<tol,['h','k','l','Sw','Vg','Xi_g']]
        print(xig)
        return xig
    def get_Sw(self,Smax=1e-2):
        Smax = min(Smax,self.Smax)
        Smax = self.df_G[['h','k','l','Sw']].loc[self.df_G.Sw<=Smax]
        print(Smax)
        return Smax

    def is_hkl(self,h,k,l):
        """ tests whether beam is part of the simulation
        Parameters
        -----------
            h,k,l
                indices

        Returns
        ---------
            True or False
        """
        return any(self.get_beam(refl=[str((h,k,l))]))

    ################################################################################
    #### display
    ################################################################################
    def show_df_G(self,n=-1,hkl=[],cols=['Sw','Fg','Fg2','Vg','Uga','xi_g'],sort='Swa'):
        formats = {'Sw':'{:>7.1e}','Swa':'{:>7.1e}',
            'Fg':'{:>6.1f}','Fg2':'{:>8.1f}',
            'Vg':'{:>6.2e}','Uga':'{:>8.1e}','xi_g':'{:>5.0f}',
            }
        df = self.df_G
        if any(hkl):
            df=df.loc[[str(tuple(h)) for h in hkl]]
        df=df.sort_values(sort)
        if type(n)==int:
            df=df[:n]
        print(df[cols].to_string(
                formatters={k: v.format for k, v in formats.items()}))

    def show_beams(self,opts='BVSkr',cond='',refl=[],name='',**kwargs):
        """ display beams"""
        # cond = self.get_cond(cond)
        hkl_idx = self.get_beam(cond=cond,refl=refl,index=False)#;print(idx)
        # hkl_idx = self.df_G.index
        # self.df_G.Ig=self.Iz[:,-1]
        self.df_G.loc[str((0,0,0)),'I']=0
        if not name : name=os.path.join(self.path,self.name+'_setup.svg')
        return EDdisp.show_frame(df_bloch=self.df_G,single_mode=False,
            opts=opts,hkl_idx=hkl_idx,name=name,
            **kwargs)

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

        if isinstance(thicks,tuple):
            self._set_beams_vs_thickness(thicks)

        ### beam selection
        if isinstance(refl,list):
            if not isinstance(beam_args,dict):beam_args={}
            beam_args['refl']= refl
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


    def show_Idyn_vs_Ikin(self,thicks=None,iZs=100,cmap='Spectral',**kwargs):
        """Display correlation plot of dynamical vs kinematical intensities

        Parameters
        -----------
        thicks
            thicknesses
        hkl
            beams to include
        """
        if thicks or not 'Iz' in self.__dict__:self._set_beams_vs_thickness(thicks)
        if isinstance(iZs,int):iZs=slice(0,None,iZs)
        # iB = self.get_Istrong(out=1)[0]#;print(iB,self.Iz[    iB,iZs],)
        # iB = np.argsort(np.sum(self.Iz,axis=1))[-2]
        z      = self.z.copy()[       iZs]
        Iz_dyn = self.Iz.copy()[    :,iZs]#/self.Iz[    iB,iZs]
        Iz_kin = self.Iz_kin.copy()[:,iZs]#/self.Iz_kin[iB,iZs]#;print(Iz_dyn[Iz_dyn>1e-3])
        nzs = z.size
        # print(Iz_kin.shape,Iz_dyn.shape)
        cs = dsp.getCs(cmap,nzs) #; print(len(cs),Iz_dyn.shape,Iz_kin.shape)
        # scat=tuple( [[I_kin,I_dyn,cs[i]] for i,(I_dyn,I_kin) in enumerate(zip(self.Iz[:,::iZ],self.Iz_kin[:,::iZ]))])
        plts=[[np.log10(I_kin),np.log10(I_dyn),[cs[i],'o'],'$z=%d\AA$' %z] for i,(z,I_dyn,I_kin) in enumerate(zip(z,Iz_dyn.T,Iz_kin.T))]
        # plts+=[ [[0,1],[0,1],[(0.5,)*3,'--'],''] ]
        return dsp.stddisp(plts,labs=['$I_{kin}$','$I_{dyn}$'],sargs={'alpha':0.5},**kwargs)

    def show_Fhkl(self,s=None,opts='m',h3D=0,**kwargs):
        """Displays structure factor over grid

        Parameters
        -----------
        opts
            see bloch_util:get_fz
        s
            slice or str('k=0' => Fhkl(k=0)) or int('l==<int>')
        """
        fz,fz_str = bloch_util.get_fz(opts)
        s,s_str = self._get_slice(s)
        tle = 'Structure factor($\AA$), showing %s in %s'  %(fz_str,s_str)
        Fhkl  = self.get_Fhkl() #.copy()
        if isinstance(s,tuple):
            Fhkl   = Fhkl[s]
            nx,ny  = np.array((np.array(Fhkl.shape)-1)/2,dtype=int)
            i,j    = np.meshgrid(np.arange(-nx,nx+1),np.arange(-ny,ny+1))
            fig,ax = dsp.stddisp(scat=[i,j,fz(Fhkl)],title=tle,**kwargs)
        else:
            h,k,l = self.get_hklF()
            fig,ax = dsp.stddisp(scat=[h,k,l,fz(self.get_Fhkl())],title=tle,rc='3d',**kwargs)
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

    def _make_img(self,
            exp:str=None,
            pred:bool=False,
            Imax:int=1,
            Nmax:int=512,
            rot:int=0,
            aperpixel:Optional[float]=None,
            fbroad=None,
            gs3:float=0.1,
            nX:int=1,
            rmax:int=0,
            thick:float=None,
            iz:int=None,
        ):
        '''Make image

        Parameters
        -----------
        exp
            path to dials or Dials (overrides Nmax,aperpixel)
        Imax
            max value for the intensities
        Nmax
            image resolution in pixel
        aperpixel
            reciprocal size for each pixel (aperture per pixel)
        rot
            in plane rotation of the image
        pred
            True to use dials predict (ignores rot)
        fbroad
            broadening function f(r2) (default np.exp(-r2/(gs3/3)**2))
        gs3
            Gaussian broadening factor for each reflections
        nX
            width factor of the Gaussian window (in pixels)
        rmax
            radius for the noise to be added
        show
            Show produced image
        kwargs
            args to be passed to the tiff viewer
        '''
        thick = self.thick
        if thick:self.set_thickness(thick)
        #dials info
        if isinstance(exp,str):exp = pt.Dials(exp)
        if exp:
            aperpixel=exp.aper
            Nmax=exp.nxy
            if not Imax>0 or Imax==1 :
                Imax=exp.Imax*5000
                print(colors.red+'setting Imax to %.1E '% Imax+colors.black)

        hkl = self.df_G.index
        pred = pred and exp
        if pred :
            hkl = [h for h in hkl if h in exp.df_pred.index]
            hkl_lost = np.setdiff1d(self.df_G.index,hkl)
            if any(hkl_lost):
                print(colors.red+"warning : reflections ignored "+colors.black)
                print(self.df_G.loc[hkl_lost].sort_values('I')['I'][-10:])
            pxy = exp.df_pred.loc[hkl,['px','py']]

        # get intensity
        px,py,I = self.df_G.loc[hkl,['px','py','I']].values.T
        Nmax=Nmax//2
        if isinstance(iz,int):
            idb=self.get_beam(refl=hkl)
            I = self.Iz[idb,iz]
            thick = self.z[iz]

        # pixel locations
        if not aperpixel:
            aperpixel = 1.1*max(px.max(),py.max())/Nmax
            print('aperpixel set to %.1E A^-1 ' %(aperpixel))
        dqx,dqy = [aperpixel]*2
        if pred:
            i,j = np.array(pxy[['px','py']].values,dtype=int).T
        else:
            if rot:
                ct,st = np.cos(np.deg2rad(rot)),np.sin(np.deg2rad(rot))
                px,py = ct*px-st*py,st*px+ct*py
            i,j = np.array([np.round(px/dqx),np.round(py/dqy)],dtype=int)+Nmax


        #### kernel (converted to pixel locations)
        if not fbroad:
            fbroad=lambda r2:np.exp(-r2/(gs3/3)**2)
            nx,ny = np.array(np.floor(gs3/np.array([dqx,dqy])),dtype=int)
        else:
            nx,ny=nX,nX
        ix,iy = np.meshgrid(range(-nx,nx+1),range(-ny,ny+1))
        x,y = ix*dqx,iy*dqy
        ## Gaussian broadening
        r2 = (x**2+y**2)
        Pb = fbroad(r2)#;dsp.stddisp(im=[x,y,Pb],pOpt='im')
        im0 = np.zeros((2*Nmax,)*2)
        for i0,j0,I0 in zip(i,j,I):
            i0x,j0y = i0+ix,j0+iy
            idx     = (i0x>=0) & (j0y>=0)  & (i0x<2*Nmax) & (j0y<2*Nmax)
            i0x,j0y = i0+ix[idx],j0+iy[idx]
            im0[i0x,j0y] += Pb[idx]/Pb[idx].sum()*I0
            # im0[i0,j0]=I0

        #### noise
        if rmax:
            h,k = np.meshgrid(range(-Nmax,Nmax),range(-Nmax,Nmax))
            r = np.sqrt(h**2+k**2);r[r==0]=1
            im0 += rmax*np.random.rand(*im0.shape)#/(rmax*r)

        return im0*Imax

    def convert2img(self,filename,template=None,**kwargs):
        im0 = self._make_img(**kwargs)#Nmax,fbroad,gs3,nX,rmax,thick,iz,rot)
        # print(im0.mean())
        fmt = template.split('.')[-1]
        if not filename:
            filename = os.path.join(self.figpath,self.name+'_%dA.%s' %(self.thick,fmt))
        if template:
            out = check_output("cp %s %s" %(template,filename),shell=True).decode()
            if out:print(out)
        bloch_util.imwrite(filename,im0)



    def convert2tiff(self,tiff_file:str=None,
        Imax:int=3e4,
        Nmax:int=512,aperpixel:Optional[float]=None,
        fbroad=None,
        gs3:float=0.1,
        nX:int=1,
        rmax:int=0,
        thick:float=None,
        iz:int=None,
        rot:int=0,
        show:bool=False,
        tif_writer_args:dict={},
        **kwargs,
    ):
        """Write intensities to a tiff file

        Parameters
        -----------
        tiff_file
            name of the file (default name if not provided)
        thick
            Thickness of the simu
        Imax
            max value for the intensities
        Nmax
            image resolution in pixel
        aperpixel
            reciprocal size for each pixel (aperture per pixel)
        fbroad
            broadening function f(r2) (default np.exp(-r2/(gs3/3)**2))
        gs3
            Gaussian broadening factor for each reflections
        nX
            width factor of the Gaussian window
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
        im0 = _make_img(Nmax,fbroad,gs3,nX,rmax,thick,iz,rot)

        if not tiff_file:
            tiff_file = os.path.join(self.path,self.name+'_%dA' %thick+'.tiff')
        I = np.array(im0*Imax,dtype='uint16')

        # ix,iy = np.meshgrid(range(2*Nmax),range(2*Nmax))
        # dsp.stddisp(im=[ix,iy,I],plots=[j,i,'bo'],xylims=[0,512,0,512],
        #     cmap='gray',caxis=[0,10],imOpt='tX',pargs={'fillstyle':'none'})

        tifffile.imwrite(tiff_file,I.T,**tif_writer_args)
        print(colors.yellow+tiff_file+colors.green+' saved'+colors.black)
        if show:
            # v=viewers.Base_Viewer(self.path,frame=1,thick=self.thick,**kwargs)
            return ut.show_tiff(tiff_file,**kwargs)
            # return v
        else:
            return tiff_file


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
        opt:  S(solve) s(save) p(plot)
        """
        solve,save,show = 'S' in opts,'s' in opts,'p' in opts
        Nmax,Smax=[a.flatten() for a in np.meshgrid(Nmax,Smax)]
        z=np.array(z)
        df=pd.DataFrame(columns=['Nmax','Smax','nbeams','Iz'])
        # simulate
        self.update_Nmax(Nmax[0])
        self._set_excitation_errors(Smax=Smax[0])
        hkl0 = self.df_G.index
        if solve:
            for i,(Nm,Sm) in enumerate(zip(Nmax,Smax)):
                 self.solve(Nmax=Nm,Smax=Sm,thicks=z,opts='v0')
                 idx='(%d,%.3f)'%(Nm,Sm)
                 idx_b = self.get_beam(refl=hkl0)
                 df.loc[idx,['Nmax','Smax','nbeams']]=[Nm,Sm,self.nbeams]
                 df.loc[idx,'Iz']=self.Iz[idx_b,:].copy()
            self.df=df
            if save:
                df.to_pickle(self.path+self.name+'_cv.pkl')

        if show:return self.show_convergence(hkl)
        # else:
        #     for i,(Nm,Sm) in enumerate(zip(Nmax,Smax)):
        #         self.update_Nmax(Nm)
        #         self._set_excitation_errors(Smax=Sm)
        #         idx='(%d,%.3f)'%(Nm,Sm)
        #         df.loc[idx,['Nmax','Smax','nbeams']]=[Nm,Sm,self.nbeams]
        #     print(df)

    def show_convergence(self,hkl,xlab='Smax',**kwargs):
        #beam selection
        valid_hkl=0
        # if isinstance(hkl,list):valid_hkl=isinstance(hkl[0],str)
        self.update_Nmax(int(self.df.Nmax.min()))
        self._set_excitation_errors(Smax=self.df.Smax.min())
        # if not valid_hkl:
        #     hkl = self.get_beam(cond='Sw<1e-2' ,opt=0)[0:2]
        #     hkl += [self.get_beam(cond='Sw>=1e-2',opt=0)[0]]
        #     print('default hkl : ',hkl)
        # Icols=hkl
        # self.df[Icols]=None
        x = self.df[xlab]
        idx_s = self.get_beam(refl=hkl)
        #Iz.shape:nbeams x nsimus
        Iz = np.array([r.Iz[idx_s,-1] for i,r in self.df.iterrows()]).T
        cs=dsp.getCs('Spectral',len(hkl))
        #plot
        plts=[ [x,I,[c,'-o'],h]  for i,h,c,I in zip(idx_s,hkl,cs,Iz)]

        # ## thickness dependence
        # cs=dsp.getCs('RdYlGn',self.df.shape[0])
        # ms = dsp.markers[:len(hkl)]
        # for idx,m in zip(idx_s,ms):
        #     plts+=[[self.z,r.Iz[idx],[c,m+'-' ],
        #          ] for c,(l,r) in zip(cs,self.df.iterrows())]
        # legElt={'$%s$' %h:'k-'+ms[i] for i,h in enumerate(hkl)}
        # legElt.update({'$%s$' %l:[c,'-'] for c,l in zip(cs,self.df.index)})
        # return dsp.stddisp(plts,labs=['$z$','$I$'],legElt=legElt,lw=2)
        return dsp.stddisp(plts,labs=['$%s$' %xlab,'$I$'],lw=2,**kwargs)

    def solve_strong_beam(self,thick=500,solve_args={},
        Imin=1e-4,dImin=0.01):
        '''Run the simulation for only the strongest beams for which I>Imin
        and iteratively removes all the reflections that affect those beams
        by less than dImin.
        Parameters :
        -----------
        solve_args : parsed to self.solve
        Imin : strong beam criteria
        dImin : perturbation beam criteria (relative error )
        thick : thickness at which the procedure is performed
        '''
        self.solve(**solve_args)
        self.set_thickness(thick=thick,v=False)
        ## gather strong beams
        hkls        = self.df_G.index
        hkl_strong  = self.df_G.loc[self.df_G.I>Imin].index
        hkl_weak    = np.setdiff1d(hkls,hkl_strong)
        # print('hkl_strong',hkl_strong.values)
        # print('hkl_weak',hkl_weak)

        self.df_hkl_weak = pd.DataFrame(index=hkl_weak,columns=['dI','dImax','dhmax'])
        solve_args['Smax']=0

        self.df_hkl_weak.loc[hkl_weak,'Iref'] = self.df_G.loc[hkl_weak,'I']
        self.df_hkl_weak.loc[hkl_weak,'Inew'] = 0
        self.df_hkl_weak.loc[hkl_weak,'in']  = True
        I0 = self.df_G.loc[hkl_strong,'I'].copy()
        Iref = I0.copy()
        for hi in hkl_weak:
            #### run without beam hi
            _hkls = np.array([eval(s) for s in hkls.drop(hi)])
            self.solve(hkl=_hkls,**solve_args)
            self.set_thickness(thick=thick,v=False)

            v = np.abs(self.df_G.loc[hkl_strong,'I']-Iref)/Iref
            # print(hi,len(_hkls),v.max(),hkl_strong[v.argmax()])
            self.df_hkl_weak.loc[hi,'dI']=','.join(v.astype(str))
            self.df_hkl_weak.loc[hi,'dImax'] = v.max()
            self.df_hkl_weak.loc[hi,'dhmax'] = hkl_strong[v.argmax()]
            if v.max()<dImin:
                hkls=hkls.drop(hi)
                Iref = self.df_G.loc[hkl_strong,'I']
                self.df_hkl_weak.loc[hi,'in'] = False


        hkl_in = self.df_hkl_weak.loc[self.df_hkl_weak['in']==True].index
        #rerun previous to include last perturbation beam if necessary
        if self.df_hkl_weak.loc[hkl_weak[-1],'in']:
            print('rerunning for %s' %hkl_weak[-1])
            _hkls = np.array([eval(s) for s in hkls])
            self.solve(hkl=_hkls,**solve_args)
            self.set_thickness(thick=thick,v=False)


        self.df_hkl_weak.loc[hkl_in,'Inew'] = self.df_G.loc[hkl_in,'I']
        self.df_hkl_strong = self.df_G.loc[hkl_strong,['I']].copy()
        # self.df_G.loc['Iref'] = 0
        self.df_hkl_strong['Iref'] = I0
        # self.df_G.loc[hkl_in    ,'Iref'] = self.df_hkl_weak.loc[hkl_in,'Iref']
        self.save()
