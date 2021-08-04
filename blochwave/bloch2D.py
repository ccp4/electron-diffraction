# from . import bloch as bl
import importlib as imp
from EDutils.display  import dfb
import numpy as np,pandas as pd,pickle5,os,glob
from crystals import Crystal
from multislice import mupy_utils as mut     ;imp.reload(mut)
from multislice import postprocess as pp     ;imp.reload(pp)
from scattering import structure_factor as sf;imp.reload(sf)
from utils import displayStandards as dsp    ;imp.reload(dsp)
from utils import physicsConstants as cst    ;imp.reload(cst)
from utils import glob_colors as colors,handler3D as h3d

################################################################################
class Bloch_B:
    def __init__(self,keV=200,K=None,u=[0,0,1],Nmax=1,
        Smax=0.2,thick=500,name='',path=None,solve=0,opts='sv'):
        ''' Bloch wave simulation class
        - cif_file : str - structure
        - name     : str - optional name (used to save object see set_name)
        - K,u,keV  : reciprocal space beam vector(see update beam)
        - Nmax     : int - max order of reflections/resolution(see update_Nmax)
        - Smax     : float - maximum excitation error to be included(see solve)
        - thick    : float - thickness of crystal(can be modified after solve)
        - solve    : bool - diagonalize the Bloch matrix and find intensities
        - opts     : str - (see solve)
        '''
        self.Nmax  = 0
        self.thick = thick
        self.update_Nmax(Nmax)
        self.set_beam(K,u,keV)
        self.set_name(name,path)
        if solve or 'S' in opts:self.solve(Smax,opts=opts)

    def set_name(self,name='',filepath=''):
        '''set name for Bloch obj(will be saved as path+name+'.pkl'):
        - name : The name of the simulation(can include the full path)
        - path : The path to the simulation folder
        By default the simulation name is '<cif><zone axis>_<keV>keV_bloch'
        '''
        if not name:
            basefile = self.get_basefile()
            u_str = ''.join(['%d' %np.round(u) for u in self.Kabc0])
            name='%s%s_%dkeV_bloch' %(basefile,u_str,np.round(self.keV))
        if not filepath:filepath=os.path.dirname(name)
        self.path = filepath                            #; print(self.filepath)
        self.name = os.path.basename(name)              #;print(self.name)

    def update_Nmax(self,Nmax):
        '''Update resolution/max order. (Updates lattice and structure factor)
        Nmax : maximum h,k,l order
        '''
        if isinstance(Nmax,int):
            if not Nmax==self.Nmax:
                self.Nmax=Nmax
                self.update_lattice()
                self.update_struct_fact()

    def set_beam(self,K=None,u=None,keV=None):
        ''' Update the beam
        - keV : float - update beam wavelength
        - u   : list - update beam direction
        - K   : list - update beam direction and wavelength
        '''
        if keV:self.k0 = 1/cst.keV2lam(keV)
        if isinstance(K,list) or isinstance(K,np.ndarray):
            self.K  = K
            self.k0 = np.linalg.norm(K)
        elif isinstance(u,list) or isinstance(u,np.ndarray):
            self.K = self.k0*np.array(u)/np.linalg.norm(u)
        self.lam = 1/self.k0
        self.keV = cst.lam2keV(self.lam)
        self.sig = cst.keV2sigma(self.keV)
        self.u   = self.K/self.k0
        self.Kuvw = self.lat_vec.dot(self.K)  #projection in reciprocal crystal basis
        self.Kabc = self.lat_vec0.dot(self.K) #projection in crystal basis
        self.Kuvw0 = self.Kuvw/(abs(self.Kuvw)[abs(self.Kuvw)>0.01].min()) #zone axis notation in reciprocal crystal basis ;print(Kabc)
        self.Kabc0 = self.Kabc/(abs(self.Kabc)[abs(self.Kabc)>0.01].min()) #zone axis notation in crystal basis ;print(Kabc)

    def set_thickness(self,thick):
        '''set thickness and update beams
        - thick : thickness
        '''
        if type(thick) in [int,float,np.int64,np.float64] :self.thick=thick
        self._set_kinematic()
        if self.solved:self._set_intensities()

    def solve(self,Smax=0.02,Nmax=None,K=None,u=None,keV=None,thick=None,
        opts='sv',Vopt0=False,v=True,name=''):
        ''' diagonalize the Blochwave matrix
        - Smax    : float- maximum excitation error to be included
        - Nmax    : int - max order of reflections/resolution(see update_Nmax)
        - K,u,keV : reciprocal space beam vector(see update beam)
        - opts  : 'str - s'(save) '0'(Vopt0) 'v'(verbose) 'H'(show H)
        - Vopt0 : bool - set Vg0=0
        - v     : bool - verbose
        '''
        self.update_Nmax(Nmax)
        self.set_beam(K,u,keV)
        # self.set_name(name,self.path)
        self._set_excitation_errors(Smax)
        self._solve_Bloch(opts,Vopt0,v)
        self._set_Vg()
        self.set_thickness(thick)
        self._set_kinematic()
        if 's' in opts:self.save()

    ################################################################################
    #### private
    ################################################################################
    def _solve_Bloch(self,opts='0',Vopt0=True,v=0):
        ''' Diagonalize the Hamiltonian
        Ug is a (2*Nmax+1)^3 tensor :
        # Ug[l] = [U(-N,-N,l) .. U(-N,0,l) .. U(-N,N,l)
        #          U( 0,-N,l) .. U( 0,0,l) .. U( 0,N,l)
        #          U(-N,-N,l) .. U( N,0,l) .. U( N,N,l)]
        '''
        Vopt0 = Vopt0 or '0' in opts
        v = v or 'v' in opts

        Gs = self.df_G[self.refl].values
        Sg  = self.df_G.Sw.values
        Ug = self._get_Ug()

        # Ug[U0_idx] = U0
        # and Ug(iG,jG) are obtained from Ug[h,k,l] where h,k,l = hlk_iG-hkl_jG
        G0 = [2*self.Nmax]*self.ndim
        Ug[tuple(G0)] = 0   #setting average potential to 0

        if v:print(colors.blue+'...assembling %dx%d matrix...' %((Sg.shape[0],)*2)+colors.black)
        H = np.diag(Sg+0J)
        for iG,Gi in enumerate(Gs) :
            U_iG = np.array([Ug[tuple(Gij+G0)] for Gij in Gi-Gs]) #;print(V_iG.shape)
            # print('Gi : ',Gi)
            # print('U_G:',Ug[tuple(U0_idx+hkl_G)])
            # print('Gi-Gj:',Gi-Gs)#[tuple(Gij+G0) for Gij in Gi-Gs])
            # print('UiG : ',U_iG)
            H[iG,:] += U_iG/(2*self.k0)  #off diagonal terms as potential

        self.H = H
        if 'H' in opts:self.show_H()

        if v:print(colors.blue+'...diagonalization...'+colors.black)
        self.gammaj,self.CjG = np.linalg.eigh(H) #;print(red+'Ek',lk,black);print(wk)
        self.invCjG = np.linalg.inv(self.CjG)
        self.solved = True


    def _set_Vg(self):
        Vg     = self._get_Ug()
        V0_idx = np.array([2*self.Nmax]*self.ndim)
        hkl    = self.df_G[self.refl].values
        Vg[tuple(V0_idx)] = 0
        Vg_G       = np.array([ Vg[tuple(hkl_G+V0_idx)] for hkl_G in hkl])

        self.set_beam_positions()
        self.df_G['Vg'] = Vg_G
        self.df_G['L']  = np.ones(Vg_G.shape)
        # self.df_G.loc[(self.df_G.h==0) & (self.df_G.k==0) & (self.df_G.l)==0]['Vg'] = 0
        self._set_zones()

    def _set_zones(self):
        refl = self.df_G[self.refl].values
        # Khkls   = np.round(hkl.dot(self.Kuvw/np.linalg.norm(self.Kuvw))*100)/100
        Khkls   = np.round(refl.dot(self.Kuvw0)*100)/100 #integer only in zone axis orientation
        Khkl,ar = np.unique(Khkls,return_inverse=True) #; print(Khkls);print(zones);print(ar)
        zones = np.argsort(Khkl)
        self.df_G['zone'] = zones[ar]

    def _set_intensities(self):
        '''get beam intensities at thickness'''
        id0 = self.is_hkl([0]*self.ndim,v=0)
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

    def set_beams_vs_thickness(self,thicks):
        ''' get Scattering matrix as function of thickness for all beams
        - thicks : tuple(ti,tf,step) or list or np.ndarray thicknesses
        '''
        self._set_thicks(thicks)
        id0 = self.is_hkl([0]*self.ndim,v=0)
        St = np.zeros((self.df_G.shape[0],self.z.size),dtype=complex)
        gammaj,CjG,invCjG = self.gammaj,self.CjG,self.invCjG
        for iT,thick in enumerate(self.z):
            # S=CjG.T.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG)
            S = CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*thick))).dot(invCjG)
            # St[:,iT]=CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG.T)[0,:]
            St[:,iT] = S[:,id0]
        self.Sz = St
        self.Iz = np.abs(self.Sz)**2

    def _set_thicks(self,thicks):
        if isinstance(thicks,tuple):thicks = np.linspace(*thicks)
        elif isinstance(thicks,float) or isinstance(thicks,int) : thicks=[thicks]
        self.z = np.array(thicks)

    ################################################################################
    #### getter
    ################################################################################
    def get_intensities(self):return self.df_G.I
    def get_kin(self):return self.df_G[self.refl+['Sw','Vg','Ig']]
    def get_zones(self):return self.df_G[self.refl+['zone']].values
    def get_Xig(self):
        self.df_G['Xi_g'] = self.k0/abs(self.df_G.Vg)
        return self.df_G.loc[self.df_G.Xi_g<1e4,self.refl+['Vg','Xi_g']]
    def get_G(self):return self.df_G[self.q].values
    def get_refl(self,idx=[]):
        if np.array(idx).size:
            return self.df_G[self.refl].iloc[idx].values
        else:
            return self.df_G[self.refl].values
    def get_Sw(self,Smax=1):
        Smax = min(Smax,self.Smax)
        print(self.df_G[self.refl+['Sw']].loc[self.df_G.Sw<=Smax])
    def is_hkl(self,Ghkl,v=1):
        Gidx = np.where(np.linalg.norm(self.get_refl() - Ghkl,axis=1)==0)[0]
        if v:print(self.df_G.iloc[Gidx])
        return Gidx[0]
    def get_Istrong(self,m={'I':1000,'Ig':1000,'Vg':100},out=0,Icols=['I']):
        cond = self.df_G.L<0
        Icols = [Icol for Icol in Icols if Icol in self.df_G.columns and Icol in m.keys()]
        for Icol in Icols :
            Imax = abs(self.df_G[Icol].max())
            cond |= (abs(self.df_G[Icol])>Imax/m[Icol])
        Istrong = self.df_G.loc[cond]
        if out:return Istrong.index.values
        else:
            if any(Istrong):print(Istrong[self.refl+['Sw']+Icols])

    def G(self,col):print(self.df_G[self.refl+[col]])

    ################################################################################
    #### display
    ################################################################################
    def show_beams_vs_thickness(self,thicks=None,strong=['I'],m={'I':1000,'Ig':1000,'Vg':10},
        cm='jet',**kwargs):
        if thicks:self.set_beams_vs_thickness(thicks)
        refl = self.get_refl().copy()
        Iz  = self.Iz.copy()
        if any(strong):
            Istrong = self.get_Istrong(out=1,Icols=strong,m=m) #;print(Istrong)
            Iz  = Iz[Istrong]
            refl = refl[Istrong]
        refl_str = [str(tuple(h)) for h in refl]
        beams=[refl,self.z,None,None,Iz]
        # beams=[hkl,self.z,np.real(self.Sz),np.imag(self.Sz),self.Iz]
        return pp.plot_beam_thickness(beams,cm=cm,**kwargs)


    def show_H(self,**kwargs):
        dsp.stddisp(im=[np.abs(self.H)],title='abs(H)', pOpt='im')

    ################################################################################
    #### misc
    ################################################################################
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
    def save(self,file=None,v=1):
        '''save this object'''
        if not file:file=os.path.join(self.path,self.name+'.pkl')
        with open(file,'wb') as out :
            pickle5.dump(self, out, pickle5.HIGHEST_PROTOCOL)
        if v:print(colors.green+"object saved\n"+colors.yellow+file+colors.black)


from wallpp import wallpaper as wallpp;imp.reload(wallpp)
class Bloch2D(Bloch_B):
    def __init__(self,file='',crys=None,eps=0.1,u=[0,1],**kwargs):
        self.file = file
        if not crys:crys=file
        # self.crys = mut.import_wallpp(crys)
        self.crys = wallpp.Wallpaper(**crys)
        self.lat_vec0 = self.crys.lattice_vectors
        self.lat_vec  = self.crys.reciprocal_vectors#/(2*np.pi)
        self.pattern  = self.crys.pattern_fract
        self.refl = ['h','k']
        self.q    = ['qx','qy']
        self.ndim=2
        self.eps=eps
        super().__init__(u=u,**kwargs)

    def update_lattice(self):
        (h,k),(qx,qz) = mut.get_lattice2D(self.lat_vec,self.Nmax)
        self.lattice = [(h,k),(qx,qz)]
    def update_struct_fact(self):
        self.hkF,self.Fhk = sf.structure_factor2D(self.pattern, 2*np.pi*self.lat_vec,hkMax=2*self.Nmax,eps=self.eps)
    def _get_Ug(self):
        #in 2D Fhk is just the Fourier transform of fv
        vg = self.Fhk.copy()
        #It needs to be converted into a scattering amplitude [A]
        #Fg[no dim] = vg[kVA^2], sig=[1/kVA], k0[A-1]
        Fg = vg*self.sig*self.k0/np.pi
        # Ug[A^-2]
        Ug =  Fg/self.crys.area
        return Ug

    def set_beam_positions(self):
        self.df_G['px'] = mut.project_beams2D(K=self.K,qxy=self.get_G())
    def _set_excitation_errors(self,Smax=0.02):
        ''' get excitation errors for Sg<Smax
        - Smax : maximum excitation error to be included
        '''
        K,K0 = self.K,self.k0
        (h,k),(qx,qy) = self.lattice

        Kx,Ky = K
        Sw = np.abs(np.sqrt((Kx-qx)**2+(Ky-qy)**2) - K0)
        if Smax:
            idx = Sw<Smax
            h,k = np.array([h[idx],k[idx]],dtype=int)
            qx,qy,Sw = qx[idx],qy[idx],Sw[idx]
        d = dict(zip(['h','k','qx','qy','Sw'],[h,k,qx,qy,Sw]))

        self.Smax   = Smax
        self.nbeams = Sw.size
        self.df_G   = pd.DataFrame.from_dict(d)
        self.solved = False

    def get_basefile(self):
        return self.file.replace('.%s' %self.file.split('.')[-1],'')

    def show_ewald(self,ax=None,**kwargs):
        fig,ax = mut.show_Ewald2D(self.K,self.lattice[1],ax=ax,legOpt=0,opt='')
        x,z = self.df_G[['qx','qy']].values.T
        scat = [x,z,80,'b']
        dsp.stddisp(ax=ax,scat=scat,**kwargs)

    def show_beams(self,opts='B',fz=None,**kwargs):
        '''Display beam values
        - opts : B(Bloch), K(kin), V(Vg), S(Sg), L(lat)
        '''
        if opts=='B' and not fz:fz=np.abs

        tle = 'u=%s, thickness=%d$\AA$'  %(str(np.round(10*self.Kuvw0)/10),self.thick)
        px = self.df_G.px.values
        pm = 1.25*abs(px).max()
        x = np.arange(-pm,pm,0.001)#np.linspace(-pm,pm,1000)
        plts = []
        dq = {'B':0.01,'K':0.01,'V':0.005,'S':0.001}
        for c,C in dfb.iterrows() :
            if c in opts:
                if not fz:
                    F = C.fz(self.df_G[C.F])
                else:
                    F = fz(self.df_G[C.F])
                y = np.zeros(x.shape)
                for p,v in zip(px,F):y += v/(1+((x-p)/dq[c])**2)
                plts += [[x,y,C.color,C.legend]]

        fig,ax = dsp.stddisp(plts,lw=2,title=tle,**kwargs)

    def show_Fhk(self,s=None,opts='m',h3D=0,**kwargs):
        '''Displays structure factor over grid
        - opts : see get_fz
        - s : slice or str('k=0' => Fhk(k=0)) or int('l==<int>')
        '''
        fz,fz_str = bl.get_fz(opts)
        # s,s_str = self._get_slice(s)
        tle = 'Structure factor($\AA$), showing %s in %s'  %(fz_str,s_str)
        h,k = self.lat
        Fhk = self.Fhk #.copy()
        fig,ax = dsp.stddisp(scat=[h,k,fz(Fhk)],title=tle,**kwargs)



class Bloch(Bloch_B):
    def __init__(self,cif_file,u=[0,0,1],**kwargs):
        self.cif_file = cif_file
        self.crys     = mut.import_crys(cif_file)
        self.lat_vec0 = np.array(self.crys.lattice_vectors)
        self.lat_vec  = np.array(self.crys.reciprocal_vectors)/(2*np.pi)
        self.pattern  = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in self.crys.atoms] )
        self.ndim = 3

        super().__init__(u=u,**kwargs)

    def update_lattice(self):
        (h,k,l),(qx,qy,qz) = mut.get_lattice(self.lat_vec,self.Nmax)
        self.lattice = [(h,k,l),(qx,qy,qz)]
    def update_struct_fact(self):
        self.hklF,self.Fhkl = sf.structure_factor3D(self.pattern, 2*np.pi*self.lat_vec,hklMax=2*self.Nmax)
    def _get_Ug(self):
        return self.Fhkl/self.crys.volume #/3
    def set_beam_positions(self):
        px,py,e0x  = mut.project_beams(K=self.K,qxyz=self.get_G(),e0=[1,0,0],v=1)
        self.e0x = e0x
        self.df_G['px'] = px
        self.df_G['py'] = py
    def _set_excitation_errors(self,Smax=0.02):
        ''' get excitation errors for Sg<Smax
        - Smax : maximum excitation error to be included
        '''
        K,K0 = self.K,self.k0
        (h,k,l),(qx,qy,qz) = self.lattice

        Kx,Ky,Kz = K
        Sw = np.abs(np.sqrt((Kx-qx)**2+(Ky-qy)**2+(Kz-qz)**2) - K0)
        if Smax:
            idx = Sw<Smax
            h,k,l = np.array([h[idx],k[idx],l[idx]],dtype=int)
            qx,qy,qz,Sw = qx[idx],qy[idx],qz[idx],Sw[idx]
        d = dict(zip(['h','k','l','qx','qy','qz','Sw'],[h,k,l,qx,qy,qz,Sw]))

        self.Smax = Smax
        self.nbeams = Sw.size
        self.df_G = pd.DataFrame.from_dict(d)
        self.solved = False

    def get_basefile(self):
        return self.cif_file.replace('.cif','')

    def show_beams(self,F='I',fopts='m',opts='xN',mag=500,cutoff=0,cmap='Greens',**kwargs):
        '''Display beam values
        - F : 'I'(intensity),'S'(scattered beams),'Vg'(potential)
        - opts  : 'x'(include projection of x axis), 'N'(normalize)
        - fopts : see get_fz
        '''
        fz,fz_str = get_fz(fopts)
        F_str = Fs_str[F]

        tle = '%s, %s, thickness=%d$\AA$'  %(F_str,fz_str,self.thick)
        px,py = self.df_G[['px','py']].values.T
        Fvals = fz(self.df_G[F])
        if 'N' in opts:Fvals/=Fvals.max()
        # mag /= Fvals.max()
        scat  = [px,py,Fvals*mag/Fvals.max(),Fvals]
        plts = [[0,0,'b+','']]
        if 'x' in opts:plts+=[ [[0,self.e0x],[0,0],'r','']]
        if not cutoff:cutoff = Fvals.max()
        # print(Fvals.max())
        fig,ax = dsp.stddisp(plts,lw=2,scat=scat,caxis=[0,cutoff],cmap=cmap,cs='S',title=tle,**kwargs)
Fs_str = {'L':'Lattice','Vg':'Potential $V_{g}$',
    'S' :'Scattered beams Bloch $S_{g}$',
    'I' :'Intensity Bloch$I_{g}$',
    'Sg':'Scattered beams Kinematic $S_{g,kin}$',
    'Ig':'Intensity kinematic $I_{g,kin}$',
    'Sw':'Excitation error $\zeta_{g}$'}



################################################################################
