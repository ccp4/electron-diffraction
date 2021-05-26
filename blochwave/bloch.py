'''Bloch wave solver'''
import importlib as imp
import numpy as np,pandas as pd,pickle5,os,glob
from crystals import Crystal
from multislice import mupy_utils as mut     ;imp.reload(mut)
from scattering import structure_factor as sf;imp.reload(sf)
from utils import displayStandards as dsp    ;imp.reload(dsp)
from utils import physicsConstants as cst    ;imp.reload(cst)
from utils import glob_colors as colors,handler3D as h3d


class Bloch:
    def __init__(self,cif_file,keV=200,K=None,u=[0,0,1],Nmax=1,
        Smax=0.2,thick=500,name='',path='',solve=0,opts='sv'):
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
        self.cif_file = cif_file
        self.crys     = mut.import_crys(cif_file)
        self.lat_vec0 = np.array(self.crys.lattice_vectors)
        self.lat_vec  = np.array(self.crys.reciprocal_vectors)/(2*np.pi)
        self.pattern  = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in self.crys.atoms] )
        self.Nmax  = 0
        self.thick = thick
        self.update_Nmax(Nmax)
        self.set_beam(K,u,keV)
        self.set_name(name,path)
        if solve or 'S' in opts:self.solve(Smax,opts=opts)

    def set_name(self,name='',path=''):
        '''set name for Bloch obj(will be saved as path+name+'.pkl'):
        - name : The name of the simulation(can include the full path)
        - path : The path to the simulation folder
        By default the simulation name is '<cif><zone axis>_<keV>keV_bloch'
        '''
        if not name:
            basefile = self.cif_file.replace('.cif','')
            Kabc = self.Kabc.copy()
            Kabc /= abs(Kabc)[abs(Kabc)>0.01].min() #;print(Kabc)
            u_str = ''.join(['%d' %np.round(u) for u in Kabc])
            name='%s%s_%dkeV_bloch' %(basefile,u_str,np.round(self.keV))
        if not path:os.path.dirname(name)
        self.path = path                            #; print(self.path)
        self.name = os.path.basename(name)

    def update_Nmax(self,Nmax):
        '''Update resolution/max order. (Updates lattice and Fhkl)
        Nmax : maximum h,k,l order
        '''
        if isinstance(Nmax,int):
            if not Nmax==self.Nmax:
                self.Nmax=Nmax
                (h,k,l),(qx,qy,qz) = mut.get_lattice(self.lat_vec,self.Nmax)
                self.lattice = [(h,k,l),(qx,qy,qz)]
                self.hklF,self.Fhkl = sf.structure_factor3D(self.pattern, 2*np.pi*self.lat_vec,hklMax=2*self.Nmax)

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
        self.u   = self.K/self.k0
        self.Kuvw = self.lat_vec.dot(self.K)  #projection in reciprocal crystal basis
        self.Kabc = self.lat_vec0.dot(self.K) #projection in crystal basis

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
        self.df_G = pd.DataFrame.from_dict(d)

    def set_thickness(self,thick):
        '''set thickness and update beams
        - thick : thickness
        '''
        if isinstance(thick,int):self.thick=thick
        self._set_intensities()

    def set_beams_vs_thickness(self,thicks):
        ''' get Scattering matrix as function of thickness for all beams
        - thicks : int or list or np.ndarray thicknesses
        '''
        if isinstance(thicks,float) or isinstance(thicks,int) : thicks=np.array([thicks])
        if isinstance(thicks,list):thicks = np.array(thicks)
        St = np.zeros((self.df_G.shape[0],Ts.size),dtype=complex)
        gammaj,CjG = self.gammaj,self.CjG
        for iT,T in enumerate(thicks):
            # S=CjG.T.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG)
            S = CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*thick))).dot(CjG.T)
            # St[:,iT]=CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG.T)[0,:]
            St[:,iT] = S[:,0]
        self.z     = thicks
        self.df_Sz = list(St)
        return St

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
        self.set_name(name,self.path)
        self._set_excitation_errors(Smax)
        self._solve_Bloch(opts,Vopt0,v)
        self._set_Vg()
        self.set_thickness(thick)
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

        hkl = self.df_G[['h','k','l']].values
        Sg  = self.df_G.Sw.values
        Ug = self.Fhkl/self.crys.volume

        # Ug[U0_idx] = U0
        # and Ug(iG,jG) are obtained from Ug[h,k,l] where h,k,l = hlk_iG-hkl_jG
        U0_idx = [2*self.Nmax]*3
        if Vopt0 : Ug[tuple(U0_idx)] = 0   #setting average potential to 0

        if v:print(colors.blue+'...assembling %dx%d matrix...' %((Sg.shape[0],)*2)+colors.black)
        H = np.diag(Sg+0J)
        for iG,hkl_G in enumerate(hkl) :
            U_iG = np.array([Ug[tuple(hkl_J+U0_idx)] for hkl_J in hkl_G-hkl]) #;print(V_iG.shape)
            # print('G : ',hkl_G)
            # print('U_G:',Ug[tuple(U0_idx+hkl_G)])
            # print('idx iG:',[tuple(hkl_J+U0_idx) for hkl_J in hkl_G-hkl])
            # print(U_iG)
            H[iG,:] += U_iG/(2*self.k0)  #off diagonal terms as potential

        self.H = H
        if 'H' in opts:self.show_H()

        if v:print(colors.blue+'...diagonalization...'+colors.black)
        self.gammaj,self.CjG = np.linalg.eigh(H) #;print(red+'Ek',lk,black);print(wk)

    def _set_Vg(self):
        Vg     = self.Fhkl/self.crys.volume
        V0_idx = [2*self.Nmax]*3
        hkl    = self.df_G[['h','k','l']].values
        Vg_G   = np.array([ Vg[tuple(hkl_G+V0_idx)] for hkl_G in hkl])
        px,py,e0x = mut.project_beams(K=self.K,qxyz=self.get_G(),e0=[1,0,0],v=1)
        self.e0x = e0x
        self.df_G['px'] = px
        self.df_G['py'] = py
        self.df_G['Vg'] = Vg_G
        self.df_G['L']  = np.ones(Vg_G.shape)

    def _set_intensities(self):
        '''get beam intensities at thickness'''
        id0 = self.is_hkl([0,0,0],v=0)
        gammaj,CjG = self.gammaj,self.CjG
        S = CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*self.thick))).dot(CjG.T)
        # S = S[:,id0]
        S = S[id0,:]
        self.df_G['S'] = S
        self.df_G['I'] = np.abs(S)**2

    ################################################################################
    #### getter
    ################################################################################
    def get_intensities(self):return self.df_G.I
    def get_hkl(self):return self.df_G[['h','k','l']].values
    def get_kin(self):return self.df_G[['h','k','l','Sw','Vg','S','I']]
    def get_G(self):return self.df_G[['qx','qy','qz']].values
    def get_Istrong(self,m=100):
        Imax = self.df_G.I.max()
        print(self.df_G[['h','k','l','Sw','Vg','I']].loc[self.df_G.I>Imax/m])
    def get_Sw(self,Smax=1):
        Smax = min(Smax,self.Smax)
        print(self.df_G[['h','k','l','Sw']].loc[self.df_G.Sw<=Smax])
    def is_hkl(self,Ghkl,v=1):
        Gidx = np.where(np.linalg.norm(self.get_hkl() - Ghkl,axis=1)==0)[0]
        if v:print(self.df_G.iloc[Gidx])
        return Gidx[0]

    ################################################################################
    #### display
    ################################################################################
    def show_beams(self,F='I',fopts='m',opts='xN',mag=500,cmap='Greens',**kwargs):
        '''Display beam values
        - F : 'I'(intensity),'S'(scattered beams),'Vg'(potential)
        - opts  : 'x'(include projection of x axis), 'N'(normalize)
        - fopts : see get_fz
        '''
        fz,fz_str = get_fz(fopts)
        F_str = {'I':'Intensity $I_{g}$','Vg':'Potential $V_{g}$','S':'Scattered beams $S_{g}$','L':'Lattice','Sw':'Excitation error $\zeta_{g}$'}[F]

        tle = '%s, %s, thickness=%d$\AA$'  %(F_str,fz_str,self.thick)
        px,py = self.df_G[['px','py']].values.T
        Fvals = fz(self.df_G[F])
        if 'N' in opts:Fvals/=Fvals.max()
        # mag /= Fvals.max()
        scat  = [px,py,Fvals*mag/Fvals.max(),Fvals]
        plts = [[0,0,'b+','']]
        if 'x' in opts:plts+=[ [[0,self.e0x],[0,0],'r','']]
        # print(Fvals.max())
        fig,ax = dsp.stddisp(plts,lw=2,scat=scat,caxis=[0,Fvals.max()],cmap=cmap,cs='S',title=tle,**kwargs)

    def show_Fhkl(self,s=None,opts='m',h3D=0,**kwargs):
        '''Displays structure factor over grid
        - opts : see get_fz
        - s : slice or str('k=0' => Fhkl(k=0)) or int('l==<int>')
        '''
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



################################################################################
#### functions
################################################################################
def load(path='',tag='',file='',):
    '''load a saved Bloch object
    filename : pickle file (.pkl)  '''
    files = [f for f in glob.glob(os.path.join(path,'*.pkl')) if tag in f]
    if tag:
        if len(files)==1:
            file = files[0]
    if file:
        with open(file,'rb') as f : obj = pickle5.load(f)
        print(colors.green+'loaded:' +colors.yellow+file+colors.black);
        return obj
    else:
        if len(files):
            print(colors.red+'mulitple simus available with tag %s:' %tag);
            print(colors.yellow,[os.path.basename(f) for f in files],colors.black)
        else:
            print(colors.red+'no simus available with in "%s" with tag %s' %(path,tag))

def get_fz(opts):
    '''mapping for function to apply
    - opts : 'r'(real) 'i'(imag) 'a'(angle) 'm'(mag) 'l'(log10(|Fhkl|+1)) '2'(mag^2) 'L'(logM)
    '''
    if not opts:opts='m'
    keys   =['r','i','a','m','l','2','L']
    fs     = [np.real,np.imag,np.angle,np.abs,logF,abs2,logM]
    fs_str = ['real part','imag part','phase','magnitude','$\log_{10}(|F|+1)$','$|F|^2$','$-\log_{10}$']
    fz     = dict(zip(keys,fs))[opts]
    fz_str = dict(zip(keys,fs_str))[opts]
    return fz,fz_str

abs2 = lambda F:np.abs(F)**2
logF = lambda F:np.log10(np.abs(F)+1)
logM = lambda F:-np.log10(np.maximum(F,1e-5))
