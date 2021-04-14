import importlib as imp
import pickle,matplotlib,os
import numpy as np, pandas as pd
import scipy.fftpack as fft
from scipy.integrate import nquad,trapz,quad
import utils.displayStandards as dsp
import utils.physicsConstants as cst
import utils.glob_colors as colors

class Multislice:
    '''python multislice
    #### Structure
    - pattern : 2D or 3D pattern (see multi2D or multi3D)
    #### Parameters
    - keV   : float - electron wavelength (keV)
    - tilt  : float 2-tuple : beam tilt (degrees)
    - dz    : float - slice thickness
    - nz    : int - number of slices along propagation axis
    - Nz    : int - number of unit cells along propagation axis (Nz takes preference over nz if set)
    - nx    : int or 2-tuple - used sampling if defined
    - Nx    : int or 2-tuple - supercell size
    #### TDS
    - TDS    : Use thermal diffuse scattering
    - nTDS   : nb configurations
    - wobble : wobble parameters (see Multislice._wobble for more info  )
    #### aliases
    - ndeg         : alias for nx
    - slice_thick  : alias for dz
    #### Save/Verbose options
    - opts : str - 'q'(save reciprocal space slices) 'x'(save real space slices)
    - v    : bool or str - verbose option
    - iZs  : int - info are saved every iZs slices
    - iZv  : int - info are displayed every iZv slices
    #### Display options
    - ppopt : TVPQXBZY (see Multislice.display for more info)
    #### Misc
    - copt : limit bandwidth on propagator option
    - sg   : change the sign of the propagator
    - eps  : scale the strength of the potential

    '''
    def __init__(self,
        pattern,
        keV=200,tilt=0,dz=1,slice_thick=None,nz=0,Nz=None,nx=2**10,Nx=1,
        TDS=False,nTDS=8,ndeg=None,wobble=0.05,
        copt=1,eps=1,sg=-1,
        iZs=1,s_opts='q',opts=None,iZv=1,v=1,
        ppopt='',name='./unknown',**kwargs):
        self.Mversion = 1.1
        self.pattern  = pattern
        self.keV  = keV
        self.dz   = self._get_thickness(dz,slice_thick)
        self.nx   = self._alias(nx,ndeg)
        self.Nx   = Nx
        self.Nz   = 0 # incremented during propagation
        self.nz   = 0 # incremented during propagation
        self.tilt = tilt
        #energy dependent
        self.lam  = cst.keV2lam(keV)
        self.sig  = cst.keV2sigma(keV)
        self.k0   = 1/self.lam
        self.eps  = eps
        self.copt = copt
        self.sg   = sg
        #TDS
        self.TDS    = TDS
        self.nTDS   = nTDS
        self.wobble = self._wobble(wobble)
        #Misc

        self._set_name(name)
        self._set_ns(v)
        nz     = self._get_nz(nz,Nz)
        s_opts = self._alias(s_opts,opts)
        self.run(nz,iZs,iZv,s_opts,v)
        if ppopt:
            self.display(ppopt,self.basename,**kwargs)


    ############################################################################
    #### parameters
    def _get_nz(self,nz,Nz):
        if Nz:nz=int(Nz*self.ns)
        return nz
    def _get_thickness(self,dz,slice_thick):
        if slice_thick:dz=slice_thick
        return dz
    def _alias(self,param,alias):
        if alias : param=alias
        return param
    def _wobble(self,wobble):
        if isinstance(wobble,int) or isinstance(wobble,float): wobble = [wobble]*5
        if isinstance(wobble,dict) :
            wobbles = [0]*5
            for k,v in wobble.items():wobbles[k]=v
            wobble = wobbles
        return np.array(wobble)

    ############################################################################
    #### Misc
    def _set_name(self,name):
        self.name = os.path.basename(name)
        self.path = os.path.realpath(os.path.dirname(name))+'/'
        self.fullname = self.path+self.name

    def _set_ns(self,v):
        self.ns = int(np.round(self.ez/self.dz))
        self.dz = self.ez/self.ns
        if v:
            print('Slice thickness and number of slices per cell')
            print('     lattice param z=%.2fA, dz=%.2fA, nzs=%d' %(self.ez,self.dz,self.ns))

    def update_z(self,nz):
        self.Nz += int(nz/self.ns)
        self.nz += nz
        self.z = np.hstack([self.z,np.arange(self.nz)*self.dz])
        self.update(nz)

    def save(self,file):
        with open(file,'wb') as out :
            pickle.dump(self, out, pickle.HIGHEST_PROTOCOL)
        print(colors.green+"object saved : \n"+colors.yellow+file+colors.black)

    ############################################################################
    #### Display
    def display(self,ppopt,name='',**kwargs):
        ''' ppopt options :
        T(Transmission)P(propagator)Q(Psi_q)X(Psi_x)B(beams)Z(Psi_xz)Y(Psi_qz)
        '''
        if 'T' in ppopt:self.Tz_show( name=name+'Tz.svg' ,**kwargs)
        if 'Z' in ppopt:self.Za_show( name=name+'Za.png' ,**kwargs)
        if 'P' in ppopt:self.Pq_show( name=name+'Pq.svg' ,**kwargs)
        if 'Q' in ppopt:self.Qz_show( name=name+'Qz.svg' ,**kwargs)
        if 'X' in ppopt:self.Xz_show( name=name+'Xz.svg' ,**kwargs)
        if 'B' in ppopt:self.Bz_show( name=name+'Bz.svg' ,**kwargs)
        if 'Z' in ppopt:self.Xxz_show(name=name+'Xxz.png',**kwargs)
        if 'Y' in ppopt:self.Qxz_show(name=name+'Xxz.png',**kwargs)


##################################################################
#### misc functions
##################################################################
def load(filename):
    '''load a saved object'''
    with open(filename,'rb') as f : multi = pickle.load(f)
    return multi

def tilts_show(tilts,mp2,iBs,iZs,**kwargs):
    ''' display beams for a sequence of tilted simulations
    tilts,mp2 : tilts array and multislice objects
    - iBs,iZs : beam and slice indices (integer allowed)
    '''
    bopt,zopt,nz=isinstance(iBs,int),isinstance(iZs,int),1
    if bopt:iBs=np.array([iBs])
    if zopt:iZs=np.array([iZs])
    if isinstance(iZs,slice):
        stop = (mp2[0].z.size+1)*(iZs.stop<0)+iZs.stop #;print(stop)
        nz = int(np.ceil((stop-iZs.start)/iZs.step))   #;print(nz)
    else:
        nz = np.array(iZs).size
    Itbz = np.zeros((tilts.size,iBs.size,nz))
    for i,t in enumerate(tilts):
        Itbz[i,:,:] = mp2[i].getB(iBs)[iZs,:].T

    if zopt and bopt:
        plts=[tilts,Itbz[:,0,iZs],'bo-']
        dsp.stddisp(plts,labs=[r'$\theta(deg)$','$I$'],**kwargs)
    if bopt and not zopt:
        csZ = dsp.getCs('Reds',nz)
        plts1 = [[tilts,Itbz[:,0,iZ],[csZ[iZ],'o-'],'%dnm' %(z0/10)] for iZ,z0 in enumerate(mp2[0].z[iZs])]
        dsp.stddisp(plts1,labs=[r'$\theta(deg)$','$I$'],**kwargs)
    if zopt and not bopt:
        cs = dsp.getCs('Spectral',iBs.size)
        # for iZ in iZs:Iz[:,iZ]/=Iz[:,iZ].max()
        plts2=[ [tilts,Itbz[:,i],[cs[i],'o-'],'iB=%d' %(iB)] for i,iB in enumerate(iBs)]
        dsp.stddisp(plts2,labs=[r'$\theta(deg)$','$I$'],**kwargs)
