import importlib as imp
import pickle,matplotlib,os
import numpy as np, pandas as pd
import scipy.fftpack as fft
# from scipy.integrate import nquad,trapz,quad
import utils.displayStandards as dsp
import utils.physicsConstants as cst
import utils.glob_colors as colors
from . import postprocess as pp     ;imp.reload(pp)
from . import pymultislice          ;imp.reload(pymultislice)

atoms = pd.DataFrame.from_dict(
    { 1:[colors.unicolor(0.9),40],
      6:[colors.unicolor(0.5),60],7:[(0,0,1),70], 8:[(1,0,0),80],
      14:[(1,0.75,0.5),100],16:[(0,1,1),100]},
    orient='index',columns=['c','s'])
Ai = {1:1.0, 6:1.5, 7:1.5, 8:1.5, 14:5.0, 16:2.0}
Zs = {1:0, 6:1, 7:2, 8:3, 14:4, 16:5}
fsplines = os.path.dirname(__file__)+'/data/splines.npy'
rx2,y_vz,b_vz,c_vz,d_vz = np.load(fsplines,allow_pickle=True)
def fv(z,y,x,za,ya,xa,Za):
    r2 = (x-xa)**2+(y-ya)**2#+(z-za)**2
    # V = np.exp(-r2/Ai[Za])
    # return  np.exp(-r2/Ai[Za])
    iz = Zs[Za]
    i  = np.argmin(abs(r2-rx2))
    z  = r2 - rx2[i];
    vz = y_vz[iz][i] + ( b_vz[iz][i] + ( c_vz[iz][i] + d_vz[iz][i] *z ) *z) *z;
    return vz
# fv = lambda z,y,x,za,ya,xa,Za:np.exp(-((np.sqrt(x**2+y**2+z**2)-np.sqrt(xa**2+ya**2+za**2))/Ai[Za]) **2)



class Multi3D(pymultislice.Multislice):
    '''multislice 3D for quick testing
    - pattern  : list of 4d-arrays - [Z,x,y,z]
    - ax,by,cz : lattice constants
    #### beam recording
    - hk    : int or 2-list of list - number/indices of beams to record
    - hkopt : str - 's'(sym) 'r'(repeat) 'o'(include origin)
    #### aliases
    - nxy : alias for nx
    - Nxy : alias for Nx
    #### other
    - **kwargs : see pymultislice.Multislice.__init__
    - s_opts : 's'(save object) 'x'(real space), 'q'(reciprocal space) 't'(transmission function)
    - temsim : bool - same implemetation as temsim if True
    Example :
    Multi3D(pattern,ax,by,cz,
        keV=200,Nx=1,tilt=0,dz=1,nz=1,
        TDS=False,nTDS=8,nx=2**10,wobble=0.05,
        copt=1,eps=1,sg=-1,
        iZs=1,opts='q',iZv=1,v=1,ppopt='')
    '''
    def __init__(self,pattern,ax,by,cz,
        nxy=None,Nxy=None,
        hk=None,hkopt='sr',temsim=True,
        **kwargs):
        self.version = 1.0
        self.ax  = ax
        self.by  = by
        self.cz  = cz
        self.ez  = self.cz
        self.nxy = nxy
        self.Nxy = Nxy
        self.h,self.k = self._set_beams(hk,hkopt)
        self.temsim   = temsim
        print(colors.red+'\t\t 3D multislice simulation '+colors.black)
        super().__init__(pattern,**kwargs)

    def _init_3D_params(self):
        self.nxy  = self._alias3D(self.nx,self.nxy)
        self.Nxy  = self._alias3D(self.Nx,self.Nxy)
        self.tilt = self._alias3D(self.tilt,None)
        self.nx,self.ny     = self.nxy
        self.Nx,self.Ny     = self.Nxy
        self.axby           = np.array([self.ax,self.by])
        self.dx,self.dy     = np.array([self.ax,self.by])*self.Nxy/self.nxy
        self.dqx,self.dqy   = 1/np.array([self.ax,self.by])*self.Nxy
        self.dxy  = np.array([self.dx,self.dy])
        self.dqxy = np.array([self.dqx,self.dqy])
        ## bandwidth limit
        self.q2max = (2/3*min(self.nxy/(2*self.axby)))**2 #prevent aliasing

        self.Psi_x = np.ones(self.nxy)#/np.sqrt(self.nx*self.ny*self.dx*self.dy) #;print((self.Psi_x**2).sum()*self.dx*self.dy)
        self.iz = 0
        self.Nhk   = self.h.size
        self.beams = np.zeros((self.Nhk,0),dtype=complex)
        self.hk = [(h0,k0) for h0,k0 in zip(self.h,self.k)]
        self.z = np.array([])

    def update(self,nz):
        self.beams = np.hstack([ self.beams,np.zeros((self.Nhk,nz),dtype=complex) ])
    def _get_z(self,iz=None):
        if not isinstance(iz,int):
            z = self.z.copy()
        else:
            z = self.z[iz]
        if self.temsim:
            z+= 0.75*self.dz
        return z

    def _alias3D(self,param,alias):
        if alias : param=alias
        if isinstance(param,int):param=[param]*2
        return np.array(param)

    def _set_beams(self,hk,hkopt):
        if isinstance(hk,int):hk=[hk]*2
        if isinstance(hk,list):
            Nh,Nk = hk
            nh,nk = 1,1
            if 'r' in hkopt : nh,nk = self.Nxy
            h,k = np.meshgrid(np.arange(0,Nh+1)*nh,np.arange(0,Nk+1)*nk)
            h,k = h.flatten(),k.flatten()
            if 'o' not in hkopt:
                h,k = h[1:],k[1:]
            if 's' in hkopt:
                hsym,ksym = np.meshgrid(np.arange(-Nh,1)*nh,np.arange(-Nk,1)*nk)
                hsym,ksym = hsym[:-1],ksym[:-1]
                h,k = np.hstack([h,hsym.flatten()]),np.hstack([k,ksym.flatten()])
        return h,k

    def save(self):
        super().save(self.fullname+'_3D.pkl')

    ##################################################################
    ###### Display
    ##################################################################
    def Vz_show(self,iz,**kwargs):
        print(colors.blue+'...Integrating projected potential...'+colors.black)
        Vz = self._projected_potential(iz)
        nx,ny = np.array(self.nxy/self.Nxy,dtype=int)
        x,y   = np.meshgrid(np.arange(nx)*self.dx, np.arange(ny)*self.dy)  #;print(x.shape,y.shape,Vz.shape)
        im    = [x,y,Vz.T]
        Za,xa,ya  = self.pattern[self.Zas[iz]:self.Zas[iz+1],:3].T
        Zatom   = atoms.loc[Za]
        scat    = [xa,ya,Zatom.s,Zatom.c]

        tle = '$V_z$ for slice z=%.2fA, iz=%d' %(self._get_z(iz),iz)
        dsp.stddisp(im=im,scat=scat,labs=['$x$','$y$'],title=tle,
            caxis=[0,Vz.max()],imOpt='c',cs='I',axPos='V',pOpt='pt',#gridOn=False,
            **kwargs)

    def Tz_show(self,iz,opts='i',save=0,load=1,**kwargs):
        '''Plot the transmission function for slice iz :
        - iz : int - slice number
        - opts : str - 'i'(show im(T)) 'r'(show Re(T))
        '''
        T = self._get_transmission_function(iz,1,self.copt,save,load).T
        nx,ny = np.array(self.nxy/self.Nxy,dtype=int)
        x,y   = np.meshgrid(np.arange(nx)*self.dx, np.arange(ny)*self.dy)  #;print(x.shape,y.shape,Vz.shape)
        if 'i' in opts:Tz = T.imag
        if 'r' in opts:Tz = T.real
        im    = [x,y,Tz]

        Za,xa,ya  = self.pattern[self.Zas[iz]:self.Zas[iz+1],:3].T
        Zatom   = atoms.loc[Za]
        scat    = [xa,ya,Zatom.s,Zatom.c]

        tle = ' $%s(T)$ for slice z=%.2fA, iz=%d' %(['Re','Im']['i' in opts],self._get_z(iz),iz)
        dsp.stddisp(im=im,scat=scat,labs=['$x$','$y$'],title=tle,
            caxis=[0,Tz.max()],pOpt='im',
            **kwargs)

    def Pq_show(self,opts='qri',**kwargs):
        P = self.Pq
        Px,Py = P[0,:],P[:,0]
        if 'q' in opts:
            x = np.arange(self.nx)
            plts = [[x,Px.real,'b','Re Px'],[x,Px.imag,'r','Im Px']]
            plts+= [[x,Py.real,'c','Re Py'],[x,Py.imag,'m','Im Py']]
            dsp.stddisp(plts,labs=['q','Px,Py'])

        if 'r' in opts:dsp.stddisp(im=[P.real],pOpt='im',title='re P',**kwargs)
        if 'i' in opts:dsp.stddisp(im=[P.imag],pOpt='im',title='im P',**kwargs)
        return Px,Py

    def Bz_show(self,iBs=slice(1,None,None),**kwargs):
        if isinstance(iBs,list):
            if isinstance(iBs[0],tuple):
                iBs = [i for i,hk0 in enumerate(self.hk) if hk0 in iBs]

        beams = self.beams[iBs,:]
        hk      = np.array(['(%d_%d,%d_%d)' %(int(h0/self.Nx),h0%self.Nx,int(k0/self.Ny),k0%self.Ny ) for h0,k0 in zip(self.h,self.k)])[iBs]
        t       = self._get_z()
        re,im   = np.real(beams),np.imag(beams)
        Ib      = np.abs(beams)**2
        beams   = [hk,t,re,im,Ib]
        return pp.plot_beam_thickness(beams,**kwargs)



    ####################################################################################################################################
    ###### main computations
    ####################################################################################################################################
    def run(self,nz,iZs,iZv,s_opts,v):
        self._init_3D_params()
        self._propagator(v)
        self._sort_atoms(v)
        self.propagate(nz,iZs,iZv,s_opts,v)
        if 's' in s_opts:self.save()

    def _propagator(self,v=1):
        if v:print(colors.blue+'...Setting propagator...'+colors.black)
        sg,copt = self.sg,self.copt
        qx,qy = np.meshgrid(fft.fftfreq(self.nx,self.dx), fft.fftfreq(self.ny,self.dy))
        q2 = qx**2+qy**2
        # if self.tilt: #isinstance(tilt,np.ndarray) or isinstance(tilt,list):
        #     kx = self.k0*np.sin(self.tilt*np.pi/180)
        #     self.Pq = np.exp(sg*1J*np.pi*self.dz*(self.q+kx)**2/self.k0)
        # else:
        self.Pq = np.exp(self.sg*1J*np.pi*self.dz*q2/self.k0)
        self.bw_mask = q2>self.q2max
        if copt:self.Pq[self.bw_mask] = 0

    def _sort_atoms(self,v=1):
        if v:print(colors.blue+'...Sorting atoms per slice...'+colors.black)
        self.pattern = self.pattern[np.argsort(self.pattern[:,3]),:]
        #replicate one time along z
        pattern2 = self.pattern.copy(); pattern2[:,3]   += self.ez
        self.pattern = np.vstack([self.pattern,pattern2])
        za  = self.pattern[:,3]
        dz  = self.dz
        if self.temsim:
            if v>1:print('TEMSIM sorting ')
            self.zis = np.arange(-0.25*dz,self.ez+dz,dz)
            self.Zas = np.cumsum(np.histogram(za,self.zis)[0]) #;print(self.Zas)
            self.Zas = np.hstack([0,self.Zas])
        else:
            self.zis = np.arange(-dz,self.ez+1.1*dz,dz)
            self.Zas = np.cumsum(np.histogram(za,self.zis)[0]) #;print(self.Zas)

    def _projected_potential(self,iz):
        dx,dy,dz,Zas = self.dx,self.dy,self.dz,self.Zas
        nx,ny = np.array(self.nxy/self.Nxy,dtype=int) #unit cell nx,ny
        zi = iz*self.dz
        Vz = np.zeros((nx,ny))
        # print(iz,Zas,Zas.shape)
        for ia in range(Zas[iz],Zas[iz+1]):
            # print(ia)
            Za,xa,ya,za = self.pattern[ia,:4]
            Za = int(Za)
            nax,nay = int(Ai[Za]/dx),int(Ai[Za]/dy)               #;print('nax,nay',nax,nay)
            ixa,iya = int(xa/dx),int(ya/dy)
            for ix in ixa+np.arange(-nax,nax+1):
                x = dx*ix
                for iy in iya+np.arange(-nay,nay+1):
                    y = dy*iy
                    Vz[ix%nx,iy%ny] += fv(zi+dz/2,y,x,za,ya,xa,Za)
        return Vz

    def _transmission_function(self,iz=None,v=0,copt=1,save=0,load=0):
        if v:print(colors.blue+'...Integrating projected potential...'+colors.black)
        Vz = self._projected_potential(iz)
        Tz = np.exp(1J*self.sig*self.eps*Vz*1e-3)
        T  = np.zeros(self.nxy,dtype=complex)
        nx,ny = Tz.shape*self.Nxy
        T[:nx,:ny] = np.tile(Tz,self.Nxy)
        if copt:               #prevent aliasing
            if v>1:print(colors.blue+'...bandwidth limit transmission...'+colors.black)
            T = fft.fft2(T)
            T[self.bw_mask] = 0
            T = fft.ifft2(T)
        if save:
            filename =self.fullname+'_T%s.npy' %str(iz).zfill(3)
            np.save(filename,T)
            print(colors.green+'transmission saved : '+colors.yellow+filename+colors.black)
        return T

    def _get_transmission_function(self,iz=None,v=0,copt=1,save=0,load=0):
        if iz>self.ns or load:
            izl = iz%self.ns
            if not izl and iz>0:izl=self.ns
            filename =self.fullname+'_T%s.npy' %str(izl).zfill(3)
            try:
                if v>1:print(colors.green+'iz=%d,is=%d, loading ' %(iz,iz%self.ns)+colors.yellow+filename+colors.black)
                T = np.load(filename)
            except:
                T = self._transmission_function(iz,v,copt,save,load)
        else:
            T = self._transmission_function(iz,v,copt,save,load)
        return T

    def propagate(self,nz,iZs,iZv,s_opts,v):
        self.update_z(nz)
        if v:print(colors.green+'...Starting propagation loop...'+colors.black)
        save = 't' in s_opts
        for i in range(nz):
            T = self._get_transmission_function(self.iz,v,self.copt,save,'l' in s_opts)
            if v>1:print(colors.blue+'...FFT...'+colors.black)
            self.Psi_q = fft.fft2(T*self.Psi_x) #periodic assumption
            self.Psi_x = fft.ifft2(self.Pq*self.Psi_q)

            msg=''
            if v and (not i%iZv or i==nz-1):
                Ix2 = np.sum(np.abs(self.Psi_x)**2)/(self.nx*self.ny) #*self.dx*self.dy
                Iq2 = np.sum(np.abs(self.Psi_q/(self.nx*self.ny))**2)#/(self.dqx*self.dqy) #parseval's theorem of the DFT
                msg+='i=%-4d,z=%-7.3f A, I=%.4f, Iq=%.4f ' %(i,self._get_z(self.iz),Ix2,Iq2)

            # if not i%iZs and s_opts :
            self.beams[:,self.iz] = self.Psi_q[self.h,self.k]/(self.nx*self.ny)
                # if 'T' in s_opts:
                #     np.save('Tz%d.npy' %i,T)
                # if 'x' in s_opts:
                #     np.save('psi2%d_x.npy' %i,np.abs(self.Psi_x)**2)
                # if 'q' in s_opts:
                #     np.save('psi2%d_q.npy' %i,np.abs(self.Psi_q)**2)
            if msg:print(colors.blue+msg+colors.black)
            self.iz+=1


if __name__=='__main__':
    print('run tests/multislice/base_3D.py to see an example')
