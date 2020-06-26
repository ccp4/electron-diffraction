import pickle
import numpy as np
import scipy.fftpack as fft
from scipy.integrate import nquad,simps,quad
import utils.displayStandards as dsp
import utils.physicsConstants as cst
import utils.glob_colors as colors
# from scipy.signal import fftconvolve
# from scipy.interpolate import interp1d
# import utils.FourierUtils as fu


class multi2D():
    '''
    - pattern,ax,bz : x,z,f
    - Nx : increase supercell size
    - dz : slice thickness
    - nz : number of slices
    - ppopt:T(Transmission)P(propagator)Q(Psi_q)X(Psi_x)B(beams)Z(Psi_xz)
    - iZs,iZv : frequency of save and verbose
    '''
    def __init__(self,
            pattern,ax,bz,
            keV=200,Nx=1,
            dz=1,nz=1,
            iZs=1,iZv=1,ppopt=''):
        self.pattern = pattern
        self.keV     = keV
        self.Nx      = Nx
        self.dz      = dz
        self.ax      = ax
        self.bz      = bz
        # self.az,self.bx =
        #energy dependent
        self.lam = cst.keV2lam(keV)
        self.sig = cst.keV2sigma(keV)
        self.k0 = 1/self.lam
        #Computations
        self._set_transmission_function()
        self._set_propagator()
        self.set_Psi0()
        self.propagate(nz,iZs,iZv)
        #display
        if 'T' in ppopt:self.Tz_show()
        if 'P' in ppopt:self.Pq_show()
        if 'Q' in ppopt:self.Qz_show()
        if 'X' in ppopt:self.Xz_show()
        if 'B' in ppopt:self.Bz_show()
        if 'Z' in ppopt:self.Xxz_show()

    def save(self,file):
        with open(file,'wb') as out :
            pickle.dump(self, out, pickle.HIGHEST_PROTOCOL)
        print(colors.green+"object saved\n"+colors.yellow+file+colors.black)

    ###################################################################
    ####### display
    def Tz_show(self,**kwargs):
        plts = [[self.x,self.Vz,'g',r'$V_z(kV\AA)$'],
                [self.x,self.T.real,'b','$re(T)$'],
                [self.x,self.T.imag,'r','$im(T)$']]
        return dsp.stddisp(plts,labs=[r'$x(\AA)$',''],title='Projected potential $V_z$, Transmission function $T$',
            **kwargs)

    def Pq_show(self,**kwargs):
        q  = fft.fftshift(self.q)
        Pq = fft.fftshift(self.Pq)
        plts = [[q,Pq.real,'b','$re$'],
                [q,Pq.imag,'r','$im$']]
        return dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$','$P_q$'],
            title='Propagator',
            **kwargs)

    def Xz_show(self,iZs=1,cmap='jet',**kwargs):
        if isinstance(iZs,int):iZs=slice(0,self.z.size,iZs)
        z  = self.z[iZs]
        Px = self.psi_xz[iZs,:]
        cs = dsp.getCs(cmap,z.size)
        plts = [ [self.x,Px[i,:],cs[i]] for i in range(z.size)]
        return dsp.stddisp(plts,labs=[r'$x(\AA)$',r'$|\Psi(x)|^2$'],
            imOpt='hc',caxis=[z.min(),z.max()],cmap=cmap,
            **kwargs)

    def Qz_show(self,iZs=1,cmap='jet',**kwargs):
        if isinstance(iZs,int):iZs=slice(0,self.z.size,iZs)
        z   = self.z[iZs]
        Pqs = self.psi_qz[iZs,:]
        Pqs[:,0]=0# do not show central beam
        # q   = self.q
        # Pqs = [Pqs[i,:] for i in range(z.size)]
        q   = fft.fftshift(self.q)
        Pqs = [fft.fftshift(Pqs[i,:]) for i in range(z.size)]
        cs  = dsp.getCs(cmap,z.size)
        plts = [[q,Pqs[i],cs[i]] for i in range(z.size)]
        return dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$',r'$|\Psi(q)|^2$'],
            imOpt='hc',caxis=[z.min(),z.max()],cmap=cmap,
            **kwargs)

    def Bz_show(self,iBs='O',tol=1e-3,cmap='jet',**kwargs):
        Ib = self.psi_qz #np.abs()**2
        if isinstance(iBs,str):
            N   = int(self.nx/2)
            iHs = fft.fftshift(np.arange(-N,N))
            if not 'a' in iBs:
                Im  = Ib[:,1:].max()    #;print(Im)
                Imax = Ib.max(axis=0)   #;print(Imax)
                iHs = iHs[Imax>Im*tol]
            iBs=iHs['O' not in iBs:]

        if isinstance(iBs,list):iBs=np.array(iBs)

        Ib   = Ib[:,iBs]
        h    = ['%d_{%d}' %(i/self.Nx,i%self.Nx) for i in iBs]
        cs   = dsp.getCs(cmap,iBs.size)
        plts = [[self.z,Ib[:,i],cs[i],'$%s$' %h[i]] for i,iB in enumerate(iBs)]
        return dsp.stddisp(plts,labs=[r'$z(\AA)$',r'$I_b$'],
        **kwargs)

    def Xxz_show(self,**kwargs):
        x,z = np.meshgrid(self.x,self.z)
        im = [x,z,self.psi_xz]
        return dsp.stddisp(im=im,labs=[r'$x(\AA)$',r'$z(\AA)$'],
        **kwargs)


    ###### Privates
    def _set_transmission_function(self):
        Nx = self.Nx
        x,z,f = self.pattern
        Vz = 1.*np.array([simps(f[i,:],z) for i in range(f.shape[0])])
        T  = np.exp(1J*self.sig*Vz)

        self.x  = np.hstack([x + self.ax*i for i in range(Nx)])
        self.Vz = np.hstack([Vz]*Nx)
        self.T  = np.hstack([T]*Nx)
        self.nx = x.size

    def _set_propagator(self):
        self.dx = self.x[1]-self.x[0]
        self.q  = fft.fftfreq(self.nx,self.dx)
        self.dq = self.q[1]-self.q[0]
        self.Pq = np.exp(1J*np.pi*self.dz*self.q**2/self.k0)
        self.nq = int(1/3*self.nx) #prevent aliasing

    def set_Psi0(self):
        Psi  = np.ones(self.T.shape,dtype=complex)
        self.Psi_x = Psi/np.sqrt(np.sum(np.abs(Psi)**2)*self.dx)
        self.psi_xz = np.zeros((0,self.nx))
        self.psi_qz = np.zeros((0,self.nx))
        self.z  = np.array([])
        self.iz = 0

    def propagate(self,nz,iZs=1,iZv=1):
        nzq,z0 = int(nz/iZs),0
        if self.z.size : z0=self.z.max()
        self.z  = np.hstack([self.z,z0+np.arange(nzq)*self.dz*iZs ])
        self.psi_xz = np.vstack([self.psi_xz,np.zeros((nzq,self.nx))])
        self.psi_qz = np.vstack([self.psi_qz,np.zeros((nzq,self.nx))])
        # self.T=fft.fftshfft(self.T)
        for i in range(nz):
            self.Psi_q = fft.fft(self.T*self.Psi_x)   #periodic assumption
            self.Psi_q[self.nq:-self.nq] = 0          #prevent aliasing
            self.Psi_x = fft.ifft(self.Pq*self.Psi_q)
            # self.Psi_x = fft.fftshift(fft.ifft(self.Pq*self.Psi_q))
            #save and print out
            msg=''
            if not i%iZv:
                Ix2 = np.sum(np.abs(self.Psi_x)**2)*self.dx
                Iq2 = np.sum(np.abs(self.Psi_q/self.nx)**2)/self.dq #parseval's theorem of the DFT
                msg+='i=%-3d I=%.4f, Iq=%.4f ' %(i,Ix2,Iq2)
            if not i%iZs:
                msg+='z=%.1f A' %(self.z[self.iz])
                self.psi_xz[self.iz,:] = np.abs(self.Psi_x)**2
                self.psi_qz[self.iz,:] = np.abs(self.Psi_q)**2
                self.iz+=1
            if msg:print(colors.green+msg+colors.black)

if __name__=='__main__':
    import wallpp.plane_group as pg
    import importlib as imp
    imp.reload(pg)
    ndeg = 2**8
    pptype,a,b,angle = 'p1',20,10,90
    pattern = np.array([[10,5,3]])

    #get the potential from wallpaper library
    p1      = pg.Wallpaper(pptype,a,b,angle,pattern,ndeg=ndeg)
    pattern = p1.get_potential_grid()
    mp1 = multi2D(pattern,a,b,keV=200,
            Nx=1,dz=b,nz=100,ppopt='Z',#XQZTP
            iZs=1,iZv=1)
