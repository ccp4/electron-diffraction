import numpy as np
from subprocess import Popen
import utils.physicsConstants as cst
import utils.displayStandards as dsp
import utils.glob_colors as colors

nearBragg_bin=dsp.get_figpath(__file__,'/../nearBragg/bin/')+'nearBragg' #;print(nearBragg_bin)
dtype = np.float128

Cs = [colors.unicolor(0.75),(1,0,0),colors.unicolor(0.25),(0,0,1),(1,0,0),(1,1,0)]
Zs = [1,3,6,7,8,16]
# Ai = np.array([1,2,2,2,3])*0.25
Ai = np.array([0.1,0.25,0.26,0.27,1.5])
# Ai =(10-0.5*np.arange(len(Zs)))/8
# self.fj = lambda ti,j:np.sqrt(np.pi)/A[i]*np.exp(-(np.pi*Ai[j]*q)**2)

class NearBragg():
    ''' Near Bragg (distances in A)
    - pattern : x,z,f (in fractional coordinates)
    - keV  : energy(keV) prevails over lam
    - lam  : wavelength(A)
    - Nx,Nz : number of unit cells in each direction
    ## detector : \n
    - z0   : float - distance of detector to crystal origin
    - npx  : float - number of pixels
    - tmax : float - max half angle (deg)
    - qmax : float - max resolution A^-1  (prevails on tmax if defined)
    OR
    - q0s  : np.ndarray - reciprocal space sampling(A^-1)
    ## misc :\n
    - method :
        - Greens2    : Double scattering
        - Greens     : single scattering exact distance calculations
        - Fresnel    : single scattering Fresnel regime
        - Fraunhofer : Single scattering Fraunhofer regime
        - Holton : call Holton code
    - fjopt : with(1) or without(0) form factor (default 1)
    '''
    def __init__(self,pattern,ax,bz,keV=200,lam=None,path='',
            Nx=1,Nz=1,eps=1,
            z0=1e10,npx=4096,tmax=5,qmax=None,q0s=None,
            method='D',fjopt=1,iZv=5):
        if keV : lam = cst.keV2lam(keV)
        self.pattern = pattern
        self.lam  = dtype(lam)

        #initialization
        self.replicate(ax,bz,Nx,Nz)
        if isinstance(q0s,np.ndarray) :
            self.set_pixels(q0s,z0)
        else :
            self.set_detector(npx,z0,tmax,qmax)
        self.I = np.zeros(self.x0s.shape,dtype=dtype)

        #form factor
        if fjopt==1 :
            self.fj = lambda ti,j:eps*np.sqrt(np.pi)/Ai[j]*np.exp(-(np.pi*Ai[j]*(np.sin(ti)/self.lam))**2)
        elif fjopt==0 :
            #no form factor for comparison with Holton code
            self.fj = lambda ti,j:np.ones(ti.shape)
        #compute
        if   method=='Greens2' or  method=='D' : self._Greens2(iZv)
        elif method=='Proba'   or  method=='P' : self._Proba(iZv)
        elif method=='Greens'  or  method=='G' : self._Greens()
        elif method=='Fresnel'                 : self._Fresnel()
        elif method=='Fraunhofer'              : self._Fraunhofer()
        elif method=='Holton'                  : self._Holton(path=path)

    ################################################################
    # init
    ################################################################
    def replicate(self,ax,bz,Nx,Nz):
        '''replicate unit cell pattern Nx,Nz times'''
        self.Nx,self.Nz = Nx,Nz
        self.ax,self.bz = ax,bz
        x,z,Za = self.pattern
        #arange
        nx = np.arange(-int(Nx/2),int(np.ceil(Nx/2)))
        nx,nz = np.meshgrid(nx,np.arange(Nz))
        nx,nz = nx.flatten(),nz.flatten()
        #replicate
        x  = np.hstack([ (x+i)*self.ax for i in nx])
        z  = np.hstack([ (z+i)*self.bz for i in nz])
        Za = np.hstack([Za]*self.Nx*self.Nz)
        #sort
        idx = np.argsort(np.array(z))
        x,z,Za = x[idx],z[idx],Za[idx]
        #store
        self.x  = np.array(x,dtype=dtype)
        self.z  = np.array(z,dtype=dtype)
        self.Za = np.array(Za,dtype=int)

    def set_detector(self,npx,z0,tmax,qmax=None):
        ''' compute pixel positions and sizes from max scattering angle and distance to sample
        '''
        if qmax:tmax=qmax*self.lam
        self.npx  = npx
        self.z0   = dtype(z0)
        self.tmax  = tmax
        self.x_max = self.tmax*self.z0
        self.px    = self.tmax*self.z0*cst.A/cst.mm/(self.npx-1)

        x0s = np.linspace(-self.x_max,self.x_max,self.npx)
        self.x0s = np.array(x0s,dtype=dtype)+(x0s[1]-x0s[0])/2
        self.q0s = np.array(x0s/self.z0s/self.lam,dtype=dtype)

    def set_pixels(self,q0s,z0):
        self.npx  = q0s.size
        self.tmax = q0s.max()*self.lam
        self.z0   = dtype(z0)
        self.q0s = np.array(q0s,dtype=dtype)
        self.x0s = np.array(self.z0*self.lam*q0s,dtype=dtype)
        self.x_max = self.x0s.max()

    ################################################################
    # display
    ################################################################
    def Pattern_show(self,**kwargs):
        '''Show diffraction pattern'''
        C = [Cs[Z] for Z in self.Za]
        Nx0,Nx1 = int(self.Nx/2),int(np.ceil(self.Nx/2))
        dsp.stddisp(scat=[self.z,self.x,C],
            ms=20,labs=[r'$z(\AA)$',r'$x(\AA)$'],
            xyTicks=[self.bz,self.ax],xylims=[0,self.Nz*self.bz,-Nx0*self.ax,Nx1*self.ax],
            **kwargs)

    # def stat_show(self,**kwargs):
    #     '''show statistics'''
    #     z = self.bz*np.arange(self.Nz+1)/10
    #     cs = dsp.getCs('Spectral',2)
    #     S = np.array([np.cumsum(Si) for Si in self.S])
    #     plts  = [[z,Si,cs[i],'%d' %i] for i,Si in enumerate(S)]
    #     plts += [[z,S.sum(axis=0),'k--']]
    #     return dsp.stddisp(plts,labs=['$z(nm)$','$Proba$'],
    #         **kwargs)

    def F_show(self,qopt=0,**kwargs):
        '''Show atomic form factors'''
        Za = np.unique(self.Za)
        t = self.getAngle()*np.pi/180
        if qopt:
            labx,x = r'$q(\AA^{-1})$',self.getQ()
        else:
            labx,x = r'$\theta (deg)$', self.getAngle()

        plts = [[x,self.fj(t,Z)**2,'b','$Z_a=%d$' %Z] for Z in Za]
        plts += [[x,self.fj(t,Z),'g--','$Z_a=%d$' %Z] for Z in Za]
        dsp.stddisp(plts,labs=[labx,r'$f_j^2(\AA)$'],**kwargs)


    def getAngle(self):return (self.x0s/self.z0)*180/np.pi
    def getQ(self):return (self.x0s/self.z0)/self.lam
    def getI(self):return self.I

    ################################################################
    # utils to run James Holton .c code
    ################################################################
    def _save_input(self,file):
        np.savetxt(file,np.vstack([self.z,np.zeros(self.x.shape),self.x,
            np.ones(self.x.shape),np.zeros(self.x.shape),np.zeros(self.x.shape)]).T)
        print(colors.yellow+file+colors.black)

    def _cmd(self,opts='s',path='',file='atoms.txt'):
        self.path = path
        self.file = file
        if 's' in opts :self._save_input(path+self.file)
        #cmd
        cmd=" %s -file %s " %(nearBragg_bin,path+self.file)
        cmd+='-lambda %f ' %(self.lam)
        cmd+='-distance %f -source_dist far ' %(self.z0*cst.A/cst.mm)
        cmd+= '-pixel %.15f -detpixels_x %d ' %(2*self.px,self.npx)
        cmd+='-detpixels_y 1 -oversample 1 '
        cmd+='-intfile %s -floatfile %s -Ifilename %s ' %(path+'intimage.img',path+'floatimage.bin',path+'I.txt')
        # cmd+='-debug -printout'
        print(colors.red+cmd+colors.black)
        if 'r' in opts:
            print(colors.green+'.........Near Bragg.......'+colors.black)
            p = Popen(cmd,shell=True);p.wait()
            self.I = np.loadtxt(path+'I.txt')
            # print(colors.green+'Near Bragg'+colors.black)

    def _Holton(self,path=''):
        '''Run James Holton code'''
        self._cmd(opts='sr',path=path,file='atoms.txt')

    ################################################################
    # Single scattering routines
    ################################################################
    def _Fraunhofer(self):
        print(colors.green+'... Running nearBragg Fraunhofer ...'+colors.black)
        for i in range(self.npx) :
            tij  = np.abs(self.x0s[i]-self.x)/(self.z0-self.z)
            R_ij = self.z + np.sqrt((self.x0s[i]-self.x)**2+(self.z0-self.z)**2)
            Rij  = self.x0s[i]*self.x/self.z0
            self.I[i] = np.abs(np.sum(
                self.fj(tij,self.Za)*np.exp(2*np.pi*1J*Rij/self.lam)/(R_ij*cst.A)))**2

    def _Fresnel(self):
        print(colors.green+'... Running nearBragg Fresnel ...'+colors.black)
        for i in range(self.npx) :
            tij  = np.abs(self.x0s[i]-self.x)/(self.z0-self.z)
            R_ij = self.z + np.sqrt((self.x0s[i]-self.x)**2+(self.z0-self.z)**2)
            Rij  = (self.x0s[i]-self.x)**2/(2*(self.z0-self.z))
            self.I[i] = np.abs(np.sum(
                self.fj(tij,self.Za)*np.exp(2*np.pi*1J*Rij/self.lam)/(R_ij*cst.A)))**2

    # def _Dual_beams(self):
    #     print(colors.green+'... Running nearBragg Dual Beams ...'+colors.black)
    #     natoms = self.z.size
    #     self.A = np.zeros((self.x0s.size),dtype=np.complex256)
    #     for i in range(natoms) :
    #         if not i%100 : print('%d/%d' %(i,natoms))
    #         #single scattering
    #         tij0  = np.abs(self.x0s-self.x[i])/(self.z0-self.z[i])
    #         Rij0 = self.z[i] + np.sqrt((self.x0s-self.x[i])**2+(self.z0-self.z[i])**2)
    #         self.A += self.fj(tij0,self.Za[i])*np.exp(2*np.pi*1J*Rij0/self.lam)/(Rij0*cst.A)
    #         #double scattering
    #         for j in range(i):
    #             tij = np.arctan((self.x[i]-self.x[j])/(self.z[i]-self.z[j]))
    #             Rij = np.sqrt((self.x[j]-self.x[i])**2+(self.z[j]-self.z[i])**2)
    #             fij = self.fj(tij,self.Za[j]) #/np.sqrt(Rij)
    #             self.A += fij*self.fj(tij0-tij,self.Za[i])*np.exp(2*np.pi*1J*(Rij0+Rij)/self.lam)/(Rij0*cst.A)
    #     self.I = np.abs(self.A)**2

    def _Greens(self):
        print(colors.green+'... Running nearBragg Greens ...'+colors.black)
        for i in range(self.npx) :
            tij  = np.abs(self.x0s[i]-self.x)/(self.z0-self.z)
            R_ij = self.z + np.sqrt((self.x0s[i]-self.x)**2+(self.z0-self.z)**2)
            self.I[i] = np.abs(np.sum(
                self.fj(tij,self.Za)*np.exp(2*np.pi*1J*R_ij/self.lam)/(R_ij*cst.A)))**2

    #################################################################
    #### Multiple scattering routines
    #################################################################
    # def _Greens2n(self,iZv=5):
    #     '''2-level dynamical scattering ignoring backward scattering'''
    #     print(colors.green+'... Running nearBragg Greens2 ...'+colors.black)
    #
    #     ## Primary scattering
    #     ##atoms x pixels
    #     xjp,zjp = self.x0s-self.x[:,None]),self.z0-self.z[:,None]
    #     tjp = np.arctan2(xjp,zjp)
    #     Rjp = np.sqrt(xjp*2+zjp**2)
    #     ej  = np.exp(2J*np.pi*k0*self.z)
    #     Fj  = np.sum(ej*self.fj(tij,self.Za)*np.exp(2J*np.pi*k0*Rjp),axis=1)
    #
    #     ## Secondary scattering
    #     # natoms x natoms
    #     xij,zij,Zi = self.x[:,None]-self.x,self.z[:,None]-self.z,self.Za[:,None]
    #     tij = np.arctan2(xij,zij)
    #     tijp = tjp[:,None] - tij
    #     Rij = np.sqrt(xij*2+zij**2)
    #     Fij = self.fj(tij,self.Za)*self.fj(tijp,Zi)*np.exp(2J*np.pi*k0*Rjp[:,None])

    def _Greens2(self,iZv=5):
        '''2-level dynamical scattering ignoring backward scattering'''
        print(colors.green+'... Running nearBragg Greens2 ...'+colors.black)
        natoms  = self.z.size
        dq      = self.q0s[1]-self.q0s[0]
        # sig_e   = sum(self.fj(self.q0s*self.lam,self.Za[0])**2)*dq
        # print('dq=%.1E Angstrom, sig_e:%.1E Angstrom' %(dq,sig_e))

        I0 = 1
        self.A = np.zeros((self.x0s.size),dtype=complex)
        Adyn = np.zeros((self.x0s.size),dtype=complex)
        # self.S[0,0]=1
        iz=1
        print(colors.yellow+'atom %d/%d slice %d, I0=%.4f' %(1,natoms, iz,I0) +colors.black)
        for i in range(natoms) :
            #### single scattering
            ti0 = (self.x0s-self.x[i])/(self.z0-self.z[i])
            #distance from atom to detector
            Ri0 = np.sqrt((self.x0s-self.x[i])**2+(self.z0-self.z[i])**2)
            # single scattering amplitude
            Akin = self.fj(ti0,self.Za[i])*np.exp(2*np.pi*1J*self.z[i]/self.lam)
            #### double scattering with backward atoms
            ibackward = self.Nx*(iz-1)
            Adyn[:]=0
            for j in range(ibackward):
                #angle between atoms
                tij = np.arctan((self.x[i]-self.x[j])/(self.z[i]-self.z[j]))
                #distance betwwen atoms
                Rij = np.sqrt((self.x[j]-self.x[i])**2+(self.z[j]-self.z[i])**2)
                #amplitude of radiation j at atom i (dimensionless)
                fij = self.fj(tij,self.Za[j])*np.exp(2*np.pi*1J*Rij/self.lam)#*dq
                #scattering amplitude from atom i due to scattering from atom j
                Adyn += fij*self.fj(ti0-tij,self.Za[i])

            # increment slice
            # print(abs(Adyn/Akin).max())
            if not (i+1)%self.Nx:
                # I0 -= self.S[1,iz]/(self.ax*self.Nx) #;print(I0)
                iz+=1
                if not iz%iZv :
                    print(colors.yellow+'atom %d/%d slice %d, I0=%.4f' %(i+1,natoms, iz,I0) +colors.black)

            # atom contribution to diffraction pattern at detector
            # print(np.abs(Akin).max(),np.abs(Adyn).max(),np.abs(Akin-Adyn).min())
            self.A += (Akin+Adyn)*np.exp(2*np.pi*1J*Ri0/self.lam)/(Ri0*cst.A) #*np.sqrt(I0)
        self.I = np.abs(self.A)**2

    # def _Greens2_old(self,iZv=5):
    #     print(colors.green+'... Running nearBragg Greens2 ...'+colors.black)
    #     natoms  = self.z.size
    #     dq      = self.q0s[1]-self.q0s[0]
    #     sig_e   = sum(self.fj(self.q0s*self.lam,self.Za[0])**2)*dq
    #     print('dq=%.1E Angstrom, sig_e:%.1E Angstrom' %(dq,sig_e))
    #
    #     I0 = 1/(self.ax*self.Nx) #1 electron/transverse area/sec
    #     self.A = np.zeros((self.x0s.size),dtype=complex)
    #     self.S = np.zeros((2,self.Nz+1))
    #     self.S[0,0]=1
    #     iz=1
    #     for i in range(natoms) :
    #         if not i%(iZv*self.Nx) : print('atom %d/%d slice %d, I0=%.4f' %(i+1,natoms, iz,I0))
    #         #### single scattering
    #         tij0 = (self.x0s-self.x[i])/(self.z0-self.z[i])
    #         R_ij = self.z[i] + np.sqrt((self.x0s-self.x[i])**2+(self.z0-self.z[i])**2)
    #         Akin = np.sqrt(I0)*self.fj(tij0,self.Za[i])*np.exp(2*np.pi*1J*R_ij/self.lam) #/(R_ij*cst.A)
    #
    #         Skin = (abs(Akin)**2).sum()*dq
    #         assert Skin<1
    #         self.S[0,iz] -= sig_e*I0
    #         self.S[1,iz] += Skin
    #         self.A += Akin
    #         if not (i+1)%self.Nx:
    #             I0 -= self.S[1,iz]/(self.ax*self.Nx) #;print(I0)
    #             iz+=1
    #     self.I = np.abs(self.A)**2

    ##
    def _Proba(self,iZv=5):
        print(colors.green+'... Running nearBragg Proba ...'+colors.black)
        natoms = self.z.size
        dq      = self.q0s[1]-self.q0s[0]
        sig_e   = sum(self.fj(self.q0s*self.lam,self.Za[0])**2)*dq
        self.sig_e = sig_e
        print('dq=%.1E Angstrom, sig_e:%.3f Angstrom' %(dq,sig_e))

        I0 = 1/(self.ax*self.Nx) #1 electron/transverse area/sec
        # self.S = np.zeros((self.Nz+1,self.x0s.size,3)) #scattering areas
        self.I = np.zeros((3,self.Nz+1))
        self.S,Iz = np.zeros((3,self.Nz+1)),np.zeros((self.Nz+1))
        self.I[0,0],self.S[0,0],Iz[0]=1,1,I0
        iz=1
        for i in range(natoms) :
            if not i%(iZv*self.Nx) : print('atom %d/%d slice %d' %(i+1,natoms, iz))
            #### single scattering
            tij0  = (self.x0s-self.x[i])/(self.z0-self.z[i])

            self.S[0,iz] -= sig_e*I0 #fraction electron scattered per unit time
            self.S[1,iz] += np.sum(self.fj(tij0,self.Za[i])**2)*dq*I0
            #### double scattering
            # print(colors.green,i,colors.black)
            ibackward = self.Nx*(iz-1)
            for j in range(ibackward):
                tij = np.arctan((self.x[i]-self.x[j])/(self.z[i]-self.z[j]))
                if abs(tij*180/np.pi)<4 : #and abs(tij)>1e-5:
                    fij = 1#self.fj(tij,self.Za[j])*dq
                    fj  = np.sum((fij*self.fj(tij0-tij,self.Za[i]))**2)*dq
                    jz  = int(j/self.Nx)+1
                    self.S[1,jz] -= fj*self.S[1,jz]/(self.ax*self.Nx)#;Iz[jz]
                    self.S[2,iz] += fj*self.S[1,jz]/(self.ax*self.Nx)

                    # iq  = np.argmin(np.abs(self.q0s-tij/self.lam))
                    # jz  = int(j/self.Nx)+1
                    # self.S[iz,:, 2] += fj*I0
                    # self.S[jz,iq,1] -= fj.sum()*dq*I0
                    # print(colors.red+'%3d' %j+colors.black+' : %.3f,%2d,%.2E,%.2E' %(tij*180/np.pi,iq,fij,fj.sum()*dq))
            if not (i+1)%self.Nx:
                I0 += self.S[0,iz]/(self.ax*self.Nx) #;print(I0)
                Iz[iz] = I0
                self.I[0,iz] = self.I[0,iz-1]+self.S[0,iz]
                self.I[1,iz] = np.sum(self.S[1,:])
                self.I[2,iz] = self.I[2,iz-1]+self.S[2,iz]
                iz+=1


        # self.I[0,:] = np.cumsum(self.S[0,:])
        # self.I[1,:] = np.cumsum(self.S[1,:])
        # self.I[2,:] = np.cumsum(self.S[2,:])
        # self.I[0,:] = 1 - np.cumsum(self.S[:,0,0])
        # self.I[1,:] = np.cumsum(self.S[:,:,1].sum(axis=1))*dq
        # self.I[2,:] = np.cumsum(self.S[:,:,2].sum(axis=1))*dq

####################################################################
# misc
####################################################################
disp_d ={
    'q':('getQ'     ,'k',r'$q(\AA^{-1})$'   ),
    'a':('getAngle' ,'k',r'$\theta(\circ)$' ),
    'I':('getI'     ,'r',r'$I$' ),
}

def display(objs,key_pair,**kwargs):
    return dsp.display_solutions(disp_d,objs,key_pair)
