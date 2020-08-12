import numpy as np
from subprocess import Popen
import utils.physicsConstants as cst
import utils.displayStandards as dsp
import utils.glob_colors as colors

nearBragg_bin=dsp.get_figpath(__file__,'/../nearBragg/bin/')+'nearBragg' #;print(nearBragg_bin)
dtype = np.float64

Cs = [colors.unicolor(0.75),(1,0,0),colors.unicolor(0.25),(0,0,1),(1,0,0),(1,1,0)]
Zs = [1,3,6,7,8,16]
Ai = np.array([1,2,2,2,3])*0.25
# Ai =(10-0.5*np.arange(len(Zs)))/8
# self.fj = lambda ti,j:np.sqrt(np.pi)/A[i]*np.exp(-(np.pi*Ai[j]*q)**2)

class NearBragg():
    ''' Near Bragg (distances in A)
    - pattern : x,z,f (in fractional coordinates)
    - keV  : energy(keV) prevails over lam
    - lam  : wavelength(A)
    - Nx,Nz : number of unit cells in each direction
    ## detector : \n
    - npx  : number of pixels
    - z0   : distance of detector to crystal origin
    - tmax : max half angle (deg)
    - qmax : max resolution A^-1 (sets tmax if defined)
    ## misc :\n
    - method : Greens,Fresnel,Fraunhofer,Holton
    - fjopt : 0,1,2 option for form factor
    '''
    def __init__(self,pattern,ax,bz,keV=200,lam=None,path='',
            Nx=1,Nz=1,npx=4096,eps=1,
            z0=1e10, tmax=5,qmax=None,
            method='Greens',fjopt=1):
        if keV : lam = cst.keV2lam(keV)
        self.pattern = pattern
        self.lam  = dtype(lam)

        #initialization
        self.replicate(ax,bz,Nx,Nz)
        self.set_detector(npx,z0,tmax,qmax)
        self.I = np.zeros(self.x0s.shape,dtype=dtype)

        #form factor
        if fjopt==0 :
            self.fj = lambda ti,j:np.ones(ti.shape)
        elif fjopt==1 :
            self.fj = lambda ti,j:np.exp(-(np.pi*Ai[j]*(np.sin(ti*2)/self.lam))**2)#*eps*np.sqrt(np.pi)/Ai[j]*
        #compute
        if method=='Greens'or  method=='G'  : self._Greens()
        elif method=='Fresnel'              : self._Fresnel()
        elif method=='Fraunhofer'           : self._Fraunhofer()
        elif method=='Holton'               : self._Holton(path=path)

    ################################################################
    # init
    ################################################################
    def replicate(self,ax,bz,Nx,Nz):
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
        if qmax:tmax=qmax*self.lam
        self.npx  = npx
        self.z0   = dtype(z0)
        self.tmax  = tmax
        self.x_max = self.tmax*self.z0
        self.px    = self.tmax*self.z0*cst.A/cst.mm/(self.npx-1)

        x0s = np.linspace(-self.x_max,self.x_max,self.npx)
        self.x0s = np.array(x0s,dtype=dtype)+(x0s[1]-x0s[0])/2

    ################################################################
    # display
    ################################################################
    def Pattern_show(self,**kwargs):
        C = [Cs[Z] for Z in self.Za]
        Nx0,Nx1 = int(self.Nx/2),int(np.ceil(self.Nx/2))
        dsp.stddisp(scat=[self.z,self.x,C],
            ms=20,labs=[r'$z(\AA)$',r'$x(\AA)$'],
            xyTicks=[self.bz,self.ax],xylims=[0,self.Nz*self.bz,-Nx0*self.ax,Nx1*self.ax],
            **kwargs)

    def F_show(self,**kwargs):
        Za = np.unique(self.Za)
        t = np.linspace(-self.tmax,self.tmax,100)
        tdeg = t*180/np.pi
        plts = [[tdeg,fj(t,Z),'b','$%d$' %Z] for Z in Za]
        dsp.stddisp(plts,labs=[r'$\theta (deg)$','$f_j$'],**kwargs)

    def getAngle(self):return (self.x0s/self.z0)*180/np.pi
    def getQ(self):return (self.x0s/self.z0)/self.lam
    def getI(self):return self.I

    ################################################################
    # run .c code
    ################################################################
    def save_input(self,file):
        np.savetxt(file,np.vstack([self.z,np.zeros(self.x.shape),self.x,
            np.ones(self.x.shape),np.zeros(self.x.shape),np.zeros(self.x.shape)]).T)
        print(colors.yellow+file+colors.black)

    def cmd(self,opts='s',path='',file='atoms.txt'):
        self.path = path
        self.file = file
        if 's' in opts :self.save_input(path+self.file)
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


    ################################################################
    # Path lengths
    ################################################################
    def _Holton(self,path=''):
        self.cmd(opts='sr',path=path,file='atoms.txt')

    def _Fraunhofer(self):
        for i in range(self.npx) :
            tij  = np.abs(self.x0s[i]-self.x)/(self.z0-self.z)
            R_ij = self.z + np.sqrt((self.x0s[i]-self.x)**2+(self.z0-self.z)**2)
            Rij  = self.x0s[i]*self.x/self.z0
            self.I[i] = np.abs(np.sum(
                self.fj(tij,self.Za)*np.exp(2*np.pi*1J*Rij/self.lam)/(R_ij*cst.A)))**2

    def _Fresnel(self):
        for i in range(self.npx) :
            tij  = np.abs(self.x0s[i]-self.x)/(self.z0-self.z)
            R_ij = self.z + np.sqrt((self.x0s[i]-self.x)**2+(self.z0-self.z)**2)
            Rij  = (self.x0s[i]-self.x)**2/(2*(self.z0-self.z))
            self.I[i] = np.abs(np.sum(
                self.fj(tij,self.Za)*np.exp(2*np.pi*1J*Rij/self.lam)/(R_ij*cst.A)))**2

    def _Greens(self):
        for i in range(self.npx) :
            tij  = np.abs(self.x0s[i]-self.x)/(self.z0-self.z)
            R_ij = self.z + np.sqrt((self.x0s[i]-self.x)**2+(self.z0-self.z)**2)
            # Rij  = R_ij
            self.I[i] = np.abs(np.sum(
            self.fj(tij,self.Za)*np.exp(2*np.pi*1J*R_ij/self.lam)/(R_ij*cst.A)))**2

    def _Dual_beams(self):
        for i in range(self.npx) :
            tij  = np.abs(self.x0s[i]-self.x)/(self.z0-self.z)
            R_ij = self.z + np.sqrt((self.x0s[i]-self.x)**2+(self.z0-self.z)**2)
            # Rij  = R_ij
            self.I[i] = np.abs(np.sum(
            self.fj(tij,self.Za)*np.exp(2*np.pi*1J*R_ij/self.lam)/(R_ij*cst.A)))**2

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
