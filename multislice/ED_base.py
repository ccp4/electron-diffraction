import pandas as pd
import numpy as np


fv = lambda X,X0 : V0*np.exp(-5*2*np.linalg.norm(X-X0,axis=1)**2)
class ED_base():
    '''Base class for simluating ED patterns (z is the beam axis)
    **\ncrystal related**:\n
    - data : x,z,Zatom
    - Nx,Nz : number of unit cells
    - ppopt:V(potential)Q(Psi_q)
    '''
    def __init__(self,data,Nx=1,Nz=1,
        ppopt=''):
        self.data=data
        self.Nx=Nx
        self.Nz=Nz

        if 'V' in ppopt:self.potential_show()
    ##################################################################
    ### display
    def potential_show(self,**kwargs):
        x,z,Za = self.data
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
