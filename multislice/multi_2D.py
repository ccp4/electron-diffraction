import pickle,matplotlib
import numpy as np, pandas as pd
import scipy.fftpack as fft
from scipy.integrate import nquad,trapz,quad
import utils.displayStandards as dsp
import utils.physicsConstants as cst
import utils.glob_colors as colors
import wallpp.lattice as lat
import wallpp.plane_group as pg     #;imp.reload(pg)
# from scipy.signal import fftconvolve
# from scipy.interpolate import interp1d
# import utils.FourierUtils as fu


class Multi2D():
    '''multislice 2D for quick testing
    - pattern : list of 2d-arrays - [x0,z0,f] where x,z=np.meshgrid(x0,z0)
    - ax,bz : lattice constants
    - tilt : beam tilt (deg)
    - Nx : increase supercell size
    - dz : slice thickness
    - nz : number of slices
    - iZs,iZv,v : frequency of save and verbose frequency,verbose option
    - ppopt:T(Transmission)P(propagator)Q(Psi_q)X(Psi_x)B(beams)Z(Psi_xz)Y(Psi_qz)
    - copt: crop option with propagator
    '''
    def __init__(self,
            pattern,ax,bz,
            keV=200,Nx=1,tilt=0,
            dz=1,nz=1,copt=1,eps=1,sg=-1,
            TDS=False,nTDS=8,ndeg=2**10,wobble=0.05,
            iZs=1,opts='q',iZv=1,v=1,ppopt=''):
        self.version = 1.0
        self.pattern = pattern
        self.keV     = keV
        self.Nx      = Nx
        self.dz      = dz
        self.ax      = ax
        self.bz      = bz
        #energy dependent
        self.lam  = cst.keV2lam(keV)
        self.sig  = cst.keV2sigma(keV)
        self.k0   = 1/self.lam
        self.eps  = eps
        self.tilt = tilt
        self.copt = copt
        #TDS
        self.TDS  = TDS
        self.nTDS = nTDS
        self.ndeg = ndeg
        self.wobble = self._wobble(wobble)

        #Computations
        if v:print(colors.red+'\t\t 2D multislice simulation '+colors.black)

        # Thermal diffuse scattering
        if self.TDS:
            self._TDS(eps,copt,sg,v)
        else:
            # standard
            self._set_transmission_function(eps,v)
            self._set_propagator(sg,copt)
            self.set_Psi0()
            if nz:self.propagate(nz,iZs,iZv,opts,v)

        display(ppopt)

    def _TDS(self,eps,copt,sg,v):
        pattern  = np.vstack([self.pattern + np.array([self.ax*i,0,0]) for i in range(self.Nx)])
        wobbles  = self.wobble[np.array(pattern[:,-1],dtype=int)][:,None]
        pattern0 = pattern
        Na = pattern0.shape[0]
        self.patterns=[]
        for iTDS in range(nTDS):
            print(colors.yellow+'configuration %d' %iTDS+ colors.black)
            pattern[:,:2] = pattern0[:,:2]+(2*np.random.rand(Na,2)-1)*wobbles#*[self.ax,self.bz]
            self.patterns+=[pattern.copy()]

            p1 = pg.Wallpaper('p1',Nx*self.ax,self.bz,90,pattern,ndeg=ndeg,gen=True)
            self.pattern = p1.get_potential_grid_p1()

            self._set_transmission_function(eps,v)
            self._set_propagator(sg,copt)
            self.set_Psi0(iTDS)
            self.propagate(nz,iZs,iZv,opts,iTDS,v)
        self.psi_qz=np.abs(self.psi_qz)**2
        self.pattern=pattern0

    def display(self,ppopt):
        if 'T' in ppopt:self.Tz_show()
        if 'Z' in ppopt:self.Za_show()
        if 'P' in ppopt:self.Pq_show()
        if 'Q' in ppopt:self.Qz_show()
        if 'X' in ppopt:self.Xz_show()
        if 'B' in ppopt:self.Bz_show()
        if 'Z' in ppopt:self.Xxz_show()
        if 'Y' in ppopt:self.Qxz_show()

    def _wobble(self,wobble):
        if isinstance(wobble,int) or isinstance(wobble,float): wobble = [wobble]*5
        if isinstance(wobble,dict) :
            wobbles = [0]*5
            for k,v in wobble.items():wobbles[k]=v
            wobble = wobbles
        return np.array(wobble)

    def save(self,file):
        with open(file,'wb') as out :
            pickle.dump(self, out, pickle.HIGHEST_PROTOCOL)
        print(colors.green+"object saved : \n"+colors.yellow+file+colors.black)

    ###################################################################
    #### display
    ###################################################################
    def V_show(self,**kwargs):
        x0,z0,fv = self.pattern
        x,z = np.meshgrid(x0,z0)

        args = {'imOpt':'cv','axPos':'V','xylims':[z0.min(),z0.max(),x0.min(),x0.max()]}
        args.update(kwargs)
        dsp.stddisp(im=[z,x,fv],labs=['$z$','$x$'],**args)

    def Za_show(self):
        atoms = pd.DataFrame.from_dict(
            { 0:[(0,0,0),15], 1:[(0.5,0.5,0.5),40],
              2:[(0,0,1),45], 3:[(1,0,0),50],
              4:[(0,1,1),80]},
            orient='index',columns=['c','s'])

        markers = list(matplotlib.markers.MarkerStyle.markers.keys())[2:-4]
        nms=len(markers)
        scat = ()
        for i,p in enumerate(self.patterns):
            Za = atoms.iloc[p[:,2]]
            # print(Za.size)
            scat += ([p[:,0],p[:,1],Za.s,Za.c,markers[i%nms]] ,)
        dsp.stddisp(scat=scat,labs=['$x$','$z$'])

    def Tz_show(self,iSz=slice(0,None,1),opts='',Vopt='VT',cmaps=['Greens','Blues','Reds'],**kwargs):
        '''Show Transmission function
        - opts or Vopt : 'V(potential) 'T'(Transmission) 'l'(slice legend)
        if 'l' is selected
        '''
        if opts:Vopt=opts
        if isinstance(iSz,int):iSz=[iSz]
        if isinstance(iSz,slice):iSz=list(np.arange(self.ns)[iSz])
        if isinstance(iSz,list):N=len(iSz)
        if isinstance(cmaps,str):cmaps=[cmaps]*3
        cs1,cs2,cs3 = dsp.getCs(cmaps[0],N),dsp.getCs(cmaps[1],N),dsp.getCs(cmaps[2],N)
        plts,legElt=[],{}
        if 'V' in Vopt:
            plts += [[self.x,self.Vz[iSz[i],:].T ,cs1[i]] for i in range(N)]
            legElt[r'$V_z(kV\AA)$']=[cs1[int(N/2)],'-']
        if 'T' in Vopt:
            plts+= [[self.x,self.T.real[iSz[i],:].T,cs2[i]] for i in range(N)]
            plts+= [[self.x,self.T.imag[iSz[i],:].T,cs3[i]] for i in range(N)]
            legElt['$re(T)$']=[cs2[int(N/2)],'-']
            legElt['$im(T)$']=[cs3[int(N/2)],'-']
        return dsp.stddisp(plts,labs=[r'$x(\AA)$',''],#title='Projected potential $V_z$, Transmission function $T$',
            legElt=legElt,**kwargs)

    def Pq_show(self,**kwargs):
        '''Show Propagator function in q space '''
        q  = fft.fftshift(self.q)
        Pq = fft.fftshift(self.Pq)
        plts = [[q,Pq.real,'b','$re$'],
                [q,Pq.imag,'r','$im$']]
        return dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$','$P_q$'],
            title='Propagator',
            **kwargs)

    def Xz_show(self,iZs=1,cmap='jet',**kwargs):
        '''Show wave propagation for slices iZs '''
        if isinstance(iZs,int):iZs=slice(0,self.z.size,iZs)
        z  = self.z[iZs]
        Px = self.psi_xz[iZs,:]
        cs = dsp.getCs(cmap,z.size)
        plts = [ [self.x,Px[i,:],cs[i]] for i in range(z.size)]
        return dsp.stddisp(plts,labs=[r'$x(\AA)$',r'$|\Psi(x)|^2$'],
            imOpt='hc',caxis=[z.min(),z.max()],cmap=cmap,
            **kwargs)

    def Qz_show(self,iZs=1,opts='',cmap='jet',**kwargs):
        '''Show wave propagation in q space for slices iZs
        - opts : O(include Origin), N(Normalize with max), S(Shape normalize) l(show label)
        '''
        # if isinstance(iZs,int):iZs=slice(0,iZs,self.z.size)
        if isinstance(iZs,int):iZs=slice(0,None,iZs)
        # if isinstance(iZs,list):iZs=slice(0,iZs,self.z.size)
        Pqs = np.zeros((self.psi_qz[iZs,:].shape))
        z   = self.z[iZs]
        Pqs = self.psi_qz[iZs,:].copy()
        if 'S' in opts:
            if 'O' not in opts:Pqs[:,0]=0# do not show central beam
            for i in range(z.size) : Pqs[i,:]/=np.sum(Pqs[i,:])
        else:
            if 'O' not in opts:Pqs[:,0]=0# do not show central beam
            if 'N' in opts:
                for i in range(z.size) : Pqs[i,:]/=Pqs[i,:].max()
        q   = fft.fftshift(self.q.copy())
        Pqs = [fft.fftshift(Pqs[i,:]) for i in range(z.size)]
        cs  = dsp.getCs(cmap,z.size)
        if 'l' in opts:
            plts = [[q,Pqs[i],cs[i] ,'$%.1fA$' %z[i]] for i in range(z.size)]
            return dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$',r'$|\Psi(q)|^2$'],
                imOpt='c',caxis=[z.min(),z.max()],cmap=cmap,axPos='V',
                **kwargs)
        else:
            plts = [[q,Pqs[i],cs[i]] for i in range(z.size)]
            return dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$',r'$|\Psi(q)|^2$'],
                imOpt='c',caxis=[z.min(),z.max()],cmap=cmap,axPos='V',
                **kwargs)

    def Bz_show(self,iBs='O',tol=1e-3,cmap='jet',sym_opt=False,plts=[],**kwargs):
        '''Show selected beam iBs as function of thickness(see getB)
        - sym_opt : special display for symmetry pairs
        '''
        # print(iBs)
        iBs,Ib = self.getB(iBs,tol,v=1)#;print(iBs)
        h    = ['%d_{%d}' %(i/self.Nx,i%self.Nx) for i in iBs]
        cs   = dsp.getCs(cmap,iBs.size)
        plts = [[self.z,Ib[:,i],[cs[i],'-'],'$%s$' %h[i]] for i,iB in enumerate(iBs)]
        if sym_opt:
            iBs,Ib = self.getB(-iBs,tol,v=1)
            h    = ['%d_{%d}' %(i/self.Nx,i%self.Nx) for i in iBs]
            plts += [[self.z,Ib[:,i],[cs[i],'--'],'$%s$' %h[i]] for i,iB in enumerate(iBs)]
        return dsp.stddisp(plts,labs=[r'$z(\AA)$',r'$I_b$'],
        **kwargs)

    def Xxz_show(self,iZs=1,iXs=1,**kwargs):
        '''Show 2D wave propagation solution'''
        if isinstance(iZs,int):iZs=slice(0,-1,iZs)
        if isinstance(iXs,int):iZs=slice(0,-1,iXs)
        x,z = np.meshgrid(self.x,self.z)
        im = [x[iZs,:],z[iZs,:],self.psi_xz[iZs,:]]
        return dsp.stddisp(im=im,labs=[r'$x(\AA)$',r'$z(\AA)$'],
        **kwargs)

    def Qxz_show(self,iZs=1,**kwargs):
        '''Show 2D wave propagation solution'''
        iZs = slice(0,-1,iZs)
        q,z = np.meshgrid(self.q,self.z[iZs])
        Pqs = self.psi_qz[iZs,:].copy()
        Pqs[:,0] = 0
        im = [q,z,Pqs]
        return dsp.stddisp(im=im,labs=[r'$q(\AA^{-1})$',r'$z(\AA)$'],
        **kwargs)
    def Ewald_show(self,deg=10,nh=20,nk=10,relrod=0,**kwargs):
        ''' Displays the Ewald sphere and the reciprocal lattice
        - deg : extent of the ewald sphere in degrees
        - nh,nk : number of reciprocal nodes
        '''
        # Ewald circle
        K  = self.k0
        dt = deg*np.pi/180
        t  = 3*np.pi/2+dt*np.linspace(-1,1,1000)
        # reciprocal lattice
        a1,a2 = lat.get_lattice_vec(lat_type='rect',a=self.ax,b=self.bz)
        b1,b2 = lat.reciprocal_lattice_2D(a1,a2)
        h,k = np.meshgrid(np.arange(-nh,nh+1),np.arange(nk))
        X = h*b1[0]+k*b2[0]
        Z = h*b1[1]+k*b2[1]
        #tilt
        rot = lambda t:np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])
        theta = self.tilt*np.pi/180
        X,Z  = rot(-theta).dot(np.stack([X.flatten(),Z.flatten()]))
        plts = []
        if relrod:
            H    = Nz*bz
            zeta = np.linspace(0,0.25/bz,100)
            Fz   = lambda i : 1/(1.1*self.ax)*np.sinc(zeta*H)**2+i/self.ax
            plts +=[[Fz(i),zeta,'b--',''] for i in range(nh)]
        plts += [[K*np.cos(t),K*np.sin(t)+K,'r','']]
        scat  = [X,Z,15,'k']
        # ax1 = np.sqrt(1**2+10**2)*ax
        # bz1 = ax1
        return dsp.stddisp(plts,scat=scat,labs=[r'$k_x(\AA^{-1})$',r'$k_z(\AA^{-1})$'],
            **kwargs)#lw=2,xylims=[0,3,0,3])#,xyTicks=[1/ax1,1/bz1])

    ##########################################
    #### get methods
    def getQ(self):return self.q
    def getI(self):return self.psi_qz[-1,:]
    def getB(self,iBs='Oa',tol=0,v=0):
        ''' get beam as function of thickness
        - iBs :
            - str - 'O'(include origin) a(all)
            - list - the list of indices of the beams
        - tol : select only beams max(I)>tol*I_max
        '''
        Ib = self.psi_qz/self.nx**2/self.dq
        if isinstance(iBs,list):np.array(iBs)
        if isinstance(iBs,str):
            N   = int(self.nx/2)
            iHs = fft.fftshift(np.arange(-N,N))
            if not 'a' in iBs:
                Im  = Ib[:,1:].max()    #;print(Im)
                Imax = Ib.max(axis=0)   #;print(Imax)
                iHs = iHs[Imax>Im*tol]
            iBs=iHs['O' not in iBs:]
        if isinstance(iBs,list):iBs=np.array(iBs)
        if v:
            return iBs,Ib[:,iBs]
        else:
            return Ib[:,iBs]

    ##################################################################
    ###### main computations
    ##################################################################
    def _set_transmission_function(self,eps=1,v=1):
        x,z,f = self.pattern
        nx = x.size
        ns = int(np.round(self.bz/self.dz))
        self.dz = self.bz/ns;
        Vz = np.zeros((ns,nx))
        iZs = np.arange(0,ns+1)*int(z.size/ns)#;print(iZs)
        if v:print(colors.blue+'...integrating projected potential...'+colors.black)
        if v:print('Slice thickness and number of slices per cell:dz=%.2fA, nzs=%d\n' %(self.dz,ns))
        for i_s in range(ns):
            s=slice(iZs[i_s],iZs[i_s+1])
            Vz[i_s,:] = eps*np.array([trapz(f[s,i],z[s]) for i in range(nx)])

        T  = np.exp(1J*self.sig*Vz)

        Nx = [self.Nx,1][self.TDS]
        self.x  = np.hstack([x + self.ax*i for i in range(Nx)])
        self.Vz = np.hstack([Vz]*Nx)
        self.T  = np.hstack([T]*Nx)             #;print(self.T.shape)
        self.nx = self.x.size                   #;print(self.nx)
        self.ns = ns

    def _set_propagator(self,sg=-1,copt=1):
        self.dx = self.x[1]-self.x[0]
        self.q  = fft.fftfreq(self.nx,self.dx)
        self.dq = self.q[1]-self.q[0]
        if self.tilt: #isinstance(tilt,np.ndarray) or isinstance(tilt,list):
            kx = self.k0*np.sin(self.tilt*np.pi/180)
            self.Pq = np.exp(sg*1J*np.pi*self.dz*(self.q+kx)**2/self.k0)
        else:
            self.Pq = np.exp(sg*1J*np.pi*self.dz*self.q**2/self.k0)
        self.nq = int(1/3*self.nx) #prevent aliasing

        if copt:self.Pq[self.nq:-self.nq] = 0

    def set_Psi0(self,iTDS=0):
        Psi  = np.ones(self.x.shape,dtype=complex)
        self.Psi_x = Psi/np.sqrt(np.sum(np.abs(Psi)**2)*self.dx)
        if self.TDS :
            if not iTDS :
                self.psi_xz = np.zeros((0,self.nx))
                self.psi_qz = np.zeros((0,self.nx),dtype=complex)
        else:
            self.psi_xz = np.zeros((0,self.nx))
            self.psi_qz = np.zeros((0,self.nx))
        self.z  = np.array([])
        self.iz = 0

    def propagate(self,nz,iZs=1,iZv=1,opts='q',iTDS=0,v=1):
        '''Propgate over nz slices and save every iZs slices'''
        nzq,z0 = int(nz/iZs),0
        if self.z.size : z0=self.z.max()
        self.z  = np.hstack([self.z,z0+self.dz+np.arange(nzq)*self.dz*iZs ])
        if 'x' in opts and not iTDS:self.psi_xz = np.vstack([self.psi_xz,np.zeros((nzq,self.nx))])
        if 'q' in opts and not iTDS:self.psi_qz = np.vstack([self.psi_qz,np.zeros((nzq,self.nx))])
        # self.T=fft.fftshfft(self.T)
        for i in range(nz):
            i_s=i%self.ns
            #print(self.T[i_s,:].shape,self.Psi_x.shape)
            self.Psi_q = fft.fft(self.T[i_s,:]*self.Psi_x) #periodic assumption
            if self.copt:self.Psi_q[self.nq:-self.nq] = 0  #prevent aliasing
            self.Psi_x = fft.ifft(self.Pq*self.Psi_q)
            # self.Psi_x = fft.fftshift(fft.ifft(self.Pq*self.Psi_q))
            #save and print out
            msg=''
            if v and (not i%iZv or i==nz-1):
                Ix2 = np.sum(np.abs(self.Psi_x)**2)*self.dx
                Iq2 = np.sum(np.abs(self.Psi_q/self.nx)**2)/self.dq #parseval's theorem of the DFT
                msg+='i=%-4d,islice=%-2d I=%.4f, Iq=%.4f ' %(i,i_s,Ix2,Iq2)
            if not i%iZs :
                if msg and v: msg+='iz=%d, z=%.1f A' %(self.iz, self.z[self.iz])
                if 'x' in opts:
                    self.psi_xz[self.iz,:] = np.abs(self.Psi_x)**2
                if 'q' in opts:
                    if self.TDS:
                        self.psi_qz[self.iz,:] += self.Psi_q
                    else:
                        self.psi_qz[self.iz,:] = np.abs(self.Psi_q)**2
                    self.iz+=1
            if msg:print(colors.green+msg+colors.black)

##################################################################
#### misc functions
##################################################################
def load(filename):
    '''load a saved Multislice object'''
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



def plot_v(x,z,fv,ax,bz,xa,za):
    print('..plotting fv...')
    dr=bz/2
    idx = np.abs(x[0]-ax/2)<dr
    xylims = np.hstack([za+dr*np.array([-1,1]),xa+dr*np.array([-1,1])])
    dsp.stddisp(im=[z[:,idx],x[:,idx],fv[:,idx]],
        labs=['$z$','$x$'],imOpt='c',axPos='V',xylims=xylims,opt='p')


##################################################################
###### Base test
##################################################################
if __name__=='__main__':
    import wallpp.plane_group as pg
    import importlib as imp
    imp.reload(pg)
    ndeg = 2**8
    pptype,a,b,angle = 'p1',20,6,90
    pattern = np.array([[10,2,3],[5,4,3]])

    #get the potential from wallpaper library
    p1      = pg.Wallpaper(pptype,a,b,angle,pattern,ndeg=ndeg)
    pattern = p1.get_potential_grid()
    mp1 = multi2D(pattern,a,b,keV=200,
            Nx=1,dz=b/2,nz=2,ppopt='',#XQZTP
            iZs=1,iZv=1)
    mp1.Tz_show(slice(0,None,1))
