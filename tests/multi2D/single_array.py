import importlib as imp
from utils import*                  ;imp.reload(dsp)
import multislice.multi_2D as MS2D  ;imp.reload(MS2D)
import scipy.fftpack as fft
plt.close('all')
path='../../multislice/docs_fig/multi2D/single_array/'



Ai = np.array([0.1,0.25,0.26,0.27,1.5])
def Va( r,Za):
    return np.exp(-(r/Ai[Za])**2)
def Vq(r,Za) :
    fv = np.zeros(r.shape)
    fv[np.abs(r)<Ai[Za]] = 0.5
    return fv

def show_potential_profile(Za=3):
    r = np.linspace(0,1,1000)
    plts = [[r,Va(r,Za),'b','atom'],[r,Vq(r,Za),'r','qdot']]
    dsp.stddisp(plts,labs=['$r$','$V(r)$'],lw=2,fonts='2',
        name=path+'potential.eps',opt='ps')

class SingleArrays:
    def __init__(self,keV=100,nx=2**14,nz=1000,ax=2000,Za=3,
            bzs=[2,5],eps=1,Vatom=True,plot_opt=False):
        self.keV=keV
        self.ax=ax
        self.bzs=np.array(bzs)
        self.eps=eps

        bz0 = bzs[0]
        self.x0,self.z0 = np.linspace(0,ax,nx),np.linspace(0,bz0,nz)
        self.xa,self.za,self.Za = ax/2,bz0/2,Za
        self.x,self.z = np.meshgrid(self.x0,self.z0)
        self.r      = np.sqrt((self.x-self.xa)**2+(self.z-self.za)**2)
        self.Vatom  = Vatom
        self.V      = [Vq,Va][Vatom]

        print('..computing fv...')
        self.fv = self.V(self.r,self.Za)
        if plot_opt:self.plot_v()

    def run(self,Nz=100,iZv=None):
        if not iZv:iZv = Nz
        mp0 = np.zeros(self.bzs.shape,dtype=object)
        pattern = [self.x0,self.z0,self.fv]
        for i,bz in enumerate(self.bzs):
            mp0[i] = MS2D.Multi2D(pattern,self.ax,bz,self.keV,
                    Nx=1,dz=bz,nz=Nz,#nz=int(Nz/bz),
                    ppopt='',iZs=1,iZv=iZv,eps=self.eps,copt=0)
        self.mp0 = mp0

    def plot_v(self):
        return MS2D.plot_v(self.x,self.z,self.fv,self.ax,self.bz0,self.xa,self.za)

    def Qz_show(self,iBz=slice(0,None),iZs=1,opts='N',name='',**kwargs):
        for i,bz in enumerate(self.bzs[iBz]):
            self.mp0[i].Qz_show(iZs=iZs,lw=2,opts=opts,title='bz=%d' %bz,
                name=path+name+'bz%d.svg' %i,**kwargs)

    def plot_bzs(self,idx=-1,**kwargs):
        bzs,mp0 = self.bzs,self.mp0

        cs = dsp.getCs('Spectral',bzs.size)
        for i,bz in enumerate(self.bzs): self.mp0[i].psi_qz[idx,0]=0
        I = [fft.fftshift(mp0[i].psi_qz[idx,:]) for i,bz in enumerate(bzs)]
        plts = [[fft.fftshift(mp0[i].q),I[i]/max(I[i]),cs[i],'bz=%d'%bz] for i,bz in enumerate(bzs)]
        tle = '''diffraction pattern with varying interlattice at
            Nz=%d space eps=%.2f, %s''' %(idx,self.eps,['Vqdot','Vatom'][self.Vatom])
        fig,ax = dsp.stddisp(plts,title=tle,**kwargs)



def show_proj_potential(s0,s1,**kwargs):
    plts=[[s0.mp0[0].x,s0.mp0[0].Vz[0,:],'b--','$\epsilon=%.2f$' %s0.eps],
          [s1.mp0[0].x,s1.mp0[0].Vz[0,:],'r','$\epsilon=%.2f$' %s1.eps],
          ]
    dsp.stddisp(plts,lw=2,labs=[r'$x(\AA)$','$V_z$'],
        name=path+'projV.eps',**kwargs)

def show_kd(s0,s1,idx,xylims=[],**kwargs):
    for i,bz in enumerate(s0.bzs): s0.mp0[i].psi_qz[idx,0]=0
    for i,bz in enumerate(s1.bzs): s1.mp0[i].psi_qz[idx,0]=0
    qa = np.linspace(0,4,300)
    qa_x,qa_y = np.meshgrid(qa,qa)
    qas = np.sqrt(qa_x**2+qa_y**2)
    t = np.linspace(0,np.pi,1000)
    st,ct,ka = np.sin(t),np.cos(t),40
    plts0 = [[ka*st,ka*(1-ct),'r','ewald']]
    for i,bz in enumerate(s0.bzs):
        I0 = s0.mp0[i].psi_qz[idx,:]
        I1 = s1.mp0[i].psi_qz[idx,:]
        plts=[[s0.mp0[i].q,I0/I0.max(),'b--','$\epsilon=%.2f$' %s0.eps],
              [s1.mp0[i].q,I1/I1.max(),'r'  ,'$\epsilon=%.2f$' %s1.eps],
              ]

        f = np.exp(-(np.pi*qas*Ai[s0.Za])**2)*np.cos(qa_x*bz/(4*np.pi))
        dsp.stddisp(plts0,im=[qa_x,qa_y,np.abs(f)],lw=2,labs=['$q_z$','$q_x$'],
            xylims=[0,xylims[1],0,xylims[1]],cmap='viridis',title='bz=%d' %bz,
            name=path+'kd%d_V.png' %(i+1),**kwargs)
        dsp.stddisp(plts,lw=2,labs=[r'$q(\AA^{-1})$','$I$'],
            name=path+'kd%d_D.eps' %(i+1),xylims=xylims,**kwargs)

# show_kd(s0,s1,1,xylims=[0,8,0,0.02],fonts='3',opt='p',caxis=[0,0.001])

# show_potential_profile(Za=3)


ax,nx = 1000,2**15
xy = ax/2 +2*np.array([-1,1])
bzs = np.array([8,16])#[[0]]

s0 = SingleArrays(bzs=bzs,eps=0.001,Vatom=True,plot_opt=False,nx=nx,ax=ax)
s0.run(Nz=101)
# s0.mp0[0].Tz_show(iSz=slice(0,None,1),Vopt='VT',lw=2,xylims=['x',xy[0],xy[1]])
s0.Qz_show(iBz=[0,-1],iZs=np.arange(0,4,1),xylims=[0,8,0,0.00001],opts='',opt='p')
# s0.plot_bzs(idx=0,lw=2,xylims=[-4,4,0,1],opt='p')
# s0.plot_bzs(idx=-1,lw=2,xylims=[-4,4,0,1],opt='p')


s1 = SingleArrays(bzs=bzs,eps=0.01,Vatom=True,plot_opt=False,nx=nx,ax=ax)
s1.run(Nz=101)
# s1.mp0[0].Tz_show(iSz=slice(0,None,1),Vopt='VT',lw=2,xylims=['x',xy[0],xy[1]])
# s1.Qz_show(iBz=[0,-1],iZs=np.arange(0,50,15),xylims=['x',-4,4],opts='N',opt='p')
# s1.plot_bzs(idx=0,lw=2,xylims=[-4,4,0,1],opt='p')
# s1.plot_bzs(idx=-1,lw=2,xylims=[-4,4,0,1],opt='p')
#
# s0.plot_bzs(idx=1,lw=2,xylims=[-4,4,0,1],opt='p')
# s1.plot_bzs(idx=1,lw=2,xylims=[-4,4,0,1],opt='p')



# show_proj_potential(s0,s1,xylims=[xy[0],xy[1],0,0.3],fonts='2',opt='ps')
