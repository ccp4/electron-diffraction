import numpy as np
from utils import*
# plt.close('all')

Natoms,dim = 3,2
a,b = 5,3       # lattice constants A
lam = 0.05      # wavelength A
Nx,Nz = 10,5   # unit cells

opts = 'IP' #'PIFn(new)'
Cs = [unicolor(0.75),(1,0,0),unicolor(0.25),(0,0,1),(1,0,0),(1,1,0)]
Zs = [1,3,6,7,8,16]


Ai = (10-0.5*np.arange(len(Zs)))/8
# Ai = dict(zip(Zs,Ai))
# fj  = lambda ti,j:np.ones(ti.shape)#exp(-Ai[j]*(np.sin(ti)/lam)**2)
fj  = lambda ti,j:np.exp(-Ai[j]*(np.sin(ti)/lam)**2)

Nx,Nz = Nx+(Nx+1)%2,Nz+(Nz+1)%2
Nx0,Nz0 = int(Nx/2),int(Nz/2)
nx,nz = np.meshgrid(np.arange(-Nx0,Nx0+1),np.arange(-Nz0,Nz0+1))
nx,nz = nx.flatten(),nz.flatten()

if 'n' in opts:
    X   = np.random.rand(Natoms,dim)
    Za  = np.random.randint(len(Zs),size=Natoms)
    np.save('atoms.npy',[X,Za,Natoms])
else:
    X,Za,Natoms = np.load('atoms.npy')

x,z = X.T
x  = np.hstack([ (x+i)*a for i in nx])
z  = np.hstack([ (z+i)*b for i in nz])
Za = np.hstack([Za]*Nx*Nz)
idx = np.argsort(np.array(z))
x,z,Za = x[idx],z[idx],Za[idx]

if 'P' in opts:
    C = [Cs[Z] for Z in Za]
    dsp.stddisp(scat=[z,x,C],ms=20,labs=[r'$z(\AA)$',r'$x(\AA)$'],
    xyTicks=[b,a],xylims=[-Nz0*b,Nz0*b,-Nx0*a,Nx0*a])

qmax = 4
tmax = qmax*lam
if 'F' in opts:
    t = np.linspace(-tmax,tmax,100)
    tdeg =t*180/np.pi
    plts =[[tdeg,fj(t,0),'b','0'],[tdeg,fj(t,5),'r','5']]
    dsp.stddisp(plts,labs=[r'$\theta (deg)$','$f$'])


npx = 4000
z0  = 10e8 #
x0s = 1*np.linspace(-tmax*z0,tmax*z0,npx)
I   = np.zeros(x0s.shape)#,dtype=complex)
Ifr = np.zeros(x0s.shape)#,dtype=complex)
# I2  = np.zeros(x0s.shape)#,dtype=complex)
#loop through pixels
### Find at the error
for i in range(npx) :
    tij = np.abs(x0s[i]-x)/(z0-z)
    Rij = np.sqrt((x0s[i]-x)**2+(z0-z)**2)
    rij = x0s[i]*x/z0           #Fraunhofer=>FFT
    rij = (x0s[i]-x)**2/(2*z0)  #Fresnel
    # rij = Rij                   #Greens
    I[i] = np.abs(np.sum(
        fj(tij,Za)*np.exp(2*np.pi*1J*Rij/lam)/Rij))**2 #print(idx)
    Ifr[i] = np.abs(np.sum(
        fj(tij,Za)*np.exp(2*np.pi*1J*rij/lam)/Rij))**2 #print(idx)
    # F=0
    # for j in range(x.size):
    #     Rij = np.sqrt((x0s[i]-x[j])**2+(z0-z[j])**2)
    #     tij = np.abs(x0s[i]-x[j])/(z0-z[j])/2
    #     rij = Rij
    #     F += fj(tij,Za[j])*np.exp(2*np.pi*1J*rij/lam)/Rij
    # I2[i] = np.abs(F)**2
I/=I.max()
Ifr/=Ifr.max()
# I2/=I2.max()

if 'I' in opts :
    t0s = (x0s/z0)*180/np.pi
    tle = r'$D=%.1f cm$, $px=%.1f\mu m$, $n_{px}=%d$' %(z0/1e8,2*x0s.max()/npx/1e4,npx)
    plts=[[t0s,I,'r','$I$'],[t0s,np.abs(fj(x0s/z0,1))**2,'g--','$f^2$'],]
    plts+=[[t0s,Ifr,'b','$I_{fr}$']]
    # plts+=[[t0s,I2,'b','$I_2$']]
    dsp.stddisp(plts,labs=[r'$\theta (deg)$','$I_{flux}$'],
        title=tle,fonts={'title':20})
