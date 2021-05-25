from utils import*
from blochwave import bloch_old as bw;imp.reload(bw)

def test_2beams3D(k0=100,N=2):
    # k0 = 100 #A^-1 => lambda=0.01
    hkl,Gs,Vg = get_pattern3D(N)
    # idxz = np.abs(Gs[:,1])<1e-3; #print(hkl[idxz][:,[0,2]])
    # Gs = Gs[idxz][:,[0,2]]

rot  = lambda t:np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])
def test_2beams(N,s,eta=1,Ts=np.linspace(1,100,100),opts=''):
    ''' Test a 2-beam configuration
    - N : number of beams such that Vg.shape= 2N+1 x 2N+1
    - s : slices of Vg beams for the G beams to perform blochwave calculation
    - np.s_[...] - actual slice
    - int - so that s=np.s_[2*N,2*N-2:2*N+s]
    - Ts : int or np.ndarray thickness
    - eta : misaligment
    - opts : T(thickness), E(Ewald), F(structure factor)
    returns :
    - St : Nbeams x Ts.size array - intentities of beams at the different thicknesses
    '''
    k0 = 5.008343065817388
    if isinstance(s,int) :
        iS = s
        if iS<=2 : s = np.s_[2*N, [2*N-2,2*N]]
        if iS>2  : s = np.s_[2*N, 2*N-2:2*N-2+iS]

        #### basic pattern
        hl,Gs,Fhl = bw.get_pattern2D(N,s)
        gx0,gz0 = Gs.T
        #orient in a 2-beam excitation setup
        gx1,gz1 = rot(eta/10).dot(np.stack([gx0.flatten(),gz0.flatten()]))
        Gs = np.array([gx1,gz1]).T

        #plot reciprocal space lattice
        Sg = (k0**2-np.linalg.norm(k0*np.array([0,1])-Gs,axis=1)**2)/(2*k0)
        if 'E' in opts:
            t = np.linspace(0,2*np.pi,500)
            plts = [[k0*np.cos(t),k0*np.sin(t)+k0,'r','']]
            scat = [np.hstack([gx0,gx1]),np.hstack([gz0,gz1]),['k']*gx0.size+['b']*gx1.size]
            txts = [[g[0],g[1],'%.2f'%s,'b'] for g,s in zip(Gs,Sg)]
            a = 1.5
            dsp.stddisp(plts,texts=txts,scat=scat,ms=20,xylims=[-a,a,-a,a],lw=2,fonts={'text':10})

        #### structure factor
        Fhl[Fhl==Fhl.max()]=0
        if 'F' in opts:
            dsp.stddisp(im=[abs(Fhl)],imOpt='cv',axPos='V',pOpt='t')

        #### SOLVE
        gammaj,CjG = bw.solve_Bloch(k0,hl,Gs,Fhl,v=0)

        #### scattering matrix
        if isinstance(Ts,float) or isinstance(Ts,int) : Ts=np.array([Ts])
        if isinstance(Ts,list):Ts=np.array(Ts)
        Ng = Gs.shape[0]
        St = np.zeros((Ng,Ts.size),dtype=complex)
        for iT,T in enumerate(Ts):
            # S=CjG.T.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG)
            S=CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG.T)
            # St[:,iT]=CjG.dot(np.diag(np.exp(2J*np.pi*gammaj*T))).dot(CjG.T)[0,:]
            St[:,iT] = S[:,0]
        #display thickness dependent intensities
        if 'T' in opts:
            It = np.abs(St)**2
            cs = dsp.getCs('jet',Ng)
            plts = [[Ts,It[iG],[cs[iG],'-o'],'%s' %hl_iG] for iG,hl_iG in enumerate(hl)]
            dsp.stddisp(plts,labs=[r'Thickness $(\AA)$','$Ig$'],gridmOn=1)
        return St,Sg


def test_2beams_rocking(etas=np.linspace(0.5,1.5,10),Ts=100,opts='wz'):
    if isinstance(Ts,float) or isinstance(Ts,int) :Ts=np.array([Ts])
    if isinstance(Ts,list):Ts=np.array(Ts)

    Fg = np.zeros((Ts.size,etas.size),dtype=complex)
    Sg = np.zeros((etas.shape))
    for iE,eta in enumerate(etas):
        St,Sg0 = test_2beams(N=2,s=2,eta=eta,Ts=Ts,opts='')
        Fg[:,iE] = St[1,:]
        Sg[iE] = Sg0[0]
        It = np.abs(Fg)**2
    if 'w' in opts:
        cs = dsp.getCs('Spectral',Ts.size)
        plts = [[Sg,It[iT],[cs[iT],'-o'],'$%.1fA$' %T] for iT,T in enumerate(Ts)]
        dsp.stddisp(plts,labs=['$s_g$','$I_g$'],lw=2)
    if 'z' in opts:
        xi_g=525
        cs = dsp.getCs('Spectral',Sg.size)
        plts = [[Ts,It[:,iS],[cs[iS],'-o'],'$w_g=%.1f$' %(Sg0*xi_g)] for iS,Sg0 in enumerate(Sg)]
        dsp.stddisp(plts,labs=['$s_g$','$I_g$'],lw=2)

# plt.close('all')
# St = test_2beams(N=2,s=2,eta=1,Ts=np.linspace(0,1000,100),opts='T')
# St = test_2beams(N=2,s=3,eta=1,Ts=np.linspace(0,1000,100),opts='ET')
# St = test_2beams(N=2,s=5,eta=1,Ts=np.linspace(0,1000,100),opts='T')
# St = test_2beams(N=2,s=5,eta=1,Ts=np.linspace(0,1000,100),opts='T')
# St = test_2beams(N=2,s=5,eta=1.5,opts='ET')

# test_2beams_rocking(etas=np.linspace(1.0,1.05,4),Ts=np.linspace(0,1000,500),opts='z')
test_2beams_rocking(etas=np.linspace(0.95,1.05,200),Ts=np.array([0.5,0.75,1,1.25,1.5])*525,opts='w')
