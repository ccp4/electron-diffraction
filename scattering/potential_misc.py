'''Some structure factors in 2D'''
from utils import*
from mpl_toolkits.mplot3d import Axes3D
from scattering_factors import*
import scipy
from scipy.interpolate import interp1d
from scipy.signal.windows import tukey
import utils.FourierUtils as fu
import wallpp.plane_group as pg
figpath=get_figpath(__file__,"/figures/")


def get_potential(fxy,a,npts=2**8,opt='p',**kwargs):
    ixy=np.arange(0,npts)*a/npts
    x,y=np.meshgrid(ixy,ixy)
    f=fxy(x,y)
    F=np.fft.fftshift(np.fft.fft2(f))
    if 'p' in opt:
        xy,kxy=fu.get_extents(a,npts)
        fu.plot_fft2d(f,F,xy,kxy,**kwargs)
    return f,F

def test_get_potential():
    k0 = 2*pi/a*np.array([1.0,0.0]); kx0,ky0=k0
    V = lambda x,y:np.cos(kx0*x+ky0*y)
    f,F = get_potential(V,a,npts=512,opt='p',figsize='12',pad=0.5)
def test_pattern_potential():
    pattern = np.array([[0.25,0.25,2], [0.5,0.5,3]])
    crystal=pg.Wallpaper('p2',a,pattern=pattern,gen=True,ndeg=npts,nh=1,nk=1,pOpt='AVp')
    x,f=crystal.get_potential_unit_cell()
    #F=np.fft.fftshift(np.fft.fft2(f))
    #xy,kxy=fu.get_extents(a,npts)
    #fu.plot_fft2d(f,F,xy,kxy,figsize='12',pad=0.5)

def test_Vx_1D(opt='R'):
    ''' opt :
        Q(plot q,Vq) w(show window function),
        R(plot r,Vr) i(show imaginary)'''
    npts,Npad = 1000,10**4
    elts,qmax = ['C'],6
    q,fq = get_elec_atomic_factors(elts,qmax=qmax,npts=npts)
    fq = fq[0]
    w = tukey(2*npts,alpha=0.5,sym=False)[npts:]
    fqw = w*fq

    #qpad,fqpad = q,fq
    qpad,fqpad = np.linspace(0,int(qmax*Npad/npts),Npad),np.zeros((Npad))
    fqpad[:npts] = fqw
    qpad,fqpad = np.concatenate((-np.flip(qpad),qpad)),np.concatenate((np.flip(fqpad),fqpad))
    if 'Q' in opt:
        plts = [[qpad,fqpad,'b','$fq_{pad}$',2]]
        if 'w' in opt : plts += [[q,fq,'b','$fq$',2],[q,fqw,'r','$fq_w$',2],[q,w,'r--','$w$',1]]
        stddisp(plts,labs=['$q(A^{-1})$','$Fq(A)$'],opt='p')

    ##perform 1D fft
    if 'R' in opt:
        r,Vr = fu.get_iFFT(qpad,fqpad)
        plts = [[r,np.abs(Vr),'g','$V$']]
        if 'i' in opt : plts += [[r,Vr.real,'b','$Re$'],[r,Vr.imag,'r','$Im$']]
        stddisp(plts,labs=['r(A)','$V$'],opt='p',lw=2)



def test_Vx_2D(KE=100,tmax=0.2,npts=100,Npad=1000,opt='QRwp'):
    '''
KE          : incident energy (keV)
tmax        : maximum scattering angle (rad)
npts,Npad   : nb points for normal and padded regions
opt         : Q(plot q,Vq) R(plot r,Vr) w(apply window) p(apply padding)'
    '''
    lam = wavelength(KE)
    ikxy = np.linspace(-1,1,npts)
    qmax = tmax/lam
    qx,qy = np.array(np.meshgrid(ikxy,ikxy))*qmax
    q = np.sqrt(qx**2+qy**2)
    fq = get_elec_atomic_factors(['C'],q)[1][0]
    fq[q>qmax]=0

    if 'w' in opt:
        wf = interp1d(ikxy*qmax,tukey(npts,alpha=0.5))
        w = np.zeros(q.shape);w[q<=qmax]=wf(q[q<=qmax])
        fqw = fq*w
        fq = fqw
    if 'p' in opt:
        qM = qmax*Npad/npts
        iqxy = np.linspace(-qM,qM,Npad)
        (qxpad,qypad),fqpad = np.meshgrid(iqxy,iqxy),np.zeros((Npad,Npad))
        qpad = np.sqrt(qxpad**2+qypad**2)
        fqpad[(np.abs(qxpad)<=qmax) & (np.abs(qypad)<=qmax)] = fq[(qx<=qmax) & (qy<=qmax)]
        qx,qy,fq,npts=qxpad,qypad,fqpad,Npad

    if 'Q' in opt:
        stddisp(im=[qx,qy,fq],labs=['$q_x$','$q_y$'],imOpt='c',pOpt='t')

    if 'R' in opt:
        Vr=np.fft.fftshift(np.fft.ifft2(fq))
        dq = qx[0,1]-qx[0,0];#print(dq)
        r = np.fft.fftshift(np.fft.fftfreq(npts,dq))
        x,y = np.meshgrid(r,r)
        im=[x,y,np.abs(Vr)];print('im_max/re_max=%.2f' %(Vr.imag.max()/Vr.real.max()))
        stddisp(im=im,labs=['$x$','$y$'],imOpt='c',pOpt='t')

    print('tmax=%.2f deg\nqmax=%.2f A^-1\nmax_res=%.2f A' %(tmax*180/pi,qmax,lam/tmax))


#if __name__ == "__main__" :
  #test_get_potential()
  #test_pattern_potential()
  #test_Vx_1D(opt='RQ')
  #test_Vx_2D(npts=50,Npad=200,opt='QRwp',KE=100,tmax=0.2)
