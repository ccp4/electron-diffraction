from utils import*
from scipy import integrate
import subprocess

figpath=get_figpath(__file__,'/figures/')


def compute_Iz(a,lam,nz,npts=1000,nA=3,filepath='dat/Iz.npy'):
    x = np.linspace(-nA*a,nA*a,npts)
    fr = lambda xp,x,nz:np.cos(pi/nz*((x-xp)/lam)**2)
    fi = lambda xp,x,nz:np.sin(pi/nz*((x-xp)/lam)**2)
    nzs = len(nz)

    # compute integrals
    Iz = np.zeros((npts,nzs))
    for i in range(nzs):#loop over distance
        for j in range(npts): #loop over observation point
            quad_r = integrate.quad(fr,-a/2,a/2,args=(x[j],nz[i]))
            quad_i = integrate.quad(fi,-a/2,a/2,args=(x[j],nz[i]))
            Iz[j,i]=np.abs(quad_r[0]+1J*quad_i[0])**2/(lam**2*nz[i])
    params={'a':a,'nz':nz,'lambda':lam}
    if filepath:
        np.save(filepath,[x,Iz,params])
        print(green+'succesfully saved : '+yellow,filepath,black)
    return x,Iz,params




#########################################################################
## Display
def plot_overlay_Iz(x,Iz,nz,a):
    nzs, mM = len(nz), Iz.max()
    cs = getCs('viridis',nzs)
    plts=[[x,Iz[:,i],cs[i],'$%.1f\lambda$' %(nz[i])] for i in range(nzs)]
    plts+=[ [[-a/2,-a/2],[0,mM],'g',''], [[a/2,a/2],[0,mM],'g',''] ]
    stddisp(plts,labs=['$x$','$I(x)$'],figsize='f',opt='p',lw=2)

def plot_Iz(x,Iz,nz,a,lam,gif=1):
    nzs, mM = len(nz), Iz.max()
    for i in range(nzs):
        tle='$a=%.1f,\lambda=%.1f, z=%.1f\lambda$' %(a,lam,nz[i] )
        plts=[[x,Iz[:,i],'b','']]
        plts+=[ [[-a/2,-a/2],[0,mM],'g',''], [[a/2,a/2],[0,mM],'g',''] ]
        stddisp(plts,labs=['$x$','$I(x)$'],figsize='f',opt='s',lw=2,
            name=figpath+'fresnel%s.svg' %(str(i).zfill(3)),title=tle,legOpt=0)
        plt.close()
    if gif : make_gif()
def make_gif():
    cmd='im2gif '+figpath+'fresnel svg'
    print(cmd+'  ...')
    out=subprocess.check_output(['/bin/bash','-i','-c',cmd]).decode()
    print(green+out+black)

###################################################################
#test#
###################################################################
def test_run_sweep(new=0):
    a,lam = 2,1
    npts,nA = 1000,4
    nz = np.linspace(0.1,2,20).tolist()+list(range(3,10))

    filepath='dat/Iz_a2_30.npy'
    if new : compute_Iz(a,lam,nz,npts=1000,filepath=filepath)
    x,Iz,params = np.load(filepath)
    #plot_overlay_Iz(x,Iz,nz,a)
    plot_Iz(x,Iz,params['nz'],params['a'],params['lambda'])

def test_FrvsFF(opt='p'):
    a,lam = 2,0.05
    nzs=[100]#[500,2000] # in lam
    for nz in nzs:
        z=nz*lam
        x,Iz,params = compute_Iz(a,lam,[nz],nA=5,npts=1000,filepath='')
        I_fr = Iz[:,0]/Iz.max()
        I_ff = np.sinc(a/lam*x/z)**2
        plts = [[x,I_fr,'b','$Fresnel$'],
                [x,I_ff,'r','$Fraunhofer$']]
        plts+=[ [[-a/2,-a/2],[0,1],'g',''], [[a/2,a/2],[0,1],'g',''] ]
        #tle = 'Fresnel vs Fraunhofer : ''
        tle = '$a=%.2f A, \lambda=%.2f A, z=%.2f A$' %(a,lam,z)
        stddisp(plts,labs=['$x$','$I(x)$'],figsize='12',opt=opt,lw=2,xylims=[-5,5,0,1],
            title=tle,fonts={'leg':15,'title':20},name=figpath+'ffvsfr_%d.png' %(nz))

#test_run_sweep(new=0)
test_FrvsFF(opt='p');#plt.close('all')
