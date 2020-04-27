from utils import*
from scattering_factors import wavelength,
figpath=get_figpath(__file__,"/figures/")

def plot_interaction_param(npts=1000,**kwargs):
    KE=np.linspace(10,1000,npts) #keV
    lam = wavelength(KE)
    sigma_e = pi/(lam*KE)                    #rad/A*keV
    sigma_0 = 2*pi*m0*lam*A*eV/h**2*(1+KE/(2*mc2))/(1+KE/mc2)
    sigma = 2*pi*m0*(1+KE/mc2)*lam*A*eV/h**2 #rad/(m*V)
    plts = [[KE,sigma*A*kV,'b',"relativistic"],
            [KE,sigma_0*A*kV,'r--','$\pi/\lambda_0 E_0$']]
    stddisp(plts,labs=['$E_0(keV)$','$\sigma(rad/kV.A)$'],**kwargs)


def test_mulslice_params():
    a0,b0=3.995,5.65
    NxNy=[[5,3],[6,4],[7,5],[8,6]]
    N = np.array([128,256,512])
    KE=200;

    lam = wavelength(KE)
    AB_0 = np.array([max(n[0]*a0,n[1]*b0) for n in NxNy])
    dtheta=lam/AB_0;print("dtheta(mrad)=",dtheta*1e3)
    tmax = (N[:,None]/2).dot(dtheta[None,:])*2/3*1e3

    tmax_book=np.array([53.6,107.1,214.3])
    print(tmax);print(tmax_book)

if __name__ == "__main__" :
    plt.close('all')
    #test_mulslice_params()
    plot_interaction_param(npts=1000,legLoc="upper right",lw=2,
        name=figpath+"scattering_param.svg",opt='s',figsize='f')
