from utils import*
from scipy.integrate import trapz,quad
from scattering_factors import*
fig_path = get_figpath(__file__,'/figures/')


def display_MottBethe_fit(qmax=3,npts=100,**kwargs):
    elts = ['C']
    q,fq_x = get_xray_atomic_factors(elts,qmax,npts)
    q,fq_e = get_elec_atomic_factors(elts,qmax,npts)

    Z,q2=6, q**2
    fq_x,fq_e=fq_x[0],fq_e[0]
    f_MB = 1/(2*pi**2*a0)*(Z-fq_x)/q2
    #f_z = 1/(2*pi**2*a0)*Z/q2
    #fe0 = Z/(3*a0)*r2

    plts = [[q,fq_e,'b','$Gauss-fit$'],
            [q,f_MB,'r*','$Mott-Bethe$'],
            #[q,f_z**2,[(1,0.7,0.7),'--'],'$Z/q^2$']
            ]
    stddisp(plts,labs=['$q(A^{-1})$','$f^e(A)$'],xylims=[0,qmax,0,3],**kwargs)

# def compute_lenz_int(Z=6,KE=100):
#     v,lam = get_v_from_KE(KE),wavelength(KE)
#     g2,k0 = 1/(1-v**2),2*pi/lam
#     theta0 = pow(Z,1/3)/(k0*a0)
#
#     t = np.linspace(0,pi,100000)
#     fq_l = lambda t:2*np.sin(t)/(t**2+theta0**2)**2
#     si=trapz(fq_l(t),dx=t[1]-t[0])
#     siq=quad(fq_l,0,pi)
#     print(si,siq[0],1/theta0**2, theta0)
#     return siq
#
# def display_Lenz():
#     KE,Zmax=100,20
#     Z,v = list(range(1,Zmax+1)),get_v_from_KE(KE)
#     sigma_Lenz = 1.87e-4*np.array(Z)**(4/3)/v**2
#     plts = [#[Z,sigma,'rs-','$th$'],
#             [Z,sigma_Lenz,'bs--','$Lenz$']]
#     stddisp(plts,labs=['$Z$','$\sigma_e(A^2)$'],opt='p',gridmOn=True,logOpt='xy')



display_MottBethe_fit(qmax=3,npts=100,lw=3,legLoc='upper right',title='Carbone'
    name=fig_path++'MottBethe.svg',opt='p')
