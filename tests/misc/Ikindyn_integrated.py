from utils import*
from scipy.integrate import trapz
plt.close('all')
path='../figures/'

Ikf = lambda t,sg,xi_g:(np.pi*t/xi_g)**2*np.sinc(t*sg)**2
Idf = lambda t,sg,xi_g:(np.pi*t/xi_g)**2*np.sinc(t/xi_g*np.sqrt(1+(sg*xi_g)**2))**2
def integrate_th_full(Ts=np.linspace(1,5000,5000),xi=500,iTs=0,**kwargs):
    '''Integrate the full rocking curve for each thickness
    - Ts : thicknesse
    - xi : Pendullosung
    - iTs : selected indices in Ts for displaying rocking curves
    '''
    nt = 50
    Ikin,Idyn = np.zeros(Ts.shape),np.zeros(Ts.shape)
    for iT,t in enumerate(Ts):
        sg = np.linspace(-1,1,20000)*nt/t
        Ikin[iT] = trapz(Ikf(t,sg,xi),sg)
        Idyn[iT] = trapz(Idf(t,sg,xi),sg)
    plts = [[Ts,Ikin,'b-','$kin$'],[Ts,Idyn,'r-','$dyn$']]
    dsp.stddisp(plts,labs=[r'$z(\AA)$','Integrated $I$'],**kwargs)

    #plot details
    if isinstance(iTs,int):
        if iTs :
            iTs=slice(0,None,iTs)
        else:
            iTs = slice(0,0)
    pltsD, cs = [], dsp.getCs('jet',Ts[iTs].size)
    for i,t in enumerate(Ts[iTs]):
        sg = np.linspace(-1,1,2000)*nt/t
        Ik = Ikf(t,sg,xi)
        Id = Idf(t,sg,xi)

        pltsD += [[sg,Id,[cs[i],'-'],r'$z=%.1f\AA$' %t]]

        plts = [[sg,Id,'r-','$dyn$']]
        # plts += [[sg,Ik,'b-o','$kin$']]
        dsp.stddisp(plts,labs=[r'$s_g$','$I$'],lw=2,title=r'$z=%f\AA$' %t)
    # dsp.stddisp(pltsD,labs=[r'$s_g$','$I$'],lw=2)#,xylims=[-1,1,0,1])

def plot_rocking_curves(Ts=np.linspace(1,5000,10),xi=500,**kwargs):
    pltsD,cs = [], dsp.getCs('viridis',Ts.size)
    for i,t in enumerate(Ts):
        sg = np.linspace(-0.005,0.005,1000)
        pltsD += [[sg,Idf(t,sg,xi),[cs[i],'-'],r'$z=%.1f\AA$' %t]]
    dsp.stddisp(pltsD,labs=[r'$s_g$','$I$'],**kwargs)


integrate_th_full(Ts=np.linspace(1,3000,100),xi=500,iTs=0,
    xylims=[0,5000,0,0.005],lw=2,name=path+'Iint_kin_dyn.svg',opt='p')
# plot_rocking_curves(Ts=np.array([0.75,1,1.25,1.5])*500,xi=500,
#     xylims=[-0.005,0.005,0,1],lw=2,name=path+'Iint_kin_dynR.svg',opt='ps')
