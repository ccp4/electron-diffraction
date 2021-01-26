from utils import*
from scipy.integrate import trapz


def integrate_th_full(Ts=np.linspace(1,5000,5000),xi=500,iTs=0):
    '''Integrate the full rocking curve for each thickness
    - Ts : thicknesse
    - xi : Pendullosung
    - iTs : selected indices in Ts for displaying rocking curves
    '''
    Ikin,Idyn = np.zeros(Ts.shape),np.zeros(Ts.shape)
    for iT,t in enumerate(Ts):
        w = np.linspace(-1,1,2000)*10*xi/t
        Ik = np.sinc(t/xi*w)**2
        Id = np.sinc(t/xi*np.sqrt(1+w**2))**2
        Ikin[iT] = trapz(Ik,w)
        Idyn[iT] = trapz(Id,w)
    plts = [[Ts,Ikin,'b-o','$kin$'],[Ts,Idyn,'r-o','$dyn$']]
    dsp.stddisp(plts,labs=[r'$z(\AA)$','Integrated $I$'],xylims=['y',0,5])

    #plot details
    if isinstance(iTs,int):
        if iTs :
            iTs=slice(0,None,iTs)
        else:
            iTs = slice(0,0)
    for t in Ts[iTs]:
        w = np.linspace(-1,1,2000)*10*xi/t
        Ik = np.sinc(t/xi*w)**2
        Id = np.sinc(t/xi*np.sqrt(1+w**2))**2
        plts = [[w,Ik,'b-o','$kin$'],[w,Id,'r-o','$dyn$']]
        dsp.stddisp(plts,labs=[r'$w$','$I$'],title=r'$z=%f\AA$' %t)




# integrate_th_full(Ts=np.linspace(1,5000,5000),xi=500)
# integrate_th_full(Ts=np.linspace(1,5000,5000),xi=500,iTs=450)
