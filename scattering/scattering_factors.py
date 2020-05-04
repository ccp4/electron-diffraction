from utils.displayStandards import*
import pandas as pd
from utils.physicsConstants import mc2,c0,A,h,emass,keV

dat_path = get_figpath(__file__,'/data/')
__all__ = ['get_xray_atomic_factors','get_elec_atomic_factors',
    'get_v_from_KE','wavelength','get_elt']


def get_xray_atomic_factors(elts,qmax=2,npts=100):
    '''elts : list of atoms by name
    The formula on the uses the semi scattering angle
    '''
    q = np.linspace(0,min(qmax,6),npts) #A
    #q2 = ((q/(4*pi))**2)[:,None]
    q2 = ((q/2)**2)[:,None]
    dfx = pd.read_pickle(dat_path+'xray.pkl').loc[elts].astype(float)
    Ai = dfx[['a%d' %i for i in range(1,5)]].values
    Bi = dfx[['b%d' %i for i in range(1,5)]].values
    C  = dfx['c'].values
    fx_q = [(ai*np.exp(-bi*q2)).sum(axis=1) + c  for ai,bi,c in zip(Ai,Bi,C)]
    return q,fx_q

def get_elec_atomic_factors(elts,q=None,qmax=2,npts=100):
    '''
    elts : list of atoms by name 'E' or atomic number Z
    opt : G(gauss fit), L(Lenz)
    '''
    if isinstance(elts[0],str) :
        Z = pd.read_pickle(dat_path+'elec.pkl')['Z'][elts].values
    else:
        Z=elts
    if not isinstance(q,np.ndarray) : q = np.linspace(0,qmax,npts)
    q2,fq_e,nelts = q**2,[],len(Z)

    fparams = np.load(dat_path+'abcd.npy',allow_pickle=True)
    for elt in range(nelts) :
        a,b=np.reshape(fparams[Z[elt],:6],(3,2)).T
        c,d=np.reshape(fparams[Z[elt],6:],(3,2)).T; #print(a,b,c,d)
        fq = np.zeros(q.shape)
        for i in range(3):
             fq += a[i]/(b[i]+q2)+c[i]*np.exp(-d[i]*q2)
        fq_e+=[fq]
    return q,fq_e

# def get_lenz_model(Z,q,E0):
#     T  = E0*(1+E0/(2*mc2))/(1+E0/(2*mc2))   #keV
#     w  = wavelength(E0)                     #A
#     k0 = 2*pi/w                             #A^-1
#     theta  = q*w                            #no dim
#     theta0 = Z**(1/3)/(k0*a0)               #no dim
#     return 2*Z/(sqrt(T)*k0)/(theta**2+theta0**2)

########################################################################
####def : misc
########################################################################
def get_v_from_KE(KE):
    ''' v = get_v_from_KE(KE)
    KE : energy (keV)
    v  : speed (normalized to c)
    '''
    v =  np.sqrt(1 - 1/(1+np.array(KE)/mc2)**2)
    return v

def wavelength(KE):
    '''w = wavelength(KE)
    KE : energy (keV)
    lam  : wavelength (Angstrom)
    '''
    #emass = 510.99906   #keV
    lam = h*c0/(np.sqrt(KE*(2*emass+KE))*keV)
    return lam/A

def get_elt(elts):
    Z = np.array(pd.read_pickle(dat_path+'elec.pkl')['Z'][elts].values,dtype=int)
    return Z

################################################################################
#### def : test
################################################################################
def _test_elect_atomic(opt='p',fmt='svg'):
    elts = ['H','C','N','O','S','P']  #+ ['Si','Cu','Au','U']
    qe,fq_e=get_elec_atomic_factors(elts,qmax=3,npts=100)
    qx,fq_x=get_xray_atomic_factors(elts,qmax=6,npts=100)
    Z = get_elt(elts) #pd.read_pickle(dat_path+'atoms.pkl')['Z'][elts].values
    cs=[unicolor(0.85),unicolor(0.4),(0,0,1),(1,0,0),(0,1,1),(1,0,1)] #+ ['k']*4
    df=pd.DataFrame(np.array([Z,cs,fq_e,fq_x]).T,columns=['Z','color','fq_e','fq_x'],index=elts)


    plts=[[qe,elt.fq_e,elt.color,'$%s$' %elt.name ] for index,elt in df.iterrows()]
    stddisp(plts,labs=['$q(A^{-1})$','$f^e(A)$'],legLoc='upper right',lw=2,opt=opt,
        figsize='f',name=figpath+'electron_atomic_scattering_factors.%s' %fmt)
    plts=[[qx,elt.fq_x,elt.color,'$%s$' %elt.name ] for index,elt in df.iterrows()]
    stddisp(plts,labs=['$q(A^{-1})$','$f^x(e)$'],legLoc='upper right',lw=2,opt=opt,
        figsize='f',name=figpath+'xray_atomic_scattering_factors.%s' %fmt)
    return df

def _test_misc(KE=200,elts=['H','C','O','Si','Fe']):
    print( '''KE=%.2f,v=%.4f,lambda=%.4f
        ''' %(KE,get_v_from_KE(KE),wavelength(KE)))
    print(green+'get_elt(',elts,') = ',red,get_elt(elts), black)
if __name__=="__main__":
    figpath=get_figpath(__file__,'/figures/');#print(figpath)
    #_test_elect_atomic(opt='p',fmt='png')
    _test_misc()
