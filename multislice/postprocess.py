from utils import *
from scattering_factors import get_elec_atomic_factors,wavelength
from structure_factor import structure_factor3D
from utils import displayStandards as dsp
import pickle

### import from TEMSIM
def plot_pattern_section(file='si_autoslic100.tif'):
    im=plt.imread(datpath+file)
    N=im.shape[0]
    stddisp([range(N),im[:,0],'b','$I_{h0}$'])

def import_beams(file,slice_thick=1,iBs=[],tol=1e-2,Obeam=False,pImax=False):
    '''
    slice_thick : slice thickness used for simulation
    iBs   : selected beam indices : default=Ibeams.max()>tol
    Obeam : include origin
    '''
    hk    = open(file).readline().rstrip().split('=  ')[1].split(' ');#print(hk)
    beams = np.loadtxt(file,skiprows=3).T
    if isinstance(slice_thick,list) :
        idx,beams = beams[0,:],beams[1:,:]
        n_slices = int(len(idx)/len(slice_thick)); #print(n_slices,len(slice_thick))
        t = np.cumsum(np.tile(slice_thick,n_slices));#print(t[:4])
    else :
        idx,beams = beams[0,:-2],beams[1:,:-2]
        t = idx*slice_thick
    #get re,im,I
    nbs = len(hk)
    re = [beams[2*i,:] for i in range(nbs)]
    im = [beams[2*i+1,:] for i in range(nbs)]
    Ib = [np.abs(r+1J*i)**2 for r,i in zip(re,im)]
    #filter
    Imax = np.array([I.max() for I in Ib]); #print(Imax,Imax.max())
    if not iBs : iBs = [i for i in range([1,0][Obeam],nbs) if Imax[i]>= tol*Imax.max()]
    if pImax : print('Imax:',Imax,'iBs:',iBs,', nbs=',nbs)
    hk,t,re,im,Ib = np.array(hk)[iBs],t, np.array(re)[iBs],np.array(im)[iBs],np.array(Ib)[iBs]
    return hk,t,re,im,Ib

#### utilities
def load_multi_obj(filename):
    '''load a saved Multislice object
    filename : pickle file (.pkl)  '''
    with open(filename,'rb') as f : multi = pickle.load(f)
    return multi

### def :display
def plot_beam_thickness(beams,rip='I',cm='Greens',**kwargs):
    ''' plot the beams as function of thickness
    - beams : beams info from import_beams or file containing them
    - rip flags : I(Intens),r(real),i(imag)
    '''
    hk,t,re,im,Ib = beams
    nbs = len(hk)
    csp,csr,csi,plts = getCs(cm,nbs),getCs('Blues',nbs),getCs('Reds',nbs),[]
    for i in range(nbs):
        if 'I' in rip : plts += [[t,Ib[i],[csp[i],'-'],'$I_{%s}$' %(hk[i])]]
        if 'r' in rip : plts += [[t,re[i],csr[i],'$re$']]
        if 'i' in rip : plts += [[t,im[i],csi[i],'$im$']]
    dsp.stddisp(plts,lw=2,labs=['$thickness(\AA)$','$I_{hk}$'],**kwargs)#,xylims=[0,t.max(),0,5])

### misc
def plot_fe(qmax=1,Nx=20,Vcell=1,**kwargs):
    q0,fq_e = get_elec_atomic_factors(['Si'],q=np.linspace(0,qmax,1000))
    vg = fq_e[0]/Vcell
    vgmax = vg.max()
    Xq =np.arange(Nx)[:,None].dot(np.array([1,1])[None,:]).T/ax
    Yq = np.tile([0,vgmax],[Nx,1]).T
    plts=[[Xq,Yq,[(0.7,1.0,0.7),'--'],''],
          [q0,vg,'g--','$v_g(Si)$']]
    stddisp(plts,labs=['$q(A^{-1})$','$v_g(q)(A^{-2})$'], xylims=[0,qmax,0,vgmax],
        legLoc='upper right',**kwargs)


###########################################################################
#### def : test
###########################################################################

if __name__ == "__main__":
    plt.close('all')
    figpath=get_figpath(__file__,'/figures/')
    datpath=get_figpath(__file__,'/dat/Silicon/Si110/')
