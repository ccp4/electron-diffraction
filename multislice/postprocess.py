from utils import *
from scattering_factors import get_elec_atomic_factors
from utils import displayStandards as dsp
datpath=get_figpath(__file__,'/dat/Silicon/Si110/')
figpath=get_figpath(__file__,'/figures/')
ax,by,cz = 5.43,5.43,1.3575
Vcell = ax**3


### not from temsim
def plot_fe(qmax=1,Nx=20,**kwargs):
    q0,fq_e = get_elec_atomic_factors(['Si'],q=np.linspace(0,qmax,1000))
    vg = fq_e[0]/Vcell
    vgmax = vg.max()
    Xq =np.arange(Nx)[:,None].dot(np.array([1,1])[None,:]).T/ax
    Yq = np.tile([0,vgmax],[Nx,1]).T
    plts=[[Xq,Yq,[(0.7,1.0,0.7),'--'],''],
          [q0,vg,'g--','$v_g(Si)$']]
    stddisp(plts,labs=['$q(A^{-1})$','$v_g(q)(A^{-2})$'], xylims=[0,qmax,0,vgmax],
        legLoc='upper right',**kwargs)

###### Imports from main_slicelib
def import_vatom(files=[],**kwargs):
    if not isinstance(files,list) : files=[files]
    plts,cs=[],getCs('jet',len(files))
    for file,c in zip(files,cs) :
        lab = file.split('_')[0]
        r,v = np.loadtxt(datpath+file).T
        plts+=[[r,v/1e3,c,lab]]
    stddisp(plts,labs=["$r(A)$","$V_a(kV)$"],xylims=[0,0.5,0,20],**kwargs)

def import_vzatom(files=[],labels=[],**kwargs):
    if not isinstance(files,list) : files=[files]
    if not labels : labels = files
    plts,cs=[],getCs('jet',len(files))
    for file,c,lab in zip(files,cs,labels) :
        #lab = file #.split('_')[0]
        r,v = np.loadtxt(datpath+file).T
        plts+=[[r,v,c,lab]]
    stddisp(plts,labs=["$r(A)$","$V_a(V.A)$"],xylims=[0,0.5,0,5000],**kwargs)

##### Import from atompot
def import_structure_factor(files,imopt='',fe=1,**kwargs):
    if not isinstance(files,list) : files=[files]
    #unit cell spacing and atomic form factor
    if not imopt:
        nx = 20
        Xq = np.arange(nx)[:,None].dot(np.array([1,1])[None,:]).T/ax
        Yq = np.tile([0,1],[nx,1]).T
        plts=[[ [0,0],[0,1],[(0.5,1.0,0.5),'--'],'$dq_{unitcell}$']]
        plts+=[[Xq,Yq,[(0.5,1.0,0.5),'--'],'']]
        if fe:
            q0,fq_e = get_elec_atomic_factors(['Si'],q=np.linspace(0,5,1000))
            fe2 = np.abs(fq_e[0]/(ax*by))**2;fe2/=fe2.max() #A^{-1}
            plts+=[[q0,fe2,'g--','$f_e^2(Si)$']]
    cs=getCs('Blues',len(files))
    for c,file in zip(cs,files):
        lab=file.split('.')[0].split('_')[1]
        dqx,dqy=1/np.array(lab.split('x'),dtype=int)/np.array([ax,by])
        re_im = np.loadtxt(datpath+file);
        #keep only half
        N=re_im.shape[1];N2=2*int(N/2)
        re,im = re_im.reshape([2,N,N])
        Fq = re+1J*im
        S = np.abs(Fq)**2;#S/=cz*A
        S/=S.max();S[S<1e-10] = 1e-10
        if imopt :
            stddisp(im=S,imOpt=imopt,xylims=[0,30,0,30],legOpt=0,**kwargs)
        else :
            q,Sq = dqx*np.arange(N2),S[0,:N2].T #np.log10(S[0,:N2].T)
            plts+=[[q,Sq,[c,'*-'],lab]]
        #plts+=[[q,np.real(Fq),'c*-',''],[q,np.imag(Fq),'m*-','']]
    if not imopt:
        stddisp(plts,labs=['$q(A^{-1})$','$I_{h0}$'],xylims=[0,2,0,Sq.max()],
                legLoc='upper right',axPos=[0.2,0.75,0.2,0.75],**kwargs)

### import from TEMSIM
def plot_pattern_section(file='si_autoslic100.tif'):
    im=plt.imread(datpath+file)
    N=im.shape[0]
    stddisp([range(N),im[:,0],'b','$I_{h0}$'])

def plot_beam_vs_thickness(file,rip='IridO',iBs=[],tol=1e-2,cm='Greens',**kwargs):
    '''rip flags :
    I(Intens),r(real),i(imag),d(display max beam I), O(0,0 beam in default iBs)
    '''
    hk = open(file).readline().rstrip().split('=  ')[1].split(' ')
    beams = np.loadtxt(file,skiprows=3).T
    nbs = len(hk)
    t,beams = beams[0,:-5]*cz,beams[1:,:-5]
    re = [beams[2*i,:] for i in range(nbs)]
    im = [beams[2*i+1,:] for i in range(nbs)]
    Ib = [np.abs(r+1J*i)**2 for r,i in zip(re,im)]
    Imax = np.array([I.max() for I in Ib]); #print(Imax,Imax.max())

    if not iBs : iBs = [i for i in range([1,0]['O' in rip],nbs) if Imax[i]>= tol*Imax.max()]
    nbs = len(iBs);
    csp,csr,csi,plts = getCs(cm,nbs),getCs('Blues',nbs),getCs('Reds',nbs),[]
    for i,idx in zip(range(nbs),iBs):
        if 'I' in rip : plts += [[t,Ib[idx],[csp[i],'-s'],'$I_{%s}$' %(hk[idx])]]
        if 'r' in rip : plts += [[t,re[idx],csr[i],'$re$']]
        if 'i' in rip : plts += [[t,im[idx],csi[i],'$im$']]
    dsp.stddisp(plts,lw=2,labs=['$thickness(A)$','$I_{hk}$'],**kwargs)#,xylims=[0,t.max(),0,5])
    if 'd' in rip : print('Imax:',Imax,'iBs:',iBs,', nbs=',nbs)


# plt.close('all')
plot_beam_vs_thickness(datpath+'si110_autoslic_beams.txt',rip='spI',cm='jet',iBs=[],tol=1e-3,
    legLoc='center left out',setPos=True,axPos=[0.1,0.1,0.7,0.8],
    name=figpath+'Si110_Ihk.svg',opt='s',figsize='f')

#plot_pattern_section(file='si_autoslic100.tif')
#import_structure_factor(['si110_2x2.txt'],imopt='c',name=figpath+'si110_S_2D.png',opt='s')
#import_structure_factor(['si110_2x2.txt','si110_10x10.txt'],opt='s',name=figpath+'si110_S_1D.svg',gridOn=0  )
    # inset={'axpos':[0.2,0.18,0.25,0.25],'xylims':[0,150,0,0.02]})
#plot_fe(qmax=3,Nx=20,opt='s',name=figpath+'Si_vg.svg',axPos=[0.2,0.125,0.75,0.8],setPos=1,lw=3,gridOn=0)
# import_vatom('Si_vatom.txt',opt='s',name=figpath+'Si_va.svg',axPos=[0.2,0.125,0.75,0.8],setPos=1,lw=3)
# import_vzatom(['Si_vzatom.txt','Si_vzatomLUT.txt'],labels=['$Si_{vz}$','$spline_{fit}$'],opt='p',name=figpath+'Si_vz.svg',axPos=[0.2,0.125,0.75,0.8],setPos=1,lw=3)
