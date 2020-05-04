from utils import*
datpath=get_figpath(__file__,'/dat/')
figpath=get_figpath(__file__,'/figures/')

# import from slicelib
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

if __name__ == "__main__":
    plt.close('all')
    # import_structure_factor(['si110_2x2.txt'],imopt='c',name=figpath+'si110_S_2D.png',opt='s')
    # import_structure_factor(['si110_2x2.txt','si110_10x10.txt'],opt='s',name=figpath+'si110_S_1D.svg',gridOn=0  )
    #     inset={'axpos':[0.2,0.18,0.25,0.25],'xylims':[0,150,0,0.02]})
    # plot_fe(qmax=3,Nx=20,opt='s',name=figpath+'Si_vg.svg',axPos=[0.2,0.125,0.75,0.8],setPos=1,lw=3,gridOn=0)
    import_vatom(datpath+'Si_vatom.txt',opt='s',name=figpath+'Si_va.svg',axPos=[0.2,0.125,0.75,0.8],setPos=1,lw=3)
    # import_vzatom(['Si_vzatom.txt','Si_vzatomLUT.txt'],labels=['$Si_{vz}$','$spline_{fit}$'],opt='p',name=figpath+'Si_vz.svg',axPos=[0.2,0.125,0.75,0.8],setPos=1,lw=3)
