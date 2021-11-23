import numpy as np

def load_Bloch(path='',tag='',file='',v=1):
    """load a saved Bloch object
    filename : pickle file (.pkl)  """
    files = [f for f in glob.glob(os.path.join(path,'*.pkl')) if tag in f]
    if tag:
        if len(files)==1:
            file = files[0]
    if file:
        with open(file,'rb') as f : obj = pickle5.load(f)
        if v:print(colors.green+'loaded:' +colors.yellow+file+colors.black);
        return obj
    else:
        if len(files):
            print(colors.red+'mulitple simus available with tag %s:' %tag);
            print(colors.yellow,[os.path.basename(f) for f in files],colors.black)
        else:
            print(colors.red+'no simus available with in "%s" with tag %s' %(path,tag))

def get_fz(opts):
    """
    mapping for function to apply

    Parameters
    ---------
    opts
        r(real) i(imag) a(angle) m(mag) l(log10(|Fhkl|+1)) 2(mag^2) L(logM)
    """
    if not opts:opts='m'
    keys   =['r','i','a','m','l','2','L']
    fs     = [np.real,np.imag,np.angle,np.abs,logF,abs2,logM]
    fs_str = ['real part','imag part','phase','magnitude','$\log_{10}(|F|+1)$','$|F|^2$','$-\log_{10}$']
    fz     = dict(zip(keys,fs))[opts]
    fz_str = dict(zip(keys,fs_str))[opts]
    return fz,fz_str

abs2 = lambda F:np.abs(F)**2
logF = lambda F:np.log10(np.abs(F)+1)
logM = lambda F:-np.log10(np.maximum(F,1e-8))


# def show_beams(self,F='I',fopts='m',opts='xN',mag=500,cutoff=0,cmap='Greens',**kwargs):
#     """Display beam values
#     - F : 'I'(intensity),'S'(scattered beams),'Vg'(potential)
#     - opts  : 'x'(include projection of x axis), 'N'(normalize)
#     - fopts : see get_fz
#     """
#     fz,fz_str = get_fz(fopts)
#     F_str = {'L':'Lattice','Vg':'Potential $V_{g}$',
#         'S' :'Scattered beams Bloch $S_{g}$',
#         'I' :'Intensity Bloch$I_{g}$',
#         'Sg':'Scattered beams Kinematic $S_{g,kin}$',
#         'Ig':'Intensity kinematic $I_{g,kin}$',
#         'Sw':'Excitation error $\zeta_{g}$'}[F]
#
#     tle = '%s, %s, thickness=%d$\AA$'  %(F_str,fz_str,self.thick)
#     px,py = self.df_G[['px','py']].values.T
#     Fvals = fz(self.df_G[F])
#     if 'N' in opts:Fvals/=Fvals.max()
#     # mag /= Fvals.max()
#     scat  = [px,py,Fvals*mag/Fvals.max(),Fvals]
#     plts = [[0,0,'b+','']]
#     if 'x' in opts:plts+=[ [[0,self.e0x],[0,0],'r','']]
#     if not cutoff:cutoff = Fvals.max()
#     # print(Fvals.max())
#     fig,ax = dsp.stddisp(plts,lw=2,scat=scat,caxis=[0,cutoff],cmap=cmap,cs='S',title=tle,**kwargs)
