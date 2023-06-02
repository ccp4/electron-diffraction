import numpy as np
import mrcfile,tifffile
from utils import glob_colors as colors
import os,glob,pickle5

def mrc2tiff(mrc_file,outpath):
    tiff_file = os.path.basename(mrc_file).replace('.mrc','.tiff')
    im = imread(mrc_file).copy()
    im -= im.min()
    im*=50/im.mean()
    # im*=1e4/im.max()
    # im = np.array(im,dtype=np.uint16)
    imwrite(os.path.join(outpath,tiff_file),im)

#### writer
def mrc_writer(filename,im0):
    mrc = mrcfile.open(filename,'r+')
    mrc.set_data(np.array(im0,dtype=type(mrc.data[0,0])))
    mrc.flush()
    mrc.close()

def tiff_writer(tiff_file,im):
    tifffile.imwrite(tiff_file,im)
img_writers = {
    'mrc' : mrc_writer,
    'tiff': tiff_writer,
}
def imwrite(filename,im):
    fmt=filename.split('.')[-1]
    img_writers[fmt](filename,im)
    print(colors.yellow+filename+colors.green+' saved'+colors.black)


#### reader
def tiff_reader(tiff_file)  :
    with open(tiff_file,'rb') as f:
        I =tifffile.imread(f)
    return I
def mrc_reader(mrc_file):
    with mrcfile.open(mrc_file) as mrc:
        return mrc.data
img_readers = {
    'mrc' : mrc_reader,
    'tiff': tiff_reader,
}
fmts = list(img_readers.keys())#['mrc','tiff']
def imread(filename):
    fmt=filename.split('.')[-1]
    return img_readers[fmt](filename)


def strong_beams(dfG,
        tol:float=1e-2,n:int=10,
        opt:str='F'):
    """keep the n strongest beams, i.e. those with intensity I>tol not including central beam

    Parameters
    ----------
    dfG
        the bloch dataframe
    tol
        strong beam criteria
    n
        number of strong beams
    opt
        include Origin(O) Friedel pairs(F)
    Returns
    ----------
    x
        list of beams
    """
    dfM=dfG.copy()
    if not 'O' in opt:dfM = dfM.drop(str((0,0,0)))
    dfM = dfM.sort_values('I',ascending=False)[:2*n]

    df = dfM.loc[dfM.I>tol]
    if not 'F' in opt:
        df=remove_friedel_pairs(df)

    if df.shape[0]>n:df = df[:n]

    # print(dfM)
    return list(df.index.values)

def remove_friedel_pairs(dfM):
    refl=[]
    for h,hkl in zip(dfM.index,dfM[['h','k','l']].values):
        if not str(tuple(-hkl)) in refl:
            refl.append(h)
    print('removing Friedel pairs')
    dfM=dfM.loc[refl]
    return dfM

def get_inp(npx,nbeams,u,keV,thicks,out=None):
    inp = """# Input file for Felix version :
# ------------------------------------

# ------------------------------------

# control flags
IWriteFLAG                = 0
IImageFLAG                = 1
IScatterFactorMethodFLAG  = 0
IBlochMethodFLAG          = 0
IMaskFLAG                 = 1
IHolzFLAG                 = 0
IAbsorbFLAG               = 0
IAnisoDebyeWallerFlag     = 0
IByteSize                 = 8
    """

    inp+="""
# radius of the beam in pixels
IPixelCount               = %d

# beam selection criteria
IMinReflectionPool        = %d
IMinStrongBeams           = %d
IMinWeakBeams             = 0
    """ %(npx,nbeams,nbeams)

    inp+="""
# crystal settings
RDebyeWallerConstant      = 0.0
RAbsorptionPer            = 7.8

# microscope settings
ROuterConvergenceAngle    = 3.03
IIncidentBeamDirection    = [%f,%f,%f]
IXDirection               = [1,1,0]
INormalDirection          = [-1,1,0]
RAcceleratingVoltage (kV) = %.1f
RAcceptanceAngle          = 0.0
    """ %(u[0],u[1],u[2],keV)

    inp+="""
# Image Output Options
RInitialThickness        = %.1f
RFinalThickness          = %.1f
RDeltaThickness          = %.1f
IReflectOut              = 1
    """ %(thicks[0],thicks[1],thicks[2])

    inp+="""
#Refinement Specific Flags
IRefineModeFLAG          = S
IWeightingFLAG           = 0
IRefineMethodFLAG        = 3
ICorrelationFLAG         = 2
IImageProcessingFLAG     = 4
RBlurRadius              = 0.0
INoofUgs                 = 40
IAtomicSites             = (1,2)
IPrint                   = 0
RSimplexLengthScale      = 20.0
RExitCriteria            = 0.00001
    """

    if out:
        with open(out,'w') as f:f.write(inp)
    else:
        return inp

def load_bloch(path='',tag=None,file='',v=1):
    """load a saved Bloch object from
    path : path to simulation folder
    tag : identifier for a bloch simu in case there are several possibilities
        OR
    filename : directly provide the pickle file (.pkl)
    """
    # first tr toget a simu from the tag
    files = [f for f in glob.glob(os.path.join(path,'*.pkl'))]
    if not isinstance(tag,type(None)):
        if isinstance(tag,int):
            file = files[tag]
        else:
            files = [f for f in glob.glob(os.path.join(path,'*.pkl')) if tag in f]
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
