from utils import*                   ;imp.reload(dsp)
from EDutils import utilities as ut  ;imp.reload(ut)
from EDutils import pets as pt       ;imp.reload(pt)
from EDutils import display as EDdisp;imp.reload(EDdisp)
from blochwave import bloch_pp as bl ;imp.reload(bl)
import os,mrcfile
from subprocess import check_output
plt.close('all')

opts='Ias' #Solve(S) I(save images) a(all) s(sum)
path='dat/simus'
tag=''

# Thicknessses to simulate (Ang)
thicks  = np.arange(50,501,50)[:1]
# selected frames to simulate (None for all)
frames  = None
# frames  = np.array([1,2,3])
# number of intermediate simulations
npts    = 3


#### simulation info
dials = pt.Dials('dat/dials')
nframes = dials.n_frames
if type(frames)==type(None):frames = np.arange(nframes)+1
alpha0  = dials.alpha[frames-1]
uvw0    = dials.uvw0[frames-1]
uvw     = ut.uvw_add_points(uvw0,npts=npts,plot=0)
# ut.uvw_add_points(uvw0,npts=npts,plot=1)
alpha   = np.linspace(alpha0[0],alpha0[-1],uvw.shape[0])

#### simulation parameters
bloch_args = {'cif_file':'biotin.cif',
    'Smax':0.1,'Nmax':8,'solve':1,'keV':200}
rock_file=path+'/rock_%s.pkl' %tag

#### kernel and simulated frames parameters
fbroad=lambda r2:np.exp(-r2**0.8/0.001)
with mrcfile.open("dat/dials/biotin_xtal1/biotin_xtal1_0000.mrc") as f:
    aper = f.extended_header["Pixel size X"][0]*1e-10 #A^-1
    nx = int(f.header.nx)
tiff_args={'fbroad':fbroad,'gs3':0.05,'nX':25,'Imax':5e6,
    'rot':30,'aperpixel':aper,'Nmax':nx,
    'tif_writer_args':{'description':"ImageCameraName: timepix"}}


if 'S' in opts:
    rock = bl.Bloch_cont(path=path,tag=tag,uvw=-uvw,Sargs=bloch_args)

##sum images
if 'I' in opts:
    rock = ut.load_pkl(file=rock_file)
    if not 'a' in opts and not 's' in opts:opts+='as'
    for i,thick in enumerate(thicks):
        figpath=os.path.join(rock.path,'tiff','%dA' %   thick)
        if not os.path.exists(figpath):os.mkdir(figpath)

        if 'a' in opts:
            print(check_output("rm -rf %s/*" %figpath,shell=True).decode())
            rock.convert2tiff(figpath=figpath,thick=thick,**tiff_args)
        if 's' in opts:
            sum_path = figpath+'/sum'
            if os.path.exists(sum_path):
                print(check_output("rm %s/*" %(sum_path),shell=True).decode())
            rock.convert2tiff(figpath=figpath,n=npts+1)#,nmax=0)
