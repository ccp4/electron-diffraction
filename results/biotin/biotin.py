from utils import*                   ;imp.reload(dsp)
from EDutils import utilities as ut  ;imp.reload(ut)
from EDutils import pets as pt       ;imp.reload(pt)
from EDutils import display as EDdisp;imp.reload(EDdisp)
from blochwave import bloch_pp as bl ;imp.reload(bl)
import os,mrcfile
from subprocess import check_output
plt.close('all')

path='dat/simus'
tag=''

# Thicknesses to simulate (Ang)
thicks  = np.arange(50,501,50)[:1]
# selected frames to simulate (None for all)
# frames  = None
frames  = np.array([1,2,3])
# number of intermediate simulations
npts    = 1



#### simulation info
dials = pt.Dials('dat/dials')
nframes = dials.n_frames
if type(frames)==type(None):frames = np.arange(nframes)+1
alpha0  = dials.alpha[frames-1]
uvw0    = dials.uvw0[frames-1]
uvw0 /= np.linalg.norm(uvw0,axis=0)
uvw = uvw0
# uvw = ut.uvw_add_points(uvw0,npts=1,plot=0)
# ut.uvw_add_points(uvw0,npts=npts,plot=1)
alpha   = np.linspace(alpha0[0],alpha0[-1],uvw.shape[0])

#### simulation parameters
bloch_args = {'cif_file':'biotin.cif',
    'Smax':0.1,'Nmax':8,'solve':1,'keV':200}
rock_file=path+'/rock_%s.pkl' %tag

### kernel and simulated frames parameters
with mrcfile.open("dat/dials/biotin_xtal1/biotin_xtal1_0000.mrc") as f:
    aper = f.extended_header["Pixel size X"][0]*1e-10 #A^-1
    nxy = int(f.header.nx)

dqx,dqy=aper,aper
fbroad=lambda r2:np.exp(-r2**0.8/0.001)
# fbroad=lambda r2:np.exp(-r2**0.8/0.0001)
nx,ny=(35,)*2
ix,iy = np.meshgrid(range(-nx,nx+1),range(-ny,ny+1))
x,y = ix*dqx,iy*dqy
r2 = (x**2+y**2)
# dsp.stddisp(im=[x,y,fbroad(r2)])


tiff_args=dict(fbroad=fbroad,nX=nx,#gs3=0.02,
    # Imax=1e3,
    pred=1,
    #aperpixel=aper,Nmax=nxy,rot=37,
    # 'tif_writer_args':{'description':"ImageCameraName: timepix"}}
    )


opts='Ias' #Solve(S) I(save images) a(all) s(sum) n(new)
if 'S' in opts:
    if 'n' in opts and os.path.exists(path):
        check_output("rm -rf %s" %path,shell=True).decode()
    rock = bl.Bloch_cont(path=path,tag=tag,uvw=-uvw,Sargs=bloch_args)

##sum images
if 'I' in opts:
    rock = ut.load_pkl(file=rock_file)
    if not 'a' in opts and not 's' in opts:opts+='as'
    for i,thick in enumerate(thicks):
        # figpath=os.path.join(rock.path,'tiff','%dA' %   thick)
        figpath=os.path.join(rock.path,'simu_dials','%dA' %   thick)

        if 'a' in opts:
            if os.path.exists(figpath) and 'n' in opts:
                check_output("rm -rf %s/*" %figpath,shell=True).decode()
            rock.make_img(fmt='mrc',figpath=figpath,
                exp=dials,
                template="dat/dials/biotin_xtal1/biotin_xtal1_0000.mrc",
                thick=thick,
                frames=np.arange(1,3),
                **tiff_args)
        if 's' in opts:
            sum_path = figpath+'/sum'
            if os.path.exists(sum_path) and 'n' in opts:
                print(check_output("rm -rf %s" %(sum_path),shell=True).decode())
            rock.sum_images(figpath=figpath,n=1,
                frames=(1,1),
                )
