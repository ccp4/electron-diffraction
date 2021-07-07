from utils import*                  ;imp.reload(dsp)
from EDutils import utilities as ut ;imp.reload(ut)
from blochwave import bloch as bl   #;imp.reload(bl)
import tifffile
plt.close('all')
opts='S' #S(Solve) s(single)

path = 'dat/bloch/'
tag = 'rx_'
nframes = 250
deg = 0.25
theta0,phi0 = 0,0 #12,27
thick = 1000

pts_file = 'dat/pets/diamond.pts'
aperpixel = 0.01
omega = deg*(nframes-1)

bloch_args = {'cif_file':'diamond','keV':200,
    'Nmax':8,'Smax':0.05,'thick':thick,
    'path':path,'solve':1,'opts':'s'}
tiff_args = {'Imax':7e4,'Nmax':256,'aperpixel':aperpixel,
    'gs3':0.04,'rmax':8e-5}

i0=5 #which simu when running single mode
frame = lambda i: '%s%s' %(tag,str(i+1).zfill(3))

def get_uvw(theta0,phi0,e1=[0,0,1],plot=0):
    u0 = ut.u_from_theta_phi(theta0,phi0)
    u1 = np.cross(e1,u0);u1/=np.linalg.norm(u1)
    uvw = ut.get_uvw(u0,u1,omega,nframes)
    if plot:ut.show_uvw(uvw,h3d=True)
    return uvw
def solve_exp(uvw):
    return [bl.Bloch(u=u,name=frame(i),**bloch_args) for i,u in enumerate(uvw)]
def load_exp():
    return [bl.load_Bloch(file=path+frame(i)+'.pkl') for i in range(nframes)]
def convert2tiff(thick,**kwargs):
    for i,b in enumerate(blochs):
        tiff_file=path+'tiff/'+frame(i)+'.tiff'
        b.set_thickness(thick)
        b.convert2tiff(tiff_file,**kwargs)

uvw = ut.get_uvw(u0=[0,0,1],u1=[0,1,0],omega=omega,nframes=nframes)
# uvw = get_uvw(theta0,phi0,plot=0)
if 's' in opts:
    if 'S' in opts:
        bloch = bl.Bloch('diamond',u=uvw[0],keV=200,
            Nmax=8,Smax=0.1,thick=thick,
            name=frame(i0),path=path,solve=1,opts='sv')
    else:
        bloch = bl.load_Bloch(file=path+'%s.pkl' %frame(i0))
    bloch.convert2tiff(path+'tiff/%s.tiff' %frame(i0),**tiff_args)
    im = tifffile.imread(path+'tiff/%s.tiff' %frame(i0))
    dsp.stddisp(im=[im],pOpt='im',caxis=[0,25],cmap='binary_r',bgcol='k')
    bloch.get_Istrong(m={'I':1e4})

    # from EDutils import viewers as vw       ;imp.reload(vw)
    # v = vw.Viewer(bloch=bloch)

else:
    if 'S' in opts:
        blochs = solve_exp(uvw)
    else:
        blochs = load_exp()
    convert2tiff(thick=thick,**tiff_args)
    lat_params = blochs[i0].crys.lattice_parameters
    ut.make_pets(pts_file,aperpixel,deg=deg,ref_cell=lat_params,tag=tag)
