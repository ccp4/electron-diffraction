from utils import*                  ;imp.reload(dsp)
from EDutils import utilities as ut ;imp.reload(ut)
from blochwave import bloch_pp as bl;imp.reload(bl)
import tifffile
plt.close('all')
path = 'dat/bloch/'

opt='p'
opts='R' # #ts:S(Setup) I(Iz), W(Excitation) zs,hkls:R(Rocking) Z(Iintegrated(Z)) s(solve)
nframes = 50
deg = 0.02

bloch_args = {'cif_file':'diamond','keV':200,'thick':200,
    'Nmax':8,'Smax':0.025,'solve':1,'opts':'sv'}

omega = np.arange(nframes)*deg
uvw = ut.get_uvw_from_theta_phi(theta=12,phi=19,omega=omega,plot=0)

cond = '(Sw<2e-2) & (Vga>1e-6) & (I>1e-3)'
# cond = '(Sw<1e-4) & (I>1e-4) '
# fz = lambda x : -np.log10(np.maximum(x,1e-10))
# cond = ''#'(Sw<1e-2) & (Vga>1e-6)'
# hkls = [[-8,-2,2],[-5,1,1],[-1,7,-1],[3,7,-1],[-2,6,0]]
refl = [[-8,-2,2],[-5,1,1],[-1,-7,1],[3,7,-1],[-2,6,0],[7,-3,-1],[5,1,-1]]
hkls = [[refl[1]]]
rock = bl.bloch_rock(tag='diamond_rand',uvw=uvw,
    omega=omega,bloch_args=bloch_args,
    thicks=(0,800,400),
    ts0=1.0,refl=refl,cond=cond,
    zs=np.arange(1,6)*50,hkls=hkls,
    Z_args={},#{'new':1},
    path=path,opts=opts,opt=opt)


# rock = ut.load_pkl('dat/bloch/rock_diamond_rand.pkl')
# b = rock.load(0)
# b.get_beam(cond=cond)




def convert2tiff(thick,**kwargs):
    for i,b in enumerate(blochs):
        tiff_file=path+'tiff/'+frame(i)+'.tiff'
        b.set_thickness(thick)
        b.convert2tiff(tiff_file,**kwargs)

if 'c' in opts:
    pts_file = 'dat/pets/diamond.pts'
    aperpixel = 0.01
    tiff_args = {'Imax':7e4,'Nmax':256,'aperpixel':aperpixel,
        'gs3':0.04,'rmax':8e-5}

    convert2tiff(thick=thick,**tiff_args)
    lat_params = blochs[i0].crys.lattice_parameters
    ut.make_pets(pts_file,aperpixel,deg=deg,ref_cell=lat_params,tag=tag)
