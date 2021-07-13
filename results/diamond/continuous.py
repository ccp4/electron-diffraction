from utils import*                  ;imp.reload(dsp)
from EDutils import utilities as ut ;imp.reload(ut)
from blochwave import bloch_pp as bl;imp.reload(bl)
import tifffile
# plt.close('all')
path = 'dat/bloch/'

opt='ps'
opts='Z' #RI' #S(Solve) s(single)
nframes = 200
deg  = 0.01

omega = np.arange(nframes)*deg
bloch_args = {'cif_file':'diamond','keV':200,
    'Nmax':8,'Smax':0.025,'solve':1,'opts':'sv'}


# uvw = ut.get_uvw(u0=[0,0,1],u1=[0,1,0],omega=omega)
# rock = bl.bloch_rock(tag='rx',uvw=uvw,
#     omega=omega,bloch_args=bloch_args,
#     thicks=(0,800,400),ts0=0,cond='(Sw<1e-3) &(Vga>1e-3)',
#     #zs = np.arange(168,800,168), hkls=[[[0,4,0]]],
#     path=path,opts='W',opt='p')



uvw = ut.get_uvw_from_theta_phi(theta=12,phi=19,omega=omega,plot=0)

# cond = '(Sw<2e-3) & (I>2e-2) & (Vga>1e-6)'
# cond = '(Sw<1e-4) & (I>1e-4) '
fz = lambda x : -np.log10(np.maximum(x,1e-10))
cond = ''#'(Sw<1e-2) & (Vga>1e-6)'
# hkls = [[-8,-2,2],[-5,1,1],[-1,7,-1],[3,7,-1],[-2,6,0]]
hkls = [[-8,-2,2],[-5,1,1],[-1,-7,1],[3,7,-1],[-2,6,0],[7,-3,-1],[5,1,-1]]
# rock = ut.load_pkl('dat/bloch/rock_diamond_r1.pkl')
# rock.integrate_rocking(cond=cond,refl=hkls,lw=2)
hkls = [[h] for h in hkls]
rock = bl.bloch_rock(tag='diamond_r1',uvw=uvw,
    omega=omega,bloch_args=bloch_args,
    thicks=(0,800,400),ts0=2.0,cond=cond,fz=fz,thick=250,
    zs=np.arange(1,6)*50,hkls=hkls,
    path=path,opts=opts,opt=opt)


# rock.set_tag('diamond_r1')

# if 's' in opts:
#     if 'S' in opts:
#         bloch = bl.Bloch('diamond',u=uvw[0],keV=200,
#             Nmax=8,Smax=0.1,thick=thick,
#             name=frame(i0),path=path,solve=1,opts='sv')
#     else:
#         bloch = bl.load_Bloch(file=path+'%s.pkl' %frame(i0))
#     bloch.convert2tiff(path+'tiff/%s.tiff' %frame(i0),**tiff_args)
#     im = tifffile.imread(path+'tiff/%s.tiff' %frame(i0))
#     dsp.stddisp(im=[im],pOpt='im',caxis=[0,25],cmap='binary_r',bgcol='k')
#     # bloch.get_Istrong(m={'I':1e4})

    # from EDutils import viewers as vw       ;imp.reload(vw)
    # v = vw.Viewer(bloch=bloch)





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
