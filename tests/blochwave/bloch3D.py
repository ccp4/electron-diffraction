from utils import*
from EDutils import utilities as ut
from blochwave import bloch         ;imp.reload(bloch)
from blochwave import bloch_pp as bl ;imp.reload(bl)
cif_file='dat/LTA.cif'

if 0:
    b0 = bloch.Bloch(cif_file,path='dat/dummy/',keV=200,
        u=[0,0,1],Nmax=4,Smax=0.05,
        opts='vst',thick=100,solve=False)
    # df_Fhkl = pd.read_pickle(b0.get_Fhkl_pkl())
    # b0.show_df_G(n=10)


if 1:
    u  = np.array([0.66865363, 0.29164837, 0.6521759 ])
    u  = u/np.linalg.norm(u)
    uvw = ut.get_uvw(u,osc=0.1,npts=3)
    name='test'

    rock = bl.Bloch_cont(path='dat/LTA/rocks/%s' %name,uvw=uvw,tag='',
        Sargs=dict(cif_file=cif_file,Smax=0.01,keV=200,Nmax=10,solve=True),
        #params = ['hkl'],
        #vals   = [ [np.array([[0,0,0],[-1,-7,4]])]*len(uvw)],
        frames=np.arange(len(uvw)) )


    # hkl_f=rock.get_full_refl(Swm=0.01)
    # rock.show_beams(hkl=hkl_f, cols=['Sw_min','Sw_max'])
    # rock=ut.load_pkl('dat/LTA/rocks/test/rock_.pkl')
    # rock.show_excitation_map(vm=0.001,nb_max=50,
        # axpos=[0.1,0.3,085,0.6],xylims=['x',0,50],
        # );plt.show()
