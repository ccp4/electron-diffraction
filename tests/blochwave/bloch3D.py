from utils import*
from EDutils import utilities as ut
from blochwave import bloch         ;imp.reload(bloch)
from blochwave import bloch_pp as bl ;imp.reload(bl)
cif_file='dat/LTA.cif'
dsp.matplotlib.use('gtk3agg')
if 1:
    args = dict(opts='')
    b0 = bloch.Bloch(cif_file,path='dat/dummy/',keV=200,
        u=[1,6,15],Nmax=8,Smax=0.01,solve=True,
        opts='svt',thick=500)

    # df_Fhkl = pd.read_pickle(b0.get_Fhkl_pkl())
    # b0.show_df_G(n=10)
    # print(b0.df_G.sort_values('Ikin')[['Swa','xi_g','Ikin']])
    if 1:
        b0.solve_strong_beam(thick=500,solve_args=args,Imin=5e-4,dImin=0.05)
        print(b0.df_beams.sort_values('Iref')[['Inew','Iref','type']].to_string())
        # print('strong   :\n',b0.df_hkl_strong.sort_values('I')[['I','Iref']])
        # print('weak in  :\n',b0.df_hkl_weak.loc[b0.df_hkl_weak['in']==True ,['dImax','dhmax','Iref','Inew']])
        # print('weak out :\n',b0.df_hkl_weak.loc[b0.df_hkl_weak['in']==False,['dImax','dhmax','Iref','Inew']])


if 0:
    u  = np.array([0.66865363, 0.29164837, 0.6521759 ])
    u  = u/np.linalg.norm(u)
    uvw = ut.get_uvw(u,osc=0.1,npts=3)
    name='test'
    if 0:
        rock = bl.Bloch_cont(path='dat/LTA/rocks/%s' %name,uvw=uvw,tag='',
            Sargs=dict(cif_file=cif_file,Smax=0.01,keV=200,Nmax=10,solve=False,opts='s0',v=0),
            # params = ['hkl'],
            # vals   = [ [np.array([[0,0,0],[-1,-7,4]])]*len(uvw)],
            frames=np.arange(len(uvw)),verbose=True)

        rock.do('solve_strong_beam',verbose=True,thick=500,solve_args=dict(opts=''),
            Imin=5e-4,dImin=0.05)

    if 1:
        rock=ut.load_pkl('dat/LTA/rocks/test/rock_.pkl')
        # rock.do('_set_beams_vs_thickness',verbose=True,thicks=(100,500,100),v=0)

        # rock.plot_rocking(refl=rock.beams.index.tolist()[:3],kin=True);plt.show()
        rock.plot_rocking(refl=rock.beams.index.tolist()[:1],zs=[100,300,500],cmap='jet',kin=True);plt.show()
        # hkl_f=rock.get_full_refl(Swm=0.01)
        # rock.show_beams(hkl=hkl_f, cols=['Sw_min','Sw_max'])
        # rock.show_excitation_map(vm=0.001,nb_max=50,
        #     axpos=[0.1,0.3,085,0.6],xylims=['x',0,50],
        #     )
        # plt.show()
