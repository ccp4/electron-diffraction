from utils import*                   ;imp.reload(dsp)
from EDutils import utilities as ut  ;imp.reload(ut)
from EDutils import pets as pt       ;imp.reload(pt)
from EDutils import display as EDdisp;imp.reload(EDdisp)
from blochwave import bloch_pp as bl ;imp.reload(bl)

plt.close('all')

opts='z' #Solve(S) r(save rocking curve)
path='dat/bloch_zenodo_fine'
tag=''

frames  = np.arange(10,120)  # frames to simulate
npts    = 10                # number of intermediate simulations
thick   = 330 #A

pets      = pt.Pets('dat/pets/glycine.pts',gen=0)
alpha0,uvw0 = pets.alpha[frames-1],pets.uvw0[frames-1]
# ut.uvw_add_points(uvw0,npts=npts,plot=1)
uvw = ut.uvw_add_points(uvw0,npts=npts,plot=0)
# uvw=uvw0
alpha = np.linspace(alpha0[0],alpha0[-1],uvw.shape[0])
# d_Sw = ut.get_excitation_errors(-pets.K0*uvw0[0],pets.lat_vec1,Nmax=7,Smax=0.02)
# hkl_pets = pets.rpl.loc[pets.rpl.eval('(F==%d)' %(frames[0]) ), ['hkl','F','I']]
# hkl = hkl_pets.hkl.values
# print(d_Sw.loc[hkl])


if 'S' in opts:
    bloch_args = {'cif_file':'dat/pets/alpha_glycine.cif','Smax':0.1,'Nmax':8,
        'solve':1,'thick':thick,'thicks':(0,thick,500),'keV':200}
    rock = bl.Bloch_cont(path=path,tag=tag,uvw=-uvw,Sargs=bloch_args)

if 1:
    rock_file=path+'/rock_%s.pkl' %tag
    rock = ut.load_pkl(file=rock_file)

    i = 11
    pets.show_exp(frames[i],init='rc')
    pets_simu = pt.Pets('dat/pets_simu/glycine.pts',gen=0)
    # pets_simu.show_exp(frames[i])
    vw=rock.show_tiff(sum_opt=True,cutoff=20,i=i+1,pargs={'xylims':[0,516,516,0]})#,pargs={'opt':'p'},h=1)

    # rock.set_beams_vs_thickness(thicks=(0,330,500))
    fbroad=lambda r2:np.exp(-r2**0.8/0.001)
    tiff_args={'fbroad':fbroad,'gs3':0.05,'nX':25,'Imax':5e6,
        'rot':-203.5,'aperpixel':0.00534,'iz':-1,'Nmax':516,
        'tif_writer_args':{'description':"ImageCameraName: timepix"}}
    ##sum images
    if 's' in opts:
        rock.convert2tiff(**tiff_args,nmax=0)
        rock.convert2tiff(n=npts,nmax=0)
        vw=rock.show_tiff(sum_opt=True,cutoff=20,i=i+1,pargs={'xylims':[0,516,516,0]})#,pargs={'opt':'p'},h=1)
    if 'p' in opts:
        pt.make_pets(pts_file='dat/pets_zenodo/glycine.pts',aperpixel=pets.aper,
            alphas=alpha0,ref_cell=[5.08760,11.80920,5.46150,90.00000,111.99200,90.00000],)
    else:
        b = rock.load(i)
        # b.set_beams_vs_thickness(thicks=(0,330,500))
        # b.convert2tiff(tiff_file='test.tiff',show=1,cutoff=50,**tiff_args)

        # vw=rock.show_tiff(cutoff=20,i=i+1,pargs={'xylims':[0,516,516,0]})#,pargs={'opt':'p'},h=1)
        # rock.convert2png(cutoff=20)


    if 'z' in opts:
        refls={
            0:[(4,1,0),(4,1,-1),(5,1,0),(5,1,1),(5,1,2),(-3,-1,2),(-3,-1,1),(3,1,-4),(3,1,-6),(3,1,-7)],
            4:[(-4,-1,0),(0,0,-2),(0,0,-4),(-3,-1,4),(-3,-1,6),(5,1,0)],
            10:[(0,0,-2),(0,0,2),(0,0,4),(0,0,6),(2,0,6),(-5,1,2),(-4,-1,-2)],
            11:[(-5,-1,8)],
        }

        refl = [str(h) for h in refls[i]]
        b = rock.load(i)
        df_pets=pets.rpl.loc[pets.rpl.eval('(F==%d) &(I>2)' %(frames[i]+1) )]
        df_jana=pets_simu.rpl.loc[pets_simu.rpl.eval('(F==%d)' %(frames[i]+1))]

        hkl_pets = df_pets.hkl
        df_pets.index=hkl_pets

        # b.show_beams_vs_thickness(refl=refl,thicks=(0,330,500),cm='jet')
        b._set_I(iZ=-1)

        b.df_G['hkl']=b.df_G.index
        hkl = b.df_G.loc[b.df_G.hkl.isin(hkl_pets),'hkl']
        hkl_jana = df_jana.loc[df_jana.hkl.isin(hkl_pets),'hkl']
        df_jana.index=df_jana.hkl

        b.df_G['Ipets']=np.nan
        b.df_G['Ijana']=np.nan
        b.df_G.loc[hkl,'Ipets']=df_pets.loc[df_pets.hkl.isin(hkl),'I']
        b.df_G.loc[hkl,'Ijana']=df_jana.loc[df_jana.hkl.isin(hkl),'I']

        print(colors.red+'frame %d' %frames[i]+colors.black)
        # print(b.df_G.loc[hkl,['I','Ipets','Ijana']])
        print(np.setdiff1d(hkl_pets,hkl))
        # print(df_pets.loc[np.setdiff1d(hkl_pets,hkl_jana),'I'])

        hkl=refl
        df=b.df_G.copy()
        df=df.drop(str((0,0,0)))
        EDdisp.show_frame(opts='SVkr',mag=500,rot=-205,df_bloch=df,hkl_idx=hkl,
            xylims=np.array([-1,1,1,-1])*2)
        # df=df.loc[hkl]
        # EDdisp.show_frame(opts='Bk',mag=1000,rot=-205,df_bloch=df,hkl_idx=hkl,
        #     xylims=np.array([-1,1,1,-1])*2)
