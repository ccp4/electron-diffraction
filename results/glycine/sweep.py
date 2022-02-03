from utils import*                  ;imp.reload(dsp)
from EDutils import utilities as ut ;imp.reload(ut)
from blochwave import bloch_pp as bl;imp.reload(bl)
from scipy.interpolate import interp1d

plt.close('all')

opts='t' #Solve(S) r(save rocking curve) t(tiff)
npts = 100
path='dat/bloch_new'

if 'S' in opts:
    bloch_args = {'cif_file':'dat/pets/alpha_glycine.cif','Smax':0.015,'Nmax':8,
        'solve':1,'thick':300,'keV':200}
    uvw = ut.get_uvw_rock(e0=[0,3,1],e1=[2,1],deg=2,npts=npts,show=0)
    rock = bl.Bloch_cont(path=path,tag='',uvw=uvw,Sargs=bloch_args)

rock_file=path+'/rock_.pkl'
rock = ut.load_pkl(file=rock_file)

if 't' in opts:
    rock.convert2tiff(thick=200,Imax=1e6,aperpixel=0.01,iz=-1)
    vw=rock.show_tiff(cutoff=20,frame=2)#,pargs={'opt':'p'},h=1)
    rock.convert2png(cutoff=20)

if 'r' in opts:
    # rock.set_beams_vs_thickness((10,200,3))
    # rock.plot_rocking(cond=lambda dfG:bl.strong_beams(dfG,tol=0.01,n=2),n=5)
    z,I = rock.get_rocking(cond=lambda dfG:bl.strong_beams(dfG,tol=0.01,n=2),n=5)
    for h,Ib in I.items() : I[h]=Ib[:,0]
    max([Ib.max() for Ib in I.values()])


    #frame info
    u = rock.df.u
    b = rock.load(10)
    Nmax=512//2
    px,py,I = b.df_G[['px','py','I']].values.T
    dq = 1.1*max(px.max(),py.max())/Nmax
    px,py = np.round(b.df_G.loc[hkl,['px','py']].values.T/dq)+Nmax
    dfpxy = pd.DataFrame(np.array([px,py]).T,columns=['px','py'],index=hkl)
    dfpxy.to_pickle(path+'pxy.pkl')


    #prepare interp dataframe
    hkl = list(I.keys())
    nbs,nframes = len(hkl),u.shape[0]
    Sw = {h:np.zeros((nframes)) for h in hkl}
    dfb = rock.beams.loc[list(I.keys())]
    for h,r in dfb.iterrows():
        Sw[h][r.Frame] = r.Sw
        Sw[h][:r.Frame[0]] = r.Sw[0]
        Sw[h][r.Frame[-1]:] = r.Sw[-1]
    # interp dataframe
    df = pd.DataFrame(columns=['Sw','I'],index=hkl)
    frame=np.arange(100)
    for h in hkl:
        df.loc[h,'Sw'] = interp1d(frame,Sw[h])
        df.loc[h,'I']= interp1d(frame,I[h])

    df.to_pickle(path+'/manim.pkl')

    df = pd.read_pickle(path+'/manim.pkl')
    # fb = [lambda t:(df.loc[h,'Sw'](t),df.loc[h,'I'](t)) for h in hkl]
    # for f in fb:print(f(0))
    fb=dict()
    for h in hkl:
        fb[h]=lambda t:(df.loc[h,'Sw'](t),df.loc[h,'I'](t))
        print(fb[h](0))
