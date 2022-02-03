from utils import*                  #;imp.reload(dsp)
from blochwave import bloch;imp.reload(bloch)
from blochwave import bloch_pp as bl;imp.reload(bl)
from multislice import pets as pt   ;imp.reload(pt)
from EDutils import utilities as ut #;imp.reload(ut)
# plt.close('all')
path = 'dat/bloch/'

frames = np.arange(10,30)
opts='B' #P(Pets) B(Bloch)c(convergence check)
opt='ps'

pets = pt.Pets('dat/pets/glycine.pts')

if 'P' in opts:
    frames1 = np.setdiff1d(frames,np.arange(20,pets.nFrames+1,20))
    pets.show_Ihkl(name=path+'figures/glycine_Ihkl.svg',opt=opt)
    pets.show_Iavg(name=path+'figures/glycine_Iavg.svg',opt=opt)
    I = pets.show_Iframes(frames=frames1,refl=refl,cond='(I>500)',opts='AH',fz=fz,
        name=path+'figures/glycine_Iframe500.svg',lw=2,opt='ps')
    I = pets.show_Iframes(frames=frames1,refl=refl,
        Imax=500,cond='(I>100) & (I<500)',opts='AH',fz=fz,
        name=path+'figures/glycine_Iframe100.svg',lw=2,opt=opt)
    I = pets.show_Iframes(frames=frames1,refl=refl,
        cond='(I>25) & (I<100)',Imax=100,opts='AH',fz=fz,
        name=path+'figures/glycine_Iframe10.svg',lw=2,opt=opt)

    # I0 = pets.get_hklI(refl=I.columns.values)


bloch_args = {'cif_file':'dat/pets/alpha_glycine.cif','Smax':0.015,'Nmax':15,
    'opts':''}

npts = 4
alpha0,uvw0 = pets.alpha[frames-1],pets.uvw0[frames-1]
uvw = ut.uvw_add_points(uvw0,npts=npts,plot=0)
alpha = np.linspace(alpha0[0],alpha0[-1],uvw.shape[0])
# ut.show_uvw(uvw)

if 'c' in opts:
    b = bloch.Bloch(cif_file='dat/pets/alpha_glycine.cif',name='glycine14',u=uvw[14],keV=200,solve=0)
    b.convergence_test(Nmax=[7,15,30],Smax=[0.011,0.02],opts='Ss',
        z=np.arange(0,1001,1),hkl=[str(s) for s in [(0,0,2),(0,0,-2),(-3,-1,4)]])

alpha=alpha[10:20]
uvw=uvw[10:20]
if 'B' in opts:
    refl = [[-3,-1,2],[3,1,-4],[4,1,0],[-3,-1,4],[0,0,-2],[0,0,2]]
    # cond = '(Sw<1e-2) & (Vga>1e-6)'
    # cond='I>0.01'
    # cond = '(Sw<1e-2) & (Vga>1e-6) & (I>0.2)'
    fz = abs #np.log10
    cond=''
    # refl = [[0,0,-2],[0,0,2]] # [[-3,-1,2]]
    # refl = [[-4,0,0]]
    # refl = [[2,1,4]]

    hkls = [[hkl] for hkl in refl]
    # rock = bl.bloch_rock(tag='glycine_ref',uvw=uvw,bloch_args=bloch_args,
    #     omega=alpha,
    #     cond=cond,refl=refl,
    #     thicks=(10,1000,50),
    #     # W_args={'opts':'tf'},#'iTs':slice(2,14)},
    #     # R_args={'opts':'f'},
    #     # S_args={'opts':'tf'},
    #     # i=4,#ts0=None,
    #     # R_args={'cmap':'Spectral'},
    #     # zs = np.arange(100,500,100),cm='hsv',hkls=hkls,
    #     path=path,opts='',opt=opt)

    rock = ut.load_rock(path,tag='glycine_ref')
    f=rock.get_frames(str(tuple([0,0,2])),cols=['Sw'])
    # df  = pets.rpl
    # df0 = df.loc[(df.F<30) & (df.I>1)].sort_values('F')[['h','k','l','F']]
    # df0['refl'] = [str(h) for h in df0[['h','k','l']].values]
    # refl = pd.unique(df0.refl)
    # cond = ''#  '(Sw<1e-2) & I>1e-2'
    # Sw = rock.Sw_vs_theta(refl=refl,cond=cond,fz=fz,
    #     iTs=slice(2,30,1),opts='tfiI')

    # cond = 'I>0.01'
    # rock.Sw_vs_theta(refl=refl,cond='',fz=-np.log10,opts='tf');
    # rock.plot_rocking(zs=np.linspace(100,800,6),refl=[[-4,1,-1]],cond='',opts='t')
    # rock.plot_integrated(cond=cond,refl=refl,new=0,lw=2)
