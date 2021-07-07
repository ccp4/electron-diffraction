from utils import*                  ;imp.reload(dsp)
from blochwave import bloch         ;imp.reload(bloch)
from EDutils import utilities as ut ;imp.reload(ut)
# from EDutils import viewers as vw   ;imp.reload(vw)
# plt.close('all')
path,figpath = 'dat','dat/figures/'
N    = '4_2'   #N-beam configuration
opts = 'sR' #s(solve), t(set_thicks), I(Iz), R(rocking curve),
opt  = 'ps'

dt = 0.05        #tilt degree around angle
nframes = 16   #nb simus
thick,nts = 1100,200

# exact 2 beam orientation
u2   = np.array([0.13846063, 0.01406432, 0.99026807])
# exact 3 beam orientation
u3   = np.array([0.23274900552973615, 0.022797452273781056, 0.9722696007768337])
# 4 beam configurations
u4_2 = np.array([0.0702957163976856, -0.00026010222764890744, 0.9975261623651619])
#tilt
omega = np.linspace(-dt,dt,nframes)
bloch_args = {'cif_file':'diamond','keV':200,
    'Nmax':8,'Smax':0.01,'thick':thick,'thicks':(0,thick,nts),
    'solve':1,'opts':'s'}


if 's' in opts:
    if   N=='2'  :u = u2
    elif N=='3'  :u = u3
    elif N=='4_2':u = u4_2
    uvw = ut.get_uvw_rock(u0=u,u1=[1,3],omega=omega)
    rock = ut.Rocking(Simu=bloch.Bloch,param='u',vals=uvw,ts=omega,
        tag='%sbeams' %N,path=path,**bloch_args)

rock = ut.load_pkl('dat/rock_%sbeams.pkl' %N)
if 't' in opts:rock.do('set_beams_vs_thickness',thicks=(0,thick,nts))

if 'I' in opts:
    b=rock.load(ts=0.015)
    b.get_Xig()
    b.show_beams_vs_thickness(thicks=(0,thick,nts),
        strong=['I'],m={'I':1e3},cm='Spectral',
        name=figpath+'%sbeam_Iz.svg' %N,opt=opt)

if 'R' in opts:
    if   N=='2':
        zs = np.arange(168,thick,168)
        bargs = {'refl':[str(tuple([0,4,0]))]}
    elif N=='3':
        zs = np.arange(50,2000,500)
        bargs = {'refl':[str(tuple([-6,-8,2]))]}
        # bargs = {'refl':[str(tuple([5,-3,-1]))]}
        # bargs = {'refl':[str(tuple(h)) for h in [[-6,-8,2],[5,-3,-1]] ]}
    elif N=='4_2':
        zs = np.arange(200,thick,250)
        # bargs = {'refl':[str(tuple(h)) for h in [[4,-8,0],[4,8,0]] ]}
        bargs = {'refl':[str(tuple(h)) for h in [[2,-6,0],[2,6,0]] ]}
    rock.plot_rocking(zs=zs,bargs=bargs,cmap='Spectral',
        lw=2,opt=opt,name=figpath+'%sbeams_diamond.svg' %N)
