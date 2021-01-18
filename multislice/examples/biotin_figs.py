import importlib as imp
import pickle
from scattering.structure_factor import structure_factor3D
import multislice.multislice as mupy; imp.reload(mupy)
import multislice.postprocess as pp ; imp.reload(pp)
import multislice.mupy_utils as mut ; imp.reload(mut)
from crystals import Crystal
from utils import*
plt.close('all')
path = '../dat/biotin/'
figpath = '../docs_fig/'
file = path+'biotin.cif'
opts = 'T '  #E(Ewald) T(tilts) B(Beams)F(structFact)


crys = Crystal.from_cif(file)
lat_params = crys.lattice_parameters[:3]
a1,a2,a3 = lat_params
b12 = np.sqrt(1/a1**2+1/a2**2)

fs = {'leg':30,'lab':40,'tick':25}
if 'E' in opts:
    #Ewald configuration [010]
    mut.ewald_sphere(lat_params,lam=cst.eV2mum(12000)*1e4,tmax=60,T=20.,nx=5 ,ny=20,opt='p',name=figpath+'E_MX.png',xylims=[-1,1.01,-0.1,2.01],xyTicks=1,legLoc='upper right',fonts=fs)
    mut.ewald_sphere(lat_params,lam=cst.keV2lam(200)     ,tmax=7 ,T=0.2,nx=20,ny=10,opt='p',name=figpath+'E_ED.png',xylims=[-3,3.01,-0.1,1.01],xyTicks=1,legLoc='upper right',fonts=fs)

if 'B' in opts:
    df = pd.read_pickle(path+'df0.pkl')
    p1 = pp.load_multi_obj(path+'biotin_m0_autoslic.pkl')
    p1.datpath = p1.datpath.replace('ronan','tarik')
    p1.get_beams(bOpt='na')

    #plot
    iBs = ['(64,0)','(0,32)','(0,64)','(32,32)','(64,64)']
    fig,ax = p1.beam_vs_thickness(iBs=iBs,cm='jet',
        fonts=fs,opt='ps',name=figpath+'biotin_m0_beams.png')
    if 'F' in opts:
        # structure factor
        lat_vec = crys.reciprocal_vectors
        pattern = np.array([ np.hstack([a.coords_cartesian,a.atomic_number]) for a in crys.atoms])
        hkl,F3d = structure_factor3D(pattern,lat_vec,hkl=None,hklMax=3,sym=0)
        #excitation error
        K = 1/cst.keV2lam(200)
        k = [2/a1,2/a2, b12, 2*b12]
        ti = np.arcsin(k/K)
        xi = K*(1-np.cos(ti))
        z = np.linspace(0,2000,2000)
        Fibs = F3d[[2,0,0,1,2],[0,1,2,1,2],0]
        Fi = [Fibs[i]*np.sin(np.pi*z*zeta_i)/(np.pi*zeta_i) for i,zeta_i in enumerate(xi)]
        Ii = np.abs(Fi)**2
        #plot kinematic
        cs = dsp.getCs('jet',len(iBs))
        plts = [[z,I/1e8,[cs[i],'--'],''] for i,I in enumerate(Ii)]
        dsp.pltPlots(ax,plts,2)
        fig.show()



if 'T' in opts:
    # df = pp.update_df_info(path+'tilts/dfm.pkl',hostpath='/data3/lii26466/multislice/biotin/tilts/')#; print(df)
    # for i,f in enumerate(df.index):
    #     p1 = pp.load_multi_obj(path+'tilts/'+f)
    #     e=p1.ssh_get('badb',file='beamstxt',hostpath='/data3/lii26466/multislice/biotin/tilts/')#;print(e)
        # p1.get_beams(bOpt='fa')

    df = pd.read_pickle(path+'tilts/df.pkl')
    tilts = np.rad2deg(np.array([t[1] for t in df.tilt.values])/    1000)
    iBs = ['(4,0)','(0,2)','(0,4)','(2,2)','(4,4)']
    It = np.zeros((tilts.size,len(iBs)))
    for i,f in enumerate(df.index):
        p1 = pp.load_multi_obj(path+'tilts/'+f)
        p1.datpath = p1.datpath.replace('ronan','tarik')
        # p1.get_beams(bOpt='fa')
        It[i,:] = p1.beam_vs_thickness(iBs=iBs,bOpt='o')[-1][:,-1]
    cs = dsp.getCs('jet',len(iBs))
    plts = [[tilts,It[:,i],[cs[i],'-o'],iB] for i,iB in enumerate(iBs)]
    dsp.stddisp(plts,lw=2,labs=[r'$\theta(deg)$','$I$'],fonts=fs,name=figpath+'biotin_tilts_rocking.png',opt='ps')
    # p1.get_beams()
    # p1.beam_vs_thickness() #iBs=['(2,17)','(8,34)'],bOpt='f')
    # p1.save_pattern()
    # p1.pattern(Iopt='Isnlg',tol=1e-4,rings=[0,0.25,0.5,1],caxis=[-6.2,0],
    #     pOpt='ptX',xylims=[-1,1,-1,1],cmap='binary',imOpt='ch',axPos=[0.2,0.13,0.75,0.75])
    # h,k,I=p1.pattern(Iopt='Isnlg',tol=1e-4,out=1)
    # np.save('biotin_m1.npy',I)
