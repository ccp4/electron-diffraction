import importlib as imp
from utils import*
import multislice.multislice as mupy        ;imp.reload(mupy)
import multislice.postprocess as pp         ;imp.reload(pp)
# import matplotlib.animation as animation
plt.close('all')

name ='../dat/Carbone_SC/'
fig_path = '../docs_fig/'

def run(Nz, *kwargs):
    return mupy.Multislice(name,tail='test',data='Carbone_SC001A.dat',
        mulslice=1,keV=200,TDS=False,
        NxNy=1024,slice_thick=1.0,repeat=[8,8,Nz],
        Nhk=10,hk_sym=1,**kwargs)
multi = run(Nz=100, opt='dsrp',fopt='f',ppopt='wB')
# multi.print_log()
# multi.beam_vs_thickness(cm='seismic')
# beams = multi.get_beams()
# multi.gif_beams()


def run_tilts(tilts,Nz,**kwargs):
    mupy.sweep_var(name,'tilt',np.deg2rad(tilts)*1000,do_prev=0,
        mulslice=1,keV=200,TDS=False,
        NxNy=1024,slice_thick=1.0,repeat=[8,8,Nz],
        Nhk=10,hk_sym=1,**kwargs)


def integrate(df_path,tilts,Nz,**kwargs):
    df = pd.read_pickle(df_path)
    iBs,iZs = np.arange(5),np.arange(100,Nz,100)
    nbs,nzs = iBs.size,iZs.size
    I,Iint = np.zeros((nts,nbs,nzs)), np.zeros((nbs,Nz))
    for it,f in enumerate(df.index) :
        multi = pp.load_multi_obj(name+f)
        hks,t,Ib = multi.get_beams()[[0,1,-1]]
        I[it,:,:] = Ib[iBs,:][:,iZs]
        Iint[iBs,:] += Ib[iBs,:]

    hks  =[hk[1:-1].split(',') for hk in hks]
    hk = ['(%d,%d)' %((int(hk[0])/8,int(hk[1])/8)) for hk in hks]
    cs = dsp.getCs('jet',nbs)
    for iz,iZ in enumerate(iZs):
        tle = 'z=$%.0f \AA$' %t[iZ]
        plts = [[tilts[:,0],I[:,iB,iz],[cs[iB],'-o'],'%s' %hk[iB]] for iB in range(nbs) ]
        dsp.stddisp(plts,labs=[r'$\theta(^{\circ})$','$I_b$'],title=tle,lw=2,
            opt='ps', name=fig_path+'SC_tilt.svg')

    plts = [[t,Iint[iB,:],[cs[iB],'-'],'%s' %hk[iB]] for iB in range(nbs) ]
    dsp.stddisp(plts,labs=[r'$z(\AA)$','$I_b$'],lw=2,**kwargs)

# nts=20
# Nz=1000
# tilts = np.vstack([np.linspace(-0.1,0.1,nts),0*np.ones((nts))]).T
# multi = run_tilts(tilts,Nz, opt='dsr',fopt='f',ppopt='',ssh='tarik-CCP4home',v=1)
# df = pp.update_df_info(name+'df.pkl',files=['beams'])
# integrate(df_path=name+'df.pkl',tilts=tilts,Nz=Nz,opt='p',name=fig_path+'SC_Iz_int.svg')
# multi = pp.load_multi_obj(name+'Carbone_SC_tilt10_mulslice.pkl')
# multi.beam_vs_thickness(cm='jet',iBs=range(5),opt='p',name=fig_path+'SC_Iz.svg')


#### load and run an individual simu
# multi = pp.load_multi_obj(name+'Carbone_SC_tilt016_mulslice.pkl')
# multi.run(ssh_alias='tarik-CCP4home')
# multi.ssh_get(ssh_alias='tarik-CCP4home','beams')
# beams = multi.get_beams()
# multi.beam_vs_thickness()
