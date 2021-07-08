import importlib as imp
import os,numpy as np
from . import bloch                 ; imp.reload(bloch)
from EDutils import utilities as ut #; imp.reload(ut)
from EDutils import display as EDdisp;imp.reload(EDdisp)



def bloch_rock(tag,uvw=None,u=None,u1=None,omega=np.linspace(-1,1,3),
    thicks=(0,1000,1000),bloch_args={},
    ts0=0,zs=[],hkls=[],cond='',fz=abs,
    path='',opts='',opt='p'):
    ''' complete a rocking curve simulation sweep with Blochwave approach
    - tag : Sweep tag
    - thicks,bloch_args : ()
    - u,u1,omega : see ut.get_uvw_rock
    - ts0,strong : beam vs thickness
    - zs,hkls : thicknesses and subset of reflections when plotting plot_rocking
    - opts : S(single) s(solve all), t(set_thickness), I(Iz), R(rocking curve),
    '''
    b_args = {'keV':200,'solve':1,'opts':'s'}
    b_args.update(bloch_args)

    pkl = os.path.join(path,'rock_%s.pkl' %tag)
    figpath = os.path.join(path,'figures')
    if 's' in opts:
        if not isinstance(uvw,np.ndarray):
            uvw  = ut.get_uvw_rock(u0=u,u1=u1,omega=omega)
        rock = ut.Rocking(Simu=bloch.Bloch,param='u',vals=uvw,ts=omega,
            tag=tag,path=path,thicks=thicks,**b_args)
    else:
        rock = ut.load_pkl(pkl)

    if 't' in opts:
        rock.do('set_beams_vs_thickness',thicks=thicks)
    # rock.do('_set_Vg')
    # cond = ''(Sw<1e-3) & (Vga>0.01)'


    b = rock.load(ts=ts0)
    refl=[]
    for hkl in hkls :refl+=hkl
    for ch in opts:
        if ch=='X' :print(b.get_Xig())
        if ch=='S' :
            figname = '%s_Sw.svg' %tag
            idx = b.get_beam(cond=cond,refl=refl)#;print(idx)
            EDdisp.show_frame(opts='SVk',df_bloch=b.df_G,single_mode=False,
                hkl_idx=idx,mag=500,xylims=2.5,
                name=os.path.join(figpath,figname),opt=opt)

        if ch=='W':
            figname = '%s_Sw_theta.svg' %tag
            rock.Sw_vs_theta(refl=refl,cond=cond, cm='hsv',fz=fz,lw=2,
                name=os.path.join(figpath,figname),opt=opt)

        if ch=='I':
            figname = '%s_Iz.svg' %tag
            refl += [[0,0,0]]
            refl = [tuple(h) for h in refl]
            b.show_beams_vs_thickness(thicks=thicks,refl=refl,cond=cond,cm='Spectral',
                name=os.path.join(figpath,figname),opt=opt)

        if ch=='R':
            for i,hkl in enumerate(hkls):
                bargs = {'refl':[tuple(h) for h in hkl]}
                figname = '%s_beams%d.svg' %(tag,i)
                rock.plot_rocking(zs=zs,bargs=bargs,cmap='Spectral',
                    lw=2,opt=opt,name=os.path.join(figpath,figname))

    return rock
