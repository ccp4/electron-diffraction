import importlib as imp
import os,numpy as np
from . import bloch                  ;imp.reload(bloch)
from EDutils import utilities as ut  ;imp.reload(ut)
from EDutils import display as EDdisp#;imp.reload(EDdisp)



def bloch_rock(tag,uvw=None,u=None,u1=None,omega=np.linspace(-1,1,3),
    thicks=(0,1000,1000),bloch_args={},rot=0,
    iTs=[],ts=None,ts0=None,i=0,thick=None,zs=[],hkls=[],cond='',fz=abs,
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

    if any([c in opts for c in 'XSWIQZ']):
        b = rock.load(i=i,ts=ts0)#;print(b.u)
        refl=[]
        for hkl in hkls :refl+=hkl
    for ch in opts:
        if   ch=='X' :print(b.get_Xig())
        elif ch=='S' :
            figname = '%s_Sw.svg' %tag
            idx = b.get_beam(cond=cond,refl=refl)#;print(idx)
            EDdisp.show_frame(opts='SVkr',df_bloch=b.df_G,single_mode=False,
                hkl_idx=idx,mag=500,rot=rot,xylims=1.5,
                name=os.path.join(figpath,figname),opt=opt)
        elif ch=='W':
            figname = '%s_theta' %tag
            rock.Sw_vs_theta(refl=refl,cond=cond,thick=thick,fz=fz,iTs=iTs,opts='t',
                cm='hsv',lw=2,
                figname=os.path.join(figpath,figname),opt=opt)
        elif ch=='I':
            figname = '%s_Iz.svg' %tag
            refl += [[0,0,0]]
            refl = [tuple(h) for h in refl]
            b.show_beams_vs_thickness(thicks=thicks,refl=refl,cond=cond,
                title=r'$\theta=%.2f$' %ts0,cm='Spectral',
                name=os.path.join(figpath,figname),opt=opt)
        elif ch=='R':
            for i,hkl in enumerate(hkls):
                bargs = {'refl':[tuple(h) for h in hkl],'cond':''}
                figname = '%s_beams%d.svg' %(tag,i)
                rock.plot_rocking(zs=zs,bargs=bargs,cmap='Spectral',
                    lw=2,opt=opt,name=os.path.join(figpath,figname))

        elif ch=='Z':
            figname = '%s_Iint.svg' %tag
            rock.plot_integrated(cond=cond,refl=refl,
                lw=2,opt=opt,name=os.path.join(figpath,figname))
        elif ch=='Q':
            figname = '%s_QQ.svg' %tag
            rock.QQplot(zs=zs,refl=refl,
                lw=2,opt=opt,name=os.path.join(figpath,figname))

    return rock
