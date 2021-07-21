import importlib as imp
import os,numpy as np
from . import bloch                   ;imp.reload(bloch)
from EDutils import display as EDdisp #;imp.reload(EDdisp)
from EDutils import rotate_exp as exp ;imp.reload(exp)
from EDutils import utilities as ut   #;imp.reload(ut)



def bloch_rock(tag,path='',opts='',
    uvw=None,u=None,u1=None,omega=np.linspace(-1,1,3),
    thicks=(0,1000,1000),bloch_args={},hkls=None,
    S_args={},W_args={},R_args={},Q_args={},I_args={},Z_args={},
    refl=[],cond='',ts0=None,i=0,iZs=None,zs=[],cm='hsv',
    opt='p'):
    ''' complete a rocking curve simulation sweep with Blochwave approach
    - tag : Sweep tag
    - uvw,u,u1,omega : beam orientation vectors : see ut.get_uvw_rock
    - bloch_args : thicknesses
    - opts : options to perform
        - generic : s(solve all), t(set_thickness),
        - ts0 dependent       : S(setup), I(Iz)
        - zs dependent        : R(rocking curve),Q(Idyn-Ikin plot)
        - refl dependent only : W(Excitation error),Z(integrated)
    - refl,cond : selected beams
    - iZs,zs : thicknesses to select (RQ)
    - i,ts0  : single simulation to select for opts(SWI)
    '''
    b_args = {'keV':200,'solve':1,'opts':'s','thick':thicks[1]}
    b_args.update(bloch_args)

    pkl = os.path.join(path,'rock_%s.pkl' %tag)
    figpath = os.path.join(path,'figures')
    if 's' in opts:
        if not isinstance(uvw,np.ndarray):
            uvw  = ut.get_uvw_rock(u0=u,u1=u1,omega=omega)
        rock = exp.Rocking(Simu=bloch.Bloch,param='u',vals=uvw,ts=omega,
            tag=tag,path=path,thicks=thicks,**b_args)
    else:
        rock = ut.load_pkl(pkl)

    if 't' in opts:
        rock.do('set_beams_vs_thickness',thicks=thicks)

    if any([c in opts for c in 'XSI']):
        i,ts0 = rock._get_ts(i,ts0)
        b = rock.load(i=i)#;print(b.u)

    for ch in opts:
        if   ch=='X' :print(b.get_Xig())
        elif ch=='S' :
            figname = '%s_Sw.svg' %tag
            idx = b.get_beam(cond=cond,refl=refl)#;print(idx)
            frame_args = {'opts':'SVkr','mag':500,'xylims':1.5,'rot':0}
            frame_args.update(S_args)
            EDdisp.show_frame(df_bloch=b.df_G,single_mode=False,
                hkl_idx=idx,name=os.path.join(figpath,figname),opt=opt,
                **frame_args)
        elif ch=='I':
            figname = '%s_Iz.svg' %tag
            refl += [[0,0,0]]
            b.show_beams_vs_thickness(thicks=thicks,refl=refl,cond=cond,
                title=r'$\theta=%.2f$' %ts0,cm='Spectral',
                name=os.path.join(figpath,figname),opt=opt,
                **I_args)
        elif ch=='W':
            figname = '%s_theta' %tag
            fz = lambda x:-np.log10(np.maximum(x,1e-10))
            sw_args = {'cm':cm,'fz':fz,'opts':'t'}#,'thick':thick}
            sw_args.update(W_args)
            rock.Sw_vs_theta(refl=refl,cond=cond,
                figname=os.path.join(figpath,figname),opt=opt,lw=2,
                **sw_args)
        elif ch=='R':
            rocking_args={'cmap':'Spectral'}
            if not isinstance(hkls,list):hkls = [refl]
            for iB,hkl in enumerate(hkls):
                figname = '%s_beams%d.svg' %(tag,iB)
                rocking_args.update(R_args)
                rock.plot_rocking(zs=zs,cond='',refl=hkl,
                    lw=2,opt=opt,name=os.path.join(figpath,figname),
                    **rocking_args)
        elif ch=='Q':
            figname = '%s_QQ.svg' %tag
            rock.QQplot(zs=zs,iZs=iZs,refl=refl,cond=cond,
                lw=2,opt=opt,name=os.path.join(figpath,figname),
                **Q_args)
        elif ch=='Z':
            figname = '%s_Iint.svg' %tag
            int_args = {'cm':cm}
            int_args.update(Z_args)
            rock.plot_integrated(cond=cond,refl=refl,
                lw=2,opt=opt,name=os.path.join(figpath,figname),
                **int_args)

    return rock




def load_pkl(file):
    with open(file,'rb') as f : obj = pickle5.load(f)
    return obj
