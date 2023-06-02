"""Bloch contiuous rotation experiment"""
import importlib as imp
import os,glob,numpy as np
from subprocess import Popen,check_output
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from utils import displayStandards as dsp#;imp.reload(dsp)
from utils import glob_colors as colors#;imp.reload(dsp)
from EDutils import display as EDdisp #;imp.reload(EDdisp)
from EDutils import rotate_exp        ;imp.reload(rotate_exp)
from EDutils import utilities as ut   #;imp.reload(ut)
from EDutils import viewers as vw     #;imp.reload(vw)
from . import bloch                   ;imp.reload(bloch)
from . import util as bu              ;imp.reload(bu)



class Bloch_cont(rotate_exp.Rocking):
    def __init__(self,**kwargs):
        """ Bloch continous rocking curve

        Parameters
        ----------
        kwargs
            to pass to :class:~EDutils.rotate_exp.Rocking
        """
        super().__init__(bloch.Bloch,**kwargs)
        self.figpath = os.path.join(self.path,'tiff')
        # self.thick=
        if not os.path.exists(self.figpath):
            Popen('mkdir %s' %self.figpath,shell=True)
        self.save()

    def set_beams_vs_thickness(self,thicks,v=1):
        self.do('_set_beams_vs_thickness',thicks=thicks,v=v)
        self.z = self.load(0).z

    def sum_images(self,n,fmt='',figpath=None,frames=()):
        ''' Sums images found in figpath by chuncks of n images
        and puts them into directory "figpath/sum"

        Parameters
        ----------
        n
            number of images to sum (each side of the central image will be n/2 images)
        figpath
            Where the images are located
        frames
            The range of frames consider in (f_init,f_end)
        '''
        #handle output style
        if not fmt:fmt = glob.glob(os.path.join(figpath,'*'))[0].split('.')[-1]
        sum_path=os.path.join(figpath,'sum')
        pad_n = int(np.ceil(np.log10(np.ceil(self.df.shape[0]/n))))
        if not os.path.exists(sum_path):
            Popen('mkdir %s' %sum_path,shell=True)

        # handles min and max frames
        nmax = self.df.shape[0]
        n_max = int(np.ceil(nmax/n))
        if not len(frames)==2:frames=(0,n_max+1)
        n_init,n_end = frames
        n_init,n_end=max(0,n_init),min(n_end,n_max)
        n2 = n//2 #int(np.ceil(n/2))

        ## check images exists
        nbounds=lambda n0:max(0,min(nmax-1,n0))
        ni,nf = nbounds(n*n_init-n2),nbounds(n*n_end+n2+1)
        filenames = np.array([figpath+'/%s.%s' %(self.load(i).name,fmt)
            for i in [ni,nf]])
        miss = [not os.path.exists(f) for f in filenames]
        if any(miss):
            print(colors.red+'Missing images : \n'+colors.black)
            print('\n'.join(filenames[miss]))
            return
        #get the size of the images
        nxy = bu.imread(filenames[0]).shape
        for i in np.arange(n_init,n_end+1):
            im = np.zeros(nxy)
            subframes = np.arange(max(0,i*n-n2),min(nmax,i*n+n2+1))
            print(colors.red,i,subframes,colors.black)
            for j in subframes:
                b = self.load(j)
                img_file=os.path.join(figpath,'%s.%s' %(b.name,fmt))
                im+=bu.imread(img_file)/n
            frame_str = str(i).zfill(pad_n)
            new_file = os.path.join(sum_path,'%s.%s' %(frame_str,fmt))
            out = check_output("cp %s %s" %(filenames[0],new_file),shell=True).decode()
            if out:print(out)
            bu.imwrite(new_file,im.T)#np.fliplr(np.flipud(im)))
            # print(colors.yellow+new_file+colors.green+' saved'+colors.black)

    def make_img(self,template=None,figpath=None,fmt='',frames=None,**kwargs):
        if not figpath:figpath=self.figpath
        nmax=self.df.shape[0]
        if template and not fmt:fmt=template.split('.')[-1]
        if type(frames)==type(None):frames=np.arange(nmax)
        if not os.path.exists(figpath):
            out=check_output('mkdir -p %s' %figpath,shell=True).decode()
            if out:print(out)
        for i in frames:
            b0 = self.load(i)
            filename=figpath+'/%s.%s' %(b0.name,fmt)
            b0.convert2img(filename,template,**kwargs)


    def convert2tiff(self,figpath=None,n=0,nmax=0,**kwargs):
        ''' Generate tiff files
        Parameters
        ------------
            figpath
                final name will be figpath+'/%s.tiff' %b.name
            kwargs
                passed to Bloch.convert2tiff
            n
                if n> 1 sum n images at a time
            nmax
                number of frames to consider
        '''
        import tifffile
        if not figpath:figpath=self.figpath
        if not nmax:nmax=self.df.shape[0]
        if n>1:
            b0=self.load(0)
            tiff_file=figpath+'/%s.tiff' %b0.name
            im=tifffile.imread(tiff_file)
            sum_path=os.path.join(figpath,'sum')
            # self.sum_path=sum_path
            pad_n = int(np.ceil(np.log10(np.ceil(self.df.shape[0]/n))))
            if not os.path.exists(sum_path):
                Popen('mkdir %s' %sum_path,shell=True)

            for i in range(0,nmax,n):
                im = np.zeros(im.shape)
                new_tiff_file = os.path.join(sum_path,'%s.tiff' %(str(i//n).zfill(pad_n)))
                # for j in range(n):
                for j in range(max(0,i+n-5),min(nmax-1,i+n+5)):
                    b = self.load(j)
                    tiff_file=figpath+'/%s.tiff' %b.name
                    im+=tifffile.imread(tiff_file)/n
                tifffile.imwrite(new_tiff_file,im)
                print(colors.yellow+new_tiff_file+colors.green+' saved'+colors.black)

        else:
            for i in range(nmax):
                b = self.load(i)
                b.convert2tiff(tiff_file=figpath+'/%s.tiff' %b.name,**kwargs)
                b.save()

    def convert2png(self,cutoff=20,n=None,**kwargs):
        import tifffile
        for i in range(self.df[:n].shape[0]):
            b = self.load(i)
            tiff_file=self.figpath+'/%s.tiff' %b.name
            im=tifffile.imread(tiff_file)
            dsp.stddisp(im=[im],
                cmap='gray',caxis=[0,cutoff],axPos='F',pOpt='p',
                opt='sc',name=tiff_file.replace('.tiff','.png'),figsize='im',**kwargs)

    def show_tiff(self,sum_opt=False,**kwargs):
        # figpath=self.figpath
        sum_path=os.path.join(self.figpath,'sum')
        if sum_opt:figpath=sum_path
        return vw.Base_Viewer(figpath,**kwargs)

    def plot_rocking(self,cond='',refl=[],**kwargs):
        if not any(refl) and not cond:cond=strong_beams
        return super().plot_rocking(cond=cond,refl=refl,**kwargs)

def strong_beams(dfG,
        tol:float=1e-2,n:int=10,
        opt:str='F'):
    """keep the n strongest beams, i.e. those with intensity I>tol not including central beam

    Parameters
    ----------
    dfG
        the bloch dataframe
    tol
        strong beam criteria
    n
        number of strong beams
    opt
        include Origin(O) Friedel pairs(F)
    Returns
    ----------
    x
        list of beams
    """
    dfM=dfG.copy()
    if not 'O' in opt:dfM = dfM.drop(str((0,0,0)))
    dfM = dfM.sort_values('I',ascending=False)[:2*n]
    if not 'F' in opt:remove_friedel_pairs(dfM)

    df = dfM.loc[dfM.I>tol]
    if df.shape[0]>n:df = df[:n]
    if df.shape[0]==0:df= dfM[:n]
    # print(df.shape,n)

    # print(dfM)
    return list(df.index.values)




# def bloch_rock(tag,path='',opts='',
#     uvw=None,u=None,u1=None,omega=np.linspace(-1,1,3),
#     thicks=(0,1000,1000),bloch_args={},hkls=None,
#     S_args={},W_args={},R_args={},Q_args={},I_args={},Z_args={},
#     refl=[],cond='',ts0=None,i=0,iZs=None,zs=[],cm='hsv',
#     opt='p'):
#     ''' complete a rocking curve simulation sweep with Blochwave approach
#     - tag : Sweep tag
#     - uvw,u,u1,omega : beam orientation vectors : see ut.get_uvw_rock
#     - bloch_args : thicknesses
#     - opts : options to perform
#         - generic : s(solve all), t(set_thickness),
#         - ts0 dependent       : S(setup), I(Iz)
#         - zs dependent        : R(rocking curve),Q(Idyn-Ikin plot)
#         - refl dependent only : W(Excitation error),Z(integrated)
#     - refl,cond : selected beams
#     - iZs,zs : thicknesses to select (RQ)
#     - i,ts0  : single simulation to select for opts(SWI)
#     '''
#     b_args = {}#'keV':200,'solve':1,'opts':'s','thick':thicks[1]}
#     b_args.update(bloch_args)
#
#     pkl = os.path.join(path,'rock_%s.pkl' %tag)
#     figpath = os.path.join(path,'figures')
#     if 's' in opts:
#         if not isinstance(uvw,np.ndarray):
#             uvw  = ut.get_uvw_rock(u0=u,u1=u1,omega=omega)
#         rock = exp.Rocking(Simu=bloch.Bloch,param='u',vals=uvw,ts=omega,
#             tag=tag,path=path,thicks=thicks,**b_args)
#     else:
#         rock = ut.load_pkl(pkl)
#
#     if 't' in opts:
#         rock.do('set_beams_vs_thickness',thicks=thicks)
#
#     if any([c in opts for c in 'XSI']):
#         i,ts0 = rock._get_ts(i,ts0)
#         b = rock.load(i=i)#;print(b.u)
#
#     for ch in opts:
#         if   ch=='X' :print(b.get_Xig())
#         elif ch=='S' :
#             figname = '%s_Sw.svg' %tag
#             idx = b.get_beam(cond=cond,refl=refl)#;print(idx)
#             frame_args = {'opts':'SVkr','mag':500,'xylims':1.5,'rot':0}
#             frame_args.update(S_args)
#             EDdisp.show_frame(df_bloch=b.df_G,single_mode=False,
#                 hkl_idx=idx,name=os.path.join(figpath,figname),opt=opt,
#                 **frame_args)
#         elif ch=='I':
#             figname = '%s_Iz.svg' %tag
#             refl += [[0,0,0]]
#             b.show_beams_vs_thickness(thicks=thicks,refl=refl,cond=cond,
#                 title=r'$\theta=%.2f$' %ts0,cm='Spectral',
#                 name=os.path.join(figpath,figname),opt=opt,
#                 **I_args)
#         elif ch=='W':
#             figname = '%s_theta' %tag
#             fz = lambda x:-np.log10(np.maximum(x,1e-10))
#             sw_args = {'cm':cm,'fz':fz,'opts':'t'}#,'thick':thick}
#             sw_args.update(W_args)
#             rock.Sw_vs_theta(refl=refl,cond=cond,
#                 figname=os.path.join(figpath,figname),opt=opt,lw=2,
#                 **sw_args)
#         elif ch=='R':
#             rocking_args={'cmap':'Spectral'}
#             if not isinstance(hkls,list):hkls = [refl]
#             for iB,hkl in enumerate(hkls):
#                 figname = '%s_beams%d.svg' %(tag,iB)
#                 rocking_args.update(R_args)
#                 rock.plot_rocking(zs=zs,cond='',refl=hkl,
#                     lw=2,opt=opt,name=os.path.join(figpath,figname),
#                     **rocking_args)
#         elif ch=='Q':
#             figname = '%s_QQ.svg' %tag
#             rock.QQplot(zs=zs,iZs=iZs,refl=refl,cond=cond,
#                 lw=2,opt=opt,name=os.path.join(figpath,figname),
#                 **Q_args)
#         elif ch=='Z':
#             figname = '%s_Iint.svg' %tag
#             int_args = {'cm':cm}
#             int_args.update(Z_args)
#             rock.plot_integrated(cond=cond,refl=refl,
#                 lw=2,opt=opt,name=os.path.join(figpath,figname),
#                 **int_args)
#
#     return rock
#
#
#
#
# def load_pkl(file):
#     with open(file,'rb') as f : obj = pickle5.load(f)
#     return obj
