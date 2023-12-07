"""Bloch contiuous rotation experiment"""
import importlib as imp
import os,glob,tifffile,numpy as np,pandas as pd
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
        # self.figpath = os.path.join(self.path,'tiffs')
        # # self.thick=
        # if not os.path.exists(self.figpath):
        #     Popen('mkdir %s' %self.figpath,shell=True)
        self.save()

    def set_beams_vs_thickness(self,thicks,v=1):
        self.do('_set_beams_vs_thickness',thicks=thicks,verbose=v)
        self.z = self.load(0).z

    def set_thickness(self,verbose=True,**kwargs):
        self.do('set_thickness',verbose,**kwargs)


    #### Frames production related
    def make_frames(self,exp,
        path='',
        im_args=dict(gs3=0.1),
        tiff_writer_args={},
        png=False,
        ):
        '''
        Parameters
        ----------
        exp : Experiment object
        '''
        i=0
        sub_frames=self.sub_frames
        rock_frames=self.rock_frames
        pad=int(np.ceil(np.log10(self.n_simus/sub_frames)))+1
        while i<self.n_simus:
            exp_frame = rock_frames[0] + i//sub_frames
            print('exp frame %d' %exp_frame)
            simus = i+np.arange(sub_frames)
            beams = np.unique(np.hstack([self.load(j).df_G.index for j in simus]))
            df = pd.DataFrame(0,index=beams,
                columns=['px','py','I'])
            df.loc[beams,['px','py']] = exp.hkl_to_pixels(beams,exp_frame+1)
            for j in simus:
                b0=self.load(j)
                df.loc[b0.df_G.index,'I'] += b0.df_G.I
            # print(df.I*1000)
            im = _make_img(df,**im_args)
            # print(im.max(),im.mean())
            #save
            tiff_file=os.path.join(path,'%s.tiff' %str(exp_frame).zfill(pad))
            tifffile.imwrite(tiff_file,im.T,**tiff_writer_args)
            print(colors.yellow+tiff_file+colors.green+' saved'+colors.black)

            i+=sub_frames


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

    def convert2png(self,path,cutoff=20,n=None,**kwargs):
        import tifffile
        tiff_files=glob.glob(path+'/*.tiff')
        png_path=os.path.join(path,'png')
        if not os.path.exists(png_path):
            cmd="if [ -d {path} ];then rm -rf {path};fi;mkdir -p {path}".format(
                path=png_path)
            check_output(cmd,shell=True).decode()
        for tiff_file in tiff_files:
            png_file=os.path.join(png_path,
                os.path.basename(tiff_file).replace('.tiff','.png'))
            im=tifffile.imread(tiff_file)
            args=dict(axPos='F',pOpt='p',opt='sc',figsize='im',cmap='gray',caxis=[0,cutoff])
            args.update(kwargs)
            dsp.stddisp(im=[im],name=png_file,**args)

    ##### visualization
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


def _make_img(df=None,
        Nmax:int=512,
        Imax=1000,
        gs3=0.1,nX=1,fbroad=None,
        show_kernel=False,
    ):
    '''Make image

    Parameters
    -----------
    df
        intensities dataframe ['px','py','I']
    Imax
        max value for the intensities
    Nmax
        image resolution in pixel
    kernel
        fbroad
            broadening function f(r2) (default np.exp(-r2/(gs3/3)**2))
        gs3
            Gaussian broadening factor for each reflections
        nX
            width factor of the Gaussian window (in pixels)
    '''
    #### kernel (converted to pixel locations)
    # fbroad,gs3,nX = [kernel[k] for k in ['fbroad','gs3','nX']]
    if not fbroad:
        # fbroad=lambda r2:np.exp(-r2/(gs3/3)**2)
        fbroad=lambda r2:np.exp(-r2/gs3)
        # nx,ny = np.array(np.floor(gs3/np.array([dqx,dqy])),dtype=int)
    # else:
    nx,ny=nX,nX
    ix,iy = np.meshgrid(range(-nx,nx+1),range(-ny,ny+1))
    # x,y = ix*dqx,iy*dqy
    r2 = (ix**2+iy**2)
    Pb = fbroad(r2)
    if show_kernel:
        dsp.stddisp(im=[ix,iy,Pb],pOpt='im');dsp.plt.show()
        return

    # get intensity
    hkl=df.index
    px,py,I = df.loc[hkl,['px','py','I']].values.T
    Nmax=Nmax//2

    # pixel locations
    i,j = np.array(df[['px','py']].values,dtype=int).T


    ## Apply Gaussian broadening
    im0 = np.zeros((2*Nmax,)*2)
    for i0,j0,I0 in zip(i,j,I):
        i0x,j0y = i0+ix,j0+iy
        idx     = (i0x>=0) & (j0y>=0)  & (i0x<2*Nmax) & (j0y<2*Nmax)
        i0x,j0y = i0+ix[idx],j0+iy[idx]
        im0[i0x,j0y] += Pb[idx]/Pb[idx].sum()*I0

    return np.array(im0/im0.max()*Imax,dtype=int)
