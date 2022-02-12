import importlib as imp
import pandas as pd, numpy as np
import os,matplotlib,cbf,tifffile,re,glob,easygui
from subprocess import check_output
from matplotlib import rc
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
# from scipy.spatial import ConvexHull, convex_hull_plot_2d
from utils import glob_colors as colors,handler3D as h3d
from utils import physicsConstants as cst
from utils import displayStandards as dsp   #; imp.reload(dsp)
# from scattering.structure_factor import structure_factor3D
# from . import rotating_crystal as rcc       #; imp.reload(rcc)
# from . import postprocess as pp             #; imp.reload(pp)
from . import pets as pt                      #;imp.reload(pt)


class Base_Viewer:
    """
    Viewer class tiff/cbf files (similar to adxv)

    Parameters
    ----------
    figpath
        folder path containing the images
    frame
        frame number to start with
    i
        image index to start with
    cutoff
        brightness cut off
    pargs
        arguments to be passed to :meth:~Base_Viewer.show_im
    thick
        thickness info
    h
        show help


    .. NOTE::
        Works also for a single image
    """
    def __init__(self,figpath:str='',
        frame:Optional[int]=None,i:Optional[int]=0,
        cutoff:int=50,pargs:dict={},
        h:bool=True,
        thick:float=5,
    ):
        self.figpath = figpath
        #### format
        d_fmt = {'cbf':self.load_cbf,'tiff':self.load_tif,'tif':self.load_tif}
        self.supported_fmts = d_fmt.keys()
        self.fmt  = self.find_format(v=h)
        self.load = d_fmt[self.fmt]

        self.figs  = np.sort(glob.glob(self.figpath+'/*.%s' %self.fmt))#;print(self.figs)
        self.nfigs = self._get_nfigs()

        if frame:
            i = min(max(frame,0),self.nfigs)-1
        else:
            frame = i+1
        self.frame  = frame       #starting frame
        self.i      = i           #starting image
        self.inc    = 1           #increment(use 'p' or 'm' to change)
        self.cutoff = cutoff
        self.thick  = thick
        self.fieldNames  = ['thick','frame','cutoff','inc']
        self.mode   = 1
        self.pargs  = pargs
        rc('text', usetex=False)

        self.fig,self.ax = dsp.stddisp(opt='')
        cid = self.fig.canvas.mpl_connect('key_press_event', self)
        self.show()

        self.msg = '''
        'up or right'  : show frame+1
        'down or left' : show frame-1
        ##
        'p' : increase increment rate
        'm' : decrease increment rate
        ##
        'pageup'   : increae cutoff brightness
        'pagedown' : decrease cutoff brightness
        'r'        : reset cutoff brightness
        ##
        'ctrl+t' : increase thickness
        'ctrl+T' : decrease thickness
        ##
        'enter' : change settings
        'S' : save image
        'h' : show help
        '''
        if not 'msgD' in self.__dict__:
            self.msgD=''
        if h:self.show_help()

    def __call__(self, event):
        # print(event.key)
        if event.key in ['up','right']:
            self.i=min(self.i+self.inc,self.nfigs-1)
            self.mode=1
        elif event.key in ['left','down']:
            self.i=max(0,self.i-self.inc)
            self.mode=-1
        self.frame=self.i+1

        if event.key=='S':
            dsp.saveFig(self.get_figname(),ax=self.ax)

        #increment rate
        if event.key=='p':
            self.inc=min(self.inc+1,100)        ;print('increment rate : %d' %self.inc)
        elif event.key=='m':
            self.inc=max(1,self.inc-1)          ;print('increment rate : %d' %self.inc)
        elif event.key=='ctrl+r':
            self.inc=1                          ;print('increment rate : %d' %self.inc)

        #brightness
        dc=25
        if event.key=='pageup':
            self.cutoff=min(self.cutoff+dc,10000)  ;print('cutoff : %d' %self.cutoff)
        elif event.key=='pagedown':
            self.cutoff=max(1,self.cutoff-dc)    ;print('cutoff : %d' %self.cutoff)
        elif event.key=='r':
            self.cutoff=50                      ;print('cutoff : %d' %self.cutoff)

        #thickness
        if event.key=='ctrl+t':
            self.thick+=5                       ;print('thickness : %d' %self.thick)
        elif event.key=='ctrl+T':
            self.thick=max(self.thick-5,5)      ;print('thickness : %d' %self.thick)

        if event.key=='h':self.show_help()
        elif event.key=='enter':self.settings()
        keys = self.call(event)

        update_keys = keys+['enter','ctrl+t','ctrl+T','pageup','pagedown','r','left','right','down','up']
        if event.key in update_keys:self.update_im()

    def settings(self):
        fieldValues = ['%d' %self.__dict__[f] for f in self.fieldNames]
        fieldNames = self.fieldNames.copy()
        self.get_fieldValues(fieldNames,fieldValues)
        dict_fv = multenterbox("Change settings","settings", fieldValues,fieldNames)
        if dict_fv:
            for f in self.fieldNames:
                self.__dict__[f] = int(dict_fv[f])
            self.set_fieldValues(dict_fv)
        self.i = self.frame-1

    def show(self):
        im=self.load(self.figs[self.i])
        self.show_im(im,**self.pargs)
        self.im = self.ax.get_images()[0]

    def update_im(self):
        fig = self.figs[self.i]
        im = self.load(fig)
        # self.im.set_data(im)
        self.ax.cla()
        dsp.stddisp(im=[im],ax=self.ax,title='image %d' %self.i,
            cmap='gray',caxis=[0,self.cutoff],pOpt='t',
            name=self.get_figname(),**self.pargs)
        # print(self.cutoff)
        self.fig.canvas.draw()

        # figname=os.path.basename(fig)
        print(colors.yellow+fig+colors.black)

    def show_help(self):
        print(colors.green+'shortcuts : '+colors.blue+self.msg+colors.magenta+self.msgD+colors.black)
    def _get_nfigs(self):
        return self.figs.size
    def get_figname(self):
        return self.figs[self.i][:-3]+'.png'

    ###################################
    ##### virtual functions
    ###################################
    def show_im(self,im,**kwargs):
        """The function used to display the frames"""
        print(self.cutoff)
        dsp.stddisp(im=[im],ax=self.ax,title='image %d' %self.i,
            cmap='gray',caxis=[0,self.cutoff],pOpt='tX',xylims=[0,516,0,516],
            name=self.get_figname(),**kwargs)

    def get_fieldValues(self,fieldNames,fieldValues):return None
    def set_fieldValues(self,dict_fv):return None
    def call(self,event):return []

    ###############################################################
    #### misc
    ###############################################################
    def find_format(self,v=1):
        fmts = np.unique([f.split('.')[-1] for f in os.listdir(self.figpath)])
        fmts = [fmt for fmt in fmts if fmt in self.supported_fmts]
        if not len(fmts):
            raise Exception('no supported format found in %s. Supported formats :' %(self.figpath),self.supported_fmts)
        fmt = fmts[0]
        if len(fmts)>1:
            print('warning multiple formats found',fmts)
            print('using %s' %fmt)
        if v:print('%s format detected' %fmt)
        return fmt

    def load_cbf(self,fig):
        try:
            content = cbf.read(fig)
        except:#UnicodeDecodeError
            self.i=self.i+self.mode
            print(colors.red+'error reading file'+colors.black)
            self.import_exp()
            return
        numpy_array_with_data = content.data
        header_metadata = content.metadata
        # print(colors.blue,header_metadata,colors.black)
        return numpy_array_with_data

    def load_tif(self,fig):
        return tifffile.imread(fig)





class Pets_Viewer(Base_Viewer):
    """Pets enhanced viewer

    Parameters
    ----------
    pets
        full path to pets folder
    init
        initial config to display : c(center), r(reflections), p(spot prediction), b(boxes)
    """
    def __init__(self,pets:str=None,
        Smax=0.025,Nmax=13,rot=203,sim=False,
        init:str='',**kwargs,
    ):
        if isinstance(pets,str):pets = pt.Pets(pets)
        if isinstance(pets,pt.Pets):
            if sim:
                figpath = os.path.join(pets.path,'tiff/simulation')
            else:
                figpath = os.path.join(pets.path,'tiff')

        self.pets = pets
        if self.pets:
            vals  = ['center','refl','pred','boxes']
            keys  = ['c','r','p','b']
            alias = ['1','2','3','4']
            self.show_d = dict(zip(keys,vals))
            self.alias  = dict(zip(alias,keys))
            self.show_opt = dict(zip(vals,[c in init for c in keys ]))
            self.Smax,self.Nmax,self.rot=Smax,Nmax,rot
        self.msgD='''
        #### pets info
        '1 or c' : show beam center
        '2 or r' : show detected reflections (in rpl)
        '3 or b' : show boxes
        '4 or p' : show lattice spot locations from info
        '''
        super().__init__(figpath=figpath,**kwargs)

    def call(self,event):
        # k = event.key
        if self.pets:
            if event.key in self.alias.keys():event.key = self.alias[event.key]
            if event.key in self.show_d.keys():
                key = self.show_d[event.key]
                self.show_opt[key] = not self.show_opt[key]

            return list(self.show_d.keys())
        else:
            return []

    def show_im(self,im,**kwargs):
        tle = "frame %d, thickness=%d $\AA$" %(self.frame,self.thick)
        print('ok')
        plts,scat,pp = [],[],[]
        sargs = {'facecolor':'none','edgecolor':(0.7,0.7,0.15),'s':50,'marker':'o'}
        if self.pets:
            # df = self.rpl
            df = self.pets.rpl
            frame = self.i+1
            if self.show_opt['refl']:
                rpl = df.loc[df.F==frame]
                plts+=[[rpl.rpx-0.5,rpl.rpy-0.5,'r+','']]
            if self.show_opt['center']  :
                cen = self.pets.cen.iloc[self.i]
                plts += [[ cen.px-0.5,cen.py-0.5,'b+','']]
            if self.show_opt['pred'] or self.show_opt['boxes']:
                px,py,I,hkl = self.pets.get_kin(frame,
                    rot=self.rot,thick=self.thick,Smax=self.Smax,Nmax=self.Nmax,pixel=True)
            if self.show_opt['pred']:
                scat  = [px,py]
            if self.show_opt['boxes']:
                npx = 15
                pp = [dsp.Rectangle((px0-npx/2,py0-npx/2),npx,npx,facecolor='none',edgecolor='r') for px0,py0 in zip(px,py)]
        print('cutoff',self.cutoff)
        dsp.stddisp(plts,ax=self.ax,patches=pp,scat=scat,im=[im],ms=20,sargs=sargs,
            xylims=[0,516,516,0],
            cmap='gray',caxis=[0,self.cutoff],title=tle,pOpt='tX',**kwargs)






# class Image_viewer:
#     '''a copy of dials.image_viewer'''
#     def __init__(self,file,sym=1,pad=None):
#         basename      = file.split('_')
#         self.im       = int(basename[-1].replace('.cbf',''))
#         self.basename = '_'.join(basename[:-1])
#         file_ids      = [int(f.split('_')[-1].replace('.cbf','')) for f in glob.glob(self.basename+'*.cbf')]
#         if pad :
#             self.pad=pad
#         else:
#             self.pad=int(np.log10(max(file_ids)))+1       #;print(self.pad)
#
#         ls = np.array(check_output("head -n25 %s" %file,shell=True).decode().split('\n'))
#         px_str      = ls[['Pixel_size' in l for l in ls]][0]
#         D_str       = ls[['Detector_distance' in l for l in ls]][0] #;print(D_str)
#         lam_str     = ls[["Wavelength" in l for l in ls]][0]
#
#         self.pxy    = np.array(re.findall("\d+e-\d+", px_str),dtype=float)
#         self.D      = 1 #float(re.findall("\d+.\d+",D_str)[0])
#         self.lam    = float(re.findall("\d+.\d+",lam_str)[0])
#
#         shape = np.array(cbf.read(file).data.shape)
#         self.image = cbf.read(file).data
#
#         if sym:
#             Nx,Ny = np.array(shape/2,dtype=int)
#             self.Nx,self.Ny = Nx,Ny
#             self.pX,self.pY = np.meshgrid(np.arange(-Nx,Nx+1),np.arange(-Ny,Ny+1))
#         else:
#             Nx,Ny = shape
#             # self.pX,self.pY = np.meshgrid(np.arange(Nx),np.arange(Ny))
#             self.pX,self.pY = np.meshgrid(np.arange(Nx),np.flipud(np.arange(Ny)))
#         self.qx,self.qy = self.px2s(self.pX,self.pY)
#
#     def s2px(self,xh):
#         '''xh : recorded coordinates on image'''
#         px,py = self.pxy
#         pX = self.lam*self.D/px*xh[:,0] + self.Nx
#         pY = -self.lam*self.D/py*xh[:,1] + self.Ny
#         return pX,pY
#
#     def px2s(self,pX,pY):
#         px,py = self.pxy
#         qx,qy = px/self.D/self.lam*pX, py/self.D/self.lam*pY
#         return qx,qy
#
#     def show_image(self,im=None,lab='q',stack=1,**kwargs):
#         rc('text', usetex=False)
#         if not im:im = self.im
#         file = "%s_%s.cbf" %(self.basename,str(im).zfill(self.pad))   #;print(filename)
#         # with open(file,'wb') as f:print(file)
#         content = cbf.read(file)
#         image = content.data
#         if stack>1:
#             for i in range(1,stack):
#                 filename = "%s_%s.cbf" %(self.basename,str(im+i).zfill(self.pad))   #;print(filename)
#                 image+=cbf.read(filename).data
#             if 'caxis' in list(kwargs.keys()):
#                 kwargs['caxis'] = list(np.array(kwargs['caxis'])*stack) #; print(kwargs['caxis'])
#         if 'q' in lab:
#             labs = ['$q_x(A^{-1})$','$q_y(A^{-1})$']
#             X,Y = self.qx,-self.qy
#         elif 'p' in lab:
#             labs = ['','']#['$p_x$','$p_y$']
#             X,Y = self.pX,self.pY
#         elif 'x' in lab:
#             labs = ['$x(mm)$','$y(mm)$']
#             px,py = self.pxy
#             X,Y = self.pX*px*1e3,self.pY*py*1e3
#         return dsp.stddisp(im=[X,Y,image],labs=labs,
#             pOpt='ptX',title="%s" %os.path.basename(file),**kwargs)





# class Frames_Viewer(Base_Viewer):
#     '''Viewer for pets'''
#     def __init__(self,pets,thick,Imag,kargs,**sargs):
#         self.pets   = pets
#         self.kargs  = kargs
#         super().__init__(thick=thick,cutoff=Imag,**sargs)
#
#     def get_fieldValues(self,fieldNames,fieldValues):
#         fieldNames+=list(self.kargs.keys())
#         fieldValues+=[str(f) for f in self.kargs.values()]
#     def set_fieldValues(self,dict_fv):
#         float_keys = np.setdiff1d(list(self.kargs.keys()),'opts')
#         for f in float_keys:
#             self.kargs[f] = eval(dict_fv[f])
#         self.kargs['opts']=dict_fv['opts']
#     def _get_nfigs(self):
#         return self.pets.nFrames
#     def get_im(self,**kwargs):
#         self.pets.show_frame(frame=self.i+1,thick=self.thick,Imag=self.cutoff,
#             **self.kargs,**kwargs)
#
#     def call(self,event):
#         chars = 'KkPhr'
#         keys = ['ctrl+K','ctrl+k', 'ctrl+H','ctrl+h','ctrl+g']
#
#         vals = [c in self.kargs['opts'] for c in chars]
#         for i,(c,k) in enumerate(zip(chars,keys)):
#             if event.key==k:
#                 vals[i] = not vals[i]
#
#         self.kargs['opts']='q'+''.join([c for c,val in zip(chars,vals) if val])
#
#         return keys




def multenterbox(msg,title,fieldValues,fieldNames):
    fieldValues = easygui.multenterbox(msg, title, fieldNames,fieldValues)
    while True:
        if fieldValues is None:
            break
        errs = list()
        for n, v in zip(fieldNames, fieldValues):
            if v.strip() == "":errs.append('"{}" is a required field.'.format(n))
        if not len(errs):
            break
        fieldValues = easygui.multenterbox("\n".join(errs), title, fieldNames, fieldValues)
    return dict(zip(fieldNames,fieldValues))
