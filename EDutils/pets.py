import importlib as imp
import os,glob, numpy as np, pandas as pd, tifffile,scipy.optimize as opt
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from subprocess import Popen,PIPE
from crystals import Crystal,Lattice
from utils import handler3D as h3D          #;imp.reload(h3D)
from utils import displayStandards as dsp   #;imp.reload(dsp)
from multislice  import mupy_utils as mut   #;imp.reload(mut)
from EDutils import viewers as vw           ;imp.reload(vw)
from EDutils import utilities as ut         #;imp.reload(ut)
from multislice.rotating_crystal import get_crystal_rotation

class Pets:
    def __init__(self,pts_file:str,
        cif_file:Optional[str]=None,gen:bool=False,dyn:bool=False):
        """ Pets importer

        parameters
        ----------
            pts_file
                full path to .pts file
            cif_file
                full path to .cif file (automatically found if None)
            gen
                force reload if True
            dyn
                load dyn_cif into .cif if True
        """
        self.name  = os.path.basename(pts_file).split('.pts')[0] #;print(self.name)
        self.path  = os.path.dirname(pts_file)
        self.out   = self.path+'/dat/'


        self.cif_file   = ut.find_cif_file(self.path,cif_file)
        self.crys       = ut.import_crys(self.cif_file)
        self.lat_vec    = np.array(self.crys.lattice_vectors)
        self.lat_vec1   = np.array(self.crys.reciprocal_vectors)/(2*np.pi)
        self.lat_params = self.crys.lattice_parameters

        gen |= not os.path.exists(self.out)
        # if os.path.exists(self.out):gen |= not os.path.exists(self.out+'rpl.txt')
        if gen:self._convert_pets()
        self._load_all(dyn)
        self.nFrames = self.uvw.shape[0]


    def _convert_pets(self):
        cmd = "%s/convert_pets.sh %s %s" %(os.path.dirname(__file__),self.path,self.name)
        p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE);p.wait();o,e=p.communicate()
        print(o.decode(),e.decode())

        A = np.loadtxt(self.out+'UB.txt')
        np.save(self.out+'UB.npy',np.reshape(A,(3,3)))

    def _add_hkl(self,df):
        hkl = self.invA.dot(df[['x','y','z']].values.T)
        df[['hx','kx','lx']]  = hkl.T
        df[['h','k','l']]     = np.array(np.round(hkl),dtype=int).T

    def _load_all(self,dyn=1):
        lam,omega,aper = np.loadtxt(self.out+'pts.txt')
        self.omega = omega
        self.aper  = aper
        self.lam   = lam
        self.K0    = 1/self.lam

        self.frames = pd.read_csv(self.out+'iml.txt',sep=',',names=['name','alpha','beta','domega','scale','calibration','ellipA','ellipP','used'])
        self.rpl = pd.read_csv(self.out+'rpl.txt',sep=',',names=['x','y','z','I','i','px','py','rpx','rpy','alpha','Im','F'])
        self.cor = pd.read_csv(self.out+'cor.txt',sep=',',names=['x','y','z','I','i','px','py','rpx','rpy','alpha','Im','F'])
        self.cen = pd.read_csv(self.out+'cenloc.txt',sep=',',names=['px','py','m','n'])
        self.xyz = pd.read_csv(self.out+'xyz.txt',sep=',',names=['x','y','z','I','u0','px','py','F','alpha','Im','u1'])
        self.hkl = pd.read_csv(self.out+'hkl.txt',sep=',',names=['h','k','l','I','i','F','L'])
        self.kin = pd.read_csv(self.out+'cif.txt',sep=',',names=['id','u','v','w','prec','alpha','beta','omega','scale'])
        self.dyn = pd.read_csv(self.out+'dyn.txt',sep=',',names=['id','u','v','w','prec','alpha','beta','omega','scale'])
        self.HKL = pd.read_csv(self.out+'HKL.txt',sep=',',names=['h','k','l','I','sig','F'])
        self.HKL_dyn = pd.read_csv(self.out+'HKL_dyn.txt',sep=',',names=['h','k','l','I','sig','F'])
        self.A   = np.load(self.out+'UB.npy')
        self.lat_params = np.loadtxt(self.out+'cell.txt')[:-1]
        self.lat = np.array(Lattice.from_parameters(*self.lat_params).lattice_vectors)
        self.invA = np.linalg.inv(self.A)

        self.cif = [self.kin,self.dyn][dyn]
        uvw   = self.cif[['u','v','w']].values
        beams = self.lat.T.dot(uvw.T).T
        self.uvw   = uvw/np.linalg.norm(uvw,axis=1)[:,None]
        self.uvw0  = beams/np.linalg.norm(beams,axis=1)[:,None]
        self.beams = self.K0*self.uvw0 #/np.linalg.norm(beams,axis=1)
        self.XYZ   = self.xyz[['x','y','z']].values.T

        self._add_hkl(self.rpl)
        self._add_hkl(self.cor)
        self._add_hkl(self.xyz)

        cx,cy = self.cen[['px','py']].iloc[self.rpl.F-1].values.T
        px,py = self.rpl[['px','py']].values.T
        qxqy  = self.aper*(np.vstack([px,-py]).T-np.vstack([cx,-cy]).T)
        self.rpl[['qx','qy']] = qxqy
        self.rpl['hkl'] = [str(tuple(h)) for h in self.rpl[['h','k','l']].values]

        self.alpha = self.cif.alpha.values
        self.hkl.index=[str(tuple(h)) for h in self.hkl[['h','k','l']].values]
        hkl = self.hkl[['h','k','l']].values
        rq = np.linalg.norm(hkl.dot(self.lat_vec1),axis=1)
        self.hkl['rq'] = rq

        # hkl = [str(tuple(h)) for h in self.xyz[['h','k','l']].values]
        # hkl0,idx,cc=np.unique(hkl,return_index=True,return_counts=True)
        # if not len(hkl)==idx.shape[0]:print('warning reflection not unique : ',hkl0[cc>1])
        # self.xyz.index=hkl

        if dyn:
            hkl,idx=np.unique([str(tuple(h)) for h in self.HKL_dyn[['h','k','l']].values],return_index=True)
            self.HKL_dyn=self.HKL_dyn.iloc[idx]
            self.HKL_dyn.index=hkl
            hkl = self.HKL_dyn[['h','k','l']].values
            self.HKL_dyn['rq'] = np.linalg.norm(hkl.dot(self.lat_vec1),axis=1)


    ###########################################################################
    #### compute :
    ###########################################################################
    def integrate_rpl(self,frames,cond='(I>10)',npx=10):
        """Perform manual integration"""
        if isinstance(frames,int):frames=np.array(frames)
        cond += ' & (F in %s) ' %str(list(frames))
        cond += ' & (px>%d) & (py>%d) & (px<%d) & (py<%d)' %tuple([npx]*2+[512-npx]*2)
        rpl = self.rpl.loc[self.rpl.eval(cond)]
        # print(rpl.shape)
        # pxy = rpl[['px','py']].values
        x0 = np.arange(-npx,npx+1)
        x,y = np.meshgrid(x0,x0)
        x1,y1 = np.meshgrid(np.arange(-npx,npx+1,0.1),np.arange(-npx,npx+1,0.1))

        tiffpath = os.path.join(self.path,'tiff')
        df=pd.DataFrame(columns=['Imax','px','py','sx','sy','noise','max_err','mean_err','min_err'])
        for i,frame in enumerate(frames):
            print(frame)
            im = tifffile.imread(os.path.join(tiffpath, '%s.tiff' %str(frame).zfill(5)))
            # dsp.stddisp(im=[im],pOpt='im',caxis=[0,50],cmap='gray')
            for j,r in rpl.loc[rpl.F==frame].iterrows():
                px,py = int(r.px),int(r.py)
                data  = im[py-npx:py+npx+1,px-npx:px+npx+1]

                p0 = (r.I,0,0,5,5,r.i)
                popt, pcov = opt.curve_fit(f_gauss2D, (x, y), data.ravel(), p0=p0)
                # print('I=%.1f,x0=%.1f,y0=%.1f,sx=%.1f,sy=%.1f,noise = %.1f' %tuple(popt))
                # dsp.stddisp(im=[data],pOpt='im',cmap='gray',caxis=[0,50])#;dsp.plt.show()
                # fig,ax = dsp.stddisp(im=[x,y,data],cmap='gray',
                #     contour=[x1,y1,gauss2D((x1,y1),*popt),8],cargs={'colors':'r'},pOpt='im',lw=2)
                # fit = gauss2D((x,y),*popt)
                # dsp.stddisp(im=[fit],pOpt='im',cmap='gray',caxis=[0,50])#;dsp.plt.show()
                # dsp.plt.show()
                err = abs(f_gauss2D((x,y),*popt)-data.ravel())
                hkl = str(tuple(r[['h','k','l']].values))+'_%d' %frame
                df.loc[hkl] = list(popt)+[err.max(),err.mean(),err.min()]
        df['Im'] = rpl.Im.values
        df['I']  = rpl.I.values
        df['i']  = rpl.i.values
        df['F']  = rpl.F.values
        df['Iint'] = df.Imax*df.sx*df.sy*np.pi*self.aper
        return df

    ###########################################################################
    #### get
    ###########################################################################
    def get_lattice(self,Nmax=5):
        return mut.get_lattice(self.lat_vec1,Nmax)

    def get_beam_dir(self,frame):
        return self.uvw[frame-1]

    def get_beam(self,frame):
        return self.beams[frame-1]

    def get_ewald(self,frame,nts=100,nps=200):
        return mut.get_ewald(self.get_beam(frame),nts,nps)

    def get_excitation_errors(self,frame,Nmax=5,Smax=0.02):
        return mut.get_excitation_errors(self.get_beam(frame),self.lat_vec1,Nmax=Nmax,Smax=Smax)

    def get_kin(self,frame,thick,Nmax=5,Smax=0.02,e0=[1,0,0],rot=0,Imag=10,pixel=False):
        """get kinematic intensities"""
        K = self.get_beam(frame)#;print(K)
        df = mut.get_kinematic_intensities(self.cif_file,K,thick,Nmax=Nmax,Smax=Smax)
        qxyz = df[['qx','qy','qz']].values
        I = df.I.values*Imag

        px,py = mut.project_beams(K,qxyz,e0)
        if rot:
            ct,st = np.cos(np.deg2rad(rot)),np.sin(np.deg2rad(rot))
            px,py = ct*px-st*py,st*px+ct*py
        hkl = df[['h','k','l']].values.T
        if pixel:
            cx,cy = self.cen.loc[frame-1,['px','py']].values
            aper  = self.aper
            px= px/aper + cx
            py=-py/aper + cy                                        #;print(px,py)

        return px,py,I,hkl

    def get_hklI(self,refl):
        return self.hkl.loc[refl,'I'].values

    # def get_beams(self,refl,cond,opts='A'):
    #     # for i,F in enumerate(self.nFrames):
    #     rpl = self.rpl.loc[self.rpl.eval(cond)]

    ###########################################################################
    #### Show Integration
    ###########################################################################
    def show_Iavg(self,fz=np.log10,**kwargs):
        hkl = self.HKL_dyn.loc[self.HKL_dyn.I>1].copy()

        hkl['q0'] = np.round(hkl['rq']*100)/100
        qs0 = np.unique(hkl['q0'])
        Iavg0 = np.zeros(qs0.shape)
        for i,q0 in enumerate(qs0):
            Iavg0[i] = hkl.loc[hkl.q0==q0,'I'].mean()
        plts = [[qs0,fz(Iavg0),'b','']]

        hkl['q1'] = np.round(hkl['rq']*10)/10
        qs1 = np.unique(hkl['q1'])
        Iavg1 = np.zeros(qs1.shape)
        for i,q1 in enumerate(qs1):
            Iavg1[i] = hkl.loc[hkl.q1==q1,'I'].mean()

        plts+= [[qs1,fz(Iavg1),'r','']]
        dsp.stddisp(plts,labs=['$q(A^{-1})$','$I_{avg}$'],lw=2,#name=name+'_Iavg.svg',
            **kwargs)

    def show_Ihkl(self,dyn:Iterable=None,**kwargs):
        '''Show integrated intensity and average integrated intensity as
         function of resolution

         parameters
         -----------
         dyn
            show dynamical refinement reflections
        '''

        hkl = self.HKL_dyn.copy()
        rq = hkl['rq']
        I  = hkl['I']
        hkl['q']  = np.round(rq*10)/10
        hkl['order'] = np.zeros(rq.shape)
        qs = np.unique(hkl['q'])
        for i,q in enumerate(qs):
            beams = hkl.q==q
            hklq  = hkl.loc[beams][['h','k','l']].values
            order = np.argsort(np.sum(hklq,axis=1))
            hkl.loc[beams,'order'] = order

        cond = 'I>1'
        rq,order,I = hkl.loc[hkl.eval(cond),['rq','order','I']].values.T
        if type(dyn)==type(None):
            cmap,c = 'jet',np.log10(I)
        else:
            hkl['color']=0
            # h_dyn = [h for h in dyn if h in hkl.index]
            hkl.loc[dyn,'color']=1
            c = hkl.loc[hkl.eval(cond),'color'].values
            cmap='RdYlGn'
        # print(rq.shape,I.shape,c.shape,hkl.shape)
        scat = [rq,order,10*I**(0.2),c]

        xylims = [0.1,2,-1,order.max()+1]
        fig,ax = dsp.stddisp(scat=scat,labs=['$q(A^{-1})$','beam order'],xyTicks=[qs,[]],xylims=xylims,
            cs='S',sargs={'cmap':cmap},imOpt='c',opt='')

        # res = np.round(100/qs)/100
        res = np.array([0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,5])
        q0,res = 1/res,np.array(res,dtype=str)
        dsp.addxAxis(ax,[],xLab=r'resolution $(A)$',c=(0.5,)*3,
            xTicks=q0,xTickLabs=res,xylims=xylims,axPos='V',#name=name+'_Ihkl.svg',
            **kwargs)

        fig.canvas.mpl_connect('button_press_event', self._on_click)

    def _on_click(self,event):
        from matplotlib.backend_bases import MouseButton
        if event.button is MouseButton.RIGHT:
            x,y = event.xdata,event.ydata
            idx = np.argmin(np.linalg.norm(self.hkl[['q','order']].values-[x,y],axis=1))
            hkl = self.hkl.iloc[idx]
            print(hkl[['rq','I']])

    def show_Iframes(self,frames,cm='hsv',refl=[],cond='',Imax=500,opts='AH',fz=abs,
        **kwargs):
        # frames = np.arange(1,self.nFrames+1)[iFs]
        refl = [str(tuple(h)) for h in refl]
        df = self.cor
        if 'A' in opts:
            df0 = df.loc[df.F<=frames[-1]]
            df0 = df0.loc[df0.eval(cond)]
            refl = list(pd.unique([str(tuple(h)) for h in df0[['h','k','l']].values]))
        if 'H' in opts:refl = [h for h in refl if h in self.hkl.index]
        # print(refl)

        nfs,nbs = np.array(frames).size,len(refl)
        I  = pd.DataFrame(np.nan*np.ones((nfs,nbs)),columns=[str(h) for h in refl])
        for i,F in enumerate(frames):
            # condF = '(F==%d) & ' %F + cond #;print(condF)
            dfF = df.loc[df.F==F]
            hkl = [str(tuple(h)) for h in dfF[['h','k','l']].values]#;print(hkl)
            # if not cond:
            idx = [i for i,refl0 in enumerate(hkl) if refl0 in refl]#;print(idx)
            hkl0 = np.array(hkl)[idx]                               #;print(hkl0)
            I.loc[i,hkl0] = dfF.iloc[idx]['I'].values

        I = I.loc[:,I.columns[I.max()<Imax]]
        refl = list(I.keys())
        Ib = fz(I.values.T)
        Imax = np.argmax(Ib,axis=1) #locate maxI

        cs,txts = dsp.getCs(cm,len(refl)),[]
        plts = [[frames,Ib[i],[cs[i],'-o'],'%s' %hkl] for i,hkl in enumerate(refl)]
        # if 't' in opts:txts = [[frames[idx],I[i,idx],'%s' %str(h),cs[i]] for i,(h,idx) in enumerate(zip(refl,Imax))]
        dsp.stddisp(plts,texts=txts,labs=['frame','$I$'],**kwargs)
        return I


    ###########################################################################
    #### display
    ###########################################################################
    def show_ewald_sphere(self,frame=None,Smax=0.01,Nmax=10,h3d=0,
        nts=100,nps=200,**kwargs):
        K = self.get_beam(frame)
        Kx,Ky,Kz = K

        (h,k,l),(qx,qy,qz) = self.get_lattice(Nmax)
        x,y,z = self.get_ewald(frame,nts,nps)
        df    = self.get_excitation_errors(frame,Smax=Smax,Nmax=Nmax)

        print('small excitation error beams',uvw)
        if frame:print('frame %d:' %frame)
        print(df[['h','k','l']])

        scat = ([qx,qy,qz,5,'b','o'],[0,0,0,'r','s'])
        scat+= ([df.qx,df.qy,df.qz,50,'g','s'],)
        surf = [x,y,z,{'alpha':0.4,'color':'r'}]
        plts = [[0,Kx],[0,Ky],[0,Kz],'r']
        # scat+=([x0,y0,z0])
        fig,ax = dsp.stddisp(plts,scat=scat,surfs=[surf],lw=3,
            labs=['x','y','z'],rc='3d',**kwargs)
        if h3d:h3d = h3D.handler_3d(fig,persp=False)

    def show_uvw(self):
        """show trajectory of orientation axis"""
        uvw = self.uvw
        ez =  np.cross(uvw[0],uvw[-1])
        tle = "orientation axis with frames"

        plts = [[[0,u[0]],[0,u[1]],[0,u[2]],'b'] for u in uvw ]
        plts +=[[ [0,ez[0]],[0,ez[1]],[0,ez[2]] , 'r' ]]
        dsp.stddisp(plots=plts,title=tle,labs=['x','y','z'],rc='3d',lw=4)

    def show_xyz(self,**kwargs):
        """show reflections in reciprocal space """
        x,y,z = self.XYZ
        tle = "reflections in reciprocal space"

        fig2,ax = dsp.stddisp(rc='3d',title=tle,scat=[x,y,z,2,'b'],labs=['x','y','z'],**kwargs)
        h3d2 = h3D.handler_3d(fig2,persp=False)

    def show_hkl(self,**kwargs):
        """Show intensities as hkl"""
        h,k,l,I = self.HKL[['h','k','l','I']].values.T
        fig2,ax = dsp.stddisp(rc='3d',scat=[h,k,l,I,'b'],labs=['h','k','l'],**kwargs)
        h3d2 = h3D.handler_3d(fig2,persp=False)

    def show_exp(self,frame=1,**kwargs):
        """show the experimental images"""
        return vw.Pets_Viewer(self,frame=frame,**kwargs)
    def show_sim(self,frame=1,**kwargs):
        """show the simulated images"""
        return vw.Pets_Viewer(self,frame=frame,sim=True,**kwargs)

    # def mp4_crystal_rotation(self,name,**kwargs):
    #     cif_file = self.cif_file
    #     uvw = [self.get_beam_dir(frame=f) for f in np.arange(self.nFrames)+1]
    #     for i,n in enumerate(uvw[:1]):
    #         frame = i+1
    #         name_x = 'figures/%s_x_%s.png' %(name,str(frame).zfill(3))
    #         name_y,name_z = name_x.replace('_x_','_y_'),name_x.replace('_x_','_z_')
    #         mut.show_cell(cif_file,name=name_x,title='%d' %frame,n=n,view=[0,0]  ,h3D=0,opt='sc',**kwargs)
    #         mut.show_cell(cif_file,name=name_y,title='%d' %frame,n=n,view=[90,0] ,h3D=0,opt='sc',**kwargs)
    #         mut.show_cell(cif_file,name=name_z,title='%d' %frame,n=n,view=[90,90],h3D=0,opt='sc',**kwargs)
    #
    #     cmd  = '/bin/bash -i -c "im2mp4 figures/%s_x_%%03d.png figures/x.mp4" && ' %name
    #     cmd += '/bin/bash -i -c "im2mp4 figures/%s_y_%%03d.png figures/y.mp4" && ' %name
    #     cmd += '/bin/bash -i -c "im2mp4 figures/%s_z_%%03d.png figures/z.mp4" '    %name
    #     p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE);p.wait()
    #     o,e = p.communicate();
    #     print(o.decode())

    # def show_sim(self,frame=1,**kwargs):
    #     mut.Viewer(self.path+'/multislice',frame=frame,**kwargs)
    #
    # def show_frames(self,thick=1000,Nmax=5,Smax=0.02,e0=[1,0,0],rot=0,Imag=10,Itol=20,opts='PKqr',**sargs):
    #     kargs = dict(zip(
    #         ['Nmax','Smax','e0','rot','Itol','opts'],
    #         [Nmax,Smax,e0,rot,Itol,opts]))
    #     return mut.Frames_Viewer(self,thick,Imag,kargs,**sargs)
    #
    # def show_frame(self,frame=0,opts='Pqr',
    #     thick=1000,Nmax=5,Smax=0.02,e0=[1,0,0],rot=0,Imag=10,Itol=20,
    #     show_hkl=True,qopt=True,rings=True,sim_pattern=None,
    #     **kwargs):
    #     ''' Show a frame with information specified by opts
    #     - opts : E(exp), P(proc), S(sim), K(kin), h(hkl),q(rec A),r(rings), k(hkl_k)
    #     '''
    #     if isinstance(opts,str):exp,proc,sim,kin,show_hkl,qopt,rings,show_hkl_kin = [c in opts for c in 'EPSKhqrk']
    #     if kin:qopt=1
    #
    #     plts,scat,txts,labs,qmax,wm = [[0,0,'b+']],(),[],['px','py'],0,5
    #     # plts=[]
    #     if qopt :
    #         wm*=self.aper
    #         labs = ['$q_x$','$q_y$']
    #     if proc:
    #         rpl0  = self.rpl.loc[self.rpl.F==frame]
    #         I,Im = rpl0[['I','Im']].values.T
    #         if qopt :
    #             px,py = rpl0[['qx','qy']].values.T
    #             # px,py = rpl0[['rpx','rpy']].values.T
    #             # cx,cy = self.cen.iloc[frame-1][['px','py']].values.T
    #             # px,py = ((np.vstack([px,-py]).T-[cx,-cy])*self.aper).T
    #         else:
    #             px,py = rpl0[['rpx','rpy']].values.T
    #             cx,cy = self.cen.iloc[frame-1][['px','py']].values.T
    #             plts = [[cx,cy,'b+']]
    #         c = (0.5,)*3
    #         scat += ([px,py,I,c,'o'],)
    #         qmax = np.ceil(max(py.max(),px.max()))
    #
    #         if show_hkl:
    #             hkl1  = rpl0[['h','k','l']]
    #             h,k,l = hkl1.values.T
    #             hx,kx,lx = rpl0[['hx','kx','lx']].values.T
    #             sx = lambda x:['%d' %round(x),'%.1f' %x][abs(x-round(x))>0.06]
    #             txts += [[x+0*wm,y+wm,'(%s,%s,%s)' %(sx(h0),sx(k0),sx(l0)),(0.5,)*3] for x,y,h0,k0,l0 in zip(px,py,hx,kx,lx)]
    #             print(h,k,l)
    #
    #     if kin:
    #         px_k,py_k,I_k,hkl_k = self.get_kin(frame,thick,Nmax,Smax,e0,rot,Imag)
    #         scat += ([px_k,py_k,I_k,'g','o'],)
    #         if show_hkl_kin:
    #             h,k,l = hkl_k
    #             Itol = I_k.max()/Itol
    #             txts += [[x,y+wm,'(%s,%s,%s)' %(h0,k0,l0),'g'] for x,y,h0,k0,l0,I in zip(px_k,py_k,h,k,l,I_k) if I>Itol]
    #             print(h,k,l)
    #         qmax = np.ceil(max(qmax,px_k.max(),py_k.max()))
    #
    #     im,bgcol = None,'k'
    #     if sim:
    #         qx,qy,I=sim_pattern
    #         qmax = np.ceil(max(qmax,qx.max(),qy.max()))
    #         im,bgcol=[qx,-qy,I],None
    #         # print(bgcol)
    #
    #     if rings and qopt:
    #         t = np.linspace(0,np.pi*2,100)
    #         ct,st = np.cos(t),np.sin(t)
    #         plts+=[[i*ct,i*st,'m','',0.5] for i in np.arange(0.25,qmax,0.25)]
    #         plts+=[[i*ct,i*st,'m','',2] for i in np.arange(1,qmax)]
    #     # if rot:
    #     #     ct,st = np.cos(np.deg2rad(rot)),np.sin(np.deg2rad(rot))
    #     #     px,py = ct*px-st*py,st*px+ct*py
    #
    #
    #     if not 'fonts' in kwargs.keys():kwargs['fonts']={'text':15}
    #     if not 'xylims' in kwargs.keys():kwargs['xylims']=[-px.max(),px.max(),-py.max(),py.max()]
    #     dsp.stddisp(plts,ms=20,scat=scat,im=im,texts=txts,bgcol=bgcol,gridOn=0,
    #         labs=labs,sargs={'alpha':0.5},**kwargs)

    ###########################################################################
    #### private compares
    ###########################################################################
    def _rotate_pattern(self,frame,qxqy):
        xyz0  = np.vstack([qx,qy,np.zeros(qx.shape)])

        alpha,beta,omega = -self.cif.iloc[frame-1][['alpha','beta','omega']]

        R = np.identity(3)
        omega_r = np.deg2rad(omega)
        ct,st  = np.cos(omega_r),np.sin(omega_r)
        Rz = np.array([[ct,st,0],[-st,ct,0],[0,0,1]]) #get_crystal_rotation(u=[0,0,-1],alpha=self.omega)
        R = Rz.dot(R)
        alpha_r = np.deg2rad(alpha)
        ct,st   = np.cos(alpha_r),np.sin(alpha_r)
        Rx = np.array([[1,0,0],[0,ct,st],[0,-st,ct]])
        R  = Rx.dot(R)
        beta_r = np.deg2rad(beta)
        ct,st   = np.cos(beta_r),np.sin(beta_r)
        Ry = np.array([[ct,0,st],[0,1,0],[-st,0,ct]])
        R  = Ry.dot(R)

        x0,y0,z0  = R.dot(xyz0)

    def _qxyz_to_pxy(self,frame,qxyz):
        alpha,beta,gamma = self.cif.iloc[frame-1][['alpha','beta','omega']]
        alpha_r,beta_r,gamma_r = -np.deg2rad([alpha,beta,gamma])
        ctx,stx = np.cos(alpha_r),np.sin(alpha_r)
        cty,sty = np.cos(beta_r) ,np.sin(beta_r)
        ctz,stz = np.cos(gamma_r),np.sin(gamma_r)

        Rx = np.array([[1,0,0],[0,ctx,stx],[0,-stx,ctx]])
        Ry = np.array([[cty,0,sty],[0,1,0],[-sty,0,cty]])
        Rz = np.array([[ctz,stz,0],[-stz,ctz,0],[0,0,1]])
        R = Ry.dot(Rx.dot(Rz))
        pxy  = R.dot(qxyz)[:2,:]
        return pxy

    def _compare_xyz_pxpy(self,frame=32,opts='oa',view=[90,90],**kwargs):
        rpl0    = self.rpl.loc[self.rpl.F==frame]
        pxc,pyc = self.cen.iloc[frame-1][['px','py']].values.T
        omega,beta,alpha = self.omega,0,rpl0.alpha.iloc[0]


        px,py = rpl0[['rpx','rpy']].values.T
        px    =  self.aper*(px-pxc)
        py    = -self.aper*(py-pyc)
        xyz0  = np.vstack([px,py,np.zeros(px.shape)])

        R = np.identity(3)
        if 'o' in opts:
            omega_r = np.deg2rad(omega)
            ct,st  = np.cos(omega_r),np.sin(omega_r)
            Rz = np.array([[ct,st,0],[-st,ct,0],[0,0,1]]) #get_crystal_rotation(u=[0,0,-1],alpha=self.omega)
            R = Rz.dot(R)
        if 'a' in opts:
            alpha_r = np.deg2rad(alpha)
            ct,st   = np.cos(alpha_r),np.sin(alpha_r)
            Rx = np.array([[1,0,0],[0,ct,st],[0,-st,ct]])
            R  = Rx.dot(R)
        if 'b' in opts:
            beta_r = np.deg2rad(beta)
            ct,st   = np.cos(beta_r),np.sin(beta_r)
            Ry = np.array([[ct,0,st],[0,1,0],[-st,0,ct]])
            R  = Ry.dot(R)

        x0,y0,z0  = R.dot(xyz0)

        xyz   = rpl0[['x','y','z']].values.T
        x,y,z = xyz

        Im = np.array(rpl0.Im.values/10,dtype=int)
        if isinstance(view,str):
            if   view=='x':labs,scat = ['y','z'],([y,z,Im,'b','o'],[y0,z0,Im,'r','o'])
            elif view=='y':labs,scat = ['x','z'],([x,z,Im,'b','o'],[x0,z0,Im,'r','o'])
            elif view=='z':labs,scat = ['x','y'],([x,y,Im,'b','o'],[x0,y0,Im,'r','o'])
            dsp.stddisp(xylims=1.5,scat=scat,labs=labs,legElt={'pxpy':'bo','xyz':'ro'},**kwargs)
        else:
            scat = ([x,y,z,Im,'b','o'],[x0,y0,z0,Im,'r','o'])
            fig,ax  = dsp.stddisp(rc='3d',view=view,xylims=1.5,scat=scat,labs=['x','y','z'],**kwargs)
            h3d = h3D.handler_3d(fig,persp=False)

    def _compare_hkl(self,frame,eps=None,Smax=0.025,Nmax=5,v=0):
        uvw  = self.get_beam_dir(frame=frame)
        df   = self.get_excitation_errors(frame,Nmax=Nmax,Smax=Smax)
        hkl0 = df[['h','k','l']]

        hkl1 = self.rpl.loc[self.rpl.F==frame] #[['hx','kx','lx']].values
        if eps:
            idx  = (abs(hkl1.hx-hkl1.h)<eps) & (abs(hkl1.kx-hkl1.k)<eps) & (abs(hkl1.lx-hkl1.l)<eps)
            hkl1 = hkl1.loc[idx]
        hkl1 = hkl1[['h','k','l']]

        rrmv = [r for r in hkl1.values if not np.where(np.linalg.norm(r-hkl0.values,axis=1)==0)[0].size]
        rrmv = np.array(rrmv,dtype=int)
        print('%d reflections in diffraction patterns ' %hkl1.shape[0])
        print('%d reflections within Smax=%.2E of the Ewald sphere ' %(hkl0.shape[0],Smax))
        print('%d reflections in diffraction pattern not in ewald construction ' %(rrmv.shape[0]))
        if rrmv.shape[0]:print('\tdetails : ',rrmv)

        ridx = []
        for r in hkl1.values:
            idx = np.where(np.linalg.norm(r-hkl0.values,axis=1)==0)[0]
            if idx.size:
                ridx+=[idx[0]]
        if v:
            print('Reflections from diffraction patterns :');print(hkl1)
            print('Reflections from excitation errors of Ewald construction :');
            print(df.iloc[ridx][['h','k','l','Sw']])
        return df.iloc[ridx][['h','k','l','Sw']]

    def _check_orientation_matrix(self):
        lab = np.identity(3)
        cry = self.A.dot(self.lattice_vectors)
        # cr2 = self.A.dot(self.reciprocal_vectors).T).T
        plts  = [ [[0,ai[0]],[0,ai[1]],[0,ai[2]],c,'$%s$' %l]  for i,(ai,c,l) in enumerate(zip(lab,['r','g','b'],['x','y','z'])) ]
        # plts += [ [[0,ai[0]],[0,ai[1]],[0,ai[2]],[c,'--'],'$%s$' %l]  for i,(ai,c,l) in enumerate(zip(cry,['r','g','b'],['a','b','c'])) ]
        x0,y0,z0 = [0.5]*3
        plts += [ [[x0,ai[0]+x0],[y0,ai[1]+y0],[z0,ai[2]+z0],[c,'--'],'$%s^{*}$' %l]  for i,(ai,c,l) in enumerate(zip(cr2,['r','g','b'],['a','b','c'])) ]
        # dsp.stddisp(plts,rc='3d',view=[0,0],name='figures/glycine_orient.png',opt='sc')





def gauss2D(X, amp, x0, y0, sx,sy,noise):
    x,y = X
    g = noise + amp*np.exp(-((x-x0)/sx)**2 - ((y-y0)/sy)**2)
    return g
f_gauss2D = lambda X,amp,x0,y0,sx,sy,noise:gauss2D(X,amp, x0, y0, sx,sy,noise).ravel()




def make_pets(pts_file:str,
    aperpixel:float,deg:float=0.0652,
    ref_cell:Sequence[float]=None,
    tag:Optional[str]=''):
    """creates a .pts file from info to process a simulated experiment with pets

    Parameters
    ----------
    pts_file
        full path to the pts file to write
    aperpixel
        aperture per pixel
    deg
        step angle between frames
    ref_cell
        cell info [ax,by,cz,alpha,beta,gamma] if known
    tag
        tag for the tiff images if there is one
    """
    center,pixel_size='AUTO',0.01,
    phi,omega= deg/2 ,0
    ax,by,cz,alpha,beta,gamma = ref_cell
    # ax,by,cz,alpha,beta,gamma = 8.1218, 9.977, 17.725, 90.0, 90.0, 90.0
    pts = '''lambda 0.025080
geometry continuous
omega  %d
phi  %.5f
virtualframes   7 5 1

aperpixel    %.6f
noiseparameters      2.5000      1.0000
saturationlimit   20000

center    256.0 256.0
centermode friedelpairs 0

beamstop no

dstarmax  3.0
dstarmaxps  3.0
i/sigma    7.00    5.00
reflectionsize  7

referencecell     %.5f    %.5f     %.5f    %.5f   %.5f    %.5f 1

#List of images
#columns: file name,alpha,beta,delta omega,frame scale,calibration correction(px/rec.angstrom),ellipt.distortion correction(amplitude),ellipt.distortion correction(phase), use for calculations(0/1)
imagelist
'''%(omega,phi,aperpixel,  ax,by,cz,alpha,beta,gamma)
    out = os.path.dirname(pts_file)
    tif_files = np.sort(glob.glob(out+'/tiff/%s*.tiff' %tag))
    alphas = np.arange(tif_files.size)*deg
    for i,tif_file in enumerate(tif_files):
        tif_file = os.path.basename(tif_file)
        pts+='%s %.4f 0.0 0.0 1.0 0 0 0  1\n' %('tiff\\'+tif_file,alphas[i])
    pts+='endimagelist\n'
    # pts_file = out+name+'.pts'
    with open(pts_file,'w') as f:
        f.write(pts)
        print(colors.green+'file saved : '+colors.yellow+pts_file+colors.black)
