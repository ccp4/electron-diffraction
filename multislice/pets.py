import importlib as imp
import os,glob, numpy as np, pandas as pd
from subprocess import Popen,PIPE
from crystals import Crystal
from utils import handler3D as h3D          ;imp.reload(h3D)
from utils import displayStandards as dsp   ;imp.reload(dsp)
from . import mupy_utils as mut             ;imp.reload(mut)
from .rotating_crystal import get_crystal_rotation

class Pets:
    def __init__(self,pts_file,cif_file=None,gen=False):
        '''
        - pts_file : str - full path to .pts file
        - cif_file : str - full path to .cif file (automatically found if None)
        - gen      : bool - regenerate if True  
        '''
        self.name  = os.path.basename(pts_file).split('.pts')[0] #;print(self.name)
        self.path  = os.path.dirname(pts_file)
        self.out   = self.path+'/dat/'

        self.cif_file   = self._get_cif_file(cif_file)
        self.crys       = Crystal.from_cif(self.cif_file)
        self.lat_vec    = np.array(self.crys.lattice_vectors)
        self.lat_vec1   = np.array(self.crys.reciprocal_vectors)/(2*np.pi)
        self.lat_params = self.crys.lattice_parameters

        gen |= not os.path.exists(self.out)
        # if os.path.exists(self.out):gen |= not os.path.exists(self.out+'rpl.txt')
        if gen:self.convert_pets()
        self.load_all()
        self.nFrames = self.uvw.shape[0]


    def convert_pets(self):

        cmd = "%s/convert_pets.sh %s %s" %(os.path.dirname(__file__),self.path,self.name)
        p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE);p.wait();o,e=p.communicate()
        print(o.decode()) #; print(e.decode())

        lam,omega,aper = np.loadtxt(self.out+'pts.txt')
        self.omega = omega
        self.aper  = aper
        self.lam   = lam
        self.K0    = 1/self.lam

        A = np.loadtxt(self.out+'UB.txt')
        np.save(self.out+'UB.npy',np.reshape(A,(3,3)))
        self.load_all()

    def load_all(self):
        self.frames = pd.read_csv(self.out+'iml.txt',sep=',',names=['name','alpha','beta','domega','scale','calibration','ellipA','ellipP','used'])
        self.rpl = pd.read_csv(self.out+'rpl.txt',sep=',',names=['x','y','z','I','i','px','py','rpx','rpy','alpha','Im','F'])
        self.cor = pd.read_csv(self.out+'cor.txt',sep=',',names=['x','y','z','I','i','px','py','rpx','rpy','alpha','Im','F'])
        self.cen = pd.read_csv(self.out+'cenloc.txt',sep=',',names=['px','py','m','n'])
        self.xyz = pd.read_csv(self.out+'xyz.txt',sep=',',names=['x','y','z','I','u0','px','py','F','alpha','Im','u1'])
        self.hkl = pd.read_csv(self.out+'hkl.txt',sep=',',names=['h','k','l','I','i','F','L'])
        self.cif = pd.read_csv(self.out+'cif.txt',sep=',',names=['id','u','v','w','prec','alpha','beta','omega','scale'])
        self.HKL = pd.read_csv(self.out+'HKL.txt',sep=',',names=['h','k','l','I','sig','F'])
        self.A   = np.load(self.out+'UB.npy')
        self.invA = np.linalg.inv(self.A)

        uvw = self.cif[['u','v','w']].values
        self.uvw = np.array([u/np.linalg.norm(u) for u in uvw])
        self.XYZ = self.xyz[['x','y','z']].values.T

        rpl_hkl = self.invA.dot(self.rpl[['x','y','z']].values.T)
        self.rpl[['hx','kx','lx']]  = rpl_hkl.T
        self.rpl[['h','k','l']]     = np.array(np.round(rpl_hkl),dtype=int).T

    ###########################################################################
    #### get
    ###########################################################################
    def get_lattice(self,Nmax=5):
        N = np.arange(-Nmax,Nmax+1)
        h,k,l = np.meshgrid(N,N,N)
        h,k,l = h.flatten(),k.flatten(),l.flatten()
        a1,b1,c1 = self.lat_vec1
        qx = h*a1[0]+k*b1[0]+l*c1[0]
        qy = h*a1[1]+k*b1[1]+l*c1[1]
        qz = h*a1[2]+k*b1[2]+l*c1[2]
        return (h,k,l),(qx,qy,qz)

    def get_beam_dir(self,frame=None):
        if frame:
            uvw = self.uvw[frame-1]
            #convert beam to xyz reciprocal basis
            lat_vec = self.lat_vec
            uvw = lat_vec.T.dot(uvw)
            uvw/=np.linalg.norm(uvw)
        return uvw

    def get_beam(self,uvw):
        K = self.K0*uvw/np.linalg.norm(uvw)
        return K

    def get_ewald(self,uvw):
        '''uvw : reciprocal space beam direction'''
        K = self.get_beam(uvw)
        Kx,Ky,Kz = K

        theta,phi = np.linspace(0,np.pi,50),np.linspace(0,2*np.pi,100)
        x = Kx+self.K0*np.outer(np.sin(theta),np.cos(phi))
        y = Ky+self.K0*np.outer(np.sin(theta),np.sin(phi))
        z = Kz+self.K0*np.outer(np.cos(theta),np.ones(phi.shape))
        return x,y,z

    def get_excitation_errors(self,uvw,Nmax=5,Smax=0.02):
        '''uvw : reciprocal space beam direction'''
        K = self.get_beam(uvw)
        Kx,Ky,Kz = K

        (h,k,l),(qx,qy,qz) = self.get_lattice(Nmax)

        Sw = np.abs(np.sqrt((Kx-qx)**2+(Ky-qy)**2+(Kz-qz)**2) - self.K0)
        idx = Sw<Smax
        h0,k0,l0 = np.array([h[idx],k[idx],l[idx]],dtype=int)
        data = np.array([h0,k0,l0,qx[idx],qy[idx],qz[idx],Sw[idx]]).T
        return pd.DataFrame(data,columns=['h','k','l','qx','qy','qz','Sw'])

    ###########################################################################
    #### display
    ###########################################################################
    def show_ewald_sphere(self,uvw=None,frame=None,Smax=0.01,Nmax=10,h3d=0,**kwargs):
        if frame:uvw = self.get_beam_dir(frame)
        K = self.get_beam(uvw)
        Kx,Ky,Kz = K

        (h,k,l),(qx,qy,qz) = self.get_lattice(Nmax)
        x,y,z = self.get_ewald(uvw)
        df    = self.get_excitation_errors(uvw,Smax=Smax,Nmax=Nmax)

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
        # rotation matrix B=MA
        # convert axis and angle, average
        uvw = self.uvw
        ez =  np.cross(uvw[0],uvw[-1])

        plts = [[[0,u[0]],[0,u[1]],[0,u[2]],'b'] for u in uvw ]
        plts +=[[ [0,ez[0]],[0,ez[1]],[0,ez[2]] , 'r' ]]
        dsp.stddisp(plots=plts,labs=['x','y','z'],rc='3d',lw=4)

    def show_xyz(self):
        x,y,z = self.XYZ
        fig2,ax = dsp.stddisp(rc='3d',scat=[x,y,z,2,'b'],labs=['x','y','z'])
        h3d2 = h3D.handler_3d(fig2,persp=False)

    def show_frame(self,frame=0,show_hkl=True,print_hkl=False,qopt=True,rings=True,
        Smax=None,Nmax=10,**kwargs):
        rpl0 = self.rpl.loc[self.rpl.F==frame]
        wm,txts = 5,[]

        px,py,I,Im = rpl0[['px','py','I','Im']].values.T

        hkl1  = rpl0[['h','k','l']]
        h,k,l = hkl1.values.T
        hx,kx,lx = rpl0[['hx','kx','lx']].values.T
        if print_hkl:
            print(h,k,l)

        cx,cy = self.cen.iloc[frame-1][['px','py']].values.T
        labs = ['px','py']
        if qopt :
            wm*=self.aper
            px,py = ((np.vstack([px,py]).T-[cx,cy])*self.aper).T
            cx,cy=0,0
            labs = ['$q_x$','$q_y$']

        plts = [[cx,cy,'b+']]
        if show_hkl:
            sx = lambda x:['%d' %round(x),'%.1f' %x][abs(x-round(x))>0.06]
            txts += [[x+0*wm,y+wm,'(%s,%s,%s)' %(sx(h0),sx(k0),sx(l0)),'g'] for x,y,h0,k0,l0 in zip(px,py,hx,kx,lx)]
        if rings and qopt:
            t = np.linspace(0,np.pi*2,100)
            ct,st = np.cos(t),np.sin(t)
            qmax = np.ceil(max(py.max(),px.max()))
            plts+=[[i*ct,i*st,'m','',0.5] for i in np.arange(0.25,qmax,0.25)]
            plts+=[[i*ct,i*st,'m','',2] for i in np.arange(1,qmax)]

        scat = ([px,py,I,'w','o'],)
        if Smax:
            uvw  = self.get_beam_dir(frame=frame)
            df   = self.get_excitation_errors(uvw,Nmax=Nmax,Smax=Smax)
            hkl0 = df[['h','k','l']]
            ridx = []
            for r in hkl0.values:
                idx = np.where(np.linalg.norm(r-hkl1.values,axis=1)==0)[0]
                if idx.size:
                    ridx+=[idx[0]]
            scat+= ([px[ridx],py[ridx],20,'g','s'],)

        if not 'fonts' in kwargs.keys():kwargs['fonts']={'text':20}
        dsp.stddisp(plts,ms=20,scat=scat,texts=txts,bgcol='k',gridOn=0,
            labs=labs,xylims=[-px.max(),px.max(),-py.max(),py.max()],**kwargs)

    def mp4_crystal_rotation(self,name,**kwargs):
        cif_file = self.cif_file
        uvw = [self.get_beam_dir(frame=f) for f in np.arange(self.nFrames)+1]
        for i,n in enumerate(uvw[:1]):
            frame = i+1
            name_x = 'figures/%s_x_%s.png' %(name,str(frame).zfill(3))
            name_y,name_z = name_x.replace('_x_','_y_'),name_x.replace('_x_','_z_')
            mut.show_cell(cif_file,name=name_x,title='%d' %frame,n=n,view=[0,0]  ,h3D=0,opt='sc',**kwargs)
            mut.show_cell(cif_file,name=name_y,title='%d' %frame,n=n,view=[90,0] ,h3D=0,opt='sc',**kwargs)
            mut.show_cell(cif_file,name=name_z,title='%d' %frame,n=n,view=[90,90],h3D=0,opt='sc',**kwargs)

        cmd  = '/bin/bash -i -c "im2mp4 figures/%s_x_%%03d.png figures/x.mp4" && ' %name
        cmd += '/bin/bash -i -c "im2mp4 figures/%s_y_%%03d.png figures/y.mp4" && ' %name
        cmd += '/bin/bash -i -c "im2mp4 figures/%s_z_%%03d.png figures/z.mp4" '    %name
        p = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE);p.wait()
        o,e = p.communicate();
        print(o.decode())


    ###########################################################################
    #### compare :
    ###########################################################################
    def compare_xyz_pxpy(self,frame=32,opts='ma'):
        rpl0    = self.rpl.loc[self.rpl.F==frame]
        pxc,pyc = self.cen.iloc[frame-1][['px','py']].values.T
        alpha   = rpl0.alpha.iloc[0]

        pxy0  = rpl0[['px','py']].values
        sxy   = self.aper*(pxy0-[pxc,pyc])
        xyz0  = np.hstack([sxy,np.zeros((sxy.shape[0],1))]).T
        # x0,y0,z0 = xyz0

        Ro = np.identity(3)
        if 'a' in opts:
            Ro    = get_crystal_rotation(u=[-1,0,0],alpha=alpha)
        if 'm' in opts:
            R180  = get_crystal_rotation(omega=65,alpha=180)
            Ro    = Ro.dot(R180)
        # Ro = get_crystal_rotation(omega=230,eo='x',alpha=-qrpl0.alpha.iloc[0])
        x0,y0,z0  = Ro.dot(xyz0)

        xyz   = rpl0[['x','y','z']].values.T
        # x,y,z  = Ro.dot(xyz)
        x,y,z = xyz

        Im = np.array(rpl0.Im.values/10,dtype=int)
        scat = ([x,y,z,Im,'b','o'],[x0,y0,z0,Im,'r','o'])
        fig,ax  = dsp.stddisp(rc='3d',view=[-90,-90],xylims=1.5,scat=scat,labs=['x','y','z'])
        h3d = h3D.handler_3d(fig,persp=False)

    def compare_hkl(self,frame,eps=None,Smax=0.025,Nmax=5,v=0):
        uvw = self.get_beam_dir(frame=frame)
        df  = self.get_excitation_errors(uvw,Nmax=Nmax,Smax=Smax)
        hkl0 = df[['h','k','l']]

        hkl1 = self.rpl.loc[self.rpl.F==frame] #[['hx','kx','lx']].values
        if eps:
            idx  = (abs(hkl1.hx-hkl1.h)<eps) & (abs(hkl1.kx-hkl1.k)<eps) & (abs(hkl1.lx-hkl1.l)<eps)
            hkl1 = hkl1.loc[idx]
        hkl1 = hkl1[['h','k','l']]

        rrmv = [r for r in hkl1.values if not np.where(np.linalg.norm(r-hkl0.values,axis=1)==0)[0].size]
        rrmv = np.array(rrmv,dtype=int)
        print('%d patterns reflections ' %hkl1.shape[0])
        print('%d lattice  reflections ' %hkl0.shape[0])
        print('reflections in diffraction pattern not in ewald construction with Smax=%.2E :' %Smax)
        print(rrmv)

        if v:
            print('from diffraction patterns :');print(hkl1)
            ridx = []
            for r in hkl1.values:
                idx = np.where(np.linalg.norm(r-hkl0.values,axis=1)==0)[0]
                if idx.size:
                    ridx+=[idx[0]]
            print('from excitation errors    :');
            print(df.iloc[ridx][['h','k','l','Sw']])
            return df.iloc[ridx][['h','k','l','Sw']]

    ###########################################################################
    #### misc :
    ###########################################################################
    def _get_cif_file(self,cif_file):
        if not cif_file:
            cif_file = self._find_files('cif')
        return cif_file


    def _find_files(self,f_type):
        files = glob.glob(self.path+'/*.%s' %f_type)
        if len(files)==1:
            file=files[0]
            return file
        else:
            msg ='''
        only 1 %s file should be found in path folder :
         %s
        but %d were found :
         %s''' %(f_type,os.path.realpath(self.path),len(files),str(files))
            raise Exception(msg)
