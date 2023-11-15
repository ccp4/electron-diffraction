# import importlib as imp
import re,os,glob,json,numpy as np, pandas as pd, tifffile,mrcfile,scipy.optimize as opt
from crystals import Crystal,Lattice
from utils import physicsConstants as cst
from . import utilities as ut               #;imp.reload(ut)
from . import readers                       #;imp.reload(readers)


class Dataset:
    def R(self,pz) :
        '''Rotation matrix at frame pz
        Parameters:
            pz : frame number
        '''
        a0 = self.info['STARTING_ANGLE']
        da = self.info['OSCILLATION_RANGE']
        f0 = self.info['DATA_RANGE'][0]
        angle = a0 + da*(pz-f0)
        R = ut.rotation_matrix(
            self.info['ROTATION_AXIS'],
            angle,
            deg=True)
        return R

    def get_image_size(self):
        if self.frame_folder:
            frame_folder=self.frame_folder
            d_frames = readers.detect_frame(os.path.join(self.path,frame_folder))
            if d_frames:
                fmt = d_frames['fmt']
                frames = glob.glob(os.path.join(self.path,frame_folder,'*.%s' %fmt))
                self.nxy = readers.read(frames[0]).shape

    def init_geom(self):
        self.aper = 1/abs(self.lam*self.F/self.dx) #;print('aper : %.4f recAngstrom' %self.aper)
        self.keV = cst.lam2keV(self.lam)

        lat = Lattice.from_parameters(*self.lat_params)
        self.Alat=np.array(lat.lattice_vectors).T
        self.Arec=np.array(lat.reciprocal_vectors).T/(2*np.pi)
        a,b,c = self.A
        a_star = np.cross(b,c)/a.dot(np.cross(b,c))
        b_star = np.cross(c,a)/b.dot(np.cross(c,a))
        c_star = np.cross(a,b)/c.dot(np.cross(a,b))
        self.UB = np.array([a_star,b_star,c_star]).T

    def get_uvw(self):
        '''(u,v,w) are the components of the orientation vector in reciprocal cartesian space
        It is such that K-g as it appears in the blochwave excitation error
        can be computed where the components of the beam g are simply calculated
        from the crystal reciprocal lattice vectors and the beam miller indices (h,k,l)
        g = lat_rec.dot(hkl)
        The beam is given is real space lab frame
        '''
        self.frames=np.arange(self.info['DATA_RANGE'][0],self.info['DATA_RANGE'][1])
        beam = self.info['INCIDENT_BEAM_DIRECTION']
        self.uvw0 = np.array([
            self.Arec.dot(np.linalg.inv(self.R(f).dot(self.UB))).dot(beam)
                for f in self.frames])
        self.n_frames = self.uvw0.shape[0]
        self.cen = pd.DataFrame([[self.orgx,self.orgy]]*self.n_frames,
            columns=['px','py'])

    def sw(self,hkl,frame):
        g_lab = self.R(frame).dot(self.UB).dot(hkl.T)
        s0 = self.info['INCIDENT_BEAM_DIRECTION']/self.lam

        K0 = 1/self.lam
        qx,qy,qz = g_lab
        Kx,Ky,Kz = s0
        Sw = (K0**2-((Kx+qx)**2+(Ky+qy)**2+(Kz+qz)**2))/(2*K0)
        return Sw

    def hkl_to_pixels(self,h,frame):
        ''' convert miller indices to pixel locations
            - h : list of miller indices
            - frame : frame location of the miller indices
        '''
        # pz = self.rpl.loc[hkl,'pz'] #frame location of reflection
        #hlist of reflections with shape :n_refls x 3
        hkl = np.array([eval(hkl) for hkl in h])

        #### #reflection in reciprocal lab frame
        r = self.R(frame).dot(self.UB).dot(hkl.T).T
        #incident beam in reciprocal lab frame
        s0 = self.info['INCIDENT_BEAM_DIRECTION']/self.lam
        # scattering vector in reciprocal lab frame
        s = r+s0                                    #;print('s : ',s)

        #### reflection in reciprocal Orthonormal crystal frame
        # r = self.Arec.dot(hkl.T).T
        # #Incident Beam vector in reciprocal Orthonormal crystal frame
        # M = self.Arec.dot(np.linalg.inv(self.R(frame).dot(self.UB)))
        # u = M.dot(self.info['INCIDENT_BEAM_DIRECTION'])
        # s0 = u/(self.lam*np.linalg.norm(u))
        # #scattering vector in reciprocal Orthonormal crystal frame
        # s1 = r+s0
        # #scattering vector in reciprocal lab frame
        # s = (np.linalg.inv(M).dot(s1.T)).T          #;print('s1 : ',s1)


        #scattering vector direction in lab frame
        norms = np.linalg.norm(s,axis=1)            #; print('norms:',norms)
        s=(s.T/np.linalg.norm(s,axis=1)).T          #; print('s:\n',s);print('norms:',np.linalg.norm(s,axis=1))

        #location of spot in real space lab frame
        x = np.array([
            np.tan(np.arctan2(s[:,0],s[:,2]))*self.F,
            np.tan(np.arctan2(s[:,1],s[:,2]))*self.F,
            self.F*np.ones(s[:,1].shape)]).T
        #print('x:\n',x)
        # spot location in detector frame
        xd = np.linalg.inv(self.ED).dot(x.T).T          #;print('xd:',xd)
        px,py = self.xd_to_px(xd)

        #package
        df_pxy = pd.DataFrame()
        df_pxy['px'] = px
        df_pxy['py'] = py
        df_pxy.index = h
        # pixel locations

        return df_pxy


    def to_shelx(self,hkl,file='',output_dir=None):
        '''converts to a .hkl file ready to use by shelx
        hkl:
        '''
        if not isinstance(output_dir,str):output_dir=self.path
        if not file:file='refl.hkl'
        file_path = os.path.join(output_dir, file)
        ut.to_shelx(hkl,file_path)
