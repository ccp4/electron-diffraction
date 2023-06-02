import importlib as imp
import os,glob,json,numpy as np, pandas as pd, tifffile,mrcfile,scipy.optimize as opt
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from subprocess import Popen,PIPE,check_output
from utils import glob_colors as colors
from crystals import Crystal,Lattice
from EDutils import utilities as ut         #;imp.reload(ut)
from multislice.rotating_crystal import get_crystal_rotation
from gemmi import cif
from . import import_ED as ED               ;imp.reload(ED)

class Dials(ED.Dataset):
    def __init__(self,path:str):
        '''Importing dials information
        It expects a .expt file and a reflections.txt file
        '''
        self.path = path
        refl_file = os.path.join(self.path,'reflections.txt')
        if not os.path.exists(refl_file):
            raise Exception('%s file not found' %refl_file)

        df = load_dials_reflections(refl_file)
        self.rpl = df[['h','k','l','I']].copy()
        self.rpl['hkl'] = [str((h,k,l)) for h,k,l in self.rpl[['h','k','l']].values]
        self.rpl[['qx','qy']] = df[['s1x','s1y']]
        self.rpl[['px','py']] = df[['o_px','o_py']]
        # self.rpl[['qy']]*=-1
        self.rpl['F'] = np.array(np.round(df['o_pz'])+1,dtype=int)
        self.index = self.rpl.hkl

        ######
        self.read_expt_file()
        self.init_geom()
        self.get_uvw()
        self.uvw0*=-1

        dyn_cif_file=os.path.join(self.path,'dials_dyn.cif_pets')
        if os.path.exists(dyn_cif_file):
            self.read_dyn_cif(dyn_cif_file)

    def xd_to_px(self,xd):
        '''
        - xd : Nx3 spot locations in the detector coordinate system
        '''
        px = xd[:,0]/xd[:,2]/self.dx  #;print('px:',px)
        py = xd[:,1]/xd[:,2]/self.dy  #;print('py:',py)
        return px,py

    def read_expt_file(self):
        expt_file=''
        for s in ['integrated','refined','indexed']:
            expt_file=os.path.join(self.path,'%s.expt' %s)
            if os.path.exists(expt_file):
                break
        try:
            with open(expt_file,'r') as f:
                expt = json.load(f)
        except:
            raise Exception('To import dials you at need a indexed.expt or refined.expt, integrated.expt ')

        print(colors.blue+'reading %s' %expt_file+colors.black)
        # expt = expt
        self.beam = -np.array(expt['beam'][0]['direction'])
        self.lam  = expt['beam'][0]['wavelength']
        detector  = expt["detector"][0]["panels"][0]
        self.orgx,self.orgy,self.F  = detector['origin']
        self.dx,self.dy             = detector['pixel_size']
        self.nxy                    = detector['image_size']
        oscillation = expt['scan'][0]['oscillation']
        self.hm     = expt['crystal'][0]['space_group_hall_symbol']

        #real space orientation matrix [a;b;c]
        self.A = np.array([
            expt['crystal'][0]['real_space_%s' %s] for s in 'abc'
            ])

        lat_c = np.linalg.norm(self.A,axis=1)
        angles = [ np.rad2deg(np.arccos(
                self.A[i].dot(self.A[j])/(lat_c[i]*lat_c[j])
            ))
            for i,j in zip([0,1,0],[1,2,2])]
        self.lat_params=[*lat_c,*angles]

        self.info={
            'ROTATION_AXIS'                 : expt['goniometer'][0]['rotation_axis'],
            'DATA_RANGE'                    : expt['scan'][0]['image_range'],
            'STARTING_ANGLE'                : oscillation[0],
            'OSCILLATION_RANGE'             : oscillation[1],
            'INCIDENT_BEAM_DIRECTION'       : self.beam,
            'DIRECTION_OF_DETECTOR_X-AXIS'  : detector['fast_axis'],
            'DIRECTION_OF_DETECTOR_Y-AXIS'  : detector['slow_axis'],
            'UNIT_CELL_CONSTANTS'           : self.lat_params,
        }
        e1,e2 = self.info['DIRECTION_OF_DETECTOR_X-AXIS'],self.info['DIRECTION_OF_DETECTOR_Y-AXIS']
        e3 = np.array([self.orgx,self.orgy,self.F])
        ED = np.array([e1,e2,e3]).T
        self.ED = ED


    def read_dyn_cif(self,dyn_cif_file):
        print(colors.blue+'reading %s' %dyn_cif_file+colors.black)
        d = load_dyn(dyn_cif_file)
        self.UB2 = d['UB']
        self.uvw = np.stack([d['u'],d['v'],d['w']]).T
        self.lat_params = d['lat_params']
        self.lat = np.array(Lattice.from_parameters(*self.lat_params).lattice_vectors)

        for k in ['alpha','n_frames'] : self.__dict__[k] = d[k]
        beams = self.lat.T.dot(self.uvw.T).T
        beams/=np.linalg.norm(beams,axis=1)[:,None]
        self.uvw0 = -beams

        # #pb dials duplicate first frame
        self.n_frames-=1
        self.alpha=self.alpha[1:]

        # self.uvw0=self.uvw0[1:]
        # print(np.cross(self.uvw0[0],self.uvw0[-1]))

        # predict_txt = os.path.join(self.path,'predict.txt')
        # if os.path.exists(predict_txt):
        #     self.df_pred = load_dials_predict(predict_txt)
        #
        # files = glob.glob(os.path.join(path,'**'),recursive=True)
        # imfile = ''
        # for ext in ['tiff','mrc']:
        #     images =[f for f in files if f.split('.')[-1]==ext]
        #     if any(images):
        #         imfile = images[0]
        #         print(colors.green+'image format %s detecte in %s' %(ext, imfile)+colors.black)
        # if imfile:
        #     if ext=='mrc':
        #         with mrcfile.open(imfile) as f:
        #             self.aper = f.extended_header["Pixel size X"][0]*1e-10 #A^-1
        #             self.nxy  = int(f.header.nx)
        #             self.Imax =  f.data.max()

def load_dials_predict(refl_txt):
    tmp = 'tmp.txt'
    cmd = "tail --lines=+25 %s | "\
        "sed -E 's/ +/ /g' | sed -E 's/^ //' | sed 's/,//g' "\
        "> %s" %(refl_txt, tmp)

    out = check_output(cmd,shell=True).decode()
    if out:print(out)
    df = pd.read_csv(tmp,
        names=[
            'h','k','l','id','panel','flag',
            'x','y','z',
            'px','py','pz',
            's1x','s1y','s1z',
            ],
        engine="python",sep=" ")

    out=check_output("rm %s" %tmp ,shell=True).decode()
    if out:print(out)
    df.index=[str(tuple(h)) for h in df[['h','k','l']].values]
    df['hkl']=df.index
    df['F'] = np.array(np.round(df['pz'])+1,dtype=int)
    return df

def load_dials_reflections(refl_txt):
    print(colors.blue+'reading %s' %refl_txt+colors.black)
    tmp = 'tmp.txt'
    cmd = "tail --lines=+39 %s | "\
        "sed -E 's/ +/ /g' | sed -E 's/^ //' | sed 's/,//g' "\
        "> %s" %(refl_txt, tmp)

    out = check_output(cmd,shell=True).decode()
    if out:print(out)
    df = pd.read_csv(tmp,
        names=[
            'h','k','l','id','i','panel','flag','I','sig',
            'c_qx','c_qy','c_qz','c_px','c_py','c_pz',
            'o_qx','o_qy','o_qz','vqx','vqy','vqz',
            'o_px','o_py','o_pz','vpx','vpy','vpz',
            's1x','s1y','s1z','Npix','rlpx','rlpy','rlpz'],
        engine="python",sep=" ")

    out=check_output("rm %s" %tmp ,shell=True).decode()
    if out:print(out)
    return df

def load_dyn_intensities(file_dyn):
    doc = cif.read_file(file_dyn )
    block = doc.sole_block()
    h = np.array([float(i) for i in list(block.find_loop('_refln_index_h'))])
    k = np.array([float(i) for i in list(block.find_loop('_refln_index_k'))])
    l = np.array([float(i) for i in list(block.find_loop('_refln_index_l'))])
    I = np.array([float(i) for i in list(block.find_loop('_refln_intensity_meas'))])
    sig = np.array([float(i) for i in list(block.find_loop('_refln_intensity_sigma'))])
    F = np.array([float(i) for i in list(block.find_loop('_refln_zone_axis_id'))])

    df = pd.DataFrame(columns=['h','k','l','I','sig','frame'])
    return df

def load_dyn(file_dyn):
    #  Now use dyn.cif_pets to set up the felix simulation input files
    doc = cif.read_file(file_dyn )
    block = doc.sole_block()

    # Extract data:
    # first the frames and their orientation for the felix.inp file
    # frame number fn
    fn = np.array([int(i) for i in list(block.find_loop('_diffrn_zone_axis_id'))])
    n_frames = len(fn)  # no. of frames
    # incident beam orientation u,v,w
    u = np.array([float(i) for i in list(block.find_loop('_diffrn_zone_axis_u'))])
    v = np.array([float(i) for i in list(block.find_loop('_diffrn_zone_axis_v'))])
    w = np.array([float(i) for i in list(block.find_loop('_diffrn_zone_axis_w'))])

    # Frame angles
    frame_alpha = np.array([float(i) for i in
                            list(block.find_loop('_diffrn_zone_axis_alpha'))])
    # Unit cell
    cell_a = float(block.find_value('_cell_length_a'))
    cell_b = float(block.find_value('_cell_length_b'))
    cell_c = float(block.find_value('_cell_length_c'))
    alpha = float(block.find_value('_cell_angle_alpha'))
    beta  = float(block.find_value('_cell_angle_beta'))
    gamma = float(block.find_value('_cell_angle_gamma'))
    UB = np.reshape(
        [float(block.find_value('_diffrn_orient_matrix_UB_%d%d' %(i,j)))
            for i,j in zip(*[x.flatten()
                for x in np.meshgrid(*(np.arange(3)+1,)*2) ])
        ],(3,3)).T
    cell=[cell_a,cell_b,cell_c]
    lat_params = [cell_a,cell_b,cell_c,alpha,beta,gamma]
    return dict(u=u,v=v,w=w,n_frames=n_frames,alpha=frame_alpha,cell=cell,lat_params=lat_params,UB=UB)
