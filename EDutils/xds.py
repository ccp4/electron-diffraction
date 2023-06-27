import importlib as imp
import os,glob,io,re, numpy as np, pandas as pd, tifffile,mrcfile,scipy.optimize as opt
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from subprocess import Popen,PIPE,check_output
from crystals import Crystal,Lattice
from utils import physicsConstants as cst
from . import utilities as ut;imp.reload(ut)
from . import import_ED as ED               ;imp.reload(ED)

class XDS(ED.Dataset):
    def __init__(self,xds_ascii:str,):
        self.read_xds_ascii(xds_ascii)

        # A is the orientation matrix in lab space (not reciprocal space)
        self.A = np.array([self.info[s] for s in
                ['UNIT_CELL_A-AXIS','UNIT_CELL_B-AXIS','UNIT_CELL_C-AXIS']
            ])
        self.lam = self.info['X-RAY_WAVELENGTH']
        self.F   = self.info['DETECTOR_DISTANCE']
        self.orgx,self.orgy = np.array([self.info['ORGX'],self.info['ORGY']])
        self.dx,self.dy = self.info['QX'],self.info['QY']
        self.nxy = np.array([self.info['NX'],self.info['NY']],dtype=int)
        self.lat_params = self.info["UNIT_CELL_CONSTANTS"]

        self.init_geom()
        self.get_uvw()

        # lat = Lattice.from_parameters(*self.info["UNIT_CELL_CONSTANTS"])
        # Alat=np.array(lat.lattice_vectors).T
        # Arec=np.array(lat.reciprocal_vectors).T/(2*np.pi)
        # self.Arec=Arec
        # self.R0 = np.linalg.inv(Alat).dot(self.A.T)
        # self.UB=Arec.dot(np.linalg.inv(Alat).dot(self.A.T))
        #
        # #the proper way
        # a,b,c = self.A
        # a_star = np.cross(b,c)/a.dot(np.cross(b,c))
        # b_star = np.cross(c,a)/b.dot(np.cross(c,a))
        # c_star = np.cross(a,b)/c.dot(np.cross(a,b))
        # self.UB2 = np.array([a_star,b_star,c_star]).T

        # print(Alat)
        # print(self.A)
        # Alab=np.linalg.inv(Alat).dot(self.A)
        # print(Alab,np.linalg.norm(Alab,axis=1))
        # print(Arec)
        # print(self.UB)


        self.get_qxy()

    def xd_to_px(self,xd):
        px = xd[:,0]/self.dx + self.orgx            #;print('px:',px)
        py = xd[:,1]/self.dy + self.orgy            #;print('py:',py)
        return px,py

    def read_xds_ascii(self,xds_ascii):
        self.path=os.path.dirname(xds_ascii)
        with open(xds_ascii,'r') as f:
            lines = f.readlines()
        info = dict()
        keys = iter(['DATA_RANGE','ROTATION_AXIS','OSCILLATION_RANGE','STARTING_ANGLE',
            'UNIT_CELL_CONSTANTS','UNIT_CELL_A-AXIS','UNIT_CELL_B-AXIS','UNIT_CELL_C-AXIS',
            'X-RAY_WAVELENGTH','INCIDENT_BEAM_DIRECTION','NX','ORGX','DETECTOR_DISTANCE',
            'DIRECTION_OF_DETECTOR_X-AXIS','DIRECTION_OF_DETECTOR_Y-AXIS',
            'END_OF_HEADER',
            ])
        lines=iter(lines)
        l,k=next(lines),next(keys)
        while '!END_OF_HEADER' not in l:
            # print(l,k)
            l = re.sub("\s{2,}"," ",l[1:].strip())
            if k in l:
                v=l.split("=")
                if len(v)>2:
                    v = ''.join(v).split(" ")
                    names,vals = v[::2],v[1::2]
                    for key,val in zip(names,vals):
                        info[key] = float(val)
                else:
                    val=np.array(v[1].split(" ")[1:],dtype=float)
                    if val.size==1:val=val[0]
                    info[k]=val
                k=next(keys)
            l = next(lines)
        # # print(info)

        self.info=info
        self.info["DATA_RANGE"]=np.array(info["DATA_RANGE"],dtype=int)

        #detector
        e1,e2 = self.info['DIRECTION_OF_DETECTOR_X-AXIS'],self.info['DIRECTION_OF_DETECTOR_Y-AXIS']
        e3 = np.cross(e1,e2)
        ED = np.array([e1,e2,e3]).T
        self.ED = ED

        lines=''.join(list(lines)[:-1])
        hkl = pd.read_csv(io.StringIO(lines),engine="python",delimiter=' +',
            names = ['h','k','l','I','sigma','px','py','pz','rlp','peak','corr','psi'],
            )
        hkl.index = [str(tuple(h)) for h in hkl[['h','k','l']].values ]
        hkl['F'] = hkl.pz.round()
        hkl['hkl'] = hkl.index
        self.rpl = hkl[['h','k','l','I','px','py','pz','F','hkl']].copy()
        self.hkl = self.rpl[['I']].copy()

    def get_qxy(self):
        ED = self.ED
        ix,iy = self.rpl.px,self.rpl.py
        orgx,orgy = self.orgx,self.orgy
        qx,qy = self.dx,self.dy
        lam = self.lam
        x = qx*(ix-orgx)*ED[0,0] + qy*(iy-orgy)*ED[0,1] + self.F*ED[0,2]
        y = qx*(ix-orgx)*ED[1,0] + qy*(iy-orgy)*ED[1,1] + self.F*ED[1,2]
        z = qx*(ix-orgx)*ED[2,0] + qy*(iy-orgy)*ED[2,1] + self.F*ED[2,2]
        s = np.array([x,y,z])
        # print('s:\n',s[:,:5].T)
        s_norms=np.linalg.norm(s,axis=0)
        # print('s norms : \n',s_norms[:5])
        s = (s/s_norms).T
        # s/=lam
        # print('s normalized : \n',s[:5,:])
        self.rpl[['s1_x','s1_y','s1_z']] = s
        # s0 = np.array([0,0,1])
        s0 = self.info['INCIDENT_BEAM_DIRECTION']
        r  = (s-s0)/lam
        # print('r:\n',r[:5,:])
        self.rpl[['r_x','r_y','r_z']] = r
        self.rpl[['qx','qy']] = self.rpl[['r_x','r_y']]
