'''Defines 2D lattices with some utilities'''
from utils import glob_colors as colors
from utils import displayStandards as dsp
import numpy as np

lat_types = ['square','rect','oblique','hex']

class Lattice2D:
    def __init__(self,lat_type='oblique',a=1,b=1,alpha=90):
        self.lat_type,self.a,self.b,self.alpha=lat_type,a,b,alpha
        self.params = {'a':a,'b':b,'alpha':alpha}
        self.lattice_vectors   =set_lattice_vec(lat_type,**self.params)
        self.reciprocal_vectors=get_reciprocal_vectors_2D(self.lattice_vectors)
        self.lat_vec = self.lattice_vectors
        self.rec_vec = self.reciprocal_vectors
        self.a1,self.a2 = self.lat_vec
        self.b1,self.b2 = self.rec_vec
        self.area=abs(np.cross(self.a1,self.a2))

    def get_vectors(self):return self.lattice_vectors
    def get_reciprocal_vectors(self):return self.reciprocal_vectors
    def get_lattice(self,nh,nk):
        ''' Get the miller indices and corresponding vectors
            - nh,nk : number of indices on each direction
        returns :
            - hk : Nx2 array of indices
            - R  : Nx2 array of vectors
        '''
        if type(nh) in [int,np.int64]:
            h,k = np.arange(nh),np.arange(nk)
            if nh<0:h = np.arange(nh,0)
            if nk<0:k = np.arange(nk,0)
        elif type(nh) in [np.ndarray,list]:
            h,k = nh,nk
        else:
            raise Exception('nh need be int or ndarray')
        h,k = np.meshgrid(h,k)
        hk  = np.vstack([h.flatten(),k.flatten()]).T
        R   = hk.dot(self.lattice_vectors) # all points of lattice
        return hk,R

    def plot_lattice(self,nh=3,nk=3,gridcolor=(0.5,)*3,hOpt=True,**kwargs) :
        '''Plot lattice grid :
        - nh,nk : number of cells along each direction
        - hOpt : additional lines for hexagonal lattice
        '''
        a1,a2=self.a1,self.a2
        hk,R = self.get_lattice(nh,nk)
        hs,ks = np.arange(nh+1),np.arange(nk+1)
        if nh<0:hs = np.arange(nh,1)
        if nk<0:ks = np.arange(nk,1)
        line_h = [h*a1 + np.array([[0,0],nk*a2]) for h in hs]
        line_k = [k*a2 + np.array([[0,0],nh*a1]) for k in ks]
        plts = [[xy[:,0],xy[:,1],gridcolor,''] for xy in line_h+line_k]
        if hOpt and self.lat_type=='hex':
            line_hex = [ np.vstack([r+a1,r+a2]) for r in R]
            plts += [[xy[:,0],xy[:,1],[gridcolor,'--'],''] for xy in line_hex]
        return dsp.stddisp(plts,**kwargs)

###########################################################################
# utils
###########################################################################
def get_reciprocal_vectors_2D(lat_vec):
    '''get reciprocal lattice vectors 2D
        - lat_vec : lattice vectors [a1;a2]
    returns :
        - rec_vec : lattice vectors [a1;a2]
    '''
    a1,a2   = lat_vec
    b1,b2   = reciprocal_lattice_2D(a1,a2)
    rec_vec = np.array([b1,b2])
    return rec_vec

def reciprocal_lattice_2D(a1,a2):
    Rot90 = np.array([[0,1],[-1,0]])
    b1 = Rot90.dot(a2)/(a1.dot(Rot90.dot(a2)))#*2*np.pi
    b2 = Rot90.dot(a1)/(a2.dot(Rot90.dot(a1)))#*2*np.pi
    return b1,b2

def set_lattice_vec(lat_type='',a=1,b=1,alpha=90,v=0) :
    ''' Computes lattice vectors from :
        - a,b,angle : lattice parameters
        - lat_type  : lattice types [hex,square,rect,oblique]
        - v : change orientation for oblique lattice
    '''
    if lat_type == 'square' :
        lattice_vec = np.array([[1,0],[0,1]])*a
    elif lat_type == 'rect' :
        lattice_vec = np.array([[a,0],[0,b]])
    elif lat_type == 'hex' :
        lattice_vec = np.array([[1,0],[np.cos(np.pi/3),np.sin(np.pi/3)]])*a
    elif lat_type == 'oblique' :
        angle=np.deg2rad(alpha)
        cp,sp = np.cos(angle),np.sin(angle)
        lattice_vec = np.array([[a,0],[b*cp,b*sp]])
        if v : lattice_vec = np.array([[a*cp,a*sp],[0,b]])
    else :
        msg ='''%s not a valid lattice type. choose either:
        %s''' %(lat_type,str(lat_types))
        raise Exception(msg)
    return lattice_vec


def params_from_lat_vec(lat_vec):
    a1,a2=lat_vec
    a = np.linalg.norm(a1)
    b = np.linalg.norm(a2)
    alpha = np.rad2deg(np.arccos(a1.dot(a2)/(a*b)))
    return {'a':a,'b':b,'alpha':alpha}

def get_miller(lattice_vec,nh,nk,ax=None):
    ''' get the miller indices and correspoding vectors in format :
            nh.nk x 2 '''
    n_h,n_k = np.arange(nh),np.arange(nk)
    miller = np.array(np.meshgrid(n_h,n_k)).T.reshape(((nh)*(nk),2))
    u_h = lattice_vec.T.dot(miller.T).T
    if ax :
        for u in u_h:
            ax.plot([0,u[0]],[0,u[1]],'s-g')
    return u_h,miller
