''' Build wallpapers from the 17 plane groups
'''
from utils import displayStandards as dsp
from utils import glob_colors as colors
from matplotlib import patches
from matplotlib import collections as coll
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
import os
import numpy as np
import pandas as pd
from . import lattice
from . import nodesTri as nt
#import lattice
#import nodesTri as nt

pp_types = ['p1','p2','pm','pg','cm','pmm','pmg','pgg','cmm',
            'p4','p4m','p4g',  'p3','p3m1','p31m','p6','p6m']; #print(pp_types)
pp_long  = ['p1','p2','p1m1','p1g1','c1m1','p2mm','p2mg','p2gg','c2mm',
            'p4','p4mm','p4gm', 'p3','p3m1','p31m','p6','p6mm']
pp_lattices = ['oblique']*2+['rect']*7+['square']*3+['hex']*5

df_wallpp = pd.DataFrame(
    np.array([pp_long,pp_lattices,]).T,
    columns = ['long_name','lat_type'],
    index = pp_types)

################################# asym unit cells
u1,u21   = np.array([[0.0,0],[1,0],[1,1],[0,1]]),  np.array([[0,0],[1/2,0],[1/2,1],[0,1]])
u22,u42  = u1/2,u21/2
u41      = np.array([[0,0,],[1/4,0],[1/4,1],[0,1]])
u22x = np.array([[0,0],[1/2,0],[1/2,1/2]])
u22y = np.array([[0,0],[1/2,0],[0,1/2]])
u3   = np.array([[0,0],[1/2,-1/3],[1,0],[1/2,1/3]])
u3m1 = np.array([[0,0],[1/2,-1/3],[1/2,1/3]])
u31m = np.array([[0,0],[1,0],[1/2,1/3]])
u6   = np.array([[0,0],[1/2,-1/3],[1/2,1/3]])
u6m  = np.array([[0,0],[1/2,0],[1/2,1/3]])
asym_cells = {'p1':u1,'p2':u21 ,'pm':u21,'pg':u21 ,'cm':u22 ,'pmm':u22 ,
             'pmg':u41,'pgg':u22,'cmm':u42,'p4':u22,'p4m':u22x,'p4g':u22y,
             'p3':u3,'p3m1':u3m1,'p31m':u31m,'p6':u6,'p6m':u6m}
df_wallpp['asym_cell'] = asym_cells.values()


################################# Atomic parameters
V0 = 1.0
Ai = np.array([0.1,0.25,0.26,0.27,1.5])
fv = lambda X,X0,Za : V0*np.exp(-(np.linalg.norm(X-X0,axis=1)/Ai[Za])**2)
#f  = np.sum(np.array([fc(X,X0) for X0 in Xasym]),axis=0)
g,marker=0.5,'o' #^'
atoms = pd.DataFrame.from_dict(                                         # covalent radius in A
    { 0:['H',(0,0,0),0.37], 1:['C',(g,g,g),0.77],
      2:['N',(0,0,1),0.75], 3:['O',(1,0,0),0.73],4:['S',(0,1,1),1.02]},
    orient='index',columns=['atom','color','size'])


class Wallpaper :
    '''Class for creating wallpapers'''
    def __init__(self,pp_type,a=1,b=1,angle=90,pattern=None,gen=True,ndeg=0,nh=3,nk=3,pOpt='',name='',path='.') :
        ''' initialize class
        pp_type   : plane group type (str)
        a,b,angle : lattice parameters (double)
        pattern   : x y f coordinates within asymetric unit cells   (see add_pattern)
        gen       : generate wallpaper applying symmetry operations (see generate_wallpp for ndeg,nh,nk)
        pOpt      : plot options                                    (see plot_unit_cells)
                  : 'w' to see the unit cell if warning is issued
            Note : gen,ndeg,nh,nk will work only if pattern is defined
        ndeg   : number of discretization points along one the direction for building the potential map
        '''
        self.pp_type  = pp_type
        self.name     = name if name else df_wallpp.loc[pp_type]['long_name']
        self.init_lattice(a,b,angle)
        self.init_unit_cell()
        self.path = path
        self.ndeg = ndeg
        if isinstance(pattern,list):pattern=np.array(pattern)
        if isinstance(pattern,np.ndarray) or isinstance(pattern,str):
            self.add_pattern(pattern,pOpt)
            if gen :
                self.generate_wallpp(nh,nk,ndeg)
                if pOpt and 'q' not in pOpt:
                    self.plot_unit_cells(pOpt,nh,nk)#,xylims=[0,1,0,0.5])

    def init_lattice(self,a,b,angle) :
        '''init lattice assigning attributes :
        lat_type    : lattice type
        lattice_vec : lattice vectors
        a,b,angle   : lattice parameters enforcing lattice constraints
        '''
        self.lat_type = df_wallpp.loc[self.pp_type]['lat_type']
        if self.lat_type == 'hex' :
            self.a,self.b,self.angle = a,a,60
            self.lattice_vec = np.array([[np.cos(np.pi/3),-np.sin(np.pi/3)],[np.cos(np.pi/3),np.sin(np.pi/3)]])*a
        elif self.lat_type == 'square' :
            self.a,self.b,self.angle = a,a,90
            self.lattice_vec = np.array([[1,0],[0,1]])*a
        elif self.lat_type == 'rect' :
            self.a,self.b,self.angle = a,b,90
            self.lattice_vec = np.array([[a,0],[0,b]])
        elif self.lat_type == 'oblique' :
            self.a,self.b,self.angle = a,b,angle
            self.lattice_vec = np.array([[a,0],[b*np.cos(np.pi/180*angle),b*np.sin(np.pi/180*angle)]])
        self.params = [self.a,self.b,np.pi/180*self.angle]

    def init_unit_cell(self) :
        '''init unit cell defining :
        - sym_ops   : symmetry operations
        - asym_cell : asymmetric unit cell (matplotlib.patches.Polygon)'''
        self.sym_ops   = wallpp_syms(self.pp_type,self.lattice_vec)
        verts    = asym_cells[self.pp_type]
        a,b,angle = self.params
        if self.lat_type == 'hex' :
            #print('HEX')
            verts[:,0] *=a
            verts[:,1] *=a*np.cos(np.pi/6)
        else :
            verts = self.lattice_vec.T.dot(verts.T).T         # 2x2.2x3 = 2x3
        self.asym_cell = patches.Polygon(verts,alpha=0.3,edgecolor='k',linewidth=1.5,label='$asym_cell$')

    def add_pattern(self,pattern,pOpt='') :
        ''' add a pattern from :
        Nx3 np.array : x y f coordinates within asymetric unit cells
        txt file     : contains  x y f
        where f : list of tags for atoms '''
        if isinstance(pattern,str) : pattern = np.loadtxt(pattern)
        self.X0 = pattern[:,[0,1]]
        self.f0 = pattern[:,2]
        idb = self.asym_cell.contains_points(self.X0)
        if (~idb).sum() :
            print(colors.red+"warning:pattern not in asymetric unit cells"+colors.black);
            pOpt += ['','W']['w' in pOpt]
            print(self.X0);print(self.asym_cell.get_verts());print(idb)
        if 'W' in pOpt:
            a,b = self.lattice_vec
            pcell = patches.Polygon(np.array([[0,0],a,a+b,b]),label='$unit_cell$')
            colls  = [coll.PatchCollection([pcell],alpha=0.1,edgecolor=(g,g,g),linewidth=1),
                      coll.PatchCollection([self.asym_cell],alpha=0.3,edgecolor='k',linewidth=1.5,label='$asym_cell$')]
            (x,y),c = self.X0.T,self.f0
            dsp.stddisp(scat=[x,y,c],colls=colls,opt=pOpt,ms=20,title=self.pp_type,equal=True)
        #self.X0 = X0[idb]

    def generate_wallpp(self,nh,nk,ndeg=0) :
        ''' Generate wallpaper applying symmetry operations :
        ndeg   : order for generating potential map              (see build_potential)
        nh,nk  : number of unit cells to be displayed
        ndeg   : number of discretization points along one the direction for building the potential map'''

        if self.pp_type=='p1' and ndeg:
            self.build_potential_p1(ndeg)
        else:
            self.Xcell,self.fcell   = self.apply_symmetries(self.X0,self.f0)
            self.Xa,self.fa         = self.repeat_pattern(self.Xcell,self.fcell,nh,nk)
            if ndeg :
                self.build_potential(ndeg)
                self.Xc,self.fc = self.apply_symmetries(self.Xp,self.fp)
                self.X,self.f   = self.repeat_pattern(self.Xc,self.fc,nh,nk)

    def build_potential_p1(self,ndeg=2**8):
        print('...building potential...')
        if isinstance(ndeg,int):ndeg=[ndeg]*2
        self.x,self.z = np.linspace(0,self.a,ndeg[0]),np.linspace(0,self.b,ndeg[1])
        x,z = np.meshgrid(self.x,self.z)
        self.Xa,self.fa=self.X0,self.f0
        self.fp=np.zeros(x.shape)
        for X0,Za in zip(self.Xa,self.fa):
            r = np.sqrt((x-X0[0])**2+(z-X0[1])**2)
            self.fp += np.exp(-(r/Ai[int(Za)])**2)
        self.Xp = np.array([x.flatten(),z.flatten()]).T

    def build_potential(self,ndeg=10) :
        ''' Build the potential map :
        ndeg : number of discretization points along one the direction for building the potential map'''
        vert_cell = self.asym_cell.get_verts()[:-1,:]   #last point is sme as first in polygon
        if vert_cell.shape[0]==3 :
            xi,eta = nt.nodesTri(ndeg).T                #get triangular elements
            N  = 0.5*np.array([-xi-eta,1+xi,1+eta]).T   #lin shape functions at xi,eta
            self.Xp  = N.dot(vert_cell)                 #transform : Nx3.3x2=Nx2
        if vert_cell.shape[0]==4 :
            # ndeg+=1
            nlin = np.arange(ndeg)/(ndeg); #np.linspace(0,1,ndeg)
            x,y = np.meshgrid(nlin,nlin)
            N = np.array([x.flatten(),y.flatten()]).T     #get rectangular grid Nx2
            vecs = vert_cell[[1,-1],:]                  #2x2 [a;b]
            self.Xp  = N.dot(vecs)                      #transform : Nx2.2x2=Nx2
        self.Xp += self.lattice_vec[0]+self.lattice_vec[1]
        #print(self.Xp)
        self.fp = np.sum(np.array([fv(self.Xp,X0,int(Za)) for X0,Za in zip(self.Xa,self.fa)]),axis=0)
        self.Xp -= self.lattice_vec[0]+self.lattice_vec[1]
        #self.fp  = fv(self.Xp,self.Xa)                  #potential

    def get_potential_grid(self):
        '''potential on a grid within unit cell
        ADAPT TO P2 etc...'''
        (x,z),f = self.Xc.T,self.fc
        Nxz   = (self.ndeg,self.ndeg)
        x,z,f = x.reshape(Nxz),z.reshape(Nxz),f.reshape(Nxz)
        return x[0,:],z[:,0],f

    def get_potential_grid_p1(self,ndeg=2**8):
        if not 'fp' in self.__dict__.keys():
            self.build_potential_p1(ndeg)
        return self.x,self.z,self.fp

    def plot_unit_cells(self,opts='AV',nh=3,nk=4,lvls=10,alpha=0.5,fg='12',**kwargs) :
        '''Display the wallpaper
        pOpt :  A(Atom), V(potential), u(unit_cell), a(asym_unit)
        '''
        fig,ax = dsp.create_fig(figsize=fg,pad=[2.,5]['t' in opts])
        colls,scat,im,cs,contours=[],[],None,None,None
        lOpt = ''
        if self.lat_type=='hex' : lOpt += 'hq'
        if 'u' in opts : lattice.plot_lattice(self.lattice_vec,max(nh,3),max(nk,3),pOpt=lOpt,ax=ax)
        if 'a' in opts : colls=[coll.PatchCollection([self.asym_cell],alpha=0.3,linewidth=1.5,edgecolor='b')]
        #plot potential
        if 'V' in opts :
            print('display potential')
            # na = int(self.Xcell.shape[0]/self.X0.shape[0])
            N = int(np.sqrt(self.Xp.shape[0]))
            (x,y),f = self.Xp.T, self.fp
            # triang = mtri.Triangulation(x, y, Delaunay(self.Xc).vertices)
            # cs=ax.tricontourf(triang,self.fc,levels=lvls)
            # cont = ax.contourf(x.reshape((N,N)),y.reshape((N,N)),f.reshape((N,N)),levels=lvls)
            im = [x.reshape((N,N)),y.reshape((N,N)),f.reshape((N,N))]#,snap=True)
            # contours = im+[lvls]
            im+=[alpha]
            # cs = ax.pcolormesh(x.reshape((N,N)),y.reshape((N,N)),f.reshape((N,N)))#,snap=True)

        #plot atoms
        if 'A' in opts :
            x,y = self.Xa.T
            s,c = atoms.iloc[self.fa][['size','color']].values.T
            scat=[x.flatten(),y.flatten(),np.int_(50*np.array(s)),c]
        dsp.stddisp(fig=fig,ax=ax,labs=[r'$x(\AA)$',r'$y(\AA)$'],
            scat=scat,colls=colls,im=im,contour=contours,
            title=self.name,name=self.path+self.name,
            **kwargs)

    def get_reciprocal_lattice_2D(self):
        a1,a2 = self.lattice_vec
        b1,b2 = lattice.reciprocal_lattice_2D(a1,a2)
        return b1,b2

    #private (should be)
    def apply_symmetries(self,X0,f0) :
        Xcell = X0
        for sym in self.sym_ops :
            Xcell = np.concatenate((Xcell,sym(X0)))
        fcell = np.tile(f0,len(self.sym_ops)+1)
        return Xcell,fcell

    def repeat_pattern(self,Xcell,fcell,nh=3,nk=4) :
        u_hk,miller = lattice.get_miller(self.lattice_vec,nh,nk)
        X = np.empty(shape=(0, 2))
        for u in u_hk :
            X = np.concatenate((X,u+Xcell),axis=0)
        f = np.tile(fcell,miller.shape[0])
        #print(X.shape,f.shape)
        return X,np.array(f,dtype=int)



################################################################
#Define symmetry operations
def mirr(x,u,A) :
    ''' Symmetry with respect to line :
    u : direction
    A : Point belonging to the line
    '''
    #print(u,np.linalg.norm(u))
    u = u/np.linalg.norm(u)
    n  = np.array([-u[1],u[0]]) #;print(n);print(x)
    xa = x-np.array(A)
    return x-2*n*xa.dot(n)[:,None]
def glax(x,u,A):
    ''' Glide axis transformation
    u : direction of the glide axis
    A : Point belonging to the line
    '''
    #u = np.linalg.norm(u)
    xM = mirr(x,u,A)
    xg = xM+0.5*u
    return xg

def wallpp_syms(pp_type,lattice_vec):
    cp,sp=np.cos(2*np.pi/3),np.sin(2*np.pi/3)
    a1,a2 = lattice_vec
    #a,b = lat_params[:2]
    O =  np.array([0,0])
    M = (a1+a2)/2; #print(M)
    Pp,Pm = a1+a2,a1-a2
    x0 = 0*np.array([1/2,1/3])
    rot = lambda t: np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])
    centr = lambda x:M+x
    invc  = lambda x:Pp-x
    invc2 = lambda x:M-x
    def mir_uA(u,A) : return lambda x:mirr(x,u,A)
    def gla_uA(u,A) : return lambda x:glax(x,u,A)
    mir2 = mir_uA(a2,a1/2)
    def rot3(i) : return lambda x:rot(2*i*np.pi/3).dot(x.T).T
    def rot3m1(i) : return lambda x:rot(2*i*np.pi/3).dot(mirr(x,Pm,O).T).T
    def rot31m(i) : return lambda x:rot(2*i*np.pi/3).dot(mirr(x,Pp,O).T).T
    def rot4(i) : return lambda x:rot(i*np.pi/2).dot((x-M).T).T+M
    def rot6(i) : return lambda x:rot(i*np.pi/3).dot(x.T).T
    def rot6m(i) : return lambda x:rot(i*np.pi/3).dot(mirr(x,a1,O).T).T
    wallpp_sym={
       'p1'  :[],
       'p2'  :[invc],
       'pm'  :[mir_uA(a2,a1/2)],
       'pg'  :[gla_uA(a1,a2/2)],
       'cm'  :[mir_uA(a1,a2/2),centr,gla_uA(a1,a2/4)],
       'pmm' :[invc,mir_uA(a1,a2/2),mir_uA(a2,a1/2)],
       'pmg' :[invc,mir_uA(a2,a1/4),gla_uA(a1,a2/2)],
       'pgg' :[invc,gla_uA(a2,a1/4),gla_uA(a1,a2/4)],
       'cmm' :[invc,mir_uA(a1,a2/2),mir_uA(a2,a1/2),centr,invc2,gla_uA(a2,a1/4),gla_uA(a1,a2/4)],
       'p4'  :[invc,rot4(1),rot4(3)],
       'p4m' :[invc,rot4(1),rot4(3),mir_uA(a1,a2/2),mir_uA(a2,a1/2),mir_uA(Pp,M),mir_uA(Pm,M)],
       'p4g' :[invc,rot4(1),rot4(3),gla_uA(a2,a1/4),gla_uA(a1,a2/4),mir_uA(Pm,M/2),gla_uA(Pp,O)],
       'p3'  :[rot3(1),rot3(2)],
       'p3m1':[rot3(1),rot3(2),rot3m1(0),rot3m1(1) ,rot3m1(2)],
       'p31m':[rot3(1),rot3(2),rot31m(0),rot31m(1) ,rot31m(2)],
       'p6'  :[rot6(i) for i in range(1,6)],
       'p6m' :[rot6(i) for i in range(1,6)]+[rot6m(i) for i in range(6)]
        }[pp_type]
    return wallpp_sym

#######################################################################
# test
#######################################################################
def test_wallpps() :
    fig_path = os.path.abspath(os.path.dirname(__file__)+'/..') +'/docs/figures/'
    a,b,angle=2,1,70
    gen,pOpt =  True,'VA' #V
    pOpt += 'aue sc'
    ndeg,nh,nk = 30,3,2
    wallpps=dict()
    #pattern asymmetric unit cell in real coordinates
    pattern=np.array([[0.9,0.3,1],[0.3,0.15,2],[0.45,0.4,3]])
    for pp_type in pp_types :#['cmm'] :#,'cmm'] :# pp_types :
        print(pp_type)
        wallpps[pp_type] = Wallpaper(pp_type,a,b,angle,pattern,
                                     gen,ndeg,nh,nk,
                                     pOpt,path=fig_path)
    #plt.show()
    return wallpps

if __name__=='__main__' :
    wallpps = test_wallpps()
    print(green+__file__.split('/',)[-1]+' success'+black)
