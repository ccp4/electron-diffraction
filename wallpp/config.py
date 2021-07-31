import importlib as imp
import pandas as pd
import numpy as np
from . import lattice as lat2D;imp.reload(lat2D)

pp_types = ['p1','p2','pm','pg','cm','pmm','pmg','pgg','cmm',
            'p4','p4m','p4g',  'p3','p3m1','p31m','p6','p6m']; #print(pp_types)

#########################################################################
#### Lattice types
#########################################################################
df_wallpp = pd.DataFrame(index=pp_types)
df_wallpp['pp_long'] = [
    'p1','p2',
    'p1m1','p1g1','c1m1','p2mm','p2mg','p2gg','c2mm',
    'p4','p4mm','p4gm',
    'p3','p3m1','p31m','p6','p6mm']
df_wallpp['lat_type'] = ['oblique']*2+['rect']*7+['square']*3+['hex']*5


#########################################################################
#### Asymmetric unit cells (fractional coordinates)
#########################################################################
u1       = np.array([[0,0],[1,0],[1,1],[0,1]])
u21      = np.array([[0,0],[0.5,0],[0.5,1],[0,1]])
u22,u42  = u1/2,u21/2
u41      = np.array([[0,0],[1/4,0],[1/4,1],[0,1]])
u22x     = np.array([[0,0],[1/2,0],[1/2,1/2]])
u22y     = np.array([[0,0],[1/2,0],[0,1/2]])
# u3   = np.array([[0,0],[1/2,-1/3],[1,0],[1/2,1/3]])
ugg  = np.array([[0.5,0],[0,0.5],[1,0.5]])
u3   = np.array([[1/3,1/3],[1,0],[2/3,2/3],[0,1]])
u3m1 = np.array([[1/3,1/3],[2/3,2/3],[0,1]])
u31m = np.array([[0,0],[1,0],[1/3,1/3]])
# u31m = np.array([[0,0],[1,0],[1/2,1/3]])
u6   = u31m
u6m  = np.array([[0,0],[1/2,0],[1/3,1/3]])
asym_cells = {
    'p1':u1,'p2':u21 ,'pm':u21,'pg':u21 ,'cm':u22 ,
    'pmm':u22 ,'pmg':u22,'pgg':ugg,'cmm':u42,
    'p4':u22,'p4m':u22x,'p4g':u22y,
    'p3':u3,'p3m1':u3m1,'p31m':u31m,'p6':u6,'p6m':u6m}
df_wallpp['asym_cell'] = asym_cells.values()

#########################################################################
#### Symmetry operations
#########################################################################
bar2   = lambda x,y:np.vstack([1-x,1-y]).T
r2     = lambda x,y:np.vstack([1-x,1-y]).T
bar41  = lambda x,y:np.vstack([y,1-x]).T
mir_v  = lambda x,y:np.vstack([1-x,y]).T
mir_h  = lambda x,y:np.vstack([x,1-y]).T
gla_v  = lambda x,y:np.vstack([1-x,y+1/2]).T
gla_h  = lambda x,y:np.vstack([x+1/2,1-y]).T
gla_v2 = lambda x,y:np.vstack([1/2-x,y+1/2]).T
gla_h2 = lambda x,y:np.vstack([x+1/2,1/2-y]).T
mir_d  = lambda x,y:np.vstack([y,x]).T
mir_e  = lambda x,y:np.vstack([1-y,1-x]).T
mir_f  = lambda x,y:np.vstack([1/2-y,1/2-x]).T
r31    = lambda x,y:np.vstack([1-x-y,x]).T
r32    = lambda x,y:np.vstack([y,1-x-y]).T
m6     = lambda x,y:np.vstack([1-x-y,y]).T
cen    = lambda x,y:np.vstack([x+1/2,y+1/2]).T

wallpp_gen={
   'p1'  :[],
   'p2'  :[bar2],
   'pm'  :[mir_v],
   'pg'  :[gla_v],
   'cm'  :[mir_v,cen],
   'pmm' :[mir_v,mir_h],
   'pmg' :[mir_v,gla_h],
   'pgg' :[gla_h2,r2],
   'cmm' :[mir_v,r2,cen],
   'p4'  :[bar41,bar2],
   'p4m' :[mir_v,bar41,bar2],
   'p4g' :[mir_f,bar41,bar2],
   'p3'  :[r31,r32],
   'p3m1':[r31,r32,mir_d],
   'p31m':[r31,r32,mir_e],
   'p6'  :[r31,r32,r2],
   'p6m' :[m6,r31,r32,r2],
    }

def replicate(nh,nk,pp_type,params,pattern):
    lat = lat2D.Lattice2D(df_wallpp.loc[pp_type,'lat_type'],**params)
    hk,R = lat.get_lattice(nh,nk)
    R = np.hstack([R,np.zeros((R.shape[0],pattern.shape[1]-2))])
    pattern1 = np.vstack([pattern +r for r in R])
    return pattern1

def generate(pp,lat_vec,pattern):
    lat1 = np.linalg.inv(lat_vec)
    pattern[:,:2] = pattern[:,:2].dot(lat1)             #fractional coordinates
    pattern1 = generate_symmetries(pp,pattern=pattern)  #generate
    pattern1[:,:2] = pattern1[:,:2].dot(lat_vec)        #real coordinates
    return pattern1

def generate_symmetries(pp,pattern) :
    '''Apply symmetries in fractional coordinates '''
    for sym in wallpp_gen[pp] :
        x,y = np.array(pattern[:,:2]).T
        new_pos = sym(x,y)
        new_pos[new_pos<0]+=1
        new_pos[new_pos>1]-=1
        # x,y = new_pos.T
        pattern = np.vstack([pattern,np.hstack([new_pos,pattern[:,2:]])])
    return pattern
    
def apply_symmetries(pp,pattern) :
    '''Apply symmetries in fractional coordinates '''
    x,y,f = np.array(pattern).T
    for sym in wallpp_sym[pp] :
        new_pos = sym(x,y)
        # new_pos[new_pos<0]+=1
        # new_pos[new_pos>1]-=1
        # x,y = new_pos.T

        pattern = np.vstack([pattern,np.hstack([new_pos,f[:,None]])])
    return pattern


#########################################################################
#### Atomic parameters
#########################################################################
# V0 = 1.0
# Ai = np.array([0.1,0.25,0.26,0.27,1.5])
# fv = lambda X,X0,Za : V0*np.exp(-(np.linalg.norm(X-X0,axis=1)/Ai[Za])**2)
# #f  = np.sum(np.array([fc(X,X0) for X0 in Xasym]),axis=0)
# g,marker=0.5,'o' #^'
# atoms = pd.DataFrame.from_dict(                                         # covalent radius in A
#     { 0:['H',(0,0,0),0.37], 1:['C',(g,g,g),0.77],
#       2:['N',(0,0,1),0.75], 3:['O',(1,0,0),0.73],4:['S',(0,1,1),1.02]},
#     orient='index',columns=['atom','color','size'])


#########################################################################
#### tests
#########################################################################
def F_pattern(npts=21,lw=0.05,x0=0.5,y0=0.5):
    x,y = np.meshgrid(np.linspace(0,1,npts),np.linspace(0,1,npts))
    F = np.zeros((npts,npts))
    F[(y>0.15) & (y<0.75) & (abs(x-0.3)<=lw)]=1
    F[(abs(y-0.75)<=lw) & (x>0.3) & (x<0.7)]=1
    F[(abs(y-0.5)<=lw) & (x>0.3) & (x<0.55)]=1
    x,y = np.meshgrid(np.linspace(0,x0,npts),np.linspace(0,y0,npts))
    pattern = np.vstack([x.flatten(),y.flatten(),F.flatten()]).T
    return pattern

if __name__=='__main__':
    from utils import*
    from . import lattice
    lat_vec = lattice.get_lattice_vec('oblique',2,3,70)
    pattern = F_pattern(npts=51,x0=1,y0=1)
    xy = pattern[:,:2]
    x = np.inv(lat_vec).dot(xy)
    dsp.stddisp(scat=[x,y,F],xylims=[0,1,0,1])#np.flipud(F.T)])

    # pattern=[[0,0, 0],[0.55,0.55, 1]]
    # pattern = apply_symmetries('cm',pattern)
    # x,y,F = pattern.T
    # dsp.stddisp(scat=[x,y,F],xylims=[0,1,0,1])#np.flipud(F.T)])
