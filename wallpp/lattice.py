'''
Defines 2D lattices with some utilities
'''
from utils import*
from utils import displayStandards as dsp
import numpy as np

pp_types = ['p1','pg','pm','cm', 'p2','pgg','pmm','cmm','pmg',
            'p4','p4m','p4g','p3','p3m1','p31m','p6','p6m']; #print(pp_types)
pp_lattices = ['oblique']+['rect']*3+['oblique']+['rect']*4+['square']*3+['hex']*5
wallpp_lat = dict(zip(pp_types,pp_lattices))

__all__ = ['get_lattice_vec','get_miller','reciprocal_lattice_2D',
    'plot_lattice']

##### 2D
def get_lattice_vec(lat_type='',wallpp_type=None,a=1,b=1,alpha=90,v=0) :
    ''' Computes lattice vectors from :
        - a,b,angle : lattice parameters
        - lat_type  : lattice type  [hex,square,rect,oblique]'''
    pi = np.pi
    if wallpp_type : lat_type = wallpp_lat[wallpp_type]
    if lat_type == 'hex' :
        lattice_vec = np.array([[1,0],[cos(pi/3),sin(pi/3)]])*a
    elif lat_type == 'square' :
        lattice_vec = np.array([[1,0],[0,1]])*a
    elif lat_type == 'rect' :
        lattice_vec = np.array([[a,0],[0,b]])
    elif lat_type == 'oblique' :
        angle=pi/180*alpha
        lattice_vec = np.array([[a,0],[b*cos(angle),b*sin(angle)]])
        if v : lattice_vec = np.array([[a*cos(angle),a*sin(angle)],[0,b]])
    else :
        printf(red,'not a valid lattice type. choose : ''', ['hex','square','rect','oblique'],
               'or select a wallpaper type : ',pp_types ,black)
        return None
    return lattice_vec

# misc
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

def reciprocal_lattice_2D(a1,a2):
    Rot90 = np.array([[0,1],[-1,0]])
    b1 = Rot90.dot(a2)/(a1.dot(Rot90.dot(a2)))#*2*np.pi
    b2 = Rot90.dot(a1)/(a2.dot(Rot90.dot(a1)))#*2*np.pi
    return b1,b2


# Display
def plot_lattice(lattice_vec,nh=3,nk=3,hOpt=True,**kwargs) :
    '''Plot lattice grid :
    - nh,nk : number of cells along each direction
    - hOpt : additional lines for hexagonal
    '''
    gridcolor=(0.75,0.75,0.75)
    a,b=lattice_vec
    n_h,n_k = np.arange(nh+1),np.arange(nk+1)
    line_h = [h*a + np.array([[0,0],nk*b]) for h in range(nh+1)]
    line_k = [k*b + np.array([[0,0],nh*a]) for k in range(nk+1)]
    plts = [[xy[:,0],xy[:,1],gridcolor,''] for xy in line_h+line_k]
    if hOpt :
        P = (a+b)
        n,m,nn,nm = a,b,nh,nk
        if nk>nh : n,m,nn,nm = b,a,nk,nh
        line_hex  = [h*n + np.array([[0,0],(nn-h)*P]) for h in range(nm,nn)]
        line_hex += [k*m + np.array([[0,0],(nm-k)*P]) for k in range(nm)]
        if not nk==nh :
            line_hex += [h*n + np.array([[0,0],nm*P]) for h in range(nn-nm)]
        plts += [[xy[:,0],xy[:,1],[gridcolor,'--'],''] for xy in line_hex]
    fig,ax=dsp.stddisp(plts,**kwargs)
                   # xylims=xylims,mg=0,changeXYlims='x' in pOpt,
                   # gridOn='g' in pOpt,ticksOn='t' in pOpt,equal='e' in pOpt,
                   # setPos='s' in pOpt,)
    return fig,ax

##########################################################################
#### test
##########################################################################
def _test_plot_all() :
    for lat in ['hex','square','rect','oblique'] :
        plot_lattice(get_lattice_vec(lat,a=1,b=3,alpha=50),nh=6,nk=3,pOpt='et')
def _test_plot_hex() :
    vec = get_lattice_vec('hex',a=1)
    plot_lattice(vec,nh=6,nk=6,xylims=[2,5,0,3],pOpt='xehp')
def _test_miller():
    vec = get_lattice_vec('oblique',a=1,b=1.5,alpha=20)
    nh,nk=1,0#3,4
    fig,ax=plot_lattice(vec,nh,nk,xylims=[0,2,0,3],pOpt='t')
    u_h,miller = get_miller(vec,nh,nk,ax);print(u_h,miller)

if __name__=='__main__' :
    #_test_plot_all()
    #_test_plot_hex()
    #_test_miller()
    print(green+'success'+black)
