from utils import*
from postprocess import*
from utils import displayStandards as dsp
from utils.handler3D import handler_3d
from crystals import Crystal
from wallpp import lattice
import multislice as mupy
import sys,numpy,os
import rotating_crystal as rcc
#sys.backtrace=0
cs = {1:unicolor(0.75),6:unicolor(0.25),7:(0,0,1),8:(1,0,0)}

#########
def show_grid(name):
    pattern,crys = APAP_xyz(name,n=[0,0,1],dopt=0)
    a1,a2,a3,alpha,beta,gamma=crys.lattice_parameters
    a,b,angle=a3,a1,beta
    Z,x,y,z = pattern[:,:4].T
    C = [cs[z] for z in Z]
    fig,ax = stddisp(scat=[z,x,C],ms=50,opt='')
    lat_vec = lattice.get_lattice_vec('oblique',a=a,b=b,alpha=angle,v=0)
    lattice.plot_lattice(lat_vec,opts='',nh=3,nk=3,
        ax=ax,labs=['$z$','$x$'],lw=2,
        pOpt='tG',pad=1,xyTicks=[1.5,5.0],xylims=[-6,15,0,21])#,xyTicksm=1.0,xylims=[0,100,0,100])

def get_unit_cell(cell_mesh,n):
    alpha,lw,c = 0.1,2,'c'
    x,y,z = cell_mesh
    planes_coords,ncells = [], x.shape[0]
    for i in range(ncells):
        planes_coords+=[[x[i,:,:],y[i,:,:],z[i,:,:]],
                        [x[:,i,:],y[:,i,:],z[:,i,:]],
                        [x[:,:,i],y[:,:,i],z[:,:,i]]]
    surfs = []
    for plane_coords in planes_coords :
        x,y,z  = plane_coords
        p_shape = x.shape
        coords = np.array([x.flatten(),y.flatten(),z.flatten()]);#print(coords.shape)
        x,y,z  = rcc.orient_crystal(coords,n_u=n,T=False);
        x,y,z  = np.reshape(x,p_shape),np.reshape(y,p_shape),np.reshape(z,p_shape)
        surfs += [[x,y,z,c,alpha,lw,c]]
    return surfs

def show_cell(file,n=[0,0,1],bopt=1,x0=[0,0,0],rep=[1,1,1],**kwargs):
    crys = Crystal.from_cif(file)
    tail = ''.join(np.array(n,dtype=str))
    n    = get_vec(n,crys,bopt)

    #unit cells,lattice vectors,atom coordinates
    cell_mesh   = crys.mesh(range(rep[0]+1),range(rep[1]+1),range(rep[2]+1))
    surfs       = get_unit_cell(cell_mesh,n)
    uvw         = rcc.orient_crystal(np.array(crys.lattice_vectors),n_u=n,T=True)
    pattern,lat_params = APAP_xyz(name,n,rep,dopt='s',tail=tail)
    E,X,Y,Z     = pattern[:,:4].T
    scat        = [X,Y,Z, [cs[int(e)] for e in E]]

    scat2=[]
    # pattern2,lat_params = APAP_xyz(name,n,rep,dopt='')
    # E2,X2,Y2,Z2 = pattern2[:,:4].T
    # scat2       = [X2,Y2,Z2, [cs[int(e)] for e in E2]];#print(E2)

    #display options
    c1,c2,c3 = (X.min()+X.max())/2,(Y.min()+Y.max())/2,(Z.min()+Z.max())/2
    w = 0.75*max(X.max()-X.min(),Y.max()-Y.min(),Z.max()-Z.min())
    xylims=[c1-w/2,c1+w/2,c2-w/2,c2+w/2,c3-w/2,c3+w/2]
    fig,ax=dsp.stddisp(scat=scat,ms=100,surfs=surfs,rc='3d',std=0)
    rcc.show_trihedron(ax,uvw=uvw,x0=[0,0,0],cs=['r','g','b'],labs=['$a$','$b$','$c$'],lw=2,rc=0.1,h=0.2)
    rcc.show_trihedron(ax,x0=x0,labs=['$x$','$y$','$z$'],lw=2,rc=0.1,h=0.2)
    dsp.stddisp(ax=ax,ms=50,scat=scat2,xylims=xylims,axPos=[0,0,1,1],pOpt='eXp',**kwargs)

    hdl = handler_3d(fig,persp=False)



def APAP_xyz(name,n=[0,0,1],rep=[1,1,1],dopt='scp',lfact=1.0,tail=''):
    path    = os.path.dirname(name) + '/'
    crys    = Crystal.from_cif(path+'paracetamol.cif')
    lat_vec = np.array(crys.lattice_vectors)
    if sum(rep)>3 : crys = crys.supercell(rep[0],rep[1],rep[2])
    pattern = np.array([[a.atomic_number]+list(lfact*a.coords_cartesian)+[a.occupancy,1.0] for a in crys.atoms])
    file   = path+'APAP%s.xyz' %tail
    pattern,lat_params = mupy.make_xyz(file,pattern,lat_vec,n,fmt='%.4f',dopt=dopt)
    return pattern,lat_params #pattern,crys # file

def get_vec(n,crys,bopt) :
    if isinstance(n,int) : n = crys.lattice_vectors[n]
    if bopt : n = np.array(n).dot(np.array(crys.lattice_vectors))
    return n

##########################################################################
if __name__ == "__main__":
    plt.close('all')
    name = 'dat/Paracetamol/paracetamol.cif'
    # p0=APAP_xyz(name,n=0,dopt=1,lfact=1.0,tail='')
    # p2=APAP_xyz(name,n=2,dopt=1,lfact=1.0,tail='')
    show_cell(name,n=[1,0,0],bopt=1,rep=[2,1,2]) #z along a
    # show_cell(name,n=[0,0,1],bopt=1) #z along c
    # show_cell(name,n=[1,0,1],bopt=1,rep=[2,2,2]) #z along a+c
    # show_cell(name,n=0,x0=[3,0,-3] ,view=[ 0,90] ,name='docs_fig/APAP100_xz.png',opt='p')
    # show_cell(name,n=0,x0=[3,0,-3] ,view=[90,90] ,name='docs_fig/APAP100_xy.png',opt='p')
    # show_cell(name,n=0,x0=[3,0,-3] ,view=[30,70] ,name='docs_fig/APAP100_3D.png',opt='p')
    # show_cell(name,n=2,x0=[-1,0,-3],view=[ 0,90] ,name='docs_fig/APAP001_xz.png',opt='p')
    # show_cell(name,n=2,x0=[-2,-2,0],view=[90,90] ,name='docs_fig/APAP001_xy.png',opt='p')
    # show_cell(name,n=2,x0=[3,-5,-3],view=[25,-65],name='docs_fig/APAP001_3D.png',opt='p')
    #show_grid(name)
