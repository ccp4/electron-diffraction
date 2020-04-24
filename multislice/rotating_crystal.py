from utils import*
import matplotlib.patches as patches # import Patch
from mpl_toolkits.mplot3d import Axes3D
#from numpy import ndarray.flatten as flatten

r0,h=1,2
npts=10

def get_plane(n=[1,0,0],u=[0,1,0],w=1,h=1,x0=[0,0,0]):
    x,y = np.meshgrid([-w/2,w/2],[-h/2,h/2])
    u1,u2 = u,np.cross(n,u)
    Xp = x*u1[0] + y*u2[0] + x0[0]
    Yp = x*u1[1] + y*u2[1] + x0[1]
    Zp = x*u1[2] + y*u2[2] + x0[2]
    return Xp,Yp,Zp

def get_cylinder(ti,tf,r0=1,h=1,npts=10,x0=[0,0,0]):
    t,z,h2 = np.linspace(ti,tf,npts),np.array([-1,1])[:,None],h/2
    Xc = r0*np.tile(np.cos(t),[2,1])  + x0[0]
    Yc = r0*np.tile(np.sin(t),[2,1])  + x0[1]
    Zc = h2*np.tile(z,[1,npts]) + x0[2]
    #Zc[0,:]=-h/2;Zc[1,:]=h/2
    return Xc,Yc,Zc

n10,n01,n11,u = [1,0],[0,1],np.array([1,-1])/sqrt(2),[0,0,1]
Xc0,Yc0,Zc0=get_cylinder(0,pi/2 ,r0,h,npts)
Xc1,Yc1,Zc1=get_cylinder(pi,3*pi/2,r0,h,npts)
Xp10,Yp10,Zp10=get_plane(n10,u,w=2*r0,h=h)
Xp01,Yp01,Zp01=get_plane(n01,u,w=2*r0,h=h)
Xp11,Yp11,Zp11=get_plane(n11,u,w=2*r0,h=h)
fig = plt.figure();ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(Xc0, Yc0, Zc0, color='b', alpha=0.2)#,edgecolor='b',linewidth=2)#
ax.plot_surface(Xc1, Yc1, Zc1, color='b', alpha=0.2 )#,edgecolor='b',linewidth=2)#
ax.plot_surface(Xp10, Yp10, Zp10, color='b', alpha=0.2,linewidth=2,edgecolor='b'   )#
ax.plot_surface(Xp01, Yp01, Zp01, color='b', alpha=0.2,linewidth=2,edgecolor='b'   )#
ax.plot_surface(Xp11, Yp11, Zp11, color='r', alpha=0.4,linewidth=2,edgecolor='r'   )#
ax.plot([0,0],[0,0],[-h/2,h/2],color ='k',linewidth=2)
plt.axis('off')
fig.show()
fig.savefig('dat/cylinder.png',transparent=True)


##############################
def get_cube_points_scatter():
    x,y,z = np.meshgrid([-3,0,3],[-1,0,1],[-1,0,1])
    x,y,z = x.flatten(),y.flatten(),z.flatten()
    idx = ~((abs(x)+abs(y)+abs(z))==2);#remove edges
    return x[idx],y[idx],z[idx]
def test_cube():
    fig = plt.figure();ax=fig.add_subplot(111,projection='3d')
    xc,yc,zc = get_cube_points_scatter()
    f = np.random.rand(xc.size)+1.*xc
    ax.scatter(xc,yc,zc,c=f,s=40,cmap='coolwarm')#color='b', alpha=0.4)#,edgecolor='b',linewidth=2)#
    #ax.plot_surface(facecolors='interp')
    fig.show()
