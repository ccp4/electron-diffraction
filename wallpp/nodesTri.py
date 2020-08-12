from utils import*
import matplotlib.tri as mtri
from scipy.spatial import Delaunay

def nodesTri(nDeg) : 
    h = 1/(nDeg-1)
    # Equally-space nodes 
    nodes = np.array([[],[]])
    for j in range(nDeg) : 
        i = np.arange(nDeg-j).T
        aux = j*np.ones(i.size)
        xy = np.array([i*h,aux*h]); #print(xy)
        nodes =  np.concatenate((nodes,xy),axis=1)
    nodes = 2*nodes.T-1
    return nodes


def test():
    nDeg = 2
    nodes = nodesTri(nDeg)
    #use
    x,y = nodes.T
    x0,y0=-0.05,-0.15
    f = (x-x0)**2 + (y-y0)**2
    fig, ax = plt.subplots()
    triang = mtri.Triangulation(x, y, Delaunay(nodes).vertices)
    ax.tricontourf(triang, f,levels=20)
    ax.triplot(triang)
    ax.plot(x,y, 'o')
    fig.tight_layout()
    plt.show()
    
if '__main__'==__name__:
    test()