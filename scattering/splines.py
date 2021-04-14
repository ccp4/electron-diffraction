from utils import*
uni = colors.unicolor
Zs = [1,6,7,8,14,16]
cs   = [uni(0.9),uni(0.5),'b','r','m','y']
legs = ['H','C','N','O','Si','S']

nRmax = 100
fspline ='data/splines.npy'

def convert():
    nZs = len(Zs)
    ss = (nZs,nRmax)
    y,b,c,d = np.zeros(ss),np.zeros(ss),np.zeros(ss),np.zeros(ss)
    for iz,iZ in enumerate(Zs):
        x,y[iz],b[iz],c[iz],d[iz] = np.loadtxt('data/spline_%d.txt' %iZ).T
    np.save(fspline,[x,y,b,c,d])
    print(colors.green +'saved : '+colors.yellow+fspline+colors.black)

def fv_spline(x0,iz):
    i = np.argmin(abs(x0-x))
    z = x0 - x[i];
    vz = y[iz][i] + ( b[iz][i] + ( c[iz][i] + d[iz][i] *z ) *z) *z;
    return vz


def show_Vz():
    r = np.linspace(0,5,1000)

    plts = [ [r,np.array([fv_spline(r0**2,iz) for r0 in r]), cs[iz],legs[iz] ] for iz,iZ in enumerate(Zs)]
    plts = plts[4]
    dsp.stddisp(plts,labs=['r','Vz'],lw=2,xylims=[0,3,0,500])


# convert()
x,y,b,c,d = np.load(fspline,allow_pickle=True)
show_Vz()
