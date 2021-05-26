from utils import*                  ;imp.reload(dsp)
from blochwave import bloch         ;imp.reload(bloch)
from EDutils import viewers as vw   ;imp.reload(vw)
plt.close('all')
path = 'dat'

b0 = bloch.load(path,'Si001')

u010 = np.array([ 0.05929628, -0.00230903,  0.99823776])
u020 = np.array([ 0.13922443, -0.00461806,  0.99025009])
u030 = np.array([ 0.39685357,  0.0069271 ,  0.9178558 ])
u040 = [ 0.26448615, -0.00923606,  0.96434526]
u = None#u040#u010
# theta,phi=0,0 #23.4077, 0.9991
theta,phi=15.3460, 358.0000

rot_wv = vw.Rotate_Viewer(b0,[theta,phi,1,1],
    u=u,Smax=0.001,Nmax=4,thick=100,
    F='Sw',opts='x',cmap='Greens',mag=100,xylims=1.5)

# b0.is_hkl([0,0,0])
# I = b0.get_intensities()
# b0.get_Istrong(m=1000)
# b0.get_Sw(Smax=0.01)

# b0.set_thickness(100)
# b0.show_beams(F='I')
# b0.show_beams(F='Vg',fopts='m',opts='',pOpt='im',cmap='jet',mag=500)
# b0.show_beams(F='S' ,fopts='m',opts='',pOpt='im',cmap='jet',mag=500)

# h,k,l = b0.get_hkl().T
# u,v,w = b0.Kuvw
# print(u*h+k*v+l*w)


def test_bloch_solve():
    b0 = bloch.Bloch('Si',path=path,keV=200)#,u=[0,0,1],Nmax=3,Smax=0.02,opts='s0')
    b0.solve(u=[0,0,1],Smax=0.01,Nmax=3, opts='s0')
    # b0.solve(u=[1,1,0],Smax=0.01,Nmax=3, opts='s0')
    # b0.solve(u=[1,1,1],Smax=0.01,Nmax=3, opts='s0')

def check_multi_spot():
    ar = b0.df_G[['px','py']].values
    u,idx,idr,count = np.unique(ar,axis=0,return_index=True,return_inverse=True,return_counts=True)

    for u0 in u[count>1]:
        # ar[idx0]
        ids = np.where(np.linalg.norm(ar-u0,axis=1)<1e-3)[0]
        print(b0.df_G.iloc[ids][['h','k','l','px','py','Vg']])
    # rep = np.setdiff1d(np.arange(133),idr[idx])

def test_set_beam(b0):
    b0.set_beam(keV=200)
    b0.set_beam(u=[1,0,1])
    b0.set_beam(K=1/cst.keV2lam(100)*np.array([1,0,0]))
def test_show_Fhkl(b0):
    b0.show_Fhkl(xylims=6,h3D=1,ms=30,cmap='Greens')
    b0.show_Fhkl(s=np.s_[:,:,3],xylims=6,ms=30,cmap='Greens')
    b0.show_Fhkl(s='h=0',opts='l',xylims=7,ms=50,cmap='Greens',pOpt='im',caxis=[0,1.5])

# test_bloch_solve()
# test_set_beam(b0)
# test_show_Fhkl(b0)
