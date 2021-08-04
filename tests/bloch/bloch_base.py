from utils import*                  ;imp.reload(dsp)
from blochwave import bloch         ;imp.reload(bloch)
plt.close('all')
path = 'dat'

def zone_axis_config():
    u0s = [[0,0,1],[1,1,1],[10,2,1],[14,5,2],[3,6,8]]
    b = {}
    for u0 in u0s:
        b0 = bloch.Bloch('diamond',Smax=0.035,u0=u0,Nmax=12,solve=0,opts='',path='dat/')
        b0.show_beams(opts='VSk',title=str(u0))
        b[str(u0)]=b0

def bloch_convergence(n=3):
    Nmaxs = np.arange(4,4*n+1,4)
    refl = [[0,0,0],[2,2,0],[4,0,0],[4,4,0]]
    refl = [str(tuple(h))  for h in refl]
    nbs = len(refl)
    I,b = np.zeros((Nmaxs.size,nbs)),[]
    for i,Nmax in enumerate(Nmaxs):
        b0 = bloch.Bloch('diamond',Smax=0.035,u=[0,0,1],Nmax=Nmax,thicks=(0,1000,10),thick=1000,opts='',path='dat/',name='diamond001_%d' %Nmax)
        I[i,:]  = b0.df_G.loc[refl,'I'].values
        b += [b0]

    cs = dsp.getCs('jet',nbs)
    plts = [[Nmaxs,I[:,iB],[cs[iB],'-o'],'%s' %h] for iB,h in enumerate(refl)]
    dsp.stddisp(plts,labs=['Nmax','I'],lw=2)
    return b

b = bloch_convergence(n=6)

# b0 = bloch.Bloch('diamond',Smax=0.1,u=[10,2,5],Nmax=8,thicks=(0,500,501),path=path)
# b0.show_beams_vs_thickness(thicks=(0,1000,1000))
# b0.QQplot(iZs=slice(10,30,2),xylims=[0.2]*2)
# b0.QQplot(iZs=slice(10,None,50),xylims=[0.2]*2)
# b0.QQplot(iZs=[100],xylims=[1]*2)

# b0 = bloch.load_Bloch('dat/bloch/N4_2beams.pkl')

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
