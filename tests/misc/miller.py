import importlib as imp
from utils import*                          ; imp.reload(dsp)
import multislice.multislice as mupy        ; imp.reload(mupy)
import multislice.postprocess as pp         ; imp.reload(pp)
import multislice.mupy_utils as mut         ; imp.reload(mut)
path ='../multislice/dat/Carbone_SC/'
figpath = '../docs_fig/'
plt.close('all')

def run(data,tail, **kwargs):
    multi = mupy.Multislice(path,tail=tail,data=data,
        mulslice=False,keV=200,TDS=False,i_slice=10,
        NxNy=1024,slice_thick=3.0,repeat=[1,1,1],
        opt='dsrp', **kwargs)
    return multi

rotx = lambda a:np.array([[1,0,0],[0,np.cos(a),np.sin(a)],[0,-np.sin(a),np.cos(a)]])
# roty = lambda a:np.array([[np.cos(a),0,np.sin(a)],[0,1,0],[-np.sin(a),0,np.cos(a)]])
# rotz = lambda a:np.array([[0,np.cos(a),np.sin(a),0],[-np.sin(a),np.cos(a),0],[0,0,1]])
#
# lat_params  = np.array([1,2,3])
# lat_rec = 1/lat_params
# abc=np.diag([1,0.5,1/3])
# Rx,Ry = rotx(np.pi/4),rotx(np.pi/3)
# xh = [0,0.5,2/3]
# h = abc.dot(xh)/lat_rec**2
#
# Rabc = Rx.dot(Ry).dot(abc)
# Rxh = Rx.dot(Ry).dot(xh)
# Rh = Rabc.T.dot(Rxh)/lat_rec**2
# print(h,Rh)

def test_001():
    lat_params  = [1,2,1]
    lat_vec     = np.diag(lat_params)
    pattern     = np.array([[6,0.5,0.5,1, 1,1],[8,0.75,0.75,0.5, 1,1]])

    name=path+'carbone001.xyz'
    n = [0,0,1] # beam along z
    u = [1,0,0] # rotation axis along x
    rep,pad = [20,20,10],4
    mut.make_xyz(name,pattern,lat_vec,lat_params,n=n,pad=pad,rep=rep,fmt='%.4f',dopt='s')

    # mut.show_grid(name,'xy')
    # mut.show_grid(name,'xz')
    # mut.show_grid(name,'yz')
    # multi = run(data=name,tail='001',Nhk=31,hk_sym=0,fopt=1,ppopt='wBP')
    multi = pp.load_multi_obj(path+'Carbone_SC_001_autoslic.pkl')
    P = multi.show_patterns(Iopt='Incs',out=False,Nmax=256,
            caxis=[0,0.01],imOpt='cv',axPos='V',cmap='gray')
            # xylims=0.4*np.array([-1,1,-1,1]),name=figpath+'900_I.png',opt='')

    #
    #
    # hs,ks = np.meshgrid(np.arange(0,3)*rep[0]*(1+2*int(pad)),np.arange(0,3)*rep[1]*(1+2*int(pad)))
    # hs,ks = hs.flatten(),ks.flatten()
    # iBs = ['(%d,%d)' %(h,k) for h,k in zip(hs,ks)]
    # multi.beam_vs_thickness(iBs=iBs,cm='Spectral',name=figpath+'C001.svg',opt='p')
    # # hs,ks = np.meshgrid(np.arange(15),np.arange(15))
    # # hs,ks = hs.flatten(),ks.flatten()
    # # iBs = ['(%d, %d)' %(h,k) for h,k in zip(hs,ks)]
    # # multi.beam_vs_thickness(iBs=iBs,cm='Spectral',name=figpath+'C001.svg',opt='p')
    #
    #
    # fig,ax = multi.pattern(Iopt='Incs',out=False,Nmax=256,
    #     #tol=1e-3,
    #     caxis=[0,0.01],
    #     # rings=[0.25,0.5],lw=2,
    #     imOpt='cv',axPos='V',cmap='gray',opt='p')
    #     # xylims=0.4*np.array([-1,1,-1,1]),name=figpath+'900_I.png',opt='')
    #
    # qx,qy,I = multi.pattern(Iopt='Incs',out=True,Nmax=256,opt='c')
    # Nmax = 256
    # s1,s2 = np.s_[Nmax,:],np.s_[:,Nmax]
    #
    # plts=[]
    # plts+=[[np.arange(I[s1].size)-Nmax,I[s1],'r']]
    # plts+=[[np.arange(I[s2].size)-Nmax,I[s2],'b']]
    # # plts = [[qx[s1],I[s1],'r','x'],[qy[s2],I[s2],'b','y']]
    # dsp.stddisp(plts,lw=2)#,xylims=['y',0,0.001])


test_001()
