from utils import*                      ;imp.reload(dsp)
import matplotlib.patches as patches,os
from wallpp import wallpaper as wallpp  ;imp.reload(wallpp)
from wallpp import lattice as lat       ;imp.reload(lat)
from wallpp import config as cg         ;imp.reload(cg)
plt.close('all')
figpath = 'figures/'
datpath = 'dat/'

def test_xygrid(new=0,pp_type='p4'):
    if new:
        pattern,lat_vec = np.load(datpath+'cygne_test2_100_%s.npy'%pp_type,allow_pickle=True)
        pattern=np.array(pattern)
        params = lat.params_from_lat_vec(lat_vec)
        pp = wallpp.Wallpaper(pp,path=datpath,pattern=pattern,**params)
        pp.save()

    pp = wallpp.load(file=os.path.join(datpath,pp_type+'_wallpp.pkl'))
    # pp.show_ref(ndeg=1000)
    # x,y = np.linspace(1.25,5,1000),np.linspace(0,3.5,1000)
    x,y = np.linspace(0,1.2,1000),np.linspace(0.5,1.4,1000);
    pp.show_xygrid(x,y,pOpt='Xet')
    # pp.plot('aug',2,3)

def test_wiki_wallpp():
    for i,pp_type in enumerate(cg.pp_types[1:]):
        pattern,lat_vec = np.load(datpath+'%s.npy' %pp_type,allow_pickle=True)
        params = lat.params_from_lat_vec(lat_vec)
        pattern=np.array(pattern).T
        pp = wallpp.Wallpaper(pp_type,pattern=pattern,**params,interp=False)
        pp.plot(opts='aug',nh=2,nk=2,title=pp_type,pOpt='tXe',
            name=figpath+'wiki_%s_%s.png' %(str(i+1).zfill(2),pp_type),opt='sc')

def test_swan_wallpp(params = {'a':1.5,'b':1,'alpha' : 80}):
    # z = np.dot(im[...,:3], [0.2989, 0.5870, 0.1140])
    # from PIL import Image
    # z = Image.open('figures/cygne100.png').convert('LA')

    im = plt.imread(figpath+'cygne_test2_100.png')
    z = np.flipud(im[:,:,2])
    nx,ny = z.shape
    x,y = np.meshgrid(np.arange(nx),np.arange(ny))
    x,y,z = x.flatten()/500-0.12,y.flatten()/500-0.15,z.flatten()
    pattern0 = np.array([x,y,z]).T

    # plt.imshow(im)
    # dsp.stddisp(scat=[x,y,z])
    # plt.show()

    wallpps=dict()
    for i,pp_type in enumerate(['pg']):#cg.pp_types):
        pattern = pattern0.copy()
        pattern[:,:2] += cg.asym_cells[pp_type].mean(axis=0)
        pp = wallpp.Wallpaper(pp_type,pattern=pattern,**params)
        name = figpath+'cygne_%s_%s.png' %(str(i).zfill(2),pp_type)
        pp.plot('uag',nh=2,nk=2,title=pp_type,pOpt='tXe',
            opt='sc',name=name)
        # wallpps[pp_type] = pp


def get_supercell(lat2D,u):
    lat_vec = lat2D.lattice_vectors
    uxy = np.array(u).dot(lat_vec)
    uxy /=np.linalg.norm(uxy)
    sign = lambda x :1 if x>=0 else -1
    nh,nk = -20*sign(uxy[1]),20*sign(uxy[0])
    hk,R = lat2D.get_lattice(nh,nk)
    idx = np.argmin(abs(R.dot(uxy)))

    d = R[idx].dot(uxy)
    m,n = hk[idx]
    R0  = R[idx]
    return R0,d,m,n

def test_mesh():
    pp = wallpp.load(file=datpath+'p3_wallpp.pkl')
    u = [3,1]
    # lat2D = lat.Lattice2D('hex',a=1)
    lat2D = pp.lat
    R,d,m,n = get_supercell(lat2D,u=u)
    # print(d,m,n)

    lat_vec=pp.lat_vec
    v=np.array([m,n]).dot(lat_vec)
    K = np.array(u).dot(lat_vec)
    O1 = -u[0]*lat_vec[0]
    O2 = O1+K
    plts = []
    plts+= [ [ [K[0],0],[K[1],0],'r','K'] ]
    plts+= [ [ [O1[0],O2[0]],[O1[1],O2[1]],'r--',''] ]
    plts+= [ [ [R[0],0],[R[1],0],'b','R'] ]

    # fig,ax = lat2D.plot_lattice(m,n,lw=2,opt='')
    # dsp.stddisp(plts,lw=2,equal=1,ax=ax)



    #mesh
    thick,Tz = 10*pp.a,np.linalg.norm(R)
    x0  = np.linspace(0,thick,1000)
    z0  = np.linspace(0,1.3*Tz,1000)

    #get im for mesh
    x,z  = np.meshgrid(x0,z0)
    xz   = np.array([x.flatten(),z.flatten()]).T
    uv   = np.array([K/np.linalg.norm(K),R/np.linalg.norm(R)])
    xz_  = xz.dot(uv)                               #oriented with respect to u
    xzf_ = xz_.dot(np.linalg.inv(lat_vec))          #fractional coords
    xzf_-= np.floor(xzf_)                           #in first unit cell


    im = np.stack([np.reshape(f(xzf_),x.shape) for f in pp.f],axis=2)
    # x1,z1 = xz_.T
    # x,z = np.reshape(x1.shape),np.reshape(x1.shape)
    # dsp.stddisp(im=[im[:,:,1]],pOpt='tXe')
    hk,Rs = pp.lat.get_lattice(np.arange(-10,11),np.arange(-10,15))
    xR,yR = (Rs.dot(np.linalg.inv(uv))).T
    scat = [xR,yR,50,'b']
    pps = [patches.Rectangle((0,0),600,Tz,alpha=0.5,color='b')]
    dsp.stddisp(scat=scat,patches=pps,im=[im],xylims=[0,x0.max(),0,z0.max()],
        im_args={'origin':'lower','extent':[0,x0.max(),0,z0.max()]})



test_wiki_wallpp()
# test_swan_wallpp()
# test_xygrid(pp_type='p3',new=0)
