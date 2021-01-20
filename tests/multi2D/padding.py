import importlib as imp
from utils import*                  ;imp.reload(dsp)
import multislice.multi_2D as ms    ;imp.reload(ms)
import wallpp.plane_group as pg     ;imp.reload(pg)
plt.close('all')
figpath='../../multislice/docs_fig/multi2D/'

ax,bz = 2,2
pattern = np.array([[1,1,1]])
ndeg = 2**11
nzx = 32

# def compare_patterns(mp1,mp2):
#     mp1.

def apply_padding(p1,npad,nx,name='pad',**kwargs):
    x0,z0 = np.linspace(-npad*ax,(nx+npad)*ax,ndeg),np.linspace(0,bz,nzx)
    x,z = np.meshgrid(x0,z0)
    Xv = np.stack([x.flatten(),z.flatten()]).T
    f = np.sum(np.array([pg.fv(Xv,X0,int(Za)) for X0,Za in zip(p1.Xa,p1.fa)]),axis=0)
    f = np.reshape(f,(nzx,ndeg))
    pattern = [x0,z0,f]

    mp0 = ms.Multi2D(pattern,ax,bz,keV=200,
            Nx=1,dz=bz,nz=0,ppopt='',
            iZs=1,iZv=10,eps=1)
    # mp0.V_show(opt='p')
    mp0.Tz_show(iSz=slice(0,None,1),Vopt='VT',lw=2,name=figpath+name+'Tz.svg',**kwargs)
    mp0.propagate(201)
    # mp0.Xxz_show()#iZs=1,iXs=2)
    # mp0.Bz_show(iBs=np.arange(1,20,4),tol=1e-2,cmap='jet',lw=2,name=figpath+name+'Bz.svg',**kwargs)
    mp0.Qz_show(iZs=50,lw=2,opts='Nl',name=figpath+name+'Qz.svg',xylims=['x',-3,3],**kwargs)
    return mp0

nx=1
# p1 = pg.Wallpaper(pp_type='p1',a=ax,b=bz,angle=90,pattern=pattern,ndeg=0,nh=nx,nk=1)
# apply_padding(p1,npad=0,nx=nx,name='pad0_nx1',opt='p')
# apply_padding(p1,npad=1,nx=nx,name='pad1_nx1',opt='p')
# apply_padding(p1,npad=2,nx=nx,name='pad2_nx1',opt='p')
# apply_padding(p1,npad=5,nx=nx,name='pad3_nx1',opt='p')

# nx=2
# p1 = pg.Wallpaper(pp_type='p1',a=ax,b=bz,angle=90,pattern=pattern,ndeg=0,nh=2,nk=1)
# apply_padding(p1,npad=0 ,nx=nx,name='pad0_nx2',opt='sp')
# apply_padding(p1,npad=5 ,nx=nx,name='pad1_nx2',opt='sp')
# apply_padding(p1,npad=10,nx=nx,name='pad2_nx2',opt='sp')
# apply_padding(p1,npad=20,nx=nx,name='pad3_nx2',opt='sp')

npad=20
# p1a = pg.Wallpaper(pp_type='p1',a=ax,b=bz,angle=90,pattern=pattern,ndeg=0,nh=1,nk=1)
# p1b = pg.Wallpaper(pp_type='p1',a=ax,b=bz,angle=90,pattern=pattern,ndeg=0,nh=2,nk=1)
# p1c = pg.Wallpaper(pp_type='p1',a=ax,b=bz,angle=90,pattern=pattern,ndeg=0,nh=5,nk=1)
p1d = pg.Wallpaper(pp_type='p1',a=ax,b=bz,angle=90,pattern=pattern,ndeg=0,nh=10,nk=1)
# msa = apply_padding(p1a,npad=npad,nx=1 ,name='padnx1',opt='p')
# msb = apply_padding(p1b,npad=npad,nx=2 ,name='padnx2',opt='p')
# msc = apply_padding(p1c,npad=npad,nx=5 ,name='padnx3',opt='p')
msd = apply_padding(p1d,npad=npad,nx=10,name='padnx4',opt='p')
