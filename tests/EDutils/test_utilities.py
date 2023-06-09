from utils import*
from EDutils import utilities as ut;imp.reload(ut)
from utils import pytest_util    ;imp.reload(pytest_util)
plt.close('all')

class C:
    def __init__(self,name,path,a,b=1,c=''):
        if not isinstance(c,str):
            raise TypeError('c must be of type str')
        self.a=a
        self.b=b
        self.c=c
        self.name=name
        self.path=path
        self.save()

    def _get_pkl(self):
        return ut.get_pkl(path=self.path,name=self.name)
    def save(self):
        print(self._get_pkl())
        ut.save(self,file=self._get_pkl())
    def __str__(self):
        return 'a=%d, b=%d, c=%s' %(self.a,self.b,self.c)

def test_sweep():

    df1 = ut.sweep_var(C,path='dat/sweep',params='a',vals=np.arange(2),tag='single')
    df2 = ut.sweep_var(C,path='dat/sweep',params=['a','b'],vals=np.random.rand(3,2),tag='double')
    df3 = ut.sweep_var(C,path='dat/sweep',params=['a','c'],vals=[[1,'test1'],[2,'test2']],tag='mixed')

    print('df1[0 ] : ',ut.load_pkl(df1.iloc[0].pkl))
    print('df2[1 ] : ',ut.load_pkl(df2.iloc[1].pkl))
    print('df3[-1] : ',ut.load_pkl(df3.iloc[-1].pkl))


@pytest_util.cmp_ref(__file__)
def test_get_uvw_cont():
    return ut.get_uvw_cont(u0=[1,0,1],u1=[1,1,1] ,nframes=100        )

@pytest_util.cmp_ref(__file__)
def test_get_uvw_rock():
    return ut.get_uvw_rock(e0=[1,0,1],e1=[1,0],npts=20,deg=1         )

@pytest_util.cmp_ref(__file__)
def test_get_uvw_CBED():
    return ut.get_uvw_CBED(u0=[1,1,1],cart=False ,npts=[10,20],deg=20)

@pytest_util.add_link(__file__)
def test_show_uvw():
    uvw=ut.get_uvw_cont(u0=[1,0,1],u1=[1,1,1] ,nframes=20)
    return ut.show_uvw(uvw,opt='',view=[50,-110])


@pytest_util.add_link(__file__)
def test_misc_orient():
    u   =ut.u_from_theta_phi(30,62)
    uvw =ut.get_uvw_from_theta_phi(30,64,nframes=6)
    fig,ax = ut.uvw_add_points(uvw,npts=1,plot=1,opt='')
    return fig,ax

@pytest_util.add_link(__file__)
def test_convert2tiff():
    import tifffile
    out,ref,dir=pytest_util.get_path(__file__)
    tiff_file = out+'convert.tiff'
    n=512
    x,y=np.meshgrid(range(-n,n),range(-n,n))
    im0=np.exp(-1e-3*x**2-1e-4*y**2)
    # dsp.stddisp(im=[im0],caxis=[0,1],pOpt='im')
    Imax=1e3
    ut.convert2tiff(tiff_file,im0,n0=n,rot=20,n=0.7*n,Imax=Imax)
    im=tifffile.imread(tiff_file)
    return dsp.stddisp(im=[im],caxis=[0,Imax],pOpt='im',opt='')



def test_import_crys():
    ut.import_crys('diamond')
    ut.import_crys('pets/alpha_glycine.cif')
    ut.import_crys()
    try:
        ut.import_crys('dummy')
    except Exception as e:
        print(colors.red,e,colors.black)


def test_remove_friedel_pairs():
    hkl = [[0,1,1],[0,1,-1],[0,-1,-1]]
    refl = [str(tuple(h)) for h in hkl]
    print(ut.remove_friedel_pairs(refl))


################################################################################
#### later
################################################################################
# def show_uvw():
#     fig,ax = dsp.create_fig(rc='3d',figsize='f')
#     ut.get_uvw_cont(u0=[1,0,1],u1=[1,1,1] ,nframes=100         )#title='cont'       ,**args)
#     ut.get_uvw_rock(e0=[1,0,1],e1=[0,1],npts=20,deg=1          )#title='rock_theta' ,**args)
#     ut.get_uvw_rock(e0=[1,0,1],e1=[1,0],npts=20,deg=1          )#title='rock_phi'   ,**args)
#     ut.get_uvw_CBED(u0=[1,1,1],cart=False ,npts=[10,20],deg=20 )#title='CBED'       ,**args)
#     plts = [[[0,u[0]],[0,u[1]],[0,u[2]],'b'] for u in uvw ]
#     return fig,ax

# test_remove_friedel_pairs()
# fig,ax=test_convert2tiff();fig.show()
# test_sweep()
# fig,ax=test_show_uvw();fig.show()
# fig,ax=test_misc_orient();fig.show()
# test_import_crys()
