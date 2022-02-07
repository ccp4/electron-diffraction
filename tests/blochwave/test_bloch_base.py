import blochwave,time,sys,pytest,os,copy
from pytest_html import extras
from utils import*                  ;imp.reload(dsp)
from utils import pytest_util       ;imp.reload(pytest_util)
from blochwave import bloch         ;imp.reload(bloch)
from blochwave import util as bloch_util;imp.reload(bloch_util)
plt.close('all')

out,ref,dir = pytest_util.get_path(__file__)
b0 = bloch.Bloch('diamond',path=out,keV=200,u=[0,0,1],Nmax=8,Smax=0.05,
    opts='svt',thick=100)

@pytest_util.cmp_ref(__file__)
def test_solve_bloch():
    g=b0.gammaj
    return g

@pytest_util.cmp_ref(__file__)
def test_beam_thickness():
    idx   = b0.get_beam()
    return b0.get_beams_vs_thickness(dict_opt=False,idx=idx)

@pytest_util.add_link(__file__)
def test_show_beam_thickness():
    return b0.show_beams_vs_thickness(thicks=(0,200,1000),
        beam_args={'cond':'(Vga>1e-4) & (Sw<1e-2)'},cm='jet',
        opt='')

# @pytest.mark.new
@pytest_util.add_link(__file__)
def test_convert2tiff():
    import tifffile
    tiff_file=out+"/I.tiff"
    v=b0.convert2tiff(tiff_file=tiff_file,
        show=False,cutoff=10,thick=200,Imax=1e7)
    im=tifffile.imread(tiff_file)
    return dsp.stddisp(im=[im],
        cmap='gray',caxis=[0,20],pOpt='t',opt='')

@pytest_util.add_link(__file__)
def test_show_beams():
    return b0.show_beams(cmap='viridis',opt='')#mag=1000,)

@pytest_util.add_link(__file__)
def test_show_Idyn_vs_Ikin():
    b0 = bloch.Bloch('diamond',path=out,keV=200,u=[2,4,1],Nmax=8,Smax=0.05,
        opts='svt',thick=100,thicks=(0,1000,1000))
    return b0.show_Idyn_vs_Ikin(iZs=[10,100,500,999],opt='')#slice(249,None,250))



@pytest.mark.lvl2
def test_set_beam():
    b0.set_beam(keV=200)
    b0.set_beam(u=[1,0,1])
    b0.set_beam(K=1/cst.keV2lam(100)*np.array([1,0,0]))

@pytest.mark.lvl2
def test_get_beam():
    b0.df_G.index = np.arange(b0.df_G.index.size)
    idx = b0.get_beam(cond={'tol':1e-2})
    idx = b0.get_beam(refl=[str((0,0,0))])
    hkl = b0.get_beam(refl=[str((0,0,0))],index=False)
    idx = b0.get_beam()

def test_beam_thickness():
    # b0.__dict__.pop('Iz')
    b0.get_beams_vs_thickness(dict_opt=True)

def test_gets():
    b0.get_Xig()
    b0.get_Sw()
    print('  0,0,0 in bloch : ' ,b0.is_hkl(0,0,0))
    print('100,0,0 in bloch : ' ,b0.is_hkl(100,0,0))

@pytest.mark.lvl2
def test_show_Fhkl():
    b0.show_Fhkl(xylims=6,h3D=1,ms=30,cmap='Greens',opt='')
    b0.show_Fhkl(s=np.s_[:,:,3],xylims=6,ms=30,cmap='Greens',opt='')
    b0.show_Fhkl(s='h=0',opts='l',xylims=7,ms=50,cmap='Greens',pOpt='im',caxis=[0,1.5],opt='')

# @pytest.mark.slow
@pytest_util.add_link(__file__)
def test_bloch_convergence():
    b=copy.copy(b0)
    refl=b.get_beam(cond={'tol':1e-10,'opt':''},index=False)#;print(idx)
    # Smax,Nmax=(0.05,0.1,0.2,0.5),(6)
    Smax,Nmax=(0.1),(6,8,10,12)
    b.convergence_test(Smax=Smax,Nmax=Nmax)#,hkl=refl)
    return b.show_convergence(hkl=refl,xlab='Nmax',opt='')

def test_load_bloch():
    b = bloch_util.load_bloch(file=out+'/diamond001_200keV_bloch.pkl')
    bloch_util.load_bloch(path='.')
    bloch_util.load_bloch(path=out)
    bloch_util.load_bloch(path=out,tag='001')

# test_load_bloch()
# print(b0.get_beam())
# test_get_beam()
# test_gets()
# test_beam_thickness()
# test_show_beams()
# test_show_Idyn_vs_Ikin()
# test_bloch_convergence()

#test_beam_thickness()
# b0 = test_bloch_solve()
# test_set_beam(b0)
# test_show_Fhkl(b0)
# test_convert2tiff()
# b = bloch_convergence(n=6)
