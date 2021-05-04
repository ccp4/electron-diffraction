import importlib as imp
from utils import*
import multislice.multislice as mupy      ;imp.reload(mupy)
import multislice.postprocess as pp       ;imp.reload(pp)
import multislice.mupy_utils as mut       ;imp.reload(mut)
import multislice.rotating_crystal as rcc ;imp.reload(rcc)
plt.close('all')

def test_base(name,**kwargs):
    multi=mupy.Multislice(name,keV=200,
        repeat=[1,1,100],NxNy=512,slice_thick=1.3575,Nhk=6,
        **kwargs)
    return multi


def show_postprocess(tag):
    multi=pp.load(name,tag)
    multi.beam_vs_thickness(tol=1e-4)
    multi.pattern(Iopt='Isnlg',tol=1e-4,rings=[0,0.25,0.5,1],caxis=[-6.2,0],
        pOpt='pt',cmap='binary',imOpt='ch',axPos=[0.2,0.13,0.75,0.75])

def test_get_tilts():
    tests = [
    'mupy.get_tilts(tx=np.arange(0,10,0.1))',
    'mupy.get_tilts(tx=np.arange(0,10,0.1),ty=0.1)',
    'mupy.get_tilts(tx=0.1,ty=np.arange(0,10,0.1))',
    'mupy.get_tilts(tx=np.arange(0,10)*0.2,ty=np.arange(0,10)*0.1)',
    'mupy.get_tilts(tx=np.arange(0,10,0.2),ty=np.arange(0,10,0.1))',
    'mupy.get_tilts(tx=0.1,ty=0)',]
    for test in tests:
        try:
            print(colors.green+test+colors.black)
            eval(test)
        except Exception as e:
            print(colors.red, e,colors.black)


if __name__ == '__main__':
    name  = '../dat/test/'
    #
    # multi = test_base(name,mulslice=False,opt='dsrfp',ppopt='w',tag='base',i_slice=20,v=1)
    # multi = test_base(name,mulslice=False,opt='dsrp',fopt='f',ppopt='uwB',v=2,ssh='tarik-CCP4home')
    # multi = test_base(name,mulslice=False,opt='dsrp ',fopt='f',ppopt='uwBP',v=2,ssh='badb')
    # multi = test_base(name,mulslice=False,fopt='f',opt='dsr',ppopt='',ssh='tarik-CCP4home',v='nctrdDR')
    # multi = test_base(name,mulslice=True,ssh=None,fopt='f',opt='dsr',v='nctrdDR')
    # multi = test_base(name,tag='TDS',mulslice=False,TDS=True,T=300,fopt='f',ppopt='wP',opt='dsrp',v='nctrdDR')
    # show_postprocess(tag,'')

    # test_get_tilts()
    # rock = mupy.Rocking(name,tx=np.arange(3)*0.05,ty=0,tag='tx',
    #     NxNy=256,repeat=[1,1,10],Nhk=3,
    #     opt = 'sr')
    # rock = pp.rock_load(name,'tx')
    # rock.update(v=1);
    # rock.plot_rocking(iBs=[(1,1),(0,1)],iZs=None,zs=[5,15,38])

    # file = '../../tests/multislice/dat/alpha_glycine.cif'
    # xyz = 'dat/alpha_gly.xyz'
    # file = 'Si'
    # xyz = name+'Si111large.xyz'
    # crys = mut.import_crys(file)
    # lat_params = [100,100,100]
    # mut.gen_xyz2(file,xyz,lat_params,n=[1,1,1],theta=0,pad=1,opts='p')
    # mut.gen_xyz('Si',n=[1,1,1],theta=0,rep=[20,20,100],pad=[2,2,0],xyz=xyz)
    # mut.show_grid(xyz,opts=['xy','xz'],popts='hv',figs='f')#,xylims=[0,1000,0,1000])
    # mut.show_grid(xyz,opts='yz')
    # mut.show_grid(xyz,opts='xz')#,xylims=[])
    # mut.show_grid3(xyz)

    multi = pp.load(name,tag='base',v=2)
    # multi.resume(opt='srfp',ppopt='w',Nz=30,i_slice=20)#Nhk=4)
    # multi.save_pattern(iz=15)
    multi.pattern(Iopt='ncsrg',Nmax=256,gs=0.3,rmax=20,tol=1e-3,xylims=2,v=1,caxis=[0,0.1])
    # qx,qy,I = multi.pattern(Iopt='cns',gs=2,rmax=20,out=1)
    # ss = [np.s_[i+int(qx.shape[0]/2),:] for i in range(4)]
    # cs = dsp.getCs('Blues',4)
    # dsp.stddisp([[qx[s],I[s],cs[i],'%d' %i] for i,s in enumerate(ss)] )
    # multi.get_beams(bOpt='fa');
    # multi.merged=0
    # multi.merge_beams()
    # multi.beam_vs_thickness(orig=1,tol=1e-3)
    # multi.pattern(iz=15,Iopt='sNc',Nmax=100,cmap='viridis',caxis=[0,0.25],xylims=5)
    # multi.patterns2gif(Iopt='scN',Nmax=100,xylims=2,caxis=[0,0.25])#,pOpt='X')
    # multi.show_patterns(Iopt='s',caxis=[0,1e10],xylims=5)
