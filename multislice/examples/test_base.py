from utils import*
import multislice.multislice as mupy
import multislice.postprocess as pp
import importlib as imp
imp.reload(mupy)
plt.close('all')

def test_base(name,**kwargs):
    multi=mupy.Multislice(name,keV=200,
        repeat=[2,2,10],NxNy=512,slice_thick=1.3575,Nhk=3,
        **kwargs)
    return multi

if __name__ == '__main__':
    name  = '../dat/test/'
    #
    # multi = test_base(name,mulslice=False,opt='dsrp',fopt='f',ppopt='uwB',v=2)
    # multi = test_base(name,mulslice=False,opt='dsrp',fopt='f',ppopt='uwB',v=2,ssh='tarik-CCP4home')
    multi = test_base(name,mulslice=False,opt='dsrp ',fopt='f',ppopt='uwB',v=2,ssh='badb',hostpath='/data3/lii26466/multislice/test/',cluster=1)
    #multi = test_base(name,mulslice=False,fopt='f',opt='dsr',ppopt='',ssh='tarik-CCP4home',v='nctrdDR')
    #multi = test_base(name,mulslice=True,ssh=None,fopt='f',opt='dsr',v='nctrdDR')
    # multi = test_base(name,tail='TDS',mulslice=False,TDS=True,T=300,fopt='f',ppopt='wP',opt='dsrp',v='nctrdDR')

    # multi=pp.load_multi_obj(name+'test_autoslic.pkl')
    # multi.pattern(Iopt='Isnlg',tol=1e-4,rings=[0,0.25,0.5,1],caxis=[-6.2,0],
    #     pOpt='pt',cmap='binary',imOpt='ch',axPos=[0.2,0.13,0.75,0.75])
