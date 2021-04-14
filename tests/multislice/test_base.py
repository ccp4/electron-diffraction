import importlib as imp
from utils import*
import multislice.pymultislice as pyms;imp.reload(pyms)
import multislice.multislice as ms    ;imp.reload(ms)
import multislice.mupy_utils as mut   ;imp.reload(mut)
import multislice.postprocess as pp   ;imp.reload(pp)
import multislice.multi_3D as MS3D    ;imp.reload(MS3D)
# from crystals import Crystal
plt.close('all')
figpath='../../multislice/docs_fig/multi3D/'
datpath='../../multislice/dat/Silicon/Si110/'
Si110 = datpath+'Si110.xyz'

# opts = 'pTl' #p(py)l(load) T(temsim)L(Load)
ndz=2
nxy=2**10
pattern,(ax,by,cz) = mut.import_xyz(Si110)
pattern = pattern[:,:4]
dz = cz/ndz

def run_py3D(load_opt=0,opts='stql',Nz=None,nz=None,v=1):
    mp_py = mut.multi3D(load_opt=load_opt,name=datpath+'Si110',
        pattern=pattern,ax=ax,by=by,cz=cz,
        keV=200,Nxy=[1,1],nxy=nxy,dz=dz,Nz=Nz,nz=nz,
        TDS=False,
        iZs=1,iZv=1,opts=opts,v=v,
        hk=3,hkopt='r',
        eps=1,copt=1,temsim=True,
        ppopt='')
    # mp_py.Bz_show(cm='jet')
    return mp_py


def run_temsim(load_opt=0,Nz=1):
    if load_opt:
        mp_ms = pp.load_multi_obj(datpath+'Si110_autoslic.pkl')
    else:
        mp_ms = ms.Multislice(datpath,mulslice=False,data=Si110,
            keV=200,slice_thick=dz,NxNy=nxy,repeat=[1,1,Nz],Nhk=5,
            v=1,opt='rs',ppopt='w',fopt='f',i_slice=1)
    # mp_ms.pattern(Iopt='Icns',cmap='gray',imOpt='ch',xylims=[0,12,0,12],
    #     name='docs_fig/Si110_autoslic_pattern.png',axPos='T',setPos=1)
    # mp_ms.beam_vs_thickness(bOpt='fa',tol=1e-4)
    return mp_ms

def compare_transmission(iz,plot=0,**kwargs):
    mp_py = pyms.load(datpath+'Si110_3D.pkl')
    Tpy = mp_py._get_transmission_function(iz,v=2,copt=1,save=0,load=1)
    Tms = np.loadtxt(datpath+'translayer.%s' %str(iz).zfill(3) )
    Tms = Tms[:,::2]+1J*Tms[:,1::2]
    deltaT = Tms-Tpy
    print('max diff re=%.2E, im=%.2E' %(abs(np.real(deltaT)).max(),abs(np.imag(deltaT)).max()) )
    if 'e' in plot:
        if 'r'in plot:dsp.stddisp(im=[deltaT.real],pOpt='im',title='re error(T) iz=%d' %iz,**kwargs)
        if 'i'in plot:dsp.stddisp(im=[deltaT.imag],pOpt='im',title='im error(T) iz=%d' %iz,**kwargs)
    if 't' in plot:
        if 'r'in plot:dsp.stddisp(im=[Tms.real],pOpt='im',title='TEMSIM Re T iz=%d' %iz,**kwargs)
        if 'r'in plot:dsp.stddisp(im=[Tpy.real],pOpt='im',title='MULT3D Re T iz=%d' %iz,**kwargs)
        if 'i'in plot:dsp.stddisp(im=[Tms.imag],pOpt='im',title='TEMSIM Im T iz=%d' %iz,**kwargs)
        if 'i'in plot:dsp.stddisp(im=[Tpy.imag],pOpt='im',title='MULT3D Im T iz=%d' %iz,**kwargs)
    return Tpy,Tms

def compare_propagator(opts='qri'):
    P = np.loadtxt(datpath+'propagator.txt')
    nx = int(P.shape[0]/2)
    Px,Py = P[:nx,0]+1J*P[:nx,1] , P[nx:,0]+1J*P[nx:,1]
    if opts:
        P = np.outer(Px,Py)
        x = np.arange(nx)
        plts = [[x,Px.real,'b','Re Px'],[x,Px.imag,'r','Im Px']]
        plts+= [[x,Py.real,'c','Re Py'],[x,Py.imag,'m','Im Py']]
        if 'q' in opts:dsp.stddisp(plts,labs=['q','Px,Py'])
        if 'r' in opts:dsp.stddisp(im=[P.real.T],pOpt='im',title='Re P temsim')
        if 'i' in opts:dsp.stddisp(im=[P.imag.T],pOpt='im',title='Im P temsim')
    return Px,Py

def compare_patterns(iz,plot=0,**kwargs):
    mp_py = pyms.load(datpath+'Si110_3D.pkl')
    I_py = mp_py.pattern(iz)
    mp_ms = pyms.load(datpath+'Si110_3D.pkl')
    # I_ms = mp_py.pattern(iz)
    qx,qy,I_ms = mp_ms.pattern(iz,file=datpath+'',out=1)

mp_py = run_py3D(load_opt=1,opts='stq',Nz=5,v=1)
mp_ms = run_temsim(load_opt=1,Nz=5)
# Px_ms,Py_ms = compare_propagator(opts='r')
# Px_py,Py_py = mp_py.Pq_show(opts='r')
# Tpy,Tms = compare_transmission(2,plot='tr')
# mp_py.Tz_show(0,'i')

hk = [hk0 for hk0 in mp_py.hk if hk0 in mp_ms.hk]
fig,ax = mp_py.Bz_show(           iBs=hk,cm='jet', linespec='-' ,opt='')#title='py  3D')
mp_ms.beam_vs_thickness(bOpt='fa',iBs=hk,cm='jet', linespec='--',ax=ax,legOpt=0)#title='temsim')
