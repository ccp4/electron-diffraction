from utils import* #.displayStandards import stddisp
import time
import multislice as mupy
from structure_factor import get_pendulossung
from postprocess import load_multi_obj as loadm

#########################################################################
#### def : runs
#########################################################################
def Si100_pendullosung_fit(filename='dat/Silicon/Si100/Si100_autoslic.pkl'):
    multi = loadm(filenames)
    hk,t,re,im,I = multi.get_beams(bOpt='f')
    I22 = I[1]/I[1].max()
    xi = get_pendulossung(name='Si',miller=[0,2,2],keV=200,opt='p'); #xi = 320
    Ikin = (pi*t/xi)**2
    Idyn = np.sin(pi*t/xi)**2
    plts=[[t,I22,'g','$I_{022}^{(autoslic)}$']]
    plts+=[[t,Ikin,'b--','$kin$'],[t,Idyn,'r--','$dyn-2$']]
    stddisp(plts,labs=['$thickness(\AA)$','$I$'],
        xylims=[0,350,0,1.2],lw=2,legLoc='upper right',axPos=1,setPos=1,
        opt=pOpt,name='docs_fig/Si100_I22.svg')#,title='$Si_{[100]}$')

def Si100_autoslic(pOpt='p',ppopt='',nz=70,**kwargs):
    name = 'dat/Silicon/Si100'
    multi = mupy.Multislice(name,mulslice=False,tail='%d' %nz,
        NxNy=512,repeat=[6,6,nz],slice_thick=1.375,Nhk=5,
        **kwargs)
    if ppopt:
        if not multi.check_simu_state() == 'done' : multi.p.wait()
        if 'b' in ppopt :
            multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
                name='docs_fig/Si100_Ihk.svg',cm='Greens')#,xylims=[])
        if 'p' in ppopt:
            name='docs_fig/Si100_autoslic_pattern_%dA.png' %(round(multi.thickness))
            multi.pattern(opt=pOpt,cmap='gray',imOpt='ch',
                name=name)
    return multi

def Si110_autoslic(opt='',fopt='w',pOpt=''):
    name = 'dat/Silicon/Si110'
    multi = mupy.Multislice(name,mulslice=False,
        slice_thick=1.375,NxNy=512,repeat=[10,10,120],Nhk=5,
        v=1,opt=opt,fopt=fopt)
    if pOpt:
        time.sleep(0.1)
        if not multi.check_simu_state(0) == 'done' : multi.p.wait()
        # multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
        #     name='docs_fig/Si110_Ihk.svg',cm='jet')#,xylims=[])
        multi.pattern(Iopt='Icnsl',cmap='gray',imOpt='ch',xylims=[0,12,0,12],
            name='docs_fig/Si110_autoslic_pattern.png',axPos='T',setPos=1,opt=pOpt)
    return multi

def Si110_mulslice(opt='',pOpt=''):
    name = 'dat/Silicon/Si110'
    multi = mupy.Multislice(name,mulslice=True,
        NxNy=512,repeat=[1,1,120],Nhk=5,
        v=1,opt=opt)
    if pOpt:
        multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
            name='docs_fig/Si110_Ihk_mulslice.svg',cm='jet')#,xylims=[])
    return multi


########################################################################
if __name__ == "__main__":
    plt.close('all')
    #multi = load_multi_obj('dat/Silicon/Si110/Si110_mulslice.pkl')
    #multi=Si100_autoslic(opt='r',ppopt='p',nz=60,fopt='f',v=1,pOpt='ps')
    #Si100_pendullosung_fit(filename='dat/Silicon/Si100/Si100_autoslic.pkl')
    #multi=Si110_mulslice(opt='r',pOpt='s')
    #multi=Si110_autoslic(opt='rs',pOpt='ps',fopt='')
