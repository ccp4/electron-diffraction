from utils import*
import importlib as imp
import multislice as mupy
from structure_factor import get_pendulossung


def sweep_var(name,param,vals,tail='',**kwargs):
    '''
    runs a set of similar simulations with one varying parameter
    - name : path to the simulation folder
    - param,vals : the parameters and values to sweep
    - kwargs : see help(Multislice)
    '''
    df = pd.DataFrame(columns=[param,state,'cputime','walltime'])
    p = []
    for i,val in (range(len(vals)),vals):
        kwargs[param]=val
        multi=mupy.Multislice(name,tail=tail+str(i),**kwargs)
        p.append( multi.run() )
        df.loc[multi.name] = [val,'init','','']
    return df,p

def stdmult(name,opt='',**kwargs):
    '''creates,saves and starts a simulation :
    - opt : 's'(save object), 'r'(run)
    '''
    multi=mupy.Multislice(name,**kwargs)
    if 's' in opt : multi.save()
    p = multi.run() if 'r' in opt else None
    return multi,p


#########################################################################
#### def : runs
#########################################################################
def Si100_autoslic(opt='',pOpt=''):
    name = 'dat/Silicon/Si100'
    multi,p = stdmult(name,mulslice=False,opt=opt,
        NxNy=512,repeat=[6,6,70],slice_thick=1.375,Nhk=5)
    if pOpt:
        if not multi.check_simu_state() == 'done' : p.wait()
        # filename = 'dat/Silicon/Si100/Si100_autoslic.pkl'
        # multi = load_multi_obj(filename)
        hk,t,re,im,I = multi.get_beams(bOpt='f')
        I22 = I[1]/I[1].max()
        xi = get_pendulossung(name='Si',miller=[0,2,2],keV=200,opt='p')
        xi=320
        Ikin = (pi*t/xi)**2
        Idyn = np.sin(pi*t/xi)**2
        #plot
        plts=[[t,I22,'g','$I_{022}^{(autoslic)}$']]
        plts+=[[t,Ikin,'b--','$kin$'],[t,Idyn,'r--','$dyn-2$']]
        stddisp(plts,labs=['$thickness(\AA)$','$I$'],
            xylims=[0,350,0,1.2],lw=2,legLoc='upper right',axPos=1,setPos=1,
            opt=pOpt,name='docs_fig/Si100_I22.svg')#,title='$Si_{[100]}$')
        multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
            name='docs_fig/Si100_Ihk.svg',cm='Greens')#,xylims=[])
    return multi

def Si110_autoslic(opt='',pOpt=''):
    name = 'dat/Silicon/Si110'
    multi,p = stdmult(name,mulslice=False,
        slice_thick=1.375,NxNy=512,repeat=[10,10,120],Nhk=5,
        v=1,opt=opt)
    if pOpt:
        if not multi.check_simu_state() == 'done' : p.wait()
        multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
            name='docs_fig/Si110_Ihk.svg',cm='jet')#,xylims=[])
    return multi

def Si110_mulslice(opt='',pOpt=''):
    name = 'dat/Silicon/Si110'
    multi,p = stdmult(name,mulslice=True,
        NxNy=512,repeat=[1,1,120],Nhk=5,v=1,opt=opt)
    if pOpt:
        multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
            name='docs_fig/Si110_Ihk_mulslice.svg',cm='jet')#,xylims=[])
    return multi

########################################################################
if __name__ == "__main__":
    plt.close('all')
    #multi = load_multi_obj('dat/Silicon/Si110/Si110_mulslice.pkl')
    #multi=Si100_autoslic(opt='',pOpt='s')
    #multi=Si110_mulslice(opt='',pOpt='s')
    #multi=Si110_autoslic(opt='',pOpt='s')
