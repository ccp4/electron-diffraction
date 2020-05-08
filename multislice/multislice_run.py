from utils import*
import time
import importlib as imp
import multislice as mupy
from structure_factor import get_pendulossung
from postprocess import load_multi_obj as loadm

def sweep_var(name,param,vals,tail='',**kwargs):
    '''
    runs a set of similar simulations with one varying parameter
    - name : path to the simulation folder
    - param,vals : the parameters and values to sweep
    - kwargs : see help(Multislice)
    '''
    df,ps = pd.DataFrame(columns=[param,'zMax(A)','normI','cputime(s)','walltime(s)']),[]
    for i,val in zip(range(np.array(vals).size),vals):
        kwargs[param]=val
        multi=mupy.Multislice(name,tail=tail+param+str(i),**kwargs)
        if not multi.check_simu_state(0) == 'done' and multi.p:
            multi.p.wait()
        #get info from log
        with open(multi.outf['log'],'r') as f : log_lines=f.readlines()
        log_lines = log_lines[-6:-4] + log_lines[-2:]
        info = [l.strip().split('=')[-1].split(' ')[1] for l in log_lines]
        df.loc[multi.name] = np.array([val]+info,dtype=float)
    df.to_csv(name+tail+'_'+param+'.pkl')
    return df

def stdmult(name,opt='',**kwargs):
    '''creates,saves and starts a simulation :
    - opt : 's'(save object), 'r'(run)
    '''
    multi=mupy.Multislice(name,**kwargs)
    if 's' in opt : multi.save()
    p = multi.run(v=1,fopt=kwargs['fopt']) if 'r' in opt else None
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
        multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
            name='docs_fig/Si100_Ihk.svg',cm='Greens')#,xylims=[])
        #### Plot beam022 and compare with 2-beam theory
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
    return multi

def Si110_autoslic(opt='',fopt='w',pOpt=''):
    name = 'dat/Silicon/Si110'
    multi,p = stdmult(name,mulslice=False,
        slice_thick=1.375,NxNy=512,repeat=[10,10,120],Nhk=5,
        v=1,opt=opt,fopt=fopt)
    if pOpt:
        time.sleep(0.1)
        if not multi.check_simu_state(0) == 'done' : p.wait()
        # multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
        #     name='docs_fig/Si110_Ihk.svg',cm='jet')#,xylims=[])
        multi.pattern(Iopt='Icnsl',cmap='gray',imOpt='ch',xylims=[0,12,0,12],
            name='docs_fig/Si110_autoslic_pattern.png',axPos='T',setPos=1,opt=pOpt)
    return multi

def Si110_mulslice(opt='',pOpt=''):
    name = 'dat/Silicon/Si110'
    multi,p = stdmult(name,mulslice=True,
        NxNy=512,repeat=[1,1,120],Nhk=5,v=1,opt=opt)
    if pOpt:
        multi.beam_vs_thickness(bOpt='f',opt=pOpt,setPos=1,axPos=1,
            name='docs_fig/Si110_Ihk_mulslice.svg',cm='jet')#,xylims=[])
    return multi

def run_wallT(pOpt='',Nmax=16):
    NxNys=2**np.arange(5,Nmax)
    df = sweep_var('dat/Silicon/Si100','NxNy',NxNys,tail='wallT_',
        mulslice=False,keV=200,repeat=[1,1,5],slice_thick=1.375,Nhk=5,
        opt='dr',v='nr',fopt='')
    Ts=df['cputime(s)'].values
    plots = [np.log2(NxNys**2),np.log10(Ts),'bs-','autoslic']
    stddisp(plots,lw=2,
        labs=['$log_2(N_{beams})$','$log_{10}$(cpuTime) (sec)'],setPos=1,axPos=1)
    return df

########################################################################
if __name__ == "__main__":
    plt.close('all')
    #multi = load_multi_obj('dat/Silicon/Si110/Si110_mulslice.pkl')
    #multi=Si100_autoslic(opt='r',pOpt='s')
    #multi=Si110_mulslice(opt='r',pOpt='s')
    multi=Si110_autoslic(opt='rs',pOpt='ps',fopt='')
    #df=run_wallT(pOpt='',Nmax=13)
