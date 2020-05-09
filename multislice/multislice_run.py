from utils import*
import time,os
import importlib as imp
import multislice as mupy
from structure_factor import get_pendulossung
from postprocess import load_multi_obj as loadm

hosts = ['','brno','tarik-CCP4home']
cpus  = dict(zip(hosts,['asus-$i_5^{(4cores)}$','brno-$i_7^{(8cores)}$','xps-$i_7^{(8cores)}$']))
ncpus = dict(zip(hosts,[4,8,8]))
info_cols = ['zmax(A)','Inorm','cputime(s)','walltime(s)']

def sweep_var(df,name,param,vals,ssh='',tail='',**kwargs):
    '''
    runs a set of similar simulations with one varying parameter
    - df : The dataframe to update
    - name : path to the simulation folder
    - param,vals : the parameters and values to sweep
    - kwargs : see help(Multislice)
    '''
    for i,val in zip(range(np.array(vals).size),vals):
        kwargs[param]=val
        multi=mupy.Multislice(name,ssh=ssh,tail=tail+param+str(i),**kwargs)
        # if not multi.check_simu_state(0) == 'done' and multi.p:
        #     multi.p.wait()
        df.loc[multi.outf['obj']] = [nan]*len(df.columns)
        df.loc[multi.outf['obj']][[param,'host','state']] = [val,ssh,'start']
    return df

def get_info(log_file):
    '''compute zmax,I,cpuTime and wallTime'''
    with open(log_file,'r') as f : log_lines=f.readlines()
    log_lines = log_lines[-6:-4] + log_lines[-2:]
    info = [l.strip().split('=')[-1].split(' ')[1] for l in log_lines]
    return np.array(info,dtype=float)

def update_df(df_name):
    df = pd.read_pickle(df_name)
    datpath = os.path.dirname(df_name)+'/'
    for name in df.index:
        print(name)
        multi = loadm(datpath+name)
        state = multi.check_simu_state(ssh_alias=df.loc[name].host,v=0)
        df.loc[name]['state'] = state
        if state=='done':
            info = get_info(multi._outf('log'))
            df.loc[name][info_cols] = info
    df.to_pickle(df_name);
    print(green+'DataFrame saved : ' +yellow+df_name+black)
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

def run_wallT(Nmax=7,Nslices=5,**kwargs):
    name    = 'dat/walltimes'
    df_name = name+'/df.pkl'
    param,NxNys ='NxNy', 2**np.arange(5,Nmax)
    base_cols = [param,'host','state']
    #prepare dataframe and run
    df = pd.DataFrame(columns=base_cols+info_cols)
    for host in hosts:
        df = sweep_var(df,name,param,NxNys,tail=host,mulslice=False,keV=200,
            repeat=[1,1,Nslices],slice_thick=1.375,Nhk=5,
            ssh=host,fopt='',**kwargs)
        df.to_pickle(df_name);print(green+'DataFrame saved : ' +yellow+df_name+black)
    return df
def postprocess_wallT(df_name,Nslices=1,opts='u',**kwargs):
    '''opts:
        - u(update), W (walltime else cputime),
        - c (per core), s(per slice) '''
    df = pd.read_pickle(df_name)
    lab = ['cputime','walltime']['W' in opts]
    tle,col,logopt = lab,lab +'(s)',''
    if 'u' in opts : df = update_df(df_name)
    if 'c' in opts : lab,tle = lab+'_pcore' ,tle+' per core'
    if 's' in opts : lab,tle = lab+'_pslice',tle+' per slice'
    if 'l' in opts : lab,logopt =lab+'_log','xy'
    plots,cs=[],['r','g','b'] #,getCs('jet',len(hosts))
    for c_i,host in zip(cs,hosts):
        df_host = df.loc[df.host==host]
        NxNys = np.array(df_host['NxNy'].values,dtype=int)
        Ts = np.array(df_host[col].values,dtype=float)
        if 'c' in opts: Ts*=ncpus[host]
        if 's' in opts: Ts/=Nslices
        # if 'l'in opts :
        #     Ts = np.log10(Ts)
        #     NxNys = np.log2(NxNys**2)
        plots += [[NxNys**2,Ts,[c_i,'s-'],cpus[host]]]
    stddisp(plots,lw=2,name='docs_fig/autoslic_%s.svg' %(lab),logOpt=logopt,
        labs=['$N_{beams}$','time (s)'],
        title=tle,setPos=1,axPos='T',**kwargs)


########################################################################
if __name__ == "__main__":
    #plt.close('all')
    #multi = load_multi_obj('dat/Silicon/Si110/Si110_mulslice.pkl')
    #multi=Si100_autoslic(opt='r',pOpt='s')
    #multi=Si110_mulslice(opt='r',pOpt='s')
    #multi=Si110_autoslic(opt='rs',pOpt='ps',fopt='')
    #run_wallT(Nmax=14,Nslices=5,opt='dsr',v='nr')
    postprocess_wallT('dat/walltimes/df.pkl',opts='cls',Nslices=5,opt='s')
