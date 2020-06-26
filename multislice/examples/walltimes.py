from utils import*
import multislice as mupy
#from postprocess import load_multi_obj as loadm
import postprocess as pp
from postprocess import hosts,cpus,ncpus,info_cols


def run_wallT(name,Nmax=7,Nslices=5,**kwargs):
    df_name = name+'/df.pkl'
    param,NxNys ='NxNy', 2**np.arange(5,Nmax)
    base_cols = [param,'host','state']
    #prepare dataframe and run
    df = pd.DataFrame(columns=base_cols+info_cols)
    for host in hosts:
        print(red+"HOST : " +host+black)
        mupy.sweep_var(name,param,NxNys,df,tail=host,mulslice=False,keV=200,
            repeat=[1,1,Nslices],slice_thick=1.375,Nhk=5,
            ssh=host,fopt='',**kwargs)
        df.to_pickle(df_name);print(green+'DataFrame saved : ' +yellow+df_name+black)
    return df

def postprocess_wallT(name,Nslices=1,opts='u',**kwargs):
    '''opts:
        - u(update), W (walltime else cputime),
        - c (per core), s(per slice) l(log_opt)'''
    df_name = name+'/df.pkl'
    df  = pd.read_pickle(df_name)
    lab = ['cputime','walltime']['W' in opts]
    tle,col,logopt = lab,lab +'(s)',''
    if 'u' in opts : df = pp.update_df_info(df_name)
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
        plots += [[NxNys**2,Ts,[c_i,'s-'],cpus[host]]]
    stddisp(plots,lw=2,name='docs_fig/autoslic_%s.svg' %(lab),logOpt=logopt,
        labs=['$N_{beams}$','time (s)'],
        title=tle,setPos=1,axPos='T',**kwargs)


########################################################################
if __name__ == "__main__":
    plt.close('all')
    name,Nxy_max,Nslices = 'dat/walltest',7,5
    #run_wallT(name,Nmax=Nxy_max,Nslices=Nslices,opt='dsr',v='nrR')
    postprocess_wallT(name,opts='ucls',Nslices=Nslices,opt='p')
