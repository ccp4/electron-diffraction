import pandas as pd
import numpy as np
import multislice as mupy
import rotating crystal as rcc



def sweep_var(name,param,vals,df=None,ssh='',tail='',do_prev=0,**kwargs):
    '''
    runs a set of similar simulations with one varying parameter
    - name          : path to the simulation folder
    - param,vals    : the parameters and values to sweep
    - df            :
        - pd.Dataframe to update(since parsed as a reference)
        - int create and save the new dataframe if 1
    - do_prev       : Used for iterative fourier transform
    - kwargs : see help(Multislice)
    '''
    do_df,save = isinstance(df,pd.core.frame.DataFrame),0
    if isinstance(df,int):
        if df : df,do_df,save = pd.DataFrame(columns=[param,'host','state']+pp.info_cols),1,1
    nvals,prev = len(vals),None
    for i,val in zip(range(nvals),vals):
        kwargs[param]=val
        if do_prev and i: prev = multi.outf['image']
        multi=Multislice(name,prev=prev,
            ssh=ssh,tail=tail+param+str(i).zfill(ceil(nvals/10)),
            **kwargs)
        if do_df:
            df.loc[multi.outf['obj']] = [nan]*len(df.columns)
            df.loc[multi.outf['obj']][[param,'host','state']] = [val,ssh,'start']
    if save :
        df.to_pickle(name+'df.pkl')
        print(green+'Dataframe saved : '+yellow+name+'df.pkl'+black)
