from utils import*
from utils import displayStandards as dsp
import multislice as mupy
from crystals import Crystal
from rotating_crystal import orient_crystal,show_unit_cell
from postprocess import*

def Li_xyz(path,n=[0,0,1],dopt=1,lfact=1.0,wobble=0,tail=''):
    file    = path+'Li%s%s.xyz' %(''.join(np.array(n,dtype=str)),tail)
    crys    = Crystal.from_database('Li')
    pattern = np.array([
        [a.atomic_number]+list(lfact*a.coords_cartesian)+[a.occupancy,wobble] for a in crys.atoms])
    # lat_params = tuple(lfact*np.array(crys.lattice_parameters[:3]))
    lat_params = lfact*np.array(crys.lattice_vectors)
    coords,lat = mupy.make_xyz(file,pattern,lat_params,n,fmt='%.4f')
    return file

#########################################################################
def Li_latsize(path,nvals=2,**kwargs):
    df = pd.DataFrame(columns=['lat_cst','host','state']+info_cols)
    vals = np.array(2**np.arange(nvals),dtype=int) # np.linspace()
    plts0,plts1,plts2,cs1,Imax0,Imax1,Imax2 = [],[],[],getCs('viridis',nvals),0,0,0
    for i,val in zip(range(nvals),vals):
        tail = str(i).zfill(ceil(nvals/10))
        data = Li_xyz(path,n=[0,0,1],dopt=1,lfact=val,tail='_'+tail)
        multi=mupy.Multislice(path,tail=tail,data=data,
            mulslice=0,keV=200,
            NxNy=64*val,slice_thick=val*1.5,Nhk=5,repeat=[1,1,int(2**10/val)],
            **kwargs)
        df.loc[multi.outf['obj']] = [nan]*len(df.columns)
        df.loc[multi.outf['obj']]['lat_cst','host','state'] = [val,'','start']

        #hk,t,re,im,I = multi.get_beams(bOpt='f',tol=1e-4)
        hk,t,re,im,I = multi.get_beams(bOpt='f',iBs=['(1,1)','(2,0)','(4,0)'])
        Imax0,Imax1,Imax2 = max(Imax0,I[0].max()),max(Imax1,I[1].max()),max(Imax2,I[2].max())
        plts0 += [[t,I[0],[cs1[i],'-'],'$I_{%s}^{(%d)}$' %(hk[0],i)]]
        plts1 += [[t,I[1],[cs1[i],'-'],'$I_{%s}^{(%d)}$' %(hk[1],i)]]
        plts2 += [[t,I[2],[cs1[i],'-'],'$I_{%s}^{(%d)}$' %(hk[2],i)]]
    name,opt='docs_fig/lattice_effect','s'
    stddisp(plts0,labs=['$thickness$ ($A$)','$I_{beams}$'],
        xylims=[0,t.max(),0,Imax0],opt=opt,lw=2,axPos=1,setPos=1,name=name+'0.svg')
    stddisp(plts1,labs=['$thickness$ ($A$)','$I_{beams}$'],
        xylims=[0,t.max(),0,Imax1],opt=opt,lw=2,axPos=1,setPos=1,name=name+'1.svg')
    stddisp(plts2,labs=['$thickness$ ($A$)','$I_{beams}$'],
        xylims=[0,t.max(),0,Imax2],opt=opt,lw=2,axPos=1,setPos=1,name=name+'2.svg')
    df.to_pickle(path+'df.pkl')
    df=update_df_info(path+'df.pkl')
    return df


#########################################################################
def Li_patterns(name,nzs=20,**kwargs):
    param = 'repeat'
    reps  = np.arange(1,nzs+1)*5
    vals  = [[10,10,nz] for nz in reps]
    mupy.sweep_var(name,param,vals,df=1,do_prev=0,
        data='Li001.xyz',mulslice=0,keV=200,tail='',
        NxNy=512,slice_thick=1.5,Nhk=0,
        **kwargs)
def Li_gif(name):
    df=update_df_info(name+'df.pkl')
    nzs = df.index.size
    xt,yt = 0,4 #text position
    zI = df[['zmax(A)','Inorm']].values;#print(zI)
    for i in range(nzs):
        istr = str(i).zfill(ceil(nzs/10))
        tle = '$z=%.0f A, I_{norm}=%.4f$' %tuple(zI[i,:])
        multi=load_multi_obj(name+'Lithium_repeat%s_autoslic.pkl' %istr)
        multi.pattern(Iopt='Insc',tol=1e-3,Nmax=64,
            cmap='gray',texts=[xt,yt,tle,'g'],xylims=[-5,5,-5,5],
            axPos=[0,0,1,1],setPos=1,gridOn=0,ticksOn=False,
            name=name+'pattern%s.png' %istr,opt='s')
    #dsp.im2gif(name+'pattern','svg')

def Li_wobble(name,nvals=np.inf,**kwargs):
    wobbles = [0.001,0.1,0.3,1,3]
    if nnvals   = min(nvals,len(wobbles))
    wobbles = np.array(wobbles)[:nvals]
    plts,cs = [],dsp.getCs('coolwarm',nvals)
    for i,val in zip(range(nvals),wobbles):
        tail = str(i).zfill(int(wobbles.size/10)+1)
        data = Li_xyz(name,n=[0,0,1],dopt=1,wobble=val,tail='_'+tail)
        multi=mupy.Multislice(name,tail=tail,data=data,
            mulslice=0,keV=100,TDS=True,n_TDS=15,
            NxNy=512,slice_thick=1.5,Nhk=5,repeat=[2,2,2],
            **kwargs)
        q,I = multi.azim_avg(out=1,opt='')
        plts+=[[q,I,cs[i],'w=%.3f' %val]]
    dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$','$I_q$'],opt='p',lw=2)

##########################################################################
if __name__ == "__main__":
    plt.close('all')
    name = 'dat/Lithium/'
    # Li_xyz(name,n=[0,0,1],dopt=1)
    # Li_patterns(name+'gif/',nzs=20, opt='drsp',fopt='',ppopt='w',v=1)
    # Li_gif(name+'gif/')
    # df = Li_latsize(name+'latsize/',nvals=3,opt='drsp',fopt='',ppopt='w',v=1)
    Li_wobble(name+'wobble/',nvals=2,opt='drsp',fopt='',ppopt='w',v=1,ssh='local_london')
