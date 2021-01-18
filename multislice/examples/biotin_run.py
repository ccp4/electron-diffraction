import sys,numpy,os
from utils import*
from crystals import Crystal
import multislice.multislice as mupy
import multislice.postprocess as pp
import multislice.rotating_crystal as rcc
import importlib as imp
imp.reload(mupy)
imp.reload(pp)
imp.reload(rcc)

plt.close('all')
path = '../dat/biotin/'
file = path+'biotin.cif'

# rcc.show_cell(file,x0=-1)

def create_xyz(opts='gps'):
    nz = np.array([1,2,1,1,1,1,1,1 ,1,1 ])#,1)
    nx = np.array([0,1,1,2,3,4,6,8 ,12,24])#,1,7,12,24])
    angles = np.arctan(nx/(4*nz))*180/np.pi
    data = ['biotin%d0%d.xyz' %(n,m) for n,m in zip(nx,nz)]
    if opts:
        # rcc.import_cif(file,path+'biotin001.xyz')
        Nx,Nz=[(1,1),(30,5)]['g' in opts]
        rcc.rotate_xyz(file,nx,nz,Nx=Nx,Nz=Nz,opt=opts)
    for dat,a in zip(data,angles) : rcc.show_grid(path+dat,title=dat+' %.1f' %a)
    print('angles:',angles)
    return data,angles

def run_main(dfname='df0.pkl',thick=10,ssh='',hostpath=''):
    cols  = ['host','state','dat','angle','tilt']
    df    = pd.DataFrame(columns=cols+pp.info_cols)
    Nx = np.array(np.round(lats[:,0].max()/lats[:,0]),dtype=int)
    Ny = np.array(np.round(Nx*lats[:,0]/lats[:,1]),dtype=int)
    Nz = np.array(np.round(thick/lats[:,2]),dtype=int)
    for iD,dat in enumerate(data):
        istr  = ('m%d' %iD).zfill(2)
        multi = mupy.Multislice(path,data=data[iD],tail=istr,
            mulslice=False,keV=200,#tilt=[tx*np.pi/180/1000,0],
            NxNy=2048,slice_thick=1.0,Nhk=5,repeat=[Nx[iD],Ny[iD],Nz[iD]],
            #TDS=True,T=300,n_TDS=15,
            opt='sr',fopt='f',v='nctr',#nctrdDR',
            ssh=ssh,hostpath=hostpath
            )
        df.loc[multi.outf['obj']] = [np.nan]*len(df.columns)
        df.loc[multi.outf['obj']][cols] = [ssh,'start',data[iD],angles[iD],0]
    df.to_pickle(path+dfname)
    print(green+'Dataframe saved : '+yellow+path+dfname+black)

def run_simus(thick=5000,nts=91,ssh=''):
    cols  = ['host','state','dat','angle','tilt']
    df    = pd.DataFrame(columns=cols+pp.info_cols)
    tilts = np.linspace(0,90,nts)
    Nx = np.array(np.round(lats[:,0].max()/lats[:,0]),dtype=int)
    Ny = np.array(np.round(Nx*lats[:,0]/lats[:,1]),dtype=int)
    Nz = np.array(np.round(thick/lats[:,2]),dtype=int)
    for i in range(nts):
        pad   = int(np.log10(nts))+1
        istr  = ('%d' %i).zfill(pad)
        iD    = np.abs(tilts[i]-angles).argmin()
        tx    = tilts[i]-angles[iD]
        multi = mupy.Multislice(path,data=data[iD],tail=istr,
            mulslice=False,keV=200,tilt=[tx*np.pi/180/1000,0],
            NxNy=512,slice_thick=1.0,Nhk=5,repeat=[Nx[iD],Ny[iD],Nz[iD]],
            #TDS=True,T=300,n_TDS=15,
            opt='sr',fopt='f',v='nctr',#nctrdDR',
            ssh=ssh,
            )
        df.loc[multi.outf['obj']] = [np.nan]*len(df.columns)
        df.loc[multi.outf['obj']][cols] = [ssh,'start',data[iD],tilts[i],tx]
    df.to_pickle(path+'df.pkl')
    print(green+'Dataframe saved : '+yellow+path+'df.pkl'+black)

        # multi.wait_simu(ssh_alias='tarik-CCP4home')
        # multi.print_log()
        # multi.postprocess(ppopt='uwP',ssh_alias='tarik-CCP4home')

def run_tilts(**kwargs):
    tx = np.linspace(0,0.1,10)#deg
    tilts = [[0,t*np.pi/180*1000] for t in tx]
    mupy.sweep_var(path+'tilts/','tilt',tilts,tail='m',df='dfm.pkl',do_prev=0,**kwargs)


def update_patterns(dfname='df.pkl',ssh='',hostpath=''):
    df = pp.update_df_info(path+dfname,hostpath=hostpath)
    for dat in df.index:
        multi=pp.load_multi_obj(path+dsp.dirname(dfname)+dat)
        multi.ssh_get(ssh,'pattern',hostpath=hostpath)
        multi.save_pattern()
    return df

def get_figs(dfname='df.pkl'):
    df = pd.read_pickle(path+dfname)
    cs = dsp.getCs('Blues',3)#,dsp.getCs('Reds'),dsp.getCs('Greens')
    plts = [[],[],[]]
    #pad = int(np.log10(df.index.size))+1
    for i,dat in enumerate(df.index):
        multi=pp.load_multi_obj(path+dsp.dirname(dfname)+dat)
        multi.pattern(tol=1e-4,Iopt='Insgl',gs=0.25,caxis=[-6.3,0],Nmax=500,
            cmap='binary',imOpt='hc',#rings=[0.1,0.2],#,1],,
            pOpt='Xt',xyTicks=1,xylims=[-3,3.01,-3,3.01],#xylims=[-2,2.01,-2,2.01],
            opt='s',name=path+'figures/%s_pattern.png' %dat.replace('.pkl',''))
        # multi.beam_vs_thickness()
    #     for j in range(len(hk)):plts[jq]+=[[t,I[j,:],cs[i],'']]
    # for j in range(3):
    #     dsp.stddisp(plts[j],title='hk=%s' %hk[j],
    #         imOpt='ch',cmap='Blues',caxis=[0,90],pOpt='tG')


# data,angles = create_xyz(opts='ps')
# lats  = np.array([rcc.show_grid(path+dat,opt='')[0] for dat in data])
# multi = run_main(thick=2000,ssh='badb',hostpath='/data3/lii26466/multislice/biotin/')

# multi = run_simus(thick=100,nts=3,ssh='')#'tarik-CCP4home')
multi = run_tilts(keV=200,mulslice=False,
    NxNy=2048,slice_thick=1.0,Nhk=5,repeat=[2,1,100],
    opt='dsr',fopt='f',v=1,
    ssh='badb',cluster=1,hostpath='/data3/lii26466/multislice/biotin/tilts/')

# df = update_patterns('tilts/df.pkl',ssh='badb',hostpath='/data3/lii26466/multislice/biotin/tilts/')
# get_figs(dfname='tilts/df.pkl')

# df = pp.update_df_info(path+'df0.pkl',hostpath='/data3/lii26466/multislice/biotin/')
# df = update_patterns('df0.pkl',ssh='badb',hostpath='/data3/lii26466/multislice/biotin/')
# get_figs(dfname='df0.pkl')
