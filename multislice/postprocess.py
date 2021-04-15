#from utils import *
import pickle5,os,glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from utils import displayStandards as dsp
from utils import glob_colors as colors
from utils.glob_colors import*
from scattering import scattering_factors as scat
# import get_elec_atomic_factors,wavelength
from scattering.structure_factor import structure_factor3D

hosts = ['','brno','tarik-CCP4home']
cpus  = dict(zip(hosts,['asus-$i_5^{(4cores)}$','brno-$i_7^{(8cores)}$','xps-$i_7^{(8cores)}$']))
ncpus = dict(zip(hosts,[4,8,8]))
info_cols = ['zmax(A)','Inorm','cputime(s)','walltime(s)']

#########################################################################
### import from TEMSIM
#########################################################################
def import_beams(file,slice_thick=1,iBs=[],tol=1e-2,Obeam=False,pImax=False):
    ''' Import beam files with some postprocessing option\n
    - slice_thick : slice thickness used for simulation
    - iBs   : selected beam indices : default=Ibeams.max()>tol
    - tol   : select beams above threshold if iBs==[]
    - Obeam : include origin
    '''
    hk    = open(file).readline().rstrip().split('=  ')[1].split(' ');#print(hk)
    #print(hk)
    beams = np.loadtxt(file,skiprows=3).T
    if isinstance(slice_thick,list) :
        idx,beams = beams[0,:],beams[1:,:]
        n_slices = int(len(idx)/len(slice_thick)); #print(n_slices,len(slice_thick))
        t = np.cumsum(np.tile(slice_thick,n_slices))+0.75*slice_thick;#print(t[:4])
    else :
        idx,beams = beams[0,:-5],beams[1:,:-5]
        t = idx*slice_thick+0.75*slice_thick
    #get re,im,I
    nbs = len(hk)
    re = [beams[2*i,:] for i in range(nbs)]
    im = [beams[2*i+1,:] for i in range(nbs)]
    Ib = [np.abs(r+1J*i)**2 for r,i in zip(re,im)]
    #filter
    Imax = np.array([I.max() for I in Ib]); #print(Imax,Imax.max())
    if any(iBs) :
        if isinstance(iBs[0],str) :
            hk = np.array(hk)
            iBs = [ (iB==hk).argmax() for iB in iBs if (iB==hk).sum()]
    else :
        iBs = [i for i in range([1,0][Obeam],nbs) if Imax[i]>= tol*Imax.max()]
    if pImax : print('Imax:',Imax,'iBs:',iBs,', nbs=',nbs)
    hk,t,re,im,Ib = np.array(hk)[iBs],t, np.array(re)[iBs],np.array(im)[iBs],np.array(Ib)[iBs]
    return hk,t,re,im,Ib

#########################################################################
#### def : DataFrame utilities
#########################################################################
def load(datpath,tail=''):
    filename = datpath
    pkls = glob.glob(datpath+'*.pkl')
    filename = pkls[0]
    if tail:
        tail_files = [pkl for pkl in pkls if tail in pkl]
        if tail_files:
            filename=tail_files[0]
        else:
            print(colors.red+'warning: tail "%s" not found in %s' %(tail,datpath)+colors.black )
    print(colors.green+'loading ' +colors.yellow+filename+colors.black )
    return load_multi_obj(filename)

def load_multi_obj(filename):
    '''load a saved Multislice object
    filename : pickle file (.pkl)  '''
    with open(filename,'rb') as f : multi = pickle5.load(f)
    return multi

def get_info(log_file):
    '''compute zmax,I,cpuTime and wallTime'''
    with open(log_file,'r') as f : log_lines=f.readlines()
    log_lines = log_lines[-2:]
    info = [l.strip().split('=')[-1].split(' ')[1] for l in log_lines]
    return np.array(info,dtype=float)

def get_info_cpu(log_file,mulslice=False):
    '''compute zmax,I,cpuTime and wallTime'''
    with open(log_file,'r') as f : log_lines=f.readlines()
    if mulslice:
        log_lines = [log_lines[-6]] + log_lines[-2:]
        info = [l.strip().split('=')[-1].split(' ')[1] for l in log_lines]
        info = [log_lines[0].strip().split(',')[0].replace('slice','')]+info #zmax
    else:
        log_lines = log_lines[-7:-5] + log_lines[-2:]
        info = [l.strip().split('=')[-1].split(' ')[1] for l in log_lines]
    return np.array(info,dtype=float)

def update_df_info(df_path,hostpath=None,files=[]):
    df = pd.read_pickle(df_path)
    datpath = os.path.dirname(df_path)+'/'
    for name in df.index:
        multi = load_multi_obj(datpath+name);#print(name)
        ssh_host = df.loc[name].host
        state    = multi.check_simu_state(ssh_alias=ssh_host,v=0,hostpath=hostpath)
        df.loc[name]['state'] = state
        info = get_info(multi._outf('log'))
        df.loc[name][info_cols[:2]] = info
        if state=='done':
            info = get_info_cpu(multi._outf('log'),mulslice=multi.is_mulslice)
            df.loc[name][info_cols] = info
            for file in files:multi.ssh_get(ssh_host,file)
    df.to_pickle(df_path);
    print(green+'DataFrame updated and saved : \n' +yellow+df_path+black)
    return df

#########################################################################
### def : display
#########################################################################
def plot_beam_thickness(beams,rip='I',linespec='-',cm='Greens',**kwargs):
    ''' plot the beams as function of thickness
    - beams : beams info from import_beams or file containing them
    - rip flags : I(Intens),r(real),i(imag)
    '''
    hk,t,re,im,Ib = beams
    nbs = len(hk)
    csp,csr,csi,plts = dsp.getCs(cm,nbs),dsp.getCs('Blues',nbs),dsp.getCs('Reds',nbs),[]
    for i in range(nbs):
        if 'I' in rip : plts += [[t,Ib[i],[csp[i],'%sx' %linespec],'$I_{%s}$' %(hk[i])]]
        if 'r' in rip : plts += [[t,re[i],csr[i],'$re$']]
        if 'i' in rip : plts += [[t,im[i],csi[i],'$im$']]
    return dsp.stddisp(plts,lw=2,labs=['$thickness(\AA)$','$I_{hk}$'],**kwargs)#,xylims=[0,t.max(),0,5])

def plot_fe(qmax=1,Nx=20,Vcell=1,**kwargs):
    q0,fq_e = scat.get_elec_atomic_factors(['Si'],q=np.linspace(0,qmax,1000))
    vg = fq_e[0]/Vcell
    vgmax = vg.max()
    Xq =np.arange(Nx)[:,None].dot(np.array([1,1])[None,:]).T/ax
    Yq = np.tile([0,vgmax],[Nx,1]).T
    plts=[[Xq,Yq,[(0.7,1.0,0.7),'--'],''],
          [q0,vg,'g--','$v_g(Si)$']]
    dsp.stddisp(plts,labs=['$q(A^{-1})$','$v_g(q)(A^{-2})$'], xylims=[0,qmax,0,vgmax],
        legLoc='upper right',**kwargs)

def plot_pattern_section(file='si_autoslic100.tif'):
    im=plt.imread(datpath+file)
    N=im.shape[0]
    std.stddisp([range(N),im[:,0],'b','$I_{h0}$'])

###########################################################################
#### def : test
###########################################################################
if __name__ == "__main__":
    plt.close('all')
