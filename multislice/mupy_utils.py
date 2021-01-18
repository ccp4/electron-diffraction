import pandas as pd, numpy as np
from utils import displayStandards as dsp
from utils import glob_colors as colors
from crystals import Crystal
import multislice as mupy
# import rotating crystal as rcc

def ewald_sphere(lat_params,lam=0.025,tmax=7,T=0.2,nx=20,ny=10,**kwargs):
    '''lam(wavelength Angstrum), tmax(angle degrees), T(Thickness mum)'''
    a1,a2,a3 = lat_params
    b1,b2,b3 = 1/a1,1/a2,1/a3
    tmax, K=np.deg2rad(tmax), 1/lam
    dqT = 10/(T*1e4)*np.array([-1,1])
    #reciprocal lattice
    h,k = np.meshgrid(b1*np.arange(-nx,nx+1),b2*np.arange(-1,ny))
    #sphere
    t = 3*np.pi/2+np.linspace(-tmax,tmax,100)
    plts = [[K*np.cos(t),K*np.sin(t)+K,'r',r'$\lambda=%.2f\AA$' %lam]]
    #rel rods
    plts += [[ [h,h],k+dqT ,'b-',''] for h,k in zip(h.flatten(),k.flatten())]

    #display
    scat=[h,k,15,'b']
    fig,ax = dsp.stddisp(plts,scat=scat,labs=['$q_x$','$q_y$'],
        lw=3,#,xylims=[-nx*b1,nx*b1,-b2,ny*b2],xyTickLabs=[[],[]],
        **kwargs)


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


def make_mulslice_datfile(dat_file,cif_file):
    ''' create a data file used by atompot(mulslice) from a cif file
    - dat_file : name of file to save
    - cif_file : the file to import
    '''
    crys = Crystal.from_cif(cif_file)
    pattern = np.array([np.hstack([a.coords_cartesian,a.atomic_number]) for a in crys.atoms] )
    deck = ' '.join(np.array(crys.lattice_parameters[:3],dtype=str))+'\n'
    deck+='0\n'
    Za = np.array(np.unique(pattern[:,-1]),dtype=int)
    for Z in Za:
        deck+=str(Z)+'\n'#; print()
        xyz = np.array(pattern[pattern[:,-1]==1][:,:3],dtype=str)
        xyz = '\n'.join([' '.join(l) for l in xyz]) #; print(xyz)
        deck+=xyz+'\n'
    deck+='\n\nComment here'
    with open(dat_file,'w') as f : f.write(deck)
    # with open(datpath,'r') as f :print(f.readlines())


def import_cif(file,xyz='',n=[0,0,1],rep=[1,1,1],dopt='s',lfact=1.0,tail=''):
    '''
    - file : cif_file
    - rep : super cell repeat
    - n : reorient of the z axis into n
    '''
    crys    = Crystal.from_cif(file)
    lat_vec = np.array(crys.lattice_vectors)
    if sum(rep)>3 : crys = crys.supercell(rep[0],rep[1],rep[2])
    pattern = np.array([[a.atomic_number]+list(lfact*a.coords_cartesian)+[a.occupancy,1.0] for a in crys.atoms])
    if xyz:make_xyz(xyz,pattern,lat_vec,n,fmt='%.4f',dopt=dopt)
    return pattern #,lat_params #pattern,crys # file

def make_xyz(name,pattern,lat_vec,n=[0,0,1],fmt='%.4f',dopt='s'):
    '''Creates the.xyz file from a given compound and orientation
    - name    : Full path to the file to save
    - pattern : Nx6 ndarray - Z,x,y,z,occ,wobble format
    - lat_vec : 3x3 ndarray - lattice vectors [a1,a2,a3]
    - n : beam direction axis
    - dopt : p(print file),s(save)
    '''
    compound = dsp.basename(name)
    ax,by,cz = np.diag(lat_vec)
    # cz = np.linalg.norm(n)
    # coords = orient_crystal(pattern[:,1:4],n_u=n)

    # if 'c' in dopt:
    #     x,y,z = coords.T
    #     ax,by = np.abs(x.min()-x.max()),np.abs(y.max()-y.min())
    #     idx,idy,idz = x<0,y<0,z<0
    #     # coords[idx,0] += ax
    #     # coords[idy,1] += by
    #     coords[idz,2] += cz
    # pattern[:,1:4] = coords
    #write to file
    if 's' in dopt :
        dir=''.join(np.array(n,dtype=str))
        header = 'one unit cell of %s\n' %(compound)
        header+= ' '.join([fmt]*3) %(ax,by,cz)
        np.savetxt(name,pattern,footer='-1',header=header,fmt='%d '+' '.join([fmt]*5),comments='')
        print(colors.green+"coords file saved : \n"+colors.yellow+name+colors.black)
        if 'p' in dopt :
            with open(name,'r') as f : print(''.join(f.readlines()))
    return pattern,[ax,by,cz]
