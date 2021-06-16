'''
A Multislice interface to build decks, launch jobs and postprocess multislice simulations\n
[TOC]\n
##Example
Create a periodic simulation from data in 'Silicon/Si110',
submit it and monitor its state :
```python
multi=Multislice('Silicon/Si110',mulslice=False,NxNy=512,repeat=[6,6,40],Nhk=2)
p = multi.run()
multi.check_simu_state()
```
Once terminated, get the beams intensities as function of thickness,
forcing to repostprocess in case it was already done in an older simulation :
```python
multi.beam_vs_thickness(bOpt='f')
```
'''
import importlib as imp
import pickle,socket,time,tifffile
import pandas as pd
import numpy as np
from math import ceil,nan
from subprocess import Popen,check_output,PIPE
from glob import glob as lsfiles
from os.path import dirname,realpath,exists,basename
from string import ascii_uppercase as SLICES #ascii_letters
from string import ascii_letters
# from matplotlib.animation import FuncAnimation
# from utils.glob_colors import*
from crystals import Crystal
import utils.glob_colors as colors
import utils.displayStandards as dsp                        ;imp.reload(dsp)
# from scattering.structure_factor import structure_factor3D
from . import postprocess as pp                             ;imp.reload(pp)
from . import mupy_utils as mut                             ;imp.reload(mut)

ssh_hosts = {
    'local_london':'BG-X550',
    'tarik-CCP4stfc':'tarik-CCP4','tarik-CCP4home':'tarik-CCP4',
    'brno':'brno',
    'badb':'badb',
    'badb.rc-harwell.ac.uk':'badb',
}

temsim_hosts={
    'BG-X550'   :'/home/ronan/Documents/git/ccp4/src/electron-diffraction/multislice/bin/',
    'tarik-CCP4':'/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/bin/',
    'brno'      :'/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/bin/',
    'badb'      :'/home/lii26466/Documents/git/ccp4/src/electron-diffraction/multislice/bin/',
    'badb.rc-harwell.ac.uk':'/home/lii26466/Documents/git/ccp4/src/electron-diffraction/multislice/bin/',
}
path_hosts={
    'BG-X550'   :'/home/ronan/Documents/git/ccp4/src/electron-diffraction/multislice/dat/',
    'tarik-CCP4':'/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/dat/',
    'brno'      :'/home/tarik/Documents/git/ccp4/src/electron-diffraction/multislice/dat/',
    'badb'      :'/data3/lii26466/multislice/',
    'badb.rc-harwell.ac.uk' :'/data3/lii26466/multislice/',
}


class Multislice:
    ''' **DATA PARAMETERS :**\n
    - `name` : path to the simulation folder
    - `tail` : trailing string in the naming convention of outputs (tag is an alias as well)
        - name pattern = *name*\_*tail*\_*program_name*
    - `data` : data file names
        - if None and mulslice simu : all *[a-zA-Z].dat and stacked alphabetically
        - if None and autoslic simu : first *.xyz found
        - Otherwise paths to filename(s) of the data files to use (must be in name)
    \n**MULTISLICE simulation parameters :**\n
    - `mulslice` : False(autoslic), True(mulslice)=>periodic
    - `keV` : float - wavelength in keV
    - `tilt` : 2-list - crystal tilt tx,ty (mrad)
    - `repeat` : 3-int-list - super cell size in the x,y,z directions
    - `NxNy` : 2-int-list or int : sampling in x and y (same sampling in x and y if only int provided)
    - `slice_thick` : float - slice thickness (A) (`dz` can also be used as an alias)
    - `i_slice` : record intensity every i_slice
    - `hk` : list of tuple - beams to record as function of depth
    - `Nhk`: int - set hk on a grid of multiples of the supercell size. Prevails over hk argument if >0. default(0)
    - `hk_sym` : True (take the Friedel pairs for the beams in hk)
    \n**Thermal parameters**:\n
    - `TDS` : bool - include thermal vibrations (as dictated by wobble parameter in .dat )
    - `T`   : Temperature in K if TDS is True
    - `nTDS` : Number of runs to average over
    \n**RUN/OUTPUT options:**\n
    - `v`  : verbose - str or bool(all if True) or int(verbose level)
        - n(naming pattern),c(cell params),t(thickness),r(run cmd)
        - d(data),D(Decks),R(full run cmd)
    - `opt` : d(save_deck), s(save_obj), r(do_run),  f(force rerun),w(ask before running again), p(do_pp)
    - `fopt` : If the simulation was previously completed :\n
        - '' : The simulation will not be run again
        - 'w' is on  (default case) : the user will be asked to confirm whether he wants to run it again
        - 'f' in fopt : The simulation is rerun
    - `ppopt` : u(update),w(wait), I(image) B(beam) P(pattern) A(azim_avg) S(reduced pattern) s(save_patterns in job)
    - `ssh` : ip address or alias of the machine on which to run the job
    - `cif_file` : name of .cif file corresponding to structure in path
    ### OBSOLETE
    - `hostpath` : path to data on host
    - `cluster` : indicate if running on a cluster using qsub
    '''
    def __init__(self,name,mulslice=False,tail='',tag=None,data=None,u=None,
                tilt=[0,0],TDS=False,T=300,n_TDS=16,
                keV=200,repeat=[1,1,1],NxNy=[512,512],slice_thick=1.0,dz=None,i_slice=1000,
                hk=[(0,0)],Nhk=0,hk_pad=None,hk_sym=0,prev=None,
                opt='sr',fopt='',ppopt='uwPs',v=1,
                ssh=None,hostpath='',cluster=False,cif_file=None):

        v = self._get_verbose_options(v)
        #attributes
        self.version     = '1.4.3b'
        self.name        = dsp.basename(name)                       #;print('basename:',self.name)
        self.datpath     = realpath(name)+'/'                       #;print('path:',self.datpath)
        self.cif_file    = self.get_cif_file(cif_file)              #;print(self.cif_file)
        self.is_mulslice = mulslice
        self.tilt        = tuple(tilt)
        self.TDS         = TDS
        self.T           = T
        self.n_TDS       = n_TDS
        self.tail        = self._set_tail(tail,tag)
        self.name        = self._set_name('n' in v)
        self.data        = self._get_datafiles(data)
        self.u           = u
        self.slices      = SLICES[:len(self.data)]
        self.outf        = self._set_filenames()
        self.keV         = keV
        self.repeat      = tuple(repeat)
        self.NxNy        = self._set_NxNy(NxNy)
        self.hk          = self._set_hk(hk,Nhk,hk_sym,hk_pad)
        self.cell_params = self._get_cell_params('c' in v)
        self.slice_thick = self._set_slice_thick(slice_thick,dz)
        self.i_slice     = i_slice
        self.thickness   = self._get_thickness('t' in v)
        self.decks       = {}
        self.p           = None

        ## make decks and run if required
        self.execute(opt,fopt,ppopt,v, ssh,cluster,hostpath,prev)

    ########################################################################
    ##### Public functions
    ########################################################################
    def make_decks(self,save=True,datpath=None,prev=None,v=False):
        '''create the decks from the information provided'''
        if self.is_mulslice:
            decks={self.outf['simu_deck']:self._mulslice_deck(prev)}
            for i in self.slices:
                decks[self.outf['slice_deck%s' %i]]=self._atompot_deck(i)
        else:
            decks={self.outf['simu_deck']:self._autoslic_deck(prev)}
        if save:
            #datapath is manually provided for remote generation
            if not datpath : datpath = self.datpath
            if v:print(colors.green+'Decks saved :'+colors.black)
            for filename,deck in decks.items():
                with open(datpath+filename,'w') as f : f.write(deck)
                if v:print(colors.yellow+datpath+filename+colors.black)
        self.decks = decks

    def save(self,v=False):
        '''save this multislice object'''
        file=self._outf('obj')
        with open(file,'wb') as out :
            pickle.dump(self, out, pickle.HIGHEST_PROTOCOL)
        if v:print(colors.green+"object saved\n"+colors.yellow+file+colors.black)

    def execute(self,opt='sr',fopt='',ppopt='w',v=1, ssh='',cluster=False,hostpath='',prev=None):
        # save_deck,save_obj,do_run,do_pp,fopt,vopt = self._get_run_options(opt,fopt,v)
        save_deck,save_obj,do_run,do_pp = [s in opt for s in 'dsrp']
        save_deck |= do_run
        if 'f' in opt: fopt += 'f'
        if 'w' in opt: fopt += 'w'
        vopt=('r' in v and 'R' not in v) + 2*('R' in v)
        self.patterns_saved,patterns_opt=0,'s' in ppopt         #patterns_saved set to 1 during job postprocess
        if save_deck : self.make_decks(save_deck,prev=prev,v=v)
        if save_obj : self.save(v=v)
        if 'd' in v : self.print_datafiles()
        if 'D' in v : self.print_decks()
        if do_run : self.p = self.run(v=vopt, fopt=fopt,ssh_alias=ssh,hostpath=hostpath,cluster=cluster,patterns_opt=patterns_opt)
        if do_pp : self.postprocess(ppopt,ssh,hostpath=hostpath)
        self._set_figpath()

    def resume(self,Nz,v=1,opt='sr',
        hk=None,Nhk=None,hk_sym=0,i_slice=None,prev=1,
        **kwargs):
        ''' resume a simulation from its last point
        - Nz : number of unit cells to run
        - hk,Nhk,hk_sym,i_slice : values which can be set
        - **kwargs : opt,fopt,ppopt,v, ssh,cluster,hostpath
        '''
        v = self._get_verbose_options(v=v)
        prev = self._set_prev(prev)
        self.repeat = self.repeat[:2]+(Nz,)
        if hk or Nhk : self.hk = self._set_hk(hk,Nhk,hk_sym)
        if i_slice : self.i_slice = i_slice
        self.thickness += self._get_thickness('t' in v)
        self.execute(v=v,prev=prev,opt='sr'+opt,**kwargs)
        self.merged=0
        self.merge_beams()

    def run(self,v=1,fopt='w',ssh_alias=None,hostpath=None,cluster=False,patterns_opt=False,):
        '''run the simulation with temsim
        - fopt : f(force rerun), w(warn ask rerun already done)
        - ssh : name of the host to run the job
        '''
        if isinstance(fopt,int):fopt='f'
        run = True
        if self.check_simu_state(v=0,ssh_alias=''):
            if v>0 : print(colors.red+"Simulation already performed in the past."+colors.black)
            if 'f' in fopt :
                msg='Force re-running'
                cmd = '''echo Deleting logfile and patterns
                if [ -f %s ];then rm %s;fi
                rm %s
                '''%(self._outf('log'),self._outf('log'),self._outf('pattern')+'*')
                p = Popen(cmd, shell=True,stderr=PIPE,stdout=PIPE)
                p.wait();p.communicate()
            else :
                if 'w' in fopt :
                    run = input(colors.red+"Re-run?(y/n) : "+colors.black)=='y'
                    msg = 're-running'
                else :
                    run,msg = False, 'not running'
            if v : print(colors.red+msg+colors.black)
        p = None
        if run:
            if ssh_alias :
                if ssh_alias=='badb':cluster=1
                cmd = self._get_job_ssh(ssh_alias,hostpath=hostpath,v=v>2,cluster=cluster,patterns_opt=patterns_opt)
                #time.sleep(1)
                #self.check_simu_state(v=0,ssh_alias=ssh_alias,hostpath=hostpath)
            else :
                self._get_job(cluster=cluster, patterns_opt=patterns_opt)
                cmd = 'bash %s' %self._outf('job')
            p = Popen(cmd,shell=True) #; print(cmd)
            if v>0 : print(colors.green+self.name+" job submitted at %s" %time.ctime()+colors.black)
            if v>1 : print(colors.magenta+cmd+colors.black)
        return p

    def postprocess(self,ppopt='wuP',ssh_alias='',tol=1e-4,figpath=None,opt='p',hostpath=''):
        '''Performs postprocessing with predefined options : \n
        - ppopt:w(wait until done) u(update) d(display) for I(image), B(beam) P(pattern)     #f(force recalculate)
        - figpath : Directory to place the figures if saving with automatic naming (default datpath)
        '''
        print(colors.blue+'...postprocessing...'+colors.black)
        if not figpath : figpath = self.datpath
        time.sleep(0.1)
        state = self.check_simu_state(ssh_alias=ssh_alias,hostpath=hostpath);print(state)
        if not state == 'done' and 'w' in ppopt:
            self.wait_simu(ssh_alias,10,hostpath)
            self.p.wait()
        # self.get_beams(iBs='a',bOpt='f')
        if ssh_alias and 'u' in ppopt:
            if 'I' in ppopt : self.ssh_get(ssh_alias,'image'     ,hostpath)
            if 'B' in ppopt : self.ssh_get(ssh_alias,'beams'     ,hostpath)
            if 'P' in ppopt : self.ssh_get(ssh_alias,'patternnpy',hostpath)
            if 'S' in ppopt : self.ssh_get(ssh_alias,'patternS'  ,hostpath)
        #convert to np.array
        if 'd' in ppopt:
            if 'I' in ppopt and opt : self.image(opt=opt,name=figpath+self.outf['imagesvg'])
            if 'B' in ppopt and opt :
                # if not exists(self._outf('beams')) or 'f' in ppopt:self.get_beams(bOpt='fa')
                self.beam_vs_thickness(bOpt='f',tol=tol,opt=opt,name=figpath+self.outf['beamssvg'])
            if 'P' in ppopt and opt :
                # if not exists(self._outf('patternnpy')) or 'f' in ppopt:self.save_pattern()
                self.pattern(Iopt='Incsl',tol=tol,imOpt='ch',cmap='gray',opt=opt,name=figpath+self.outf['patternsvg'])

    ###################################################################
    ################ OUTPUT FUNCTIONS
    ###################################################################
    def get_cif_file(self,cif_file):
        if not cif_file:cif_file=self.name+'.cif'
        return self.datpath+cif_file
    def get_structure_factor(self,**sf_args):
        '''computes structure factor 3D
        - sf_args : see (mut.get_structure_factor3D)
        returns :
        - (qx,qy,qz),Fhkl
        '''
        return mut.get_structure_factor(self.cif_file,**sf_args)
        # crys = Crystal.from_cif(self.cif_file)
        # pattern = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in crys.atoms] )
        # lat_vec = np.array(crys.reciprocal_vectors)
        # (h,k,l),Fhkl = structure_factor3D(pattern, lat_vec, **sf_args)
        # qx = h/crys.lattice_parameters[0]
        # qy = k/crys.lattice_parameters[1]
        # qz = l/crys.lattice_parameters[2]
        # return (qx,qy,qz),Fhkl

    def image(self,opt='I',cmap='jet',**kwargs):
        '''Displays the 2D image out of simulation
        - opt : I(intensity)
        '''
        im=plt.imread(self._outf('image'))
        if 'I' in opt :
            N = int(im.shape[1]/2)
            real,imag = im[:,:N],im[:,N:];#print(real.max())
            im = real**2+imag**2
            im =im/im.max()
        dsp.stddisp(labs=['$x$','$y$'],im=im,legOpt=0,imOpt='c',**kwargs)

    def gif_beams(self,step=5,cmap='gray'):
        '''make a gif of the saved beams pattern'''
        hk,ts,re,im,Ib = self.get_beams(iBs=[],tol=1e-5,bOpt='fa')
        ax,by = self.cell_params[:2]
        Nx,Ny,Nz = self.repeat
        h,k = np.array(self.hk).T
        h,k = h/ax/Nx,k/by/Ny
        # nzf=int(np.log10(Nz/step))+1
        Ib = np.log10(Ib)
        Imax=Ib.max()
        gif_path=self.datpath+'gif/'#+self.name
        p=Popen("if [ ! -d %s ];then mkdir %s;fi" %(gif_path,gif_path),shell=True);p.wait()#.decode()
        ##### display
        # fig,ax = dsp.create_fig(figsize='f')
        # ims = []
        for i,it in enumerate(ts[::step]):
            tle = '$z=%.0f \AA$' %it
            # im = ax.scatter(h,k,Ib[:,i],cmap=cmap,vmin=-1,vmax=,animated=True)
            # ims = [im]
            dsp.stddisp(scat=[h,k,Ib[:,i]],labs=['$q_x(\AA^{-1})$','$q_y(\AA^{-1})$'],
                caxis=[0,Imax],imOpt='ch',#xylims=[-5,5,-5,5],
                texts=[0,4,tle,'g'],
                name=gif_path+'%s.png' %str(i).zfill(nzf),opt='sc')

        # dsp.standardDisplay(labs=['$q_x(\AA^{-1})$','$q_y(\AA^{-1})$'],
        #     ,caxis=[0,Imax],imOpt='ch',#xylims=[-5,5,-5,5],
        #     texts=[0,4,tle,'g'],
        #     name=gif_path+'%s.png' %str(i).zfill(nzf),opt='sc')

        # ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
        #     repeat_delay=1000)
        # ani.save('%s/beams.mp4' %self.datpath)
        # fig.show()

        # cmd = "im2gif %s && eog %s.gif" %(gif_path,gif_path)
        # cmd="convert -delay 20 -loop 0 %s*.png %s" %(gif_path,self.datpath+self.name+'.gif')
        # print("creating the gif")
        # print(check_output(cmd,shell=True).decode())

    def beam_vs_thickness(self,bOpt='',orig=0,iBs=[],tol=1e-2,**kwargs):
        ''' plot beam vs thickness
        - bOpt : O(include Origin),p(print Imax),n or f(force new)
        - iBs : beams (recorded beams are found in self.hk)
            - list of str ['(h0,k0)','(h1,k1)'..]
            - indices
        - kwargs : see help(plot_beam_thickness)
        '''
        hk,t,re,im,Ib = np.load(self._outf('beams'),allow_pickle=True)
        # hk,t,re,im,Ib = self.get_beams(iBs=[],tol=1e-5,bOpt='fa')
        # print(hk)
        if orig:bOpt+='O'
        if any(iBs) :
            if isinstance(iBs[0],str) :
                iBs = [ (iB==hk).argmax() for iB in iBs if (iB==hk).sum()]
                print(iBs)
            if isinstance(iBs[0],tuple):
                iBs = [i for i,hk0 in enumerate(self.hk) if hk0 in iBs]
        else :
            Imax = np.array([I.max() for I in Ib]) #;print(Imax)
            iBs = [i for i in range('O' not in bOpt,len(hk)) if Imax[i]>= tol*Imax.max()]
        beams = np.array(hk)[iBs],t, np.array(re)[iBs],np.array(im)[iBs],np.array(Ib)[iBs]

        if 'o' in bOpt:
            return beams
        else:
            pp.plot_beam_thickness(beams,**kwargs)

    def merge_beams(self):
        if not self.merged:
            beamnpy = self._outf('beams')
            p = Popen('cp %s %s' %(beamnpy,beamnpy+'%s' %str(int(self.thickness)).zfill(4)),stderr=PIPE,stdout=PIPE,shell=True)
            p.wait();p.communicate()
            beams_files = np.sort(lsfiles(self._outf('beams')+'*'))[1:]     #; print(beams_files)
            hk,t,re,im,I  = np.load(beams_files[0],allow_pickle=True)
            z=t.max()
            for file in beams_files[1:]:
                hk_i,t_i,re_i,im_i,I_i = np.load(file,allow_pickle=True)
                t_i +=z                                                     #;print(file,z)
                t,I = np.hstack([t,t_i]),np.hstack([I,I_i])
                z   += t_i.max()
            # self.merged=1
            np.save(beamnpy,[hk,t,re,im,I])
            print(colors.green+'beams file merged : \n'+colors.yellow+beamnpy+colors.black)
        else:
            print(colors.pruple + 'beams already merged' + colors.black)

    def get_beams(self,iBs=[],tol=1e-2,bOpt=''):
        ''' get the beams as recorded during:\n
        - iBs : selected beam indices : default=Ibeams.max()>tol
        - bOpt : O(include Origin),p(print Imax),n or f(force new)
        hk,t,re,im,Ib = beams
        '''
        if 'a' in bOpt : iBs = range(len(self.hk))
        new = 'n' in bOpt or 'f' in bOpt #;print(bOpt)
        if not exists(self._outf('beams')) or new:
            beams = pp.import_beams(self._outf('beamstxt'),self.slice_thick,iBs,tol,'O' in bOpt,'p' in bOpt)
            np.save(self._outf('beams'),beams)
            print(colors.green+'beams file saved : \n'+colors.yellow+self._outf('beams')+colors.black)
        beams = np.load(self._outf('beams'),allow_pickle=True)
        return beams

    def save_patterns(self,force=0,save_opt=1,v=0):
        if force : self.patterns_saved=0
        if not self.patterns_saved:
            print(colors.blue+'...saving patterns to npy...'+colors.black)
            i_files = np.sort([f.replace(self._outf('pattern'),'') for f in lsfiles(self._outf('pattern')+'*')])
            for i in i_files :
                self.save_pattern(i=i,v=v)
            self.patterns_saved=1
            if save_opt:self.save(v=1)
    def save_pattern(self,iz=None,i='',v=1):
        if iz:
            patterns = np.sort(lsfiles(self._outf('pattern')+'*'))
            txt_file = patterns[min(patterns.size-1,iz+1)]
            i=str(iz).zfill(3)
        else:
            txt_file=self._outf('pattern')+i
        if v:print('loading pattern %s' %txt_file)
        im = np.loadtxt(txt_file)
        real,imag = im[:,0:-1:2],im[:,1::2];#print(real.max())
        im = real**2+imag**2
        npy_file=self._outf('pattern').replace('.txt','')+i.replace('.','')+'.npy'
        np.save(npy_file,im,allow_pickle=True)
        if v:print(colors.green+'file saved : ' +colors.yellow+npy_file+colors.black)

    def show_patterns(self,**kwargs):
        patterns = lsfiles(self._outf('pattern').replace('.txt','')+'0*.npy')
        return pp.Multi_Pattern_viewer(self,patterns,figpath=self.datpath,**kwargs)

    # def _get_patterns(self):return
    
    def patterns2gif(self,name=None,v=0,**kwargs):
        if not name:
            self._set_figpath()
            name=self.figpath+self.name+'_pattern.gif'
        self.save_patterns(v=0)

        name = name.replace('.gif','')
        patterns = np.sort(lsfiles(self._outf('pattern').replace('.txt','')+'*.npy'))[1:]
        print(colors.blue+'...saving patterns to png...'+colors.black)
        for iz in range(len(patterns)):
            figname='%s%s.png' %(name,str(iz).zfill(4))
            self.pattern(iz=iz,name=figname,opt='sc',v=v,**kwargs)
        print(colors.blue+'...saving patterns to gif...'+colors.black)
        cmd='im2gif '+name +' png'
        out=check_output(['/bin/bash','-i','-c',cmd]).decode()
        print(colors.green+out+colors.black)

    def set_thicks(self):
        self.dzs = self.i_slice*self.slice_thick
        self.zs = np.arange(0,self.thickness,self.dzs)

    def get_iz(self,thick,v=0):
        self.dzs = self.i_slice*self.slice_thick
        self.zs = np.arange(self.dzs,self.thickness+self.slice_thick,self.dzs)

        iz = np.argmin(abs(self.zs-thick))
        if v==1: print('thick=%.2f, actual thick=%.2f' %(thick,self.zs[iz]))
        if v==2:return iz,self.zs[iz]
        return iz

    def pattern(self,iz=None,file=None,rmax=10,Iopt='Ncs',out=0,tol=1e-6,Nmax=None,gs=3,Imax=3e4,
        rot=0,rings=[],v=1,cmap='binary',pOpt='im',title='',name=None,**kwargs):
        '''Displays the 2D diffraction pattern out of simulation
        - Iopt : t(tiff), c(crop I[r]<tol), n(normalize), s(fftshift), l(logscale), q(quarter only) r(rand) g(good)
        - Nmax : crop display beyond Nmax
        - rings : list or array - of resolutions for rings to display
        - gs    : broadening range
        - rmax  : indicative of noise level fall off
        - rot   : rotate pattern anticlock by rot
        - kwargs : see stddisp
        returns : [qx,qy,I]
        '''
        if isinstance(iz,int):
            npy_files = self._outf('patternnpy').replace('.npy','[0-9]*.npy')
            patterns = np.sort(lsfiles(npy_files))#;print(patterns)
            izs = np.array([p.split('pattern')[-1].replace('.npy','') for p in patterns],dtype=int)
            # print(patterns)
            idx = np.where(izs-iz==0)[0]#;print(idx)
            if idx.size:
                file = patterns[idx[0]]#;print(file)
                zi = self.i_slice*self.slice_thick*(iz+1)#;print(zi)
        if not file:
            file = self._outf('patternnpy')
            if not exists(file):file=self._outf('patternS')#;print('ok')
            zi = self.thickness
        if not title:title = 'z=%d A' %(zi)

        if v:print('loading %s at z=%.1fA' %(file,zi))
        im = np.load(file)
        if v>1:print('original image shape',im.shape)
        ax,by = self.cell_params[:2]
        Nh,Nk = self.repeat[:2];
        Nx,Ny = np.array(np.array(self.NxNy)/2,dtype=int)
        NMAX = 1024
        if Nmax : Nmax = min(Nmax,Nx,Ny)
        # if 'I' in Iopt :
            # real,imag = im[:,0:-1:2],im[:,1::2];#print(real.max())
            # im = real**2+imag**2

            #print('normalizing')
        if 'N' in Iopt:
            im/=(4*Nx*Ny)**2
            if v>1:print(im.max(),im.sum())
        elif 'n' in Iopt : #normalize:
            im00 = im[0,0];im[0,0] = 0
            mMax = im.max();im /= mMax
            if 'l' in Iopt :
                im[0,0] = im00/mMax
            else:
                im[0,0] = 1

        if 'c' in Iopt and not Nmax:
            #print('croping')
            idx = im[:Nx,:Ny]>tol#*im.max()
            h,k = np.meshgrid(np.arange(Nx),np.arange(Ny))
            r = np.sqrt(h**2+k**2)
            Nmax = ceil(r[idx].max()); #print('Pattern Nmax=%d > %E ' %(Nmax,tol))
            Nmax = min(Nmax,NMAX)#;print(Nx,Ny,Nmax)
        if not Nmax : Nmax = min(NMAX,Nx,Ny)

        if 's' in Iopt : im = np.fft.fftshift(im)   #;print('fftshift')

        #print('preparing')
        if 'q' in Iopt:
            im0 = im[Nx:Nx+Nmax,Ny:Ny+Nmax];del(im)
            h,k = np.meshgrid(np.arange(Nmax),np.arange(Nmax))
        else:
            im0 = im[Nx-Nmax:Nx+Nmax,Ny-Nmax:Ny+Nmax];del(im)#;print(im0.shape)
            h,k = np.meshgrid(np.arange(-Nmax,Nmax),np.arange(-Nmax,Nmax))

        if 'g' in Iopt:
            # print(np.where(im0>10*tol))
            i,j = np.where(im0>10*tol) #index of spots
            # print(i,j)
            im00 = im0.copy()
        if 'g' in Iopt:
            gs3 = gs
            dqx,dqy = Nh/ax,Nk/by
            nx,ny = int(np.floor(gs3/dqx)),int(np.floor(gs3/dqy)) ;
            if v>1: print('Gaussian window function size : ', nx,ny)
            ix,iy = np.meshgrid(range(-nx,nx+1),range(-ny,ny+1))
            x,y = ix*dqx,iy*dqy
            Pb = np.exp(-(x**2+y**2)/(gs3/3)**2)#;dsp.stddisp(im=[x,y,Pb],pOpt='im')
            # Pb = 1/(x**2+y**2+0.01)**10
            # print((Pb/Pb.sum()).max())
            for i0,j0 in zip(i,j):
                i0x,j0y = i0+ix,j0+iy #;print(i0x,j0y)
                idx         =  (i0x>=0) & (j0y>=0)  & (i0x<2*Nmax) & (j0y<2*Nmax)
                i0x,j0y     = i0+ix[idx],j0+iy[idx]
                # print(idx)#i0x,j0y)
                im0[i0,j0] -= im00[i0,j0]
                im0[i0x,j0y] += Pb[idx]/Pb[idx].sum()*im00[i0,j0]
                # im0[i0,j0] = 100#im00
                # im0[i0x,j0y] = np.maximum(im0[i0x,j0y],tol*1e-2)

        if 'r' in Iopt :
            r = np.sqrt(h**2+k**2);r[r==0]=1
            im0 += np.random.rand(im0.shape[0],im0.shape[1])/(rmax*r)

        if 'l' in Iopt : #logscale the data
            im0[im0<tol] = tol*1e-2
            im0 = np.log10(im0)

        qx,qy = h/Nh/ax,k/Nk/by
        if rot:
            alpha = np.deg2rad(rot)
            ct,st = np.cos(alpha),np.sin(alpha)
            qx,qy = ct*qx-st*qy,st*qx+ct*qy
        N = [1,4]['q' in Iopt]
        if out : return qx,qy,np.array(im0*Imax,dtype='uint16')

        self._set_figpath()
        if not name:name=self.figpath+basename(file).replace('.npy','.png')
        if 't' in Iopt:
            I = np.array(im0*Imax,dtype='uint16')
            tiff_file = name.replace('.png','.tif')
            tifffile.imwrite(tiff_file,I)
            print(colors.yellow+tiff_file+colors.green+' saved'+colors.black)

        t = np.linspace(0,2*np.pi/N,100)
        ct,st = np.cos(t),np.sin(t)
        plts = [[r*ct,r*st,'g--',''] for r in rings]
        if v:print('displaying pattern:',im0.shape)

        return dsp.stddisp(plts,labs=[r'$q_x(\AA^{-1})$','$q_y(\AA^{-1})$'],im=[qx,qy,im0],
            cmap=cmap,pOpt=pOpt,title=title,name=name,**kwargs)

    def azim_avg(self,tol=1e-6,Iopt='Incsl',out=0,**kwargs):
        ''' Display the average azimuthal diffraction pattern intensities
        - out : get data only
        '''
        qx,qy,I = self.pattern(out=1,opt='',tol=tol,Iopt=Iopt)
        #brute force
        print('brute force averaging')
        qr = np.round(np.sqrt(qx**2+qy**2)*1e10)/1e10
        qr0 = np.unique(qr)#,return_index=1,return_inverse=1)
        I0 = np.array([ I[np.where(qr==q)].mean() for q in qr0])
        if out:
            return qr0,I0
        else:
            plts = [[qr0,I0,'b-','']]
            dsp.stddisp(plts,labs=[r'$q(\AA^{-1})$','$I_q$'],**kwargs)

    ########################################################################
    ###### print functions
    ########################################################################
    def print_datafiles(self,data=None):
        '''show data files *.dat or .xyz depending on mulslice'''
        if not data : data = self.data
        #if not isinstance(data,list) : data = [data]
        self._print_header('%s FILES' %(['.xyz','*.dat'][self.is_mulslice]))
        for dat in data :
            print(colors.green+dat+colors.black)
            with open(self.datpath+dat,'r') as f: print(''.join(f.readlines()))
            print("\n")

    def print_decks(self,decks=None):
        '''show decks *.in'''
        self._print_header('*.in FILES')
        if not decks:decks=list(self.decks.keys())
        if not isinstance(decks,list) : decks = [decks]
        for deck in decks:
            print(colors.green+deck+colors.black)
            try:
                with open(self.datpath+deck,'r') as f: print(''.join(f.readlines()))
            except FileNotFoundError:
                print(self.decks[deck])
                print(colors.red+"WARNING:deck file %s not on disk" %(self.datpath+deck)+colors.black)
    def print_job(self):
        '''print bash script job file '''
        with open(self._outf('job'),'r') as f:
            self._print_header('.sh FILE')
            print(''.join(f.readlines()))
    def print_log(self,head_opt=0):
        '''print the logfile of running multislice .log'''
        try:
            with open(self._outf('log'),'r') as f:
                self._print_header('.log FILE')
                log = f.readlines()
            if head_opt:
                for i,l in enumerate(log):
                    if  'Sorting atoms' in l :break
                print(''.join(log[:i+1]))
            else:
                print(''.join(log))
        except FileNotFoundError:
            print(colors.red+'logfile not created yet, run a simu first'+colors.black)
            return 0
    def log_info(self,v=1):
        with open(self._outf('log'),'r') as f : lines=f.readlines()
        if self.is_mulslice:
            log_lines = [log_lines[-6]] + log_lines[-2:]
            info = [l.strip().split('=')[-1].split(' ')[1] for l in log_lines]
            info = [log_lines[0].strip().split(',')[0].replace('slice','')]+info #zmax
        else:
            state = self.check_simu_state()
            zmax,Imax,cpuT,wallT='0','1','0','0'
            if not state=='init':
                zl = [i for i,l in enumerate(lines) if 'z=' in l ][-1]
                zmax,Imax =[l.strip().split('=')[-1].split(' ')[1] for l in lines[zl:zl+2]]
            if state=='done':
                wl = [i for i,l in enumerate(lines) if 'wall time' in l ][0]
                cpuT,wallT = [l.strip().split('=')[-1].split(' ')[1] for l in lines[wl-1:wl+1]]
        info = np.array([zmax,Imax,cpuT,wallT],dtype=float)
        if v:
            # print(info,log_lines)
            print('zmax=%.1f,Imax=%.4f,cpuT=%.1f,wallT=%.1f' %tuple(info))
        return info

    def wait_simu(self,ssh_alias='',t=1,hostpath=''):
        state = 0
        while not state=='done':
            state=self.check_simu_state(ssh_alias,v=0,hostpath=hostpath)
            print(state)
            time.sleep(t)

    def check_simu_state(self,ssh_alias=None,v=False,hostpath=''):
        '''see completion state of a simulation '''
        if ssh_alias :
            e = self.ssh_get(ssh_alias,'log',hostpath=hostpath)
            if e:
                print(colors.red+e+colors.black)
                if 'No such file' in e:return 'not started'
        try:
            with open(self._outf('log'),'r') as f :
                lines = f.readlines()  #print(self._outf('log')) #f.readlines())
                if len(lines)>=2 :
                    l1,l2=lines[-2:];#print(l1,l2)
                else :
                    return 'empty'
            # get state
            if self.is_mulslice:
                if 'slice ' in l2 :
                    slice       = int(l2.split(',')[0].split(' ')[-1]);
                    tot_slice   = len(self.data)*self.repeat[2]
                    state="%d%%" %(100*slice/tot_slice)
                elif 'elapsed time' in l2 : state='done'
                else : state='init'
            else:
                zl = [i for i,l in enumerate(lines) if 'z=' in l ]
                if not sum(['Sorting atoms' in l for l in lines]) : state='init'
                elif 'END' in l2 : state='done'
                elif sum(['wall time' in l for l in lines]) : state='processing'
                elif len(zl) :
                    l1=lines[zl[-1]]
                    state="%d%%" %int(100*float(l1[3:8])/self.thickness)
                else:
                    state='undefined'
            #print state info
            if v :
                if state=='init': print(colors.green+"Initializing..."+colors.black)
                elif state=='done' : print(colors.green+"Done : \n"+colors.black, l1,l2)
                else : print(colors.green+state.ljust(5)+black )
            return state
        except FileNotFoundError:
            if v : print(colors.red+'logfile not created yet, run a simu first'+colors.black)
            return 0

    ########################################################################
    ######## Private functions
    ########################################################################
    def _get_verbose_options(self,v):
        if isinstance(v,bool) : v=['','nctrdD'][v]
        if isinstance(v,int) : v=''.join(['nctr','DR','d'][:min(v,3)])
        return v
    # def _get_run_options(self,opt,fopt):
    #     save_deck,save_obj,do_run,do_pp = [s in opt for s in 'dsrp']
    #     save_deck |= do_run
    #     if 'f' in opt: fopt += 'f'
    #     if 'w' in opt: fopt += 'w'
    #     vopt=('r' in v and 'R' not in v) + 2*('R' in v)
    #     return save_deck,save_obj,do_run,do_pp,fopt,vopt

    def _get_datafiles(self,data):
        dat_type = ['xyz','dat'][self.is_mulslice]
        err_msg=colors.red+'no *.' +dat_type+' files found in :\n'+colors.yellow+self.datpath+colors.black
        if not data : data = lsfiles(self.datpath+'*.'+dat_type)
        if not data : raise Exception(err_msg)
        if not isinstance(data,list) : data=[data]; #print(data[0].split('.')[-1])
        if not data[0].split('.')[-1]==dat_type : print(err_msg)# raise Exception(err_msg)
        data = [dsp.basename(f) for f in data ]
        return data

    def _set_filenames(self):
        outf = {'obj'       : self.name+'.pkl',
                'simu_deck' : self.name+'.in',
                'log'       : self.name+'.log',
                'image'     : self.name+'.tif',
                'imagesvg'  : self.name+'.svg',
                'pattern'   : self.name+'_pattern.txt',
                'patternnpy': self.name+'_pattern.npy',
                'patternS'  : self.name+'_patternS.npy',
                'patternsvg': self.name+'_pattern.svg',
                'beamstxt'  : self.name+'_beams.txt',
                'beams'     : self.name+'_beams.npy',
                'beamssvg'  : self.name+'_beams.svg',
                'job'       : self.name+'.sh',
                }
        if self.is_mulslice:
            for dat,i in zip(self.data,self.slices):
                outf['slice_data%s' %i] = dat
                outf['slice_imag%s' %i] = dat.replace('.dat','.tif')
                outf['slice_deck%s' %i] = dat.replace('.dat',self.tail+'.in')
        else :
            outf['data'] = self.data[0]
        return outf

    def _set_figpath(self):
        self.figpath = self.datpath+'figures/'
        if not exists(self.figpath):
            Popen('mkdir %s ' %(self.figpath),shell=True)

    def _get_cell_params(self,v=False):
        if self.is_mulslice :
            cz=[]
            for i in list(self.slices):
                with open(self._outf('slice_data%s' %i)) as f:
                    ax,by,cz_i = np.array(f.readline().rstrip().split(' '),dtype=float)
                    cz.append(cz_i)
            cell_params = [ax,by,cz]
        else :
            with open(self._outf('data')) as f:
                f.readline()
                ax,by,cz = np.array(f.readline().rstrip().split(' '),dtype=float)
                cell_params = [ax,by,[cz]]
        if v : print('ax=%.3fA, by=%.3f, cz=' %(ax,by), cz)
        return cell_params

    def _set_tail(self,tail,tag=None):
        if isinstance(tag,str):tail=tag
        tail = '_'+tail if tail else ''
        return tail
    def _set_name(self,v=False):
        name = self.name+self.tail+['_autoslic','_mulslice'][self.is_mulslice]
        if v : print('Simu name pattern = '+name)
        return name
    def _set_NxNy(self,NxNy):
        if isinstance(NxNy,int) or  isinstance(NxNy,np.int64) :
            NxNy = [NxNy,NxNy]
        return tuple(NxNy)
    def _set_hk(self,hk,Nhk,sym=0,hk_pad=None):
        Nh,Nk,Nl=self.repeat
        if hk_pad:
            if isinstance(hk_pad,int):hk_pad=[hk_pad]*2
            Nh,Nk=hk_pad
        if Nhk:
            h,k = np.meshgrid(range(-sym*Nhk,Nhk),range(-sym*Nhk,Nhk))
            hk=[(Nh*h,Nk*k) for h,k in zip(h.flatten(),k.flatten())];#print(hk)
        return hk
    def _set_slice_thick(self,slice_thick,dz=None):
        if dz:slice_thick=dz
        if self.is_mulslice : slice_thick = self.cell_params[2]
        return slice_thick
    def _get_thickness(self,v=True):
        cz = np.array(self.cell_params[2]).sum()
        thickness = self.repeat[2]*cz
        if v : print('simulated thickness = %.3f A, nslices=%d' %(thickness,thickness/self.slice_thick))
        return thickness
    def _set_prev(self,prev):
        zstr = str(int(self.thickness)).zfill(4)
        if isinstance(prev,bool) or isinstance(prev,int):
            if prev:
                prev = self.outf['image']+zstr
        if prev:
            files = ['log', 'patternnpy','patternS','image','beams']
            cmd = ''
            for ext in files:
                file = self._outf(ext)
                old_file = file + zstr
                cmd += "if [ -f %s ]; then cp %s %s;fi\n" %(file,file,old_file)
            # files = list(self.outf.keys());files.remove('obj');files.remove('data');files.remove('beamstxt')
            # for ext in files:
            #     file = self._outf(ext)
            #     cmd += "if [ -f %s ]; then rm %s;fi\n" %(file,file)
            p = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
            p.wait();p.communicate()
        return prev

    ########################################################################
    #### Job
    ########################################################################
    def _get_job(self,temsim=None,cluster=False,datpath=None,patterns_opt=False):
        logfile     = self.outf['log'] #self._outf('log')
        simu_deck   = self.outf['simu_deck'] #self._outf('simu_deck')
        prog        = ['autoslic','mulslice'][self.is_mulslice]
        if not temsim : temsim = temsim_hosts[socket.gethostname()]
        header = self._version_header()
        #write the cmd
        job = "#!/bin/bash\n"
        if cluster:job += "#$ -j y\n#$ -cwd\n#$ -V -w e\n"
        job += '\ncd %s \n' %self.datpath
        job += 'printf "%s" > %s \n\n' %(header,logfile) #overwrite logfile
        if self.is_mulslice:
            for i in self.slices :
                deck = self.outf['slice_deck%s' %i]
                job += 'cat %s | %s >> %s \n' %(deck,temsim+'atompot',logfile)
        job += 'cat %s | %s >> %s\n' %(simu_deck,temsim+prog,logfile)

        #### postprocess
        pyexe='python3'
        if cluster:pyexe='/home/lii26466/anaconda3/bin/python3'
        # import numpy as np;
        # beams = pp.import_beams('%s',%s);
        # np.save('%s',beams)
        # ''' %(self._outf('beamstxt'),self.slice_thick,self._outf('beams'))
        job+='printf "\n\tPOSTPROCESSING\n"  >> %s \n\n' %(logfile)
        pycode='''import multislice.postprocess as pp;import numpy as np;
        multi = pp.load_multi_obj('%s');
        datpath = multi.datpath;
        multi.datpath='./';
        multi.get_beams(bOpt='fa');
        multi.save_pattern();
        qx,qy,It1 = multi.pattern(Iopt='Ncs',out=True,Nmax=260);
        np.save(multi.outf['patternS'],[qx,qy,It1]);
        ''' %(self.outf['obj'])
        if patterns_opt:
            pycode+='''multi.save_patterns(save_opt=0);'''
            pycode+='''multi.datpath = datpath;multi.save();'''
        job +='%s -c "%s" >>%s 2>&1 \n' %(pyexe, pycode.replace('\n',''),logfile)
        job+='printf "END" >>%s\n' %(logfile)

        #save job
        if not datpath:datpath=self.datpath
        with open(datpath+self.outf['job'],'w') as f : f.write(job)
        print(colors.red+colors.yellow+datpath+self.outf['job']+colors.black)

    def _get_job_ssh(self,ssh_alias,hostpath=None,v=False,cluster=False,patterns_opt=False):
        #save local datpath
        datpath  = self.datpath
        hostpath = self._get_hostpath(ssh_alias,hostpath)
        self.datpath = hostpath

        #save updated deck and job
        host    = ssh_hosts[ssh_alias]
        temsim  = temsim_hosts[host]
        self._get_job(temsim,cluster,datpath,patterns_opt=patterns_opt)
        self.make_decks(save=True,datpath=datpath)

        #create directory on remote if does not exist
        cmd = 'ssh %s "if [ ! -d %s ]; then mkdir -p %s;echo %s:%s created;fi"' %(ssh_alias,hostpath,hostpath,ssh_alias,hostpath)
        print(colors.blue+check_output(cmd,shell=True).decode()+colors.black)
        cmd = 'ssh %s "if [ -f %s ]; then  echo 1;fi"' %(ssh_alias,hostpath+self.data[0])
        # dat_exist = 0
        dat_exist = check_output(cmd,shell=True).decode()#;print(hostpath+self.data[0],dat_exist)
        #copy files over to remote
        cmd  = 'cd %s && scp ' %(datpath)
        if not dat_exist:cmd += '%s ' %(' '.join(self.data))
        cmd += '%s ' %(' '.join(list(self.decks.keys())))
        cmd += '%s ' %(self.outf['job'])
        cmd += '%s ' %(self.outf['obj'])
        cmd += '%s:%s' %(ssh_alias,hostpath)
        # check_output(cmd,shell=True).decode()
        p = Popen(cmd,shell=True);p.wait()

        # submit job command (qsub for clusters)
        if cluster:
            cmd = 'source /etc/bashrc && qsub -S /bin/bash '
        else :
            cmd = 'bash '
        cmd += '%s' %(hostpath+self.outf['job'])
        cmd = 'ssh %s "%s"' %(ssh_alias,cmd)
        #restore
        self.datpath = datpath
        self.make_decks(save=True)
        return cmd

    def ssh_get(self,ssh_alias,file,hostpath=None,dest_path=None):
        if not dest_path : dest_path = self.datpath
        hostpath = self._get_hostpath(ssh_alias,hostpath)
        cmd = 'scp %s:%s %s' %(ssh_alias,hostpath+self.outf[file],dest_path)
        p = Popen(cmd,shell=True,stderr=PIPE)#stdout=PIPE

        p.wait()
        o,e = p.communicate()
        # print(o.decode())
        return e.decode()

    ########################################################################
    #### decks
    ########################################################################
    def _atompot_deck(self,i):
        # TODO (prevent mulslice from adding the extension automatically for this case )
        dat_file = self._outf('slice_data%s' %i)
        img_file = self._outf('slice_imag%s' %i).replace('.tif','')
        deck  = "%s\n" %dat_file                #.dat file
        deck += "%s\n" %img_file                #.tif file
        deck += "n\n"                           #record fft proj
        deck += "%d %d\n" %(self.NxNy)          #sampling
        deck += "%d %d 1\n" %(self.repeat[:2])  #super cell
        deck += "n\n"                           #thermal displacement
        return deck

    def _mulslice_deck(self,prev=None):
        prev_run = "%s\n" %['y','n'][prev==None]
        if prev : prev_run += "%s\n" %(self.datpath+prev)
        Nz,sl = self.repeat[2],ascii_letters[:len(self.data)]
        deck = "%d(%s)\n" %(Nz,sl)              #Stack sequence
        for i in self.slices :                  #atom pots
            deck+="%s\n" %self._outf('slice_imag%s' %i)#
        deck += "%s\n" %(self._outf('image'))   #image file
        deck += "n\n"                           #partial coherence
        deck += prev_run                        #start previous run
        deck += "%.2f\n" %self.keV              #wavelength
        deck += "0 0\n"                         #crystal tilt
        deck += "%.4f %.4f \n" %self.tilt       #beam tilt
        deck += "y\n"                           #record beams
        deck += "%s\n" %(self._outf('beamstxt'))#     filename
        deck += "%d\n" %len(self.hk)            #     nbeams
        for hk in self.hk :                     #     beam indices
            deck += "%d %d\n" %(hk)
        #deck += "n\n"                           #record cross sections
        #deck +=
        return deck

    def _autoslic_deck(self,prev=None):
        prev_run = "%s\n" %['y','n'][prev==None]
        if prev : prev_run += "%s\n" %(self.datpath+prev)
        #Deck
        deck  = "%s\n" %self._outf('data')      #.xyz file
        deck += "%d %d %d\n" %self.repeat       #super cell
        deck += "%s\n" %self._outf('image')     #image file
        deck += "n\n"                           #partial coherence
        deck += prev_run                        #start previous run
        if not prev:
            deck += "%.4f\n" %self.keV              #wavelength
            deck += "%d %d\n" %self.NxNy            #sampling
        deck += "%.4f %.4f \n" %self.tilt       #beam tilt
        deck += "%f\n" %(self.slice_thick)      #slice thickness
        deck += "y\n"                           #record beams
        deck += "%s\n" %self._outf('beamstxt')  #filename
        deck += "%d\n" %len(self.hk)            #   nbeams
        for hk in self.hk :                     #   - beams -
            deck += "%d %d\n" %(hk)             #   beam indices
        deck += "%s\n" %(['n','y'][self.TDS])   #thermal vibration?
        if self.TDS:                            #   - Params -
            deck+='%.3f\n' %self.T              #   temperature
            deck+="%d\n" %self.n_TDS            #   nb configurations
        deck += "n\n"                           #intensity cross section
        deck += "y\n"                           #diffraction pattern
        deck += "%s\n" %self._outf('pattern')   #   filename
        deck += "%d" %self.i_slice;
        return deck

    ########################################################################
    #### misc private
    def _version_header(self):
        header = '''
------------------------------------------------------------------------------
------------------------------------------------------------------------------
This header has been produced by multislice library
version : %s
date : %s
author : Tarik Ronan Drevon
e-mail : tarik.drevon@stfc.ac.uk
------------------------------------------------------------------------------
------------------------------------------------------------------------------
''' %(self.version,time.ctime())
        return header
    def _print_header(self,msg,w=70):
        head='#'*w
        print(colors.blue+head)
        print('\t\t\t' +msg+ ' :')
        print(head+colors.black)

    def _get_hostpath(self,ssh_alias,hostpath):
        if not hostpath :
            local_hostpath  = path_hosts[socket.gethostname()]
            remote_hostpath = path_hosts[ssh_hosts[ssh_alias]]
            hostpath  = self.datpath.replace(local_hostpath,remote_hostpath)
            # simu_folder = self.datpath.split('/')[-2]
            # hostpath += simu_folder+'/'
            # hostpath = self.datpath.replace('ronan','tarik').replace('CCP4','git/ccp4')
        return hostpath

    def _outf(self,file):
        return self.datpath+self.outf[file]

#########################################################################
##### utilities
#########################################################################
def get_tilts(tx=np.arange(0,10,0.05),ty=0):
    def convert_tilts(tilts,ntilts):
        if isinstance(tilts,float) or isinstance(tilts,int):
            tilts = [tilts]*ntilts
        return np.deg2rad(tilts)*1000

    if isinstance(tx,list) : tx=np.array(tx)
    if isinstance(ty,list) : ty=np.array(ty)
    nts = 0
    if isinstance(tx,np.ndarray):
        nts=tx.size
    if isinstance(ty,np.ndarray):
        if nts:
            if not nts==ty.size:
                raise Exception('tx and ty have different sizes')
        nts = ty.size
    if not nts:
        raise Exception('at list one array of tilts must be provided ')
    txs = convert_tilts(tx,nts)
    tys = convert_tilts(ty,nts)
    tilts = [[tx,ty] for tx,ty in zip(txs,tys)]
    return tilts

class Rocking:
    def __init__(self,name,tx=np.arange(0,10,0.05),ty=0,tag='',**kwargs):
        ''' simulate rocking curve
        - tx : tilt parameters around x(degrees)
            - float - constant value
            - list or np.ndarray : actual list of tilts
        - ty : tilt parameters around y(same behaviour as tx)
        - tag : tag will then be '<tag>_tilt<nb>'
        '''
        if tag:tag+='_'
        self.path = name
        self.tag  = tag
        self.df_path = self.path+self.tag+'tilts.pkl'
        self.tx,self.ty = tx,ty
        self.tilts = get_tilts(tx,ty)
        self.df = sweep_var(name,param='tilt',vals=self.tilts,tail=self.tag,df=self.tag+'tilts.pkl',
            **kwargs)
        self.save(v=1)

    def save(self,v=0):
        '''save this object'''
        file = self.path+self.tag+'rock.pkl'
        with open(file,'wb') as out :
            pickle.dump(self, out, pickle.HIGHEST_PROTOCOL)
        if v:print(colors.green+"object saved\n"+colors.yellow+file+colors.black)


    #############################################################################
    # utilities
    #############################################################################
    def load(self,i):
        return pp.load_multi_obj(self.path+self.df.index[i])

    def update(self,files=[],v=1):
        df = pp.update_df_info(self.df_path,files)
        if v:
            np.set_printoptions(precision=3)
            pd.set_option("precision", 2)
            print(df[['tilt','state','zmax(A)','Inorm']])

    def plot_rocking(self,iBs,iZs=-1,zs=None,**kwargs):
        '''plot rocking curve for beams iBs at thickness zs
        - iBs : list - beam indices [1,2,5] or ids [(1,0),(1,1)]
        - iZs : int or list - slice indices (last slice by default )
        - zs : float or list - selected thickness (in A) to show(takes preference over iZs if set)
        '''
        nbs,nzs,iZs = self._init_rocking(iBs,iZs,zs)

        I = np.zeros((len(self.tilts),nbs,nzs))
        for i,pkl in enumerate(self.df.index):
            multi = pp.load_multi_obj(self.path+pkl)
            hk,z,re,im,Ib = multi.beam_vs_thickness(bOpt='o',iBs=iBs)           #;print(Ib.shape)
            I[i] = np.array(Ib[:,iZs])

        ts = self.tx
        cs,ms,plts = dsp.getCs('jet',nbs), dsp.markers,[]
        legElt = { '%s' %hk[i]:[cs[i],'-'] for i,iB in enumerate(iBs)}
        for iz,iZ in enumerate(z[iZs]):
            legElt.update({'$z=%d A$' %(iZ):['k',ms[iz]+'-']})
            plts += [[ts,I[:,i,iz],[cs[i],ms[iz]+'-'],''] for i,iB in enumerate(iBs)]
        dsp.stddisp(plts,labs=[r'$\theta$(deg)','$I$'],legElt=legElt,**kwargs)

    #############################################################################
    #### misc
    #############################################################################
    def _init_rocking(self,iBs,iZs,zs):
        multi = pp.load_multi_obj(self.path+self.df.index[0])
        hk,z,re,im,Ib = multi.beam_vs_thickness(bOpt='o',iBs=iBs)
        if isinstance(zs,float) or isinstance(zs,float):zs = [zs]
        if isinstance(zs,list) or isinstance(zs,np.ndarray):
            iZs = [np.argmin(abs(z-z0)) for z0 in zs]
        if isinstance(iZs,int):iZs=[iZs]
        nbs = np.array(hk).size
        nzs = z[iZs].size
        return nbs,nzs,iZs

def sweep_var(name,param,vals,df=1,ssh='',tail='',do_prev=0,**kwargs):
    '''runs a set of similar simulations with one varying parameter
    - name          : path to the simulation folder
    - param,vals    : the parameters and values to sweep
    - df            :
        - pd.Dataframe to update(since parsed as a reference)
        - int create and save the new dataframe if 1
    - do_prev       : Used for iterative fourier transform
    - kwargs : see help(Multislice)
    '''
    do_df,save = isinstance(df,pd.core.frame.DataFrame),0
    dfname = 'df.pkl'
    if isinstance(df,str):dfname,df=df,1
    if isinstance(df,int):
        if df : df,do_df,save = pd.DataFrame(columns=[param,'host','state']+pp.info_cols),1,1
    nvals,prev = len(vals),None
    for i,val in zip(range(nvals),vals):
        print(colors.red+param+':',val,colors.black)
        kwargs[param]=val
        if do_prev and i: prev = multi.outf['image']
        multi=Multislice(name,prev=prev,
            ssh=ssh,tail=tail+param+str(i).zfill(ceil(nvals/10)),
            **kwargs)
        if do_df:
            df.loc[multi.outf['obj']] = [nan]*len(df.columns)
            df.loc[multi.outf['obj']][[param,'host','state']] = [val,ssh,'start']
    if save :
        df.to_pickle(name+dfname)
        print(colors.green+'Dataframe saved : '+colors.yellow+name+dfname+colors.black)
    return df
