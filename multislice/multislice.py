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
import pickle,socket,time
import pandas as pd
import numpy as np
from matplotlib import rc
from math import ceil,nan
from subprocess import Popen,check_output,PIPE
from glob import glob as lsfiles
from os.path import dirname,realpath,exists
from string import ascii_uppercase as SLICES #ascii_letters
from string import ascii_letters
# from matplotlib.animation import FuncAnimation
# from utils.glob_colors import*
from crystals import Crystal
import utils.glob_colors as colors
import utils.displayStandards as dsp
from scattering.structure_factor import structure_factor3D
from . import postprocess as pp;imp.reload(pp)

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
    'badb.rc-harwell.ac.uk':'/data3/lii26466/multislice/',
}


class Multislice:
    ''' **DATA PARAMETERS :**\n
    - `name` : path to the simulation folder
    - `tail` : trailing string in the naming convention of outputs
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
    - `slice_thick` : float - slice thickness (A)
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
    - `ppopt` : u(update),w(wait), I(image) B(beam) P(pattern) A(azim_avg)
    - `ssh` : ip address or alias of the machine on which to run the job
    - `hostpath` : path to data on host
    - `cif_file` : name of .cif file corresponding to structure in path
    '''
    def __init__(self,name,mulslice=False,tail='',data=None,
                tilt=[0,0],TDS=False,T=300,n_TDS=16,
                keV=200,repeat=[2,2,1],NxNy=[512,512],slice_thick=1.0,i_slice=100,
                hk=[(0,0)],Nhk=0,hk_sym=0,prev=None,
                opt='sr',fopt='',ppopt='uwP',v=1,
                ssh=None,hostpath='',cluster=False,cif_file=None):
        #output options
        if isinstance(v,bool) : v=['','nctrdD'][v]
        if isinstance(v,int) : v=''.join(['nctr','DR','d'][:min(v,3)])
        save_deck,save_obj,do_run,do_pp = [s in opt for s in 'dsrp']
        save_deck |= do_run
        if 'f' in opt: fopt += 'f'
        if 'w' in opt: fopt += 'w'
        vopt=('r' in v and 'R' not in v) + 2*('R' in v)

        #attributes
        self.version     = '1.4.1'
        self.name        = dsp.basename(name)                       #;print('basename:',self.name)
        self.datpath     = realpath(name)+'/'                       #;print('path:',self.datpath)
        self.cif_file    = self.get_cif_file(cif_file)              #;print(self.cif_file)
        self.is_mulslice = mulslice
        self.tilt        = tuple(tilt)
        self.TDS         = TDS
        self.T           = T
        self.n_TDS       = n_TDS
        self.tail        = '_'+tail if tail else ''
        self.name        = self._set_name('n' in v)
        self.data        = self._get_datafiles(data)
        self.slices      = SLICES[:len(self.data)]
        self.outf        = self._set_filenames()
        self.keV         = keV
        self.repeat      = tuple(repeat)
        self.NxNy        = self._set_NxNy(NxNy)
        self.hk          = self._set_hk(hk,Nhk,hk_sym)
        self.cell_params = self._get_cell_params('c' in v)
        self.slice_thick = self._set_slice_thick(slice_thick)
        self.i_slice     = i_slice
        self.thickness   = self._get_thickness('t' in v)
        self.decks       = {}
        self.p           = None

        ## make decks and run if required
        if save_deck : self.make_decks(save_deck,prev=prev,v=v)
        if save_obj : self.save(v=v)
        if 'd' in v : self.print_datafiles()
        if 'D' in v : self.print_decks()
        if do_run : self.p = self.run(v=vopt, fopt=fopt,ssh_alias=ssh,hostpath=hostpath,cluster=cluster)
        if do_pp : self.postprocess(ppopt,ssh,hostpath=hostpath)

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

    def run(self,v=1,fopt='w',ssh_alias=None,hostpath=None,cluster=False):
        '''run the simulation with temsim
        - fopt : f(force rerun), w(warn ask rerun already done)
        - ssh : name of the host to run the job
        '''
        if isinstance(fopt,int):fopt='f'
        run = True
        if self.check_simu_state(v=0,ssh_alias='') == 'done' :
            if v>0 : print(colors.red+"Simulation already performed in the past."+colors.black)
            if 'f' in fopt :
                msg='Force re-running.'
            else :
                if 'w' in fopt :
                    run = input(colors.red+"Re-run?(y/n) : "+colors.black)=='y'
                    msg = 're-running'
                else :
                    run,msg = False, 'not running'
            if v>1 : print(colors.red+msg+colors.black)
        p = None
        if run:
            if ssh_alias :
                if ssh_alias=='badb':cluster=1
                cmd = self._get_job_ssh(ssh_alias,hostpath=hostpath,v=v>2,cluster=cluster)
                #time.sleep(1)
                #self.check_simu_state(v=0,ssh_alias=ssh_alias,hostpath=hostpath)
            else :
                self._get_job(cluster=cluster)
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
        if not figpath : figpath = self.datpath
        time.sleep(0.1)
        state = self.check_simu_state(ssh_alias=ssh_alias,hostpath=hostpath);print(state)
        if not state == 'done' and 'w' in ppopt:
            self.wait_simu(ssh_alias,10,hostpath)
            self.p.wait()
        # self.get_beams(iBs='a',bOpt='f')
        print(colors.blue+'...postprocessing...'+colors.black)
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
        - sf_args : see (structure_factor3D)
        returns :
        - (qx,qy,qz),Fhkl
        '''
        crys = Crystal.from_cif(self.cif_file)
        pattern = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in crys.atoms] )
        lat_vec = np.array(crys.reciprocal_vectors)
        (h,k,l),Fhkl = structure_factor3D(pattern, lat_vec, **sf_args)
        qx = h/crys.lattice_parameters[0]
        qy = k/crys.lattice_parameters[1]
        qz = l/crys.lattice_parameters[2]
        return (qx,qy,qz),Fhkl

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

    def beam_vs_thickness(self,bOpt='',iBs=[],tol=1e-2,**kwargs):
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

    def save_patterns(self):
        i_files = [f.replace(self._outf('pattern'),'') for f in lsfiles(self._outf('pattern')+'*')]
        for i in i_files :
            self.save_pattern(i)
    def show_patterns(self,**kwargs):
        patterns = lsfiles(self._outf('pattern').replace('.txt','')+'*.npy')
        return Pattern_viewer(self,patterns,figpath=self.datpath,**kwargs)

    def save_pattern(self,i=''):
        print('loading pattern %s' %i)
        txt_file=self._outf('pattern')+i
        im = np.loadtxt(txt_file)
        print('saving')
        real,imag = im[:,0:-1:2],im[:,1::2];#print(real.max())
        im = real**2+imag**2
        #im = np.fft.fftshift(im)
        npy_file=self._outf('pattern').replace('.txt','')+i.replace('.','')+'.npy'
        np.save(npy_file,im,allow_pickle=True)
        print(colors.green+'file saved : ' +colors.yellow+npy_file+colors.black)



    def pattern(self,file='',Iopt='Incsl',out=0,tol=1e-6,qmax=None,Nmax=None,gs=3,rings=[],v=1,**kwargs):
        '''Displays the 2D diffraction pattern out of simulation
        - Iopt : I(intensity), c(crop I[r]<tol), n(normalize), s(fftshift), l(logscale), q(quarter only) r(rand) g(good)
        - Nmax : crop display beyond Nmax
        - rings : list or array - of resolutions for rings to display
        - gs     : broadening parameter
        - kwargs : see stddisp
        returns : [qx,qy,I]
        '''
        # if not exists(self._outf('patternnpy')):self.save_pattern()
        if not file : file = self._outf('patternnpy')
        im = np.load(file)
        if v>1:print('original image shape',im.shape)
        ax,by = self.cell_params[:2]
        Nh,Nk = self.repeat[:2];
        Nx,Ny = np.array(np.array(self.NxNy)/2,dtype=int)
        if Nmax : Nmax = min(Nmax,Nx,Ny)
        # if 'I' in Iopt :
            # real,imag = im[:,0:-1:2],im[:,1::2];#print(real.max())
            # im = real**2+imag**2
        if 'n' in Iopt : #normalize
            #print('normalizing')
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
            Nmax = min(Nmax,256)#;print(Nx,Ny,Nmax)
        if not Nmax : Nmax = min(256,Nx,Ny)

        if 's' in Iopt : im = np.fft.fftshift(im)   #;print('fftshift')

        #print('preparing')
        if 'q' in Iopt:
            im0 = im[Nx:Nx+Nmax,Ny:Ny+Nmax];del(im)
            h,k = np.meshgrid(np.arange(Nmax),np.arange(Nmax))
        else:
            im0 = im[Nx-Nmax:Nx+Nmax,Ny-Nmax:Ny+Nmax];del(im)#;print(im0.shape)
            h,k = np.meshgrid(np.arange(-Nmax,Nmax),np.arange(-Nmax,Nmax))

        if 'r' in Iopt :
            r = np.sqrt(h**2+k**2);r[r==0]=1
            im0 += np.random.rand(im0.shape[0],im0.shape[1])/(10*r)
        if 'l' in Iopt : #logscale the data
            if 'g' in Iopt:
                i,j = np.where(im0>10*tol) #index of spots
            im0[im0<tol] = tol*1e-2
            if 'g' in Iopt:
                x,y = np.meshgrid(range(-5,6),range(-5,6))
                Pb = np.exp(-gs*(x**2+y**2)**2)
                # Pb = 1/(x**2+y**2+0.01)**10
                for i0,j0 in zip(i,j):
                    im00 = im0[i0,j0]
                    im0[i0+x,j0+y] += im0[i0,j0]*Pb
                    im0[i0,j0] = im00
                    im0[i0+x,j0+y] = np.maximum(im0[i0+x,j0+y],tol*1e-2)
            im0 = np.log10(im0)

        if out : return h/Nh/ax, k/Nk/by, im0
        N = [1,4]['q' in Iopt]
        t = np.linspace(0,2*np.pi/N,100)
        ct,st = np.cos(t),np.sin(t)
        plts = [[r*ct,r*st,'g--',''] for r in rings]
        if v:print('displaying pattern:',im0.shape)
        return dsp.stddisp(plts,labs=[r'$q_x(\AA^{-1})$','$q_y(\AA^{-1})$'],im=[h/Nh/ax,k/Nk/by,im0],**kwargs)

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
                if not sum(['Sorting' in l for l in lines]) : state='init'
                elif 'z=' == l1[:2] : state="%d%%" %int(100*float(l1[3:8])/self.thickness)
                elif 'END' in l2 : state='done'
                elif sum(['wall time' in l for l in lines]) : state='processing'
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
    def _set_name(self,v=False):
        name = self.name+self.tail+['_autoslic','_mulslice'][self.is_mulslice]
        if v : print('Simu name pattern = '+name)
        return name
    def _set_NxNy(self,NxNy):
        if isinstance(NxNy,int) or  isinstance(NxNy,np.int64) :
            NxNy = [NxNy,NxNy]
        return tuple(NxNy)
    def _set_hk(self,hk,Nhk,sym=0):
        if Nhk:
            Nh,Nk,Nl=self.repeat
            h,k = np.meshgrid(range(-sym*Nhk,Nhk),range(-sym*Nhk,Nhk))
            hk=[(Nh*h,Nk*k) for h,k in zip(h.flatten(),k.flatten())];#print(hk)
        return hk
    def _set_slice_thick(self,slice_thick):
        if self.is_mulslice : slice_thick = self.cell_params[2]
        return slice_thick
    def _get_thickness(self,v=True):
        cz = np.array(self.cell_params[2]).sum()
        thickness = self.repeat[2]*cz
        if v : print('thickness = %.3f A'%thickness)
        return thickness

    ########################################################################
    #### Job
    ########################################################################
    def _get_job(self,temsim=None,cluster=False,datpath=None):
        logfile     = self.outf['log'] #self._outf('log')
        simu_deck   = self.outf['simu_deck'] #self._outf('simu_deck')
        prog        = ['autoslic','mulslice'][self.is_mulslice]
        if not temsim : temsim = temsim_hosts[socket.gethostname()]
        header = self._version_header()
        #write the cmd
        job = "#!/bin/bash\n"
        if cluster:job += "#$ -j y\n#$ -cwd\n#$ -V -w e\n"
        job += '\ncd %s \n' %self.datpath
        job += 'printf "%s" > %s \n' %(header,logfile) #overwrite logfile
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
        job+='printf "\n\n\t\t POSTPROCESSING\n"  >> %s' %(logfile)
        pycode='''import multislice.postprocess as pp;import numpy as np;
        multi = pp.load_multi_obj('%s');
        multi.datpath='./';
        multi.get_beams(bOpt='fa');
        multi.save_pattern();
        qx,qy,It1 = multi.pattern(Iopt='Ics',out=True,Nmax=260);
        np.save(multi.outf['patternS'],[qx,qy,It1]);
        ''' %(self.outf['obj'])
        job +='%s -c "%s" >>%s 2>&1 ' %(pyexe, pycode.replace('\n',''),logfile)
        job+='END\n'

        #save job
        if not datpath:datpath=self.datpath
        with open(datpath+self.outf['job'],'w') as f : f.write(job)
        print(colors.red+colors.yellow+datpath+self.outf['job']+colors.black)

    def _get_job_ssh(self,ssh_alias,hostpath=None,v=False,cluster=False):
        #save local datpath
        datpath  = self.datpath
        hostpath = self._get_hostpath(ssh_alias,hostpath)
        self.datpath = hostpath

        #save updated deck and job
        host    = ssh_hosts[ssh_alias]
        temsim  = temsim_hosts[host]
        self._get_job(temsim,cluster,datpath)
        self.make_decks(save=True,datpath=datpath)

        #create directory on remote if does not exist
        cmd = 'ssh %s "if [ ! -d %s ]; then mkdir %s;echo %s:%s created;fi"' %(ssh_alias,hostpath,hostpath,ssh_alias,hostpath)
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
            hostpath = path_hosts[ssh_hosts[ssh_alias]]
            simu_folder = self.datpath.split('/')[-2]
            hostpath += simu_folder+'/'
            # hostpath = self.datpath.replace('ronan','tarik').replace('CCP4','git/ccp4')
        return hostpath

    def _outf(self,file):
        return self.datpath+self.outf[file]

#########################################################################
##### utilities
#########################################################################
def sweep_var(name,param,vals,df=1,ssh='',tail='',do_prev=0,**kwargs):
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



class Pattern_viewer():
    def __init__(self,multi,patterns,figpath,i=0,**args):
        ''' View cbf files
        - exp_path : path to images
        - figpath : place to save the figures
        - i : starting image
        '''
        self.multi = multi
        self.figpath = figpath
        self.patterns  = np.sort(patterns);#print(patterns)
        self.args=args
        self.nfigs = self.patterns.size
        self.fig,self.ax = dsp.stddisp()
        cid = self.fig.canvas.mpl_connect('key_press_event', self)
        self.i=i     #starting image
        self.inc=1   #increment(use 'p' or 'm' to change)
        self.mode=1
        rc('text', usetex=False)
        self.update()

    def update(self):
        self.ax.cla()
        print("%d/%d" %(self.i,self.nfigs))
        self.multi.pattern(fig=self.fig,ax=self.ax,file=self.patterns[self.i],
            title='pattern %d' %self.i,**self.args,opt='')
        self.fig.canvas.draw()

    def __call__(self, event):
        # print(event.key)
        if event.key in ['up','right']:
            self.i=min(self.i+self.inc,self.nfigs-1)
            self.mode=1
            self.update()
        elif event.key in ['left','down']:
            self.i=max(0,self.i-self.inc)
            self.mode=-1
            self.update()

        if event.key=='s':
            dsp.saveFig(self.figpath+'pattern%s.png' %str(self.i).zfill(3),ax=self.ax)

        if event.key=='p':
            self.inc=min(self.inc+1,100);print(self.inc)
        if event.key=='m':
            self.inc=max(1,self.inc-1);print(self.inc)

# def coords2grid(coords,cz):
#     '''arrange coordinates so they fit on periodic grid
#     Rotates the crystal around z so x and y are positive
#     - z < 0 => z_cz
#     '''
#     x,y,z = coords.T
#     Rz = lambda t : np.array([cos(t),sint(t),0],[-sint(t),cos(t),0],[0,0,1]])
#     if x[0]<0 & y[0]>0 : coords = Rz(-pi/2).dot(coords.T).T
#     if x[0]>0 & y[0]<0 : coords = Rz(pi/2).dot(coords.T).T
#     if x[0]<0 & y[0]<0 : coords = Rz(pi).dot(coords.T).T
#     idz = coords[:,2]<0
#     coords[idz] = coords[idz]+cz
#     return coords
########################################################################
#def : test
########################################################################
def test_base(name,**kwargs):
    multi=Multislice(name,keV=200,
        repeat=[2,2,5],NxNy=128,slice_thick=1.3575,Nhk=3,
        **kwargs)
    return multi

if __name__ == '__main__':
    name  = 'dat/test/'
    #multi = test_base(name,mulslice=False,opt='dsrp',fopt='',ppopt='I',ssh='',v=1)
    #multi = test_base(name,mulslice=False,fopt='f',opt='dsr',ppopt='',ssh='tarik-CCP4home',v='nctrdDR')
    # multi = test_base(name,mulslice=True,ssh=None,fopt='f',opt='dsr',v='nctrdDR')
    # multi = test_base(name,tail='TDS',mulslice=False,TDS=True,T=300,fopt='f',ppopt='wP',opt='dsrp',v='nctrdDR')
