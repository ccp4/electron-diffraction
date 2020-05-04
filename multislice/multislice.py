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
from utils import*
from subprocess import Popen #,check_output
from glob import glob as lsfiles
from os.path import dirname,realpath,exists
from postprocess import plot_beam_thickness,import_beams
import pickle

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
    - `mulslice` : True(mulslice)=>periodic, False(autoslic)
    - `keV` : float - wavelength in keV
    - `repeat` : 3-int-list - super cell size in the x,y,z directions
    - `NxNy` : 2-int-list or int : sampling in x and y (same sampling in x and y if only int provided)
    - `slice_thick` : float - slice thickness (A)
    - `hk` : list of tuple - beams to record as function of depth
    - `Nhk`: int - set hk on a grid of multiples of the supercell size. Prevails over hk argument if >0. default(0)
    \n**OUTPUT options:**\n
    - `v`  : verbose - str or bool(all if True) or int(verbose level)
        - n(naming pattern),c(cell params),t(thickness),r(run cmd)
        - d(data),D(Decks),R(full run cmd)
    - `save_deck` : save the decks to disk
    '''
    def __init__(self,name,mulslice,tail='',data=None,
                keV=200,repeat=[2,2,1],NxNy=[512,512],slice_thick=1.0,
                hk=[(0,0)],Nhk=0,
                save_deck=True,save_obj=False,run=False,opt=None,
                v=True):
        if isinstance(v,bool) : v=['','nctrdD'][v]
        if isinstance(v,int) : v=''.join(['nct','D','dR'][:min(v,3)])
        if opt : save_deck,save_obj,run = [s in opt for s in ['dsr']]
        self.name        = basename(name)                           #;print('basename:',self.name)
        self.datpath     = realpath(name)+'/'                       #;print('path:',self.datpath)
        self.is_mulslice = mulslice
        self.tail        = '_'+tail if tail else ''
        self.data        = self._get_datafiles(data,'d' in v)       #;print(self.data)
        self.cell_params = self._get_cell_params('c' in v)
        self.keV         = keV
        self.repeat      = tuple(repeat)
        self.NxNy        = self._set_NxNy(NxNy)
        self.hk          = self._set_hk(hk,Nhk)
        self.slice_thick = self._set_slice_thick(slice_thick)
        self.thickness   = self._get_thickness('t' in v)
        self.exe         = get_figpath(__file__,'/')+'run_temsim.sh'#;print(self.exe)
        self.name        = self._set_name('n' in v)
        self.outf        = {
            'log'     : self.datpath+self.name+'.log',
            'pattern' : self.datpath+self.name+'.tif',
            'beamstxt': self.datpath+self.name+'_beams.txt',
            'beams'   : self.datpath+self.name+'_beams.npy',
            }
        ## start doing things here
        self.decks  = self.make_decks(save_deck, 'D' in v)
        if save_obj : self.save()
        if run : self.run(v=[1,2]['r' in v,'R' in v])
    ########################################################################
    ##### Public functions
    ########################################################################
    def make_decks(self,save=True,v=True):
        '''create the decks from the information provided'''
        if self.is_mulslice:
            decks={dat.replace('.dat',self.tail+'.in') : self._atompot_deck(dat) for dat in self.data} #atompot first
            decks[self.name+'.in']=self._mulslice_deck()
        else:
            decks={self.name+'.in':self._autoslic_deck()}
        if save:
            print(green+'Decks saved :'+black)
            for filename,deck in decks.items():
                with open(self.datpath+filename,'w') as f : f.write(deck)
                print(yellow+self.datpath+filename+black)
        if v :
            self.decks=decks
            self.print_datafiles()
        return decks

    def run(self,v=1):
        '''run the simulation with temsim'''
        temsim='/home/ronan/Documents/CCP4/src/multislice/kirkland/temsim/'
        #very ugly sorry
        if self.is_mulslice:
            atompot_decks = [self.datpath+d for d in self.decks.keys() if 'mulslice' not in d ]
            simu_deck = [self.datpath+d for d in self.decks.keys() if 'mulslice' in d ][0]
            cmd='> %s && ' %(self.outf['log'])
            for deck in atompot_decks :
                cmd += '''cat %s | %s >> %s &&
                    ''' %(deck,temsim+'atompot',self.outf['log'])
            cmd += '''cat %s | %s >> %s
                ''' %(simu_deck,temsim+'mulslice',self.outf['log'])
        else:
            simu_deck = self.datpath+list(self.decks.keys())[0]
            cmd = '''cat %s | %s > %s
                ''' %(simu_deck,temsim+'autoslic',self.outf['log'])
        # cmd = '''
        #     %s -d %s -s %s -v > %s
        #     '''%(self.exe,deck_list,prog,self.outf['log'])
        if v>0:print(green+self.name+"job submitted"+black)
        if v>1:print(red+cmd+black)
        p = Popen(cmd,shell=True)
        return p

    def get_beams(self,iBs=[],tol=1e-2,bOpt=''):
        ''' get the beams as recorded during:\n
        - iBs : selected beam indices : default=Ibeams.max()>tol
        - bOpt : O(include Origin),p(print Imax),n or f(force new)
        hk,t,re,im,Ib = beams
        '''
        if 'a' in bOpt : iBs = range(len(self.hk))
        new = 'n' in bOpt or 'f' in bOpt; #print(bOpt)
        if not exists(self.outf['beams']) or new:
            beams = import_beams(self.outf['beamstxt'],self.slice_thick,iBs,tol,'O' in bOpt,'p' in bOpt)
            np.save(self.outf['beams'],beams)
            print(green+'beams file saved\n'+yellow+self.outf['beams']+black)
        beams = np.load(self.outf['beams'])
        return beams

    def beam_vs_thickness(self,bOpt='',iBs=[],tol=1e-2,**kwargs):
        '''
        - bOpt : O(include Origin),p(print Imax),n or f(force new)
        - kwargs : see help(plot_beam_thickness)
        '''
        beams = self.get_beams(iBs,tol,bOpt)
        plot_beam_thickness(beams,**kwargs)

    def pattern(self,opt='I',cmap='jet',**kwargs):
        '''Displays the 2D pattern out of simulation
        - opt : I(intensity)
        '''
        im=plt.imread(self.outf['pattern'])
        if 'I' in opt :
            N = int(im.shape[1]/2)
            real,imag = im[:,:N],im[:,N:];#print(real.max())
            im = real**2+imag**2
            im =im/im.max()
        stddisp(labs=['$x$','$y$'],im=im,legOpt=0,imOpt='c',**kwargs)

    def save(self):
        '''save this multislice object'''
        file=self.datpath+self.name +'.pkl'
        with open(file,'wb') as out :
            pickle.dump(self, out, pickle.HIGHEST_PROTOCOL)
        print(green+"object saved\n"+yellow+file+black)

    ########################################################################
    ###### print functions
    def print_datafiles(self,data=None):
        '''show data files *.dat or .xyz depending on mulslice'''
        if not data : data = self.data
        if not isinstance(data,list) : data = [data]
        self._print_header('%s FILES' %(['.xyz','*.dat'][self.is_mulslice]))
        for dat in data :
            print(green+dat+black)
            with open(self.datpath+dat,'r') as f: print(''.join(f.readlines()))
            print("\n")
    def print_decks(self,decks=None):
        '''show decks *.in'''
        self._print_header('*.in FILES')
        if not decks:decks=list(self.decks.keys())
        if not isinstance(decks,list) : decks = [decks]
        for deck in decks:
            print(green+deck+black)
            try:
                with open(self.datpath+deck,'r') as f: print(''.join(f.readlines()))
            except FileNotFoundError:
                print(self.decks[deck])
                print(red+"WARNING:deck file %s not on disk" %(self.datpath+deck)+black)
    def print_log(self):
        '''print the logfile of running multislice .log'''
        try:
            with open(self.outf['log'],'r') as f:
                self._print_header('.log FILE')
                print(''.join(f.readlines()))
        except FileNotFoundError:
            print(red+'logfile not created yet, run a simu first'+black)
            return 0

    def check_simu_state(self,v=True):
        '''see completion state of a simulation '''
        try:
            with open(self.outf['log'],'r') as f :
                l1,l2=f.readlines()[-2:];#print(l1,l2)
            # get state
            if self.is_mulslice:
                if 'slice' in l2 :
                    slice = int(l2.split(',')[0].split(' ')[-1]);
                    tot_slice = (len(self.decks.keys())-1)*self.repeat[2];
                    state="%d%%" %(100*slice/tot_slice)
                elif 'elapsed time' in l2 : state='done'
                else : state='init'
            else:
                if 'z=' == l1[:2] : state="%d%%" %int(100*float(l1[3:8])/self.thickness)
                elif 'wall time' in l2 : state='done'
                else : state='init'
            #print state info
            if v :
                if state=='init': print(green+"Initializing..."+black)
                elif state=='done' : print(green+"Done : \n"+black, l1,l2)
                else : print(green+state.ljust(5)+black )
            return state
        except FileNotFoundError:
            if v : print(red+'logfile not created yet, run a simu first'+black)
            return 0

    ########################################################################
    ######## Private functions
    ########################################################################
    def _get_datafiles(self,data,v=True):
        dat_type = ['xyz','dat'][self.is_mulslice]
        err_msg=red+'no *.' +dat_type+' files found in :\n'+yellow+self.datpath+black
        if not data : data = lsfiles(self.datpath+'*.'+dat_type)
        if not data : raise Exception(err_msg)
        if not isinstance(data,list) : data=[data]
        #print(data[0].split('.')[-1])
        if not data[0].split('.')[-1]==dat_type : printf(err_msg)# raise Exception(err_msg)
        if self.is_mulslice :
            data = [basename(f) for f in data ]
        else :
            data = basename(data[0])
        if v :
            self.data = data
            self.print_datafiles()
        return data
    def _get_cell_params(self,v=False):
        if self.is_mulslice :
            cz=[]
            for dat in list(self.data):
                with open(self.datpath+dat) as f:
                    ax,by,cz_i = np.array(f.readline().rstrip().split(' '),dtype=float)
                    cz.append(cz_i)
        else :
            with open(self.datpath+self.data) as f:
                f.readline()
                ax,by,cz = np.array(f.readline().rstrip().split(' '),dtype=float)
        cell_params = [ax,by,[cz]]
        if v : print('ax=%.3fA, by=%.3f, cz=' %(ax,by), cz)
        return cell_params

    def _set_name(self,v=False):
        name = self.name+self.tail+['_autoslic','_mulslice'][self.is_mulslice]
        if v : print('Simu name pattern = '+self.name)
        return name
    def _set_NxNy(self,NxNy):
        if isinstance(NxNy,int) : NxNy = [NxNy,NxNy]
        return tuple(NxNy)
    def _set_hk(self,hk,Nhk):
        if Nhk:
            Nh,Nk,Nl=self.repeat
            h,k = np.meshgrid(range(Nhk),range(Nhk))
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
    #### decks
    def _atompot_deck(self,dat):
        datfile = self.datpath+dat
        datimg = datfile.replace('.dat',self.tail)
        deck  = "%s\n" %(datfile)               #.dat file
        deck += "%s\n" %(datimg)                #.tif file
        deck += "n\n"                           #record fft proj
        deck += "%d %d\n" %(self.NxNy)          #sampling
        deck += "%d %d 1\n" %(self.repeat[:2])  #super cell
        deck += "n\n"                           #thermal displacement
        return deck

    def _mulslice_deck(self):
        deck = "%d(%s)\n" %(self.repeat[2],'abcdefghijklmnopqrstuv'[:len(self.data)]) #
        for dat in self.data :                  #potential image files
            deck+="%s\n" %(self.datpath+self.tail+dat.replace('.dat','.tif'))
        deck += "%s\n" %(self.outf['pattern'])  #pattern file
        deck += "n\n"                           #partial coherence
        deck += "n\n"                           #start previous run
        deck += "%.2f\n" %self.keV              #wavelength
        deck += "0 0\n"                         #crystal tilt
        deck += "0 0\n"                         #beam tilt
        deck += "y\n"                           #record beams
        deck += "%s\n" %(self.outf['beamstxt']) #     filename
        deck += "%d\n" %len(self.hk)            #     nbeams
        for hk in self.hk :                     #     beam indices
            deck += "%d %d\n" %(hk)
        #deck += "n\n"                           #record cross sections
        #deck +=
        return deck

    def _autoslic_deck(self):
        deck  = "%s\n" %(self.datpath+self.data)#.xyz file
        deck += "%d %d %d\n" %self.repeat       #super cell
        deck += "%s\n" %(self.outf['pattern'])  #pattern file
        deck += "n\n"                           #partial coherence
        deck += "n\n"                           #start previous run
        deck += "%.4f\n" %self.keV              #wavelength
        deck += "%d %d\n" %self.NxNy            #sampling
        deck += "0 0\n"                         #crystal tilt
        deck += "%f\n" %(self.slice_thick)      #slice thickness
        deck += "y\n"                           #record beams
        deck += "%s\n" %(self.outf['beamstxt']) #     filename
        deck += "%d\n" %len(self.hk)            #     nbeams
        for hk in self.hk :                     #     beam indices
            deck += "%d %d\n" %(hk)
        deck += "n\n"                           #thermal vibration
        deck += "n\n"                           #intensity cross section
        return deck

    ########################################################################
    #### misc private
    def _print_header(self,msg,w=70):
        head='#'*w
        print(yellow+head)
        print('\t\t\t' +msg+ ' :')
        print(head+black)

########################################################################
#def : test
########################################################################
def test(name,Nxy=2048,Nhk=3):
    multi=Multislice(name,keV=200,Nhk=Nhk,tail='1',
        repeat=[6,6,40],NxNy=[Nxy,Nxy],slice_thick=1.3575,
        mulslice=False,v='',save=True)
    #multi.beam_vs_thickness()
    return multi
if __name__ == '__main__':
    name = 'dat/test/'
    # multi  = test(name,Nxy= 1024,Nhk=3)
    #p=multi.run()
    #multi=Multislice(name,mulslice=False,v='ndD',tail='Test',save=False);
    #multi=Multislice(name,mulslice=True,v='n',tail='test',save=False);
