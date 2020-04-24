from utils import*
from subprocess import check_output
from os.path import dirname,realpath
from glob import glob as lsfiles

class Multislice:
    '''
- `name`     : path to the simulation folder
- `mulslice` : True(mulslice)=>periodic, False(autoslic)
- `tail`     : trailing string in the naming convention of outputs
    - name pattern = *name*\_*tail*\_*program_name*
- `data`     : data file names
    - if None and mulslice simu : all *[a-zA-Z].dat and stacked alphabetically
    - if None and autoslic simu : first *.xyz found
    - Otherwise paths to filename(s) of the data files to use (must be in name)
- `keV,repeat,NxNy,slice_thick` : multislice simulation parameters
- `hk` : beams to record as function of depth
- `v`  : verbose/display option : n(naming pattern),d(data),D(Decks),r(runtime), all if True
- `save` : save the decks to disk

###Example
```python
multi=Multislice(name,mulslice=False,v='ndD',tail='Test',save=True)
multi.run()
```
    '''
    def __init__(self,name,mulslice,tail='',data=None,
                keV=200,repeat=[2,2,1],NxNy=[512,512],slice_thick=1.0,
                hk=[(0,0)],
                v=True,save=True):
        self.name        = basename(name)                     #;print('basename:',self.name)
        self.datpath     = realpath(name)+'/'                 #;print('path:',self.datpath)
        self.is_mulslice = mulslice
        self.tail        = '_'+tail if tail else ''
        self.data        = self._get_datafiles(data)          #;print(self.data)
        self.keV         = keV
        self.repeat      = repeat
        self.NxNy        = NxNy
        self.hk          = hk
        self.slice_thick = slice_thick
        self.exe         = 'run_temsim.sh' #temsim launcher
        ###### start
        if isinstance(v,bool) : v='ndD' #{1:'D',2:'dD'}[v]
        self._set_name('n' in v)
        self.log    = self.datpath+self.name+'.log'
        self.decks  = self.make_decks(save)
        if v :
            if 'd' in v : self.print_datafiles()
            if 'D' in v : self.print_decks()
    ##### Public functions
    def run(self,quiet=True):
        '''run the simulation with temsim'''
        cmd = '%s %s > %s' %(self.exe,self.decks.keys(),self.log)
        p = Popen('echo '+cmd,shell=True)
        return p
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
            print(self.decks[deck])
    def print_log(self):
        '''print the logfile of running multislice .log'''
        try:
            with open(self.log,'r') as f :
                self._print_header('.log FILE')
                ''.join(f.readlines())
        except FileNotFoundError:
            print(red+'logfile not created yet, run a simu first'+black)

    ######## Private functions
    def _get_datafiles(self,data):
        dat_type = ['xyz','dat'][self.is_mulslice]
        err_msg=red+'no *.' +dat_type+' files found in :\n'+yellow+self.datpath+black
        if not data : data = lsfiles(self.datpath+'*.'+dat_type)
        if not data : raise Exception(err_msg)
        if not isinstance(data,list) : data=[data]
        #print(data[0].split('.')[-1])
        if not data[0].split('.')[-1]==dat_type : raise Exception(err_msg)
        if self.is_mulslice :
            data = [basename(f) for f in data ]
        else :
            data = basename(data[0])
        return data

    def _set_name(self,check_name=False):
        self.name += self.tail+['_autoslic','_mulslice'][self.is_mulslice]
        if check_name :
            print(green+'simu name pattern : '+yellow+self.name+black)

    def make_decks(self,save=True):
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
        return decks

    def _atompot_deck(self,dat):
        datfile = self.datpath+dat
        deck  = "%s\n" %(datfile)
        deck += "%s\n" %(datfile.replace('.dat',self.tail+'.tif'))
        deck += "%d %d\n" %(self.NxNy[0],self.NxNy[1])
        deck += "%d %d 1\n" %(self.repeat[0],self.repeat[1])#,self.repeat[2])
        deck += "n\n" #thermal displacement
        return deck

    def _mulslice_deck(self):
        deck = "%d(%s)\n" %(self.repeat[2],'abcdefghijklmnopqrstuv'[:len(self.data)])
        for dat in self.data :
            deck+="%s\n" %(self.datpath+dat.replace('.dat','.tif'))
        deck += "%s\n" %(self.datpath+self.name+'.tif') #output
        deck += "n\n" #partial coherence
        deck += "n\n" #start previous run
        deck += "%.2f\n" %self.keV
        deck += "0 0\n" #crystal tilt
        deck += "0 0\n" #beam tilt
        deck += "n\n"   #record cross sections
        #deck +=
        return deck

    def _autoslic_deck(self):
        deck  = "%s\n" %(self.datpath+self.data)
        deck += "%d %d %d\n" %(self.repeat[0],self.repeat[1],self.repeat[2])
        deck += "%s\n" %(self.datpath+self.name+'.tif') #output
        deck += "n\n" #partial coherence
        deck += "n\n" #start previous run
        deck += "%.2f\n" %self.keV
        deck += "%d %d\n" %(self.NxNy[0],self.NxNy[1])
        deck += "0 0\n" #crystal tilt
        deck += "%f\n" %(self.slice_thick)
        deck += "y\n"   ######## record beams
        deck += "%s\n" %(self.datpath+self.name+'_beams.txt')
        deck += "%d\n" %len(self.hk)
        for hk in self.hk :
            deck += "%d %d\n" %(hk[0],hk[1])
        deck += "n\n"   #thermal vibration
        deck += "n\n"   #intensity cross section
        return deck

    #### misc private
    def _print_header(self,msg,w=70):
        head='#'*w
        print(yellow+head)
        print('\t\t\t' +msg+ ' :')
        print(head+black)

if __name__ == '__main__':
    name='../decks/test/Si110'
    #multi=Multislice(name,mulslice=False,v='ndD',tail='Test',save=False);
    #multi=Multislice(name,mulslice=True,v='n',tail='test',save=False);
