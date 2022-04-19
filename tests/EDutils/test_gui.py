from utils import*                  #;imp.reload(dsp)
import time
from EDutils import gui as vw   ;imp.reload(vw)
from utils import pytest_util    ;imp.reload(pytest_util)
from pynput.keyboard import Key, Controller
plt.close('all')
out,ref,dir=pytest_util.get_path(__file__)

# v=vw.Gui(cif_file='pets/alpha_glycine.cif',pets_opts='BKVr',xylims=1)
# v=vw.Gui(cif_file='diamond',path=out,pets_opts='Bkr',xylims=1,help=0)
v=vw.Gui(cif_file='diamond',path=out,config=out+'/config.pkl')#pets_opts='Bkr',xylims=1,help=0)

# tests = range(10,12)#[]
# tests = [0]
@pytest_util.add_link(__file__)
def test_viewer() :
    v=vw.Gui(cif_file='pets/alpha_glycine.cif',pets_opts='BKVr',xylims=1,
        pargs={'opt':''})

    return v.fig,v.ax

dt=0.5
def press(key):
    time.sleep(dt)
    print('pressing %s' %key)
    keyboard.press(  key)
    keyboard.release(key)

def test_keyboard():
    keyboard = Controller()
    from threading import Thread
    thread = Thread(target=test_q)
    thread.start()
    plt.show()

def test_q():
    press(key='3')
    press(key=' ')
    press(key='2')
    press(key='3')
    press(key='4')
    press(key='5')

    for k in v.keymap.df_keys.key.values:
        press(key='5')


# test_keyboard()

# fig,ax=test_viewer();fig.show()
# if -1 in tests:v = vw.Viewer()  #nothing should raise exception
# if 0 in tests:v = vw.Gui(path='dat/glycine',init_opts='',pets_opts='SKr')           #path alone (contains cif_file)
# if 1 in tests:v = vw.Viewer(path='dat/pets',init_opts='',rot=206)    #path alone (contains pts_file)
# if 2 in tests:v = vw.Viewer(cif_file='dat/glycine/alpha_glycine.cif')       #cif_file alone
# if 3 in tests:v = vw.Viewer(pets='dat/pets/glycine.pts',init_opts='i',F='Ig',rot=206)              #pets alone
# if 4 in tests:v = vw.Viewer(bloch='dat/bloch/Si110_200keV_bloch.pkl') #bloch alone
#
# if 5 in tests:v = vw.Viewer(path='dat/Si',cif_file='Si',u=[0,0,1])                  #path and cif
# if 6 in tests:v = vw.Viewer(path='dat/Si',bloch='dat/bloch/Si110_200keV_bloch.pkl') #path and bloch
# if 7 in tests:v = vw.Viewer(path='dat/glycine',pets='dat/pets/glycine.pts')         #path and pets
#
# #path and cif and bloch
# if 8 in tests:v = vw.Viewer(path='dat/Si',cif_file='Si',bloch='dat/bloch/Si110_200keV_bloch.pkl')
#
# #path and cif and pets
# if 9 in tests:v = vw.Viewer(path='dat/glycine',cif_file='dat/glycine/alpha_glycine.cif',pets='dat/pets/glycine.pts')
#
# #path and pets and bloch
# if 10 in tests:v = vw.Viewer(path='dat/glycine',
#     pets='dat/pets/glycine.pts',bloch='dat/glycine/alpha_glycine-103_200keV_bloch.pkl')
#
# #path and cif_file pets and bloch
# if 11 in tests:v = vw.Viewer(path='dat/glycine',cif_file='dat/glycine/alpha_glycine.cif',
#     pets='dat/pets/glycine.pts',bloch='dat/glycine/0001_bloch.pkl')
