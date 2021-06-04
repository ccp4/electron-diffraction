from utils import*                  ;imp.reload(dsp)
from EDutils import viewers as vw   ;imp.reload(vw)
plt.close('all')
tests = [1]

if -1 in tests:v = vw.Viewer()  #nothing should raise exception
if 0 in tests:v = vw.Viewer(path='dat/glycine')           #path alone (contains cif_file)
if 1 in tests:v = vw.Viewer(path='dat/pets',init_opts='u',rot=206)    #path alone (contains pts_file)
if 2 in tests:v = vw.Viewer(cif_file='dat/glycine/glycine.cif')       #cif_file alone
if 3 in tests:v = vw.Viewer(pets='dat/pets/glycine.pts')              #pets alone
if 4 in tests:v = vw.Viewer(bloch='dat/bloch/Si110_200keV_bloch.pkl') #bloch alone

if 5 in tests:v = vw.Viewer(path='dat/Si',cif_file='Si',u=[0,0,1])                  #path and cif
if 6 in tests:v = vw.Viewer(path='dat/Si',bloch='dat/bloch/Si110_200keV_bloch.pkl') #path and bloch
if 7 in tests:v = vw.Viewer(path='dat/glycine',pets='dat/pets/glycine.pts')         #path and pets

#path and cif and bloch
if 8 in tests:v = vw.Viewer(path='dat/Si',cif_file='Si',bloch='dat/bloch/Si110_200keV_bloch.pkl')

#path and cif and pets
if 9 in tests:v = vw.Viewer(path='dat/glycine',cif_file='dat/glycine/glycine.cif',pets='dat/pets/glycine.pts')

#path and pets and bloch
if 10 in tests:v = vw.Viewer(path='dat/glycine',
    pets='dat/pets/glycine.pts',bloch='dat/glycine/bloch/glycine001.pkl')

#path and cif_file pets and bloch
if 11 in tests:v = vw.Viewer(path='dat/glycine',cif_file='dat/glycine/glycine.cif',
    pets='dat/pets/glycine.pts',bloch='dat/glycine/bloch/glycine001.pkl')
