from utils import*                  ;imp.reload(dsp)
from EDutils import viewers as vw   ;imp.reload(vw)
plt.close('all')


# tests = range(10,12)#[]
tests = [0]

if -1 in tests:v = vw.Viewer()  #nothing should raise exception
if 0 in tests:v = vw.Viewer(path='dat/glycine',init_opts='',pets_opts='SKr',rots={'K':206,'P':0,'M':98})           #path alone (contains cif_file)
if 1 in tests:v = vw.Viewer(path='dat/pets',init_opts='',rot=206)    #path alone (contains pts_file)
if 2 in tests:v = vw.Viewer(cif_file='dat/glycine/alpha_glycine.cif')       #cif_file alone
if 3 in tests:v = vw.Viewer(pets='dat/pets/glycine.pts',init_opts='i',F='Ig',rot=206)              #pets alone
if 4 in tests:v = vw.Viewer(bloch='dat/bloch/Si110_200keV_bloch.pkl') #bloch alone

if 5 in tests:v = vw.Viewer(path='dat/Si',cif_file='Si',u=[0,0,1])                  #path and cif
if 6 in tests:v = vw.Viewer(path='dat/Si',bloch='dat/bloch/Si110_200keV_bloch.pkl') #path and bloch
if 7 in tests:v = vw.Viewer(path='dat/glycine',pets='dat/pets/glycine.pts')         #path and pets

#path and cif and bloch
if 8 in tests:v = vw.Viewer(path='dat/Si',cif_file='Si',bloch='dat/bloch/Si110_200keV_bloch.pkl')

#path and cif and pets
if 9 in tests:v = vw.Viewer(path='dat/glycine',cif_file='dat/glycine/alpha_glycine.cif',pets='dat/pets/glycine.pts')

#path and pets and bloch
if 10 in tests:v = vw.Viewer(path='dat/glycine',
    pets='dat/pets/glycine.pts',bloch='dat/glycine/alpha_glycine-103_200keV_bloch.pkl')

#path and cif_file pets and bloch
if 11 in tests:v = vw.Viewer(path='dat/glycine',cif_file='dat/glycine/alpha_glycine.cif',
    pets='dat/pets/glycine.pts',bloch='dat/glycine/0001_bloch.pkl')
