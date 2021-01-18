import importlib as imp
from utils import*
from crystals import Crystal
from scattering import structure_factor as sf;imp.reload(sf)

crys = Crystal.from_database('Si')
pattern = np.array([np.hstack([a.coords_fractional,0*a.atomic_number]) for a in crys.atoms] )
lat_vec = np.array(crys.reciprocal_vectors)
hkl,Fhkl=sf.structure_factor3D(pattern,lat_vec,hkl=None,hklMax=2,sym=0,v='')
h,k,l = hkl
I = np.abs(Fhkl)**2

id0 = ~((h%2==0) & (k%2==0) & (l%2==0)) & ~((h%2==1) & (k%2==1) &(l%2==1))

print(np.abs(I[id0]))                   #should be 0
print(np.abs(I[~id0 & ((h+k+l)%4==0)])) #should be 64
print(np.abs(I[~id0 & ((h+k+l)%2==1)])) #should be 32
print(np.abs(I[~id0 & ((h+k+l)%4==2)])) #should be 0
