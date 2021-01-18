import importlib as imp
import pickle
from scattering.structure_factor import structure_factor3D
import multislice.multislice as mupy; imp.reload(mupy)
import multislice.postprocess as pp ; imp.reload(pp)
import multislice.mupy_utils as mut ; imp.reload(mut)
from crystals import Crystal
from utils import*
plt.close('all')
path = '../dat/epicza/'
figpath = '../docs_fig/'
file = path+'epicza.cif'
opts = 'E '  #E(Ewald) T(tilts) B(Beams)F(structFact)


if 'E' in opts:
    #Ewald configuration [010]
    mut.ewald_sphere(lat_params,lam=cst.eV2mum(12000)*1e4,tmax=60,T=20.,nx=5 ,ny=20,opt='p',name=figpath+'E_MX.png',xylims=[-1,1.01,-0.1,2.01],xyTicks=1,legLoc='upper right',fonts=fs)
    mut.ewald_sphere(lat_params,lam=cst.keV2lam(200)     ,tmax=7 ,T=0.2,nx=20,ny=10,opt='p',name=figpath+'E_ED.png',xylims=[-3,3.01,-0.1,1.01],xyTicks=1,legLoc='upper right',fonts=fs)
