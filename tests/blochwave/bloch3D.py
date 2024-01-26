from utils import*
from blochwave import bloch         ;imp.reload(bloch)
b0 = bloch.Bloch('dat/LTA.cif',path='dat/dummy/',keV=200,
    u=[0,0,1],Nmax=4,Smax=0.05,
    opts='vt',thick=100)
df_Fhkl = pd.read_pickle(b0.get_Fhkl_pkl())
