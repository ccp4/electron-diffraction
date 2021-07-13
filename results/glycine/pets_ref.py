import pandas as pd,numpy as np
import crystals
from multislice import pets as pt#   ;imp.reload(pt)
from multislice import mupy_utils as mut#   ;imp.reload(pt)

file = 'dat/pets/refinements/out.txt'
FoFc = 'dat/pets/refinements/FoFc.txt'
#pets = pt.Pets('dat/pets/glycine.pts')

df = pd.read_csv(FoFc,index_col=None,header=0,delimiter='  *',engine='python')


# Zs = {'C':6,'O':8,'N':7,'H':1}
# df = pd.read_csv(file,index_col=0,names=['x','y','z'],delimiter=' ')
# df['Z'] = [c[0] for c in df.index]
# df['Za'] = [Zs[Z] for Z in df.Z]
#
# lat_vec,lat_params = pets.lat,pets.lat_params[:3]
# pattern = np.hstack([df[['Za','x','y','z']].values,np.ones((df.shape[0],2))])
# mut.make_xyz('dat/test.xyz',pattern,lat_vec,lat_params)
# mut.show_grid3('dat/test.xyz',ms=50)
#

# equiv= ['x,y,z',
#         '-x,-y,-z',
#         '-x+1/2,y+1/2,-z+1/2',
#         'x+1/2,-y+1/2,z+1/2',]
# symmetry_operators = list(map(crystals.CIFParser.sym_ops_from_equiv, equiv))
# atoms = [crystals.Atom(c.Z, coords = [c.x,c.y,c.z]) for n,c in df.iterrows()]
#
# unitcell = crystals.crystal.symmetry_expansion(atoms, symmetry_operators)
# glycine = crystals.Crystal(unitcell, lat_vec)
