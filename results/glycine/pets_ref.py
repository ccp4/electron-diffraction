import pandas as pd,numpy as np,crystals
from utils import*
from multislice import pets as pt          ;imp.reload(pt)
from multislice import mupy_utils as mut   #;imp.reload(mut)
# I(integrate I over frames and compare with hkl.I)
# i(integrate on frame manually) c(Fo vs Fc_dyngo), o(check Fo=sqrt(Io) ok)
opts='c'
file = 'dat/pets/refinements/out.txt'
FoFc = 'dat/pets/refinements/FoFc.txt'
pets = pt.Pets('dat/pets/glycine.pts',dyn=0)


#### sum/mean I over frames and compare with hkl.I
if 'I' in opts:
    # refl = str(tuple([0,0,2]))
    # I,F = pets.rpl.loc[pets.rpl.hkl==refl,['I','F']].values.T
    # dsp.stddisp([F,I,'b-',''],lw=2,labs=['F','I'])
    refl = pets.hkl.index
    sumI = [pets.rpl.loc[pets.rpl.hkl==r,'I'].sum() for r in refl]
    hklI = pets.hkl.I.values
    dsp.stddisp([sumI,hklI,'bo',''],labs=['sum(I)','hkl(I)'])


### frame integrate beams (manual fit) and compare to rpl.I
if 'i' in opts:
    df_Iint = pets.integrate_rpl(frames=np.arange(10,20))
    plts = [[df_Iint.Iint,df_Iint.Imax,'ro']]
    dsp.stddisp(plts,labs=['$I_{int}$','$I_{max}$'])
    plts = [[df_Iint.I,df_Iint.Im,'bo']]
    dsp.stddisp(plts,labs=['$I$','$Im$'])



## import dynamical refinement and keep unique reflections
df = pd.read_csv(FoFc,index_col=None,header=0,delimiter='  *',engine='python')
hkl = [str(tuple(h)) for h in df[['h','k','l']].values]
hkl,idx,c = np.unique(hkl,return_index=True,return_counts=True)
df = df.iloc[idx]
df.index = hkl
## Fo vs Fc from dyngo dynamical refinement
if 'c' in opts:
    plts=[]
    plts+=[[df.Fo,df.Fc,'bo']]
    dsp.stddisp(plts,labs=['$F_o$','$F_c$'])
    print('Rfactor:%.1f',abs(np.sum(df.Fo**2)-np.sum(df.Fc**2))/sum(df.Fo**2))

## checked that df.Fo=sqrt(hkl_dyn.Io)
if 'o' in opts:
    pets.HKL_dyn['hkl'] = [str(tuple(h)) for h in pets.HKL_dyn[['h','k','l']].values]
    idx = pets.HKL_dyn
    pets.HKL_dyn.index=pets.HKL_dyn['hkl']
    I = np.array([pets.HKL_dyn.loc[[h],'I'][0] for h in df.index])
    Fo  = np.sqrt(I)
    print('Fo-sqrt(I)=',(abs(df['Fo']-Fo)/Fo).max())
    # plts=[[df['Fo'],Fo,'rs']]
    # dsp.stddisp(plts,labs=['$F_o$','$F_c$'])






# if 0:
#     Zs = {'C':6,'O':8,'N':7,'H':1}
#     df = pd.read_csv(file,index_col=0,names=['x','y','z'],delimiter=' ')
#     df['Z'] = [c[0] for c in df.index]
#     df['Za'] = [Zs[Z] for Z in df.Z]
#
#     lat_vec,lat_params = pets.lat,pets.lat_params[:3]
#     pattern = np.hstack([df[['Za','x','y','z']].values,np.ones((df.shape[0],2))])
#     mut.make_xyz('dat/test.xyz',pattern,lat_vec,lat_params)
#     mut.show_grid3('dat/test.xyz',ms=50)
#
#     equiv= ['x,y,z',
#             '-x,-y,-z',
#             '-x+1/2,y+1/2,-z+1/2',
#             'x+1/2,-y+1/2,z+1/2',]
#     symmetry_operators = list(map(crystals.CIFParser.sym_ops_from_equiv, equiv))
#     atoms = [crystals.Atom(c.Z, coords = [c.x,c.y,c.z]) for n,c in df.iterrows()]
#
#     unitcell = crystals.crystal.symmetry_expansion(atoms, symmetry_operators)
#     glycine = crystals.Crystal(unitcell, lat_vec)
