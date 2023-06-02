from utils import *
import importlib as imp
from EDutils import dials_utils as dials;imp.reload(dials)

d = dials.Dials('dat/dials/biotin')
F=43

df = d.rpl.loc[d.rpl.F==F][:5].copy()
h = df['hkl'].tolist()
df.index=h
df_pxy = d.hkl_to_pixels(h,F)
print(df_pxy)
print(df.loc[h, ["px", "py"]])
