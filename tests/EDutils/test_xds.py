from utils import *
import importlib as imp
from EDutils import xds;imp.reload(xds)

x = xds.XDS('dat/xds/XDS_ASCII.HKL')
# print(x.info)
# print(x.hkl[['qx','qy','px','py']])
if 0:
    hkl = x.rpl[:1000]
    print(hkl[['rx','ry','rz','r_x','r_y','r_z']])
    dists = np.linalg.norm(hkl[['rx','ry','rz']].values-hkl[['r_x','r_y','r_z']].values,axis=1)
    norms = np.linalg.norm(hkl[['rx','ry','rz']].values,axis=1)
    print('min max mean')
    print(norms.min(),norms.max(),norms.mean())
    print(dists.min(),dists.max(),dists.mean())
    # print('average accuracy : %.2f%%' %(100*dists.mean()/norms.mean()))
    print('average accuracy : %.2f%%' %(100*(dists/norms).mean()))

if 0:
    h = x.rpl.loc[x.rpl.F==100][:5].index.tolist()
    df_pxy = x.hkl_to_pixels(h,100)
    print(df_pxy)
    print(x.rpl.loc[h, ["px", "py"]])

if 0:
    hkl=x.rpl[['h','k','l','I']].copy()
    hkl['sig']=0.01
    x.to_shelx(hkl)
