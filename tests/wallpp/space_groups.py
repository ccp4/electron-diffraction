import importlib as imp
import httplib2, re
import pandas as pd,numpy as np,matplotlib.pyplot as plt
from utils import displayStandards as dsp;imp.reload(dsp)
plt.close('all')
df_file = 'dat/space_groups.pkl'

def get_space_groups():
    http = httplib2.Http()
    df = pd.DataFrame(columns=['space group','ops'])
    for i in range(1,231):
      headers, body = http.request("http://img.chem.ucl.ac.uk/sgp/large/%saz3.htm" %str(i).zfill(3))
      lines = body.decode().strip().split('\r\n')
      tle = [l for l in lines if '<TITLE>Space Group' in l][0]
      spg = tle.split(':')[1].split(';')[0]
      df.loc[i,'space group'] = spg

      print(spg)
      l1 = [i for i,l in enumerate(lines) if 'Symmetry operators' in l][0]
      l2 = [i for i,l in enumerate(lines) if '</PRE></FONT>' in l][0]
      ops = lines[l1+4:l2]
      df.loc[i,'ops'] = ops
      df.to_pickle(df_file)

# get_space_groups()

df = pd.read_pickle(df_file)
x,y,z = 0.1,0.2,0.7
dx={}
for i in df.index:
    c=df.loc[i]
    spg = c['space group']                              #;print(spg)
    ops = [l for l in c.ops if l.strip() != '']         #;print(ops)
    X = np.array([np.array(eval(op)) for op in ops])    #;print(X)
    X[X<0]+=1
    X[X>1]-=1
    dx[i] = X

    x0,y0,z0 = X.T
    c0 = dsp.getCs('Spectral',x0.size)
    fig,(ax0,ax1,ax2)=dsp.create_fig('f',rc=[1,3],pad=5)#;plt.subplots(1,3)
    dsp.stddisp(scat=[x0,y0,c0],ms=80,xylims=[0,1,0,1],xyTickLabs=[],pOpt='tXG',labs=['$x$','$y$'],ax=ax0,opt='')
    dsp.stddisp(scat=[y0,z0,c0],ms=80,xylims=[0,1,0,1],xyTickLabs=[],pOpt='tXG',labs=['$y$','$z$'],ax=ax2,opt='')
    dsp.stddisp(scat=[x0,z0,c0],ms=80,xylims=[0,1,0,1],xyTickLabs=[],pOpt='tXG',labs=['$x$','$z$'],ax=ax1,
        title='%d : $%s$' %(i,spg),name='figures/%s.png' %str(i).zfill(3),
        opt='sc')

    # dsp.stddisp(scat=[x0,y0,z0,'b'],ms=50,rc='3d',title='$%s$' %spg,
    #     xylims=1,#[0,1,0,1,0,1],
    #     name='figures/%s.png' %str(i).zfill(3),opt='sc')
