import numpy as np
import os

f = np.load(os.path.dirname(__file__)+'/data/abcd.npy')

def Va(r,Z,ip=np.arange(3)) :
    al,ag=150.4121417,266.5985798;
    V=np.zeros(r.shape)
    for i in ip:
        ai,bi = f[Z][2*i:2*i+2]
        ci,di = f[Z][6+2*i:6+2*i+2]
        # print(ai,bi)
        # print(ci,di)
        fi = ai/r*np.exp(-2*np.pi*r*np.sqrt(bi))
        gi = ci*di**(-3/2)*np.exp(-np.pi**2*r**2/di)
        V+=al*fi + ag*gi
    return V
