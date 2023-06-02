import blochwave,time,sys,pytest,os,copy
from subprocess import check_output
from pytest_html import extras
from utils import glob_colors as colors
from utils import*                      ;imp.reload(dsp)
from utils import pytest_util           ;imp.reload(pytest_util)
from blochwave import bloch             ;imp.reload(bloch)
from blochwave import util as bloch_util;imp.reload(bloch_util)
from EDutils import pets as pets_imp    ;imp.reload(pets_imp)
from EDutils import utilities as ut     ;imp.reload(ut)

opts='r'
thick=330
pets = ut.load_pkl('dat/pets/pets.pkl')

if 'r' in opts:
    pets.make_eldyn(F=3,thickness=thick,phisteps=50,thr=4,path='dyngo',
        # bargs={'Nmax':4},
        v=0)
    pets.run_dyngo()
    # pets.dyngo_path=os.path.join(pets.path,'dyngo')
    # b0 = pets.run_blochwave(F=F,Nmax=4,Smax=0.02)

df = pets.get_edout(F=F,opts='')
