from utils import*
from postprocess import*
from utils import displayStandards as dsp
from crystals import Crystal
import multislice as mupy
import sys,numpy,os
import rotating_crystal as rcc
import importlib as imp
imp.reload(mupy)
imp.reload(rcc)


path='dat/biotin/'
file = path+'biotin.cif'
xyz  = path+'001.xyz'

# rcc.show_cell(file,x0=-1)
pattern = rcc.import_cif(file,xyz,rep=[1,1,1],dopt='s')


def run_simu(**kwargs):
    Ns = 10
    multi=mupy.Multislice(path,#tail=tail,data=
        mulslice=False,keV=200,
        NxNy=2048,slice_thick=1.0,Nhk=5,repeat=[2*Ns,Ns,5],
        #TDS=True,T=300,n_TDS=15,
        opt='sr',fopt='',v='nctr',#nctrdDR',
        # ssh='tarik-CCP4home',
        )
    # multi.wait_simu(ssh_alias='tarik-CCP4home')
    # multi.print_log()
    # multi.postprocess(ppopt='uwP',ssh_alias='tarik-CCP4home')
    multi.save_pattern()
    # multi.pattern(tol=1e-4,Iopt='Inscq',#caxis=[-9,0],#Nmax=100,
    #     cmap='binary',imOpt='hc',pOpt='Xt',rings=[0.1,0.2,1])
    return multi
multi = run_simu()
