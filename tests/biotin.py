from utils import*
import multislice.multislice as mupy
import multislice.postprocess as pp
import importlib as imp
imp.reload(mupy)
path = '../multislice/dat/biotin/'

p1 = pp.load_multi_obj(path+'biotin_m1_autoslic.pkl')
