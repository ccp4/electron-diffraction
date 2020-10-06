import multislice.multislice as mupy
import multislice.postprocess as pp
import importlib as imp
imp.reload(mupy)
path = '../dat/biotin/'

# p1 = pp.load_multi_obj('dat/biotin/biotin_1_autoslic.pkl')
p1 = pp.load_multi_obj(path+'biotin_m2_autoslic.pkl')
p1.beam_vs_thickness(iBs=['(2,17)','(8,34)'],bOpt='f')
# p1.save_pattern()
# p1.pattern(Iopt='Isnlg',tol=1e-4,rings=[0,0.25,0.5,1],caxis=[-6.2,0],
#     pOpt='ptX',xylims=[-1,1,-1,1],cmap='binary',imOpt='ch',axPos=[0.2,0.13,0.75,0.75])
