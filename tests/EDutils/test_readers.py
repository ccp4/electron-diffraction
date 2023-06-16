from utils import*
import glob,os
from EDutils import readers ;imp.reload(readers)
paths=[
    '/home/tarik/Documents/data/ipraglifozin/exp',
    '/home/tarik/Documents/data/epicza/EPICZA_ED_Dataset_1/',
    '/home/tarik/Documents/git/ccp4/src/electron-diffraction/results/biotin/dat/dials/biotin_xtal1',
    '/home/tarik/Documents/git/ccp4/src/ED_processing/PETS2/Glycine/tiff',
    '/home/tarik/Documents/git/ccp4/src/ED_processing/PETS2/Glycine/tiff/processed',
    '/home/tarik//Documents/data/RR-3/RR-3/2/XDS',
]

def test_reader(**kwargs):
    for path in paths:
        im=readers.read(glob.glob(os.path.join(path,'*.*'))[0])
    dsp.stddisp(im=[im],**kwargs)

def test_save(**kwargs):
    for path in paths:
        frames = readers.detect_frame(path)
        file=glob.glob(os.path.join(path,'*.%s' %frames['fmt']))[0]
        readers.save_image(file,**kwargs)

def test_detect_frame():
    for path in paths:
        print(readers.detect_frame(path))
# test_reader()
# test_detect_frame()
test_save(dir='dat/images',zmax=10,cmap='hot')
