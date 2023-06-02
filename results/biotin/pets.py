from utils import*                   ;imp.reload(dsp)
from EDutils import pets as pt       ;imp.reload(pt)
from blochwave import util as bu     ;imp.reload(bu)
import os,mrcfile,glob
from subprocess import check_output

opts = 'pt' #p(pets) t(tiff)
path='dat/pets_sim'

def mrc2pets(mrcs,out):
    mrc_files = glob.glob(mrcs)

    if not os.path.exists(out):
        check_output('mkdir -p %s' % out,shell=True).decode()
    if 't' in opts:
        for mrc_file in mrc_files :
            bu.mrc2tiff(mrc_file,out)
    return mrc_files


mrc_files = mrc2pets(mrcs=path+'/mrc/*.mrc', out=path+'/tiff')
n_frames = len(mrc_files)

mrc = mrcfile.open("dat/dials/biotin_xtal1/biotin_xtal1_0000.mrc")
aper = mrc.extended_header["Pixel size X"][0]*1e-10 #A^-1
deg  = (mrc.extended_header['End tilt angle']-mrc.extended_header['Start tilt angle'])/n_frames
deg = deg[0]
mrc.close()

pt.make_pets(path+'/biotin.pts',aperpixel=aper,deg=deg,
    dmM=(0.2,2),pxy=(1024,)*2,
    ref_cell=(5,10,20,90,90,90))
