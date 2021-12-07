from utils import*
import cbf,tifffile, glob,os

exp = '../simus/ireloh2/exp/'
out = 'ireloh2/'


def convert_files(exp,out):
    cbf_files = glob.glob(exp+'*.cbf')
    print(colors.green+'file saved'+colors.black)
    for file in cbf_files :
        try:
            data = cbf.read(file).data
            base_file = os.path.basename(file).replace('.cbf','')
            tif_file = out+'tiff/'+base_file+'.tiff'
            tifffile.imwrite(tif_file,data)
            print(colors.yellow+tif_file+colors.black)
        except:
            print('error with file %s' %file)
def make_pets(pts_file,aperpixel,deg=0.0652,ref_cell=None):
    center,pixel_size='AUTO',0.01,
    phi,omega= deg/2 ,0
    ax,by,cz,alpha,beta,gamma = ref_cell
    # ax,by,cz,alpha,beta,gamma = 8.1218, 9.977, 17.725, 90.0, 90.0, 90.0
    pts = '''lambda 0.025080
geometry continuous
omega  %d
phi  %.5f
virtualframes   7 5 1

aperpixel    %.6f
noiseparameters      2.5000      1.0000
saturationlimit   20000

center    %s
beamstop no

dstarmax  1.800
dstarmaxps  2.0
i/sigma    7.00    5.00
reflectionsize  10

referencecell     %.5f    %.5f     %.5f    %.5f   %.5f    %.5f 1

#List of images
#columns: file name,alpha,beta,delta omega,frame scale,calibration correction(px/rec.angstrom),ellipt.distortion correction(amplitude),ellipt.distortion correction(phase), use for calculations(0/1)
imagelist
'''%(omega,phi,aperpixel,center,  ax,by,cz,alpha,beta,gamma)
    tif_files = np.sort(glob.glob(out+'/tiff/*.tiff'))
    alphas = np.arange(tif_files.size)*deg
    for i,tif_file in enumerate(tif_files):
        tif_file = os.path.basename(tif_file)
        pts+='%s %.4f 0.0 0.0 1.0 0 0 0  1\n' %('tiff\\'+tif_file,alphas[i])
    pts+='endimagelist\n'
    # pts_file = out+name+'.pts'
    with open(pts_file,'w') as f:
        f.write(pts)
        print(colors.green+'file saved : '+colors.yellow+pts_file+colors.black)

# convert_files(exp,out)
make_pets('ireloh')
