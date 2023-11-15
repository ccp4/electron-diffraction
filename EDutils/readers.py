import glob,os,re
import tifffile,mrcfile,cbf,numpy as np
from utils import glob_colors as colors
from utils import displayStandards as dsp
import cv2

nb_colors=100
cmaps = ['hot','viridis','Spectral','Greys']

cs  = {cmap : (np.array(dsp.getCs(cmap,nb_colors))*255).astype(np.uint8)
    for cmap in cmaps}

def read(file):
    fmt=file.split('.')[-1]#;print(fmt)
    return img_readers[fmt](file)

def tiff_reader(tiff_file)  :
    with open(tiff_file,'rb') as f:
        I = tifffile.imread(f)
    return I

def cbf_reader(cbf_file):
    content = cbf.read(cbf_file)
    numpy_array_with_data = content.data
    # header_metadata = content.metadata
    # print(colors.blue,header_metadata,colors.black)
    return numpy_array_with_data

def mrc_reader(mrc_file):
    with mrcfile.open(mrc_file) as mrc:
        return mrc.data

def smv_reader(file,head=False):
    with open(file, "rb") as f:
        buf= f.read()

    header = [l.strip('\n;') for l in buf[:512].decode().split('\n')]
    header = {k: v for k, v in map(lambda s: s.split("="), header[1:-3])}
    if head: print(header)
    h_bytes = int(header['HEADER_BYTES'])
    if 'TYPE' in header:
        assert(header['TYPE']=='unsigned_short')
    nx,ny=int(header['SIZE1']),int(header['SIZE2'])
    img = np.reshape(
        np.frombuffer(buf[h_bytes:], dtype=np.uint16),
        (nx,ny))
    return img


img_readers = {
    'mrc' : mrc_reader,
    'tiff': tiff_reader,
    'tif' : tiff_reader,
    'cbf' : cbf_reader,
    'smv' : smv_reader,
    'img' : smv_reader,
}
fmts = list(img_readers.keys())



def detect_frame(path):
    '''returns a structure with frames info'''
    files = glob.glob(os.path.join(path,'*.*'))
    if not len(files):return
    #detect format
    fmt,i = '',0
    while fmt not in fmts and i<len(files):
        fmt = files[i].split('.')[-1]
        i+=1

    if fmt:
        frames    = np.sort(glob.glob(os.path.join(path,'*.%s' %fmt))).tolist()
        nb_frames = len(frames)
        min_pad   = int(np.ceil(np.log10(nb_frames)))
        frame0 = os.path.basename(frames[0])
        frame1 = os.path.basename(frames[-1])
        #get the frame string which contains padding
        #if there are other occurences of numbers in the frame, take the one that changes from frame to frame
        min_frame_str = np.array(re.findall(r"[0-9]{%d,}" %(min_pad),frame0),dtype=str)
        max_frame_str = np.array(re.findall(r"[0-9]{%d,}" %(min_pad),frame1),dtype=str)

        if len(min_frame_str)>1:
            idx =~(min_frame_str==max_frame_str)
            min_frame_str=min_frame_str[idx]
            max_frame_str=max_frame_str[idx]
        min_frame_str = min_frame_str[0]
        max_frame_str = max_frame_str[0]

        prefix,suffix = frame0.replace('.%s' %fmt,'').split(min_frame_str)
        pad       = len(min_frame_str)
        min_frame = int(min_frame_str)
        max_frame = int(max_frame_str)

        # print('frames info : ', frame0,frame_str,pad,prefix,' ',min_frame)
        # max_frame    = int(os.path.basename(frames[-1]).replace('.'+fmt,'').replace(prefix,''))
        frames_dict = {
            'nb_frames' : nb_frames,
            'min_frame' : min_frame,
            'max_frame' : max_frame,
            'pad'       : pad,'fmt':fmt,
            'template'  : '%s%%s%s.%s' %(prefix,suffix,fmt),
        }
        # print(frames_dict)
    return frames_dict


def save_image(file,zmax=1000,cmap='hot',dir='',v=False):
    if v:print('reading...')
    im  = read(file)    #;print(im.shape)

    if v:print('processing...')
    idx = np.maximum(0,np.floor((np.minimum(im,zmax-1)/zmax)*nb_colors)).astype(np.int)

    fmt = file.split('.')[-1]
    out = os.path.join(dir,os.path.basename(file).replace(fmt,'jpg'))
    if v:print('writing...')
    success=cv2.imwrite(out, cs[cmap][idx,:])

    if success:
        print(colors.green,'file saved : ')
        print(colors.yellow,out,colors.black)
    else:
        print(colors.red,'saving failed',colors.black)
