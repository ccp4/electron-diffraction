from utils import*
from wallpp import gui         ;imp.reload(gui)
from wallpp import config as cg;imp.reload(cg)
import sys,subprocess
plt.close('all')

# wp = gui.Wallpaper_gui(image='figures/cygne2_100.png',config=1,
#     path='dat',pg='p3')

def test_wiki_wallpp():
    for pp_type in cg.pp_types[1:]:
        image ='figures/%s.png' %pp_type
        wp = gui.Wallpaper_gui(image=image,path='dat',pg=pp_type,
            config=1,opt=0)

def convert_figures():
    files = sys.argv[1:]#'figures/p4g.png'
    df = cg.df_wallpp
    for f in files:
        pp_type = f.split('-')[1]                           #;print(pp_type)
        # pp_type = f.split('.png')[0]                           #;print(pp_type)
        if not pp_type in df.index:                         #
            pp_type = df.loc[df['pp_long'==pp_type]].name   #;print(pp_type)
        i = np.arange(17)[np.array(cg.pp_types)==pp_type][0]
        # old   = '%s.png' %pp_type                           #;print(image)
        image = 'wallpp_%s_%s.png' %(str(i).zfill(2),pp_type)    #;print(image)
        cmd = "convert %s -resize 400x400 %s" %(f,image)    #;print(cmd)
        # cmd = "mv %s %s" %(old,image)    #;print(cmd)
        subprocess.Popen(cmd,shell=True)

convert_figures()
def F(npts=21,lw=0.05,x0=0.5,y0=0.5):
    x,y = np.meshgrid(np.linspace(0,1,npts),np.linspace(0,1,npts))
    F = np.zeros((npts,npts))
    F[(y>0.15) & (y<0.75) & (abs(x-0.3)<=lw)]=1
    F[(abs(y-0.75)<=lw) & (x>0.3) & (x<0.7)]=1
    F[(abs(y-0.5)<=lw) & (x>0.3) & (x<0.55)]=1
    x,y = np.meshgrid(np.linspace(0,x0,npts),np.linspace(0,y0,npts))
    pattern = np.vstack([x.flatten(),y.flatten(),F.flatten()]).T
    np.save('dat/F.npy',pattern)
