import importlib as imp
# from subprocess import Popen,PIPE
from utils import*                      ;imp.reload(dsp)
import multislice.mupy_utils as mut     ;imp.reload(mut)
import multislice.pets as pp            ;imp.reload(pp)
plt.close('all')
pts_path = '../glycine.pts'
cif_file = 'glycine.cif'


pets = pp.Pets(pts_path,gen=0)#,lam=0.02508,aper=0.005340,omega=230,gen=1)
frame = 18

# pets.show_exp(frame,v=0)
# pets.show_sim()

# pets.compare_xyz_pxpy(frame=frame,opts=''  ,view='x',name='xyz_vs_pxpy_raw_x.png',opt='sc')
# pets.compare_xyz_pxpy(frame=frame,opts=''  ,view='z',name='xyz_vs_pxpy_raw_z.png',opt='sc')
# pets.compare_xyz_pxpy(frame=frame,opts='oa',view='x')#,name='xyz_vs_pxpy_prc_x.png',opt='sc')
# pets.compare_xyz_pxpy(frame=frame,opts='oa',view='z')#,name='xyz_vs_pxpy_prc_z.png',opt='sc')
# pets.compare_xyz_pxpy(frame=frame,opts='oa')

# pets.show_uvw()
# pets.show_xyz()
# pets.show_xyz(view=[0,0],opt='sc',name='pets_glycine_xyz.png')
# pets.show_hkl(view=[0,0],opt='sc',name='pets_glycine_hkl.png')
# pets.show_frame(frame=frame,name='pattern_%d.png' %frame, opt='psc')
# pets.show_ewald_sphere(frame=19,Nmax=5,Smax=0.025,h3d=1)
# pets.show_hkl()
# K = pets.get_beam(frame)
# df = mut.get_kinematic_pattern(cif_file,K,Nmax=5,Smax=0.02)
pets.show_frames(frame=frame,thick=1000,Nmax=7,Smax=0.015,rot=206,Imag=10,opts='PKh',
    v=0,pargs={'xylims':1.5})

# df = pets.compare_hkl(frame,Nmax=8,Smax=0.025,eps=None,v=0)

# frames = range(1,65) #[5,10,15,25,30]
# for f in frames:
#     pets.show_frame(frame=f,Smax=0.01,name='glycine_exp%d.png' %f,opt='sc')

# xylims = [-5,10,-5,10,-12,1]
# mut.show_cell(cif_file,title='%d' %frame,xylims=xylims,n=uvw,
#     axPos=[0,0,1,0.95],view=[0,0],h3D=1,opt='p')
# pets.mp4_crystal_rotation('glycine',xylims=xylims)
