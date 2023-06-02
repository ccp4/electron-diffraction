from utils import*                   ;imp.reload(dsp)
from EDutils import pets as pets_imp ;imp.reload(pets_imp)
from utils import pytest_util        ;imp.reload(pytest_util)
plt.close('all')
pts_path = 'dat/pets/glycine/glycine.pts'
#make sure pets folder exists in dat
pt = pets_imp.Pets(pts_path,gen=False)
F = 98
df = pt.rpl.loc[pt.rpl.F==F][:5].copy()
h = df['hkl'].tolist()
df.index=h
df_pxy = pt.hkl_to_pixels(h,F)
print(df_pxy)
print(df.loc[h, ["rpx", "rpy"]])


def test_import_pets():
    pets = pets_imp.Pets(pts_path,gen=True)#,lam=0.02508,aper=0.005340,omega=230,gen=1)

@pytest_util.add_link(__file__)
def test_show_exp():
    pets = pets_imp.Pets(pts_path)
    vw=pets.show_exp(frame=65,h=True,pargs={'opt':''})
    return vw.fig,vw.ax

# pets.compare_xyz_pxpy(frame=frame,opts='oab')
# test_import_pets()
# test_show_exp()
# pets.show_uvw()
# pets.show_xyz()
# pets.show_hkl()


# pets.show_sim()
# pets.compare_xyz_pxpy(frame=frame,opts='oab')

# pets.show_xyz(view=[0,0],opt='sc',name='pets_glycine_xyz.png')
# pets.show_hkl(view=[0,0],opt='sc',name='pets_glycine_hkl.png')
# pets.show_frame(frame=frame,name='pattern_%d.png' %frame, opt='psc')
# pets.show_ewald_sphere(frame=19,Nmax=5,Smax=0.025,h3d=1)
# K = pets.get_beam(frame)
# df = mut.get_kinematic_pattern(cif_file,K,Nmax=5,Smax=0.02)
# pets.show_frames(frame=frame,thick=1000,Nmax=7,Smax=0.015,rot=206,Imag=10,opts='PKh',
#     v=0,pargs={'xylims':1.5})

# df = pets.compare_hkl(frame,Nmax=8,Smax=0.025,eps=None,v=0)

# frames = range(1,65) #[5,10,15,25,30]
# for f in frames:
#     pets.show_frame(frame=f,Smax=0.01,name='glycine_exp%d.png' %f,opt='sc')

# xylims = [-5,10,-5,10,-12,1]
# mut.show_cell(cif_file,title='%d' %frame,xylims=xylims,n=uvw,
#     axPos=[0,0,1,0.95],view=[0,0],h3D=1,opt='p')
# pets.mp4_crystal_rotation('glycine',xylims=xylims)
