import importlib as imp
from utils import*
import multislice.pymultislice as ms;imp.reload(ms)
import multislice.multi_3D as MS3D  ;imp.reload(MS3D)
import multislice.mupy_utils as mut ;imp.reload(mut)
from crystals import Crystal
plt.close('all')
path='../../multislice/docs_fig/multi3D/'
new = 1
name = 'dat/test.pkl'


#cif file for alpha-glycine from :
#https://pubs.acs.org/doi/suppl/10.1021/cg101059c/suppl_file/cg101059c_si_002.pdf
cif_file = 'dat/alpha_glycine.cif'
# xyz_file = cif_file.replace('.cif','.xyz')
crys = Crystal.from_cif(cif_file)
lat_vec = np.array(crys.lattice_vectors)
pattern = mut.import_cif(cif_file)#,xyz=xyz_file)#

pattern[:,3]-=pattern[:,3].min()

lat_params = lat_vec.max(axis=0)-lat_vec.min(axis=0)
# mut.plot_grid(pattern,lat_params,'xy',equal=True,ms=20)
# mut.plot_grid(pattern,lat_params,'xz',equal=True,ms=20)
# mut.plot_grid(pattern,lat_params,'yz',equal=True,ms=20)


if new:
    ax,by,cz = lat_params
    mp = MS3D.Multi3D(pattern,ax,by,cz,
        keV=200,Nxy=[3,2],nxy=2**10,Nz=1,dz=1,
        TDS=False,
        iZs=1,iZv=1,opts='q',v=1,
        hk=3,hkopt='sro',
        eps=1,copt=0,
        ppopt='B')#XQZTP)
    mp.save(name)

else:
    mp = ms.load(name)

mp.Bz_show(cm='jet')
# mp.Vz_show(iz=0,cmap='viridis')
# mp.Vz_show(iz=1,cmap='viridis')
# mp.Vz_show(iz=2,cmap='viridis')
