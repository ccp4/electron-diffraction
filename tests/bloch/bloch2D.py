from utils import*
from blochwave import bloch2D as bl          ;imp.reload(bl)
from multislice import mupy_utils as mut     ;imp.reload(mut)
import multislice.multi_2D as ms             ;imp.reload(ms)
plt.close('all')


file = 'dat/p1c.py'
thick=100
# p1 = mut.import_wallpp(file)
# p1.plot_unit_cells(opts='uAa',nh=3,nk=3)

b0  = bl.Bloch2D(file,keV=200,u=[1,10],thick=thick,eps=0.5,
    Nmax=10,Smax=0.1,solve=1)

b0.show_ewald()
b0.show_beams(opts='B',fz=np.abs)
b0.show_beams(opts='SVB')
b0.show_beams_vs_thickness(thicks=(0,thick,1000))
b0.G('Sw')
