from utils import *
from wallpp import lattice#;imp.reload(lattice)
# plt.close('all')
lat_types = lattice.lat_types
lats = {lat:lattice.Lattice2D(lat,a=2,b=3,alpha=50) for lat in lat_types }


def plot_lattices() :
    for lat_type,lat in lats.items():
        lat.plot_lattice(nh=6,nk=5, lw=2, equal=1)

def test_lattice_exception():
    try:
        lattice.Lattice2D('rubish')
    except Exception as e:
        print(colors.red,e,colors.black)

plot_lattices()
test_lattice_exception()
