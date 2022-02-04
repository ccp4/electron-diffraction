from utils import*
from EDutils import viewers as vw;imp.reload(vw)
from utils import pytest_util    ;imp.reload(pytest_util)

# @pytest_util.add_link(__file__)
def test_show_exp():
    v= vw.Base_Viewer('../bloch/out/test_continuous/tiff',i=10,pargs={'opt':'p'})
    return v.fig,v.ax
test_show_exp()
