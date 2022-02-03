import importlib as imp
from utils import*
from crystals import Crystal
from scattering import scattering_factors as sf;imp.reload(sf)

q=np.linspace(0,2,1000)
fe = sf.get_fe(np.array([6]),q=q)

im,xylims,s='data/sf_kirkland.png',[0,2,0,20],q
im,xylims,s='data/sf_gruza.png',[0,1.5,0,6],q/2 #
fig = dsp.image_bg(im,labs=[r'$q=\sin\theta/\lambda (\AA^{-1})$','$f_e(\AA)$'],
    plots=[s,fe,'b-'],lw=2,xylims=xylims)
dsp.fix_pos(fig)
