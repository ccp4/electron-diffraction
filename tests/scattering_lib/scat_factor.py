import importlib as imp
from utils import*
from crystals import Crystal
from scattering import scattering_factors as sf;imp.reload(sf)
opts='S'

q=np.linspace(0,2,1000)

# df = pd.read_csv(sf.dat_path+'ab_shelx.csv',index=0)

if 'S' in opts:
    fe_EJK = sf.get_fe(np.array([8,14]),q=q,fit='kirkland')
    fe_SLX = sf.get_fe(np.array(['O','Si']),q=q,fit='shelxl')
    cs = ['b','r']

    plts =[[q,fe_EJK[:,i],[cs[i],'-' ],''] for i in range(2)]
    plts+=[[q,fe_SLX[:,i],[cs[i],'--'],''] for i in range(2)]
    legElt = {'O':'r','Si':'b','shelx':'k--','ejk':'k-'}

    dsp.stddisp(plts,legElt=legElt,
        labs=[r'$q(\AA^{-1})$','$f_e(\AA)$'],
        lw=2,xylims=[0,1.5,0,6])

## Comparing Kirkland and Gruza
if 'G' in opts:
    fe = sf.get_fe(np.array([8,14]),q=q)
    im,xylims,s='data/sf_kirkland.png',[0,2,0,20],q
    im,xylims,s='data/sf_gruza.png',[0,1.5,0,6],q/2 #
    fig = dsp.image_bg(im,
        labs=[r'$q=\sin\theta/\lambda (\AA^{-1})$','$f_e(\AA)$'],
        plots=[s,fe,'b-'],lw=2,xylims=xylims)
    dsp.fix_pos(fig)

plt.show()
