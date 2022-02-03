from utils import*                  #;imp.reload(dsp)
from blochwave import bloch_pp as bl;imp.reload(bl)
from blochwave import bloch         ;imp.reload(bloch)
from multislice import pets as pt   ;imp.reload(pt)
from EDutils import utilities as ut ;imp.reload(ut)
plt.close('all')
path = 'dat/bloch/optim/'

opts = 'F' #S(single), F(Full optim)
frame = 14

pets = pt.Pets('dat/pets/glycine.pts')
tag = 'glycine%s_' %str(frame).zfill(3)
bloch_args = {'cif_file':'dat/pets/alpha_glycine.cif',
    'keV':200,'Smax':0.025,'Nmax':9,'solve':1,'opts':'s'}

rpl = pets.rpl.loc[pets.rpl.F==frame]
rpl.index = rpl['hkl'].values.T
refl = list(rpl.index.values)

def get_cor(b,rpl,hkl):
    Irpl = np.array(rpl.loc[hkl,'I'],dtype=float)
    Iz = b.get_beams_vs_thickness(hkl,dict_opt=0)
    cor = np.abs([np.corrcoef(Irpl, I0)[0,1] for I0 in Iz.T])
    iZ  = np.argmax(cor)
    Iz0 = Iz[:,iZ]
    return hkl,cor,iZ,Iz0,Irpl

# single frame optimization
if 'S' in opts:
    rock = ut.load_rock(path,tag)

    b = bloch.Bloch(u=pets.uvw0[frame],path=path,name='test0',**bloch_args)
    # b.show_beams(refl=refl,opts='SVk',rot=203,xylims=2)
    hkl = b.get_beam(refl=refl,cond='',opt=0)
    # print(np.setdiff1d(refl,hkl))
    # b.set_beams_vs_thickness(thicks=(0,1000,1000));b.save()
    hkl,cor,iZ,Iz0,Irpl = get_cor(b,rpl,hkl)
    dsp.stddisp([Irpl,Iz0*Irpl.mean()/Iz0.mean(),'bo'],labs=['$I_o$','$I_c$'])

if 'F' in opts:
    if 0:
        u0 = pets.uvw0[frame]
        uvw = ut.get_uvw_CBED(u0=u0,deg=1,npts=5,plot=0,h3d=1)
        omega = np.arange(uvw.shape[0])
        # hkl = rock._get_refl(refl=refl,cond='')[0].tolist()
        rock = bl.bloch_rock(tag=tag,uvw=uvw,bloch_args=bloch_args,omega=omega,
            cond='',refl=refl,
            thicks=(0,500,100),
            W_args={'opts':'tf'},#'iTs':slice(2,14)},
            R_args={'opts':'f'},
            # S_args={'opts':'tf'},
            # ts0=0,
            # R_args={'cmap':'Spectral'},
            zs = np.arange(100,500,100),cm='hsv',hkls=[],
            path=path,opts='',opt='p')

    rock = ut.load_rock(path,tag)
    max_cor,max_cors = 0,np.zeros(rock.df.index.shape)
    for i,name in enumerate(rock.df.index):
        b = rock.load(i)
        hkl = b.get_beam(refl=refl,cond='',opt=0)
        if len(refl)-len(hkl)<10:
            print(i,name)
            hkl,cor,iZ,Iz,Irpl = get_cor(b,rpl,hkl)
            max_cors[i] = cor[iZ]
            if cor[iZ]>max_cor:
                max_cor = cor[iZ]
                hkl0,cor0,iZ0,Iz0,Irpl0,i0 = hkl,cor,iZ,Iz,Irpl,i

    # dsp.stddisp([Irpl0,Iz0*Irpl0.mean()/Iz0.mean(),'bo'],labs=['$I_o$','$I_c$'])
