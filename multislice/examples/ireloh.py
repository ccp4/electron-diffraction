import importlib as imp
from utils import*                  ; imp.reload(dsp)
import multislice.multislice as mupy; imp.reload(mupy)
import multislice.postprocess as pp ; imp.reload(pp)
import multislice.mupy_utils as mut ; imp.reload(mut)
import multislice.rotating_crystal as rcc
from crystals import Crystal
plt.close('all')
path = '../dat/ireloh/'
figpath = '../docs_fig/ireloh/'
cif_file = path+'ireloh.cif'
opts = 'FT'  # F(structure factor) E(Ewald) T(tilts) B(Beams)
hklM = 20

# mut.show_image('/home/tarik/Documents/data/ireloh/IRELOH_ED_Dataset_1/n14_a004_0484.cbf')
cbf_v = mut.Viewer_cbf(exp_path='/home/tarik/Documents/data/ireloh/IRELOH_ED_Dataset_1/',figpath=figpath)

# name = path+'ireloh001.xyz'
name = path+'ireloh484.xyz'
crys = Crystal.from_cif(cif_file)
# mut.import_cif(cif_file,name,n=[0,0,1],rep=[1,1,1])#dopt='s')

# multi = mupy.Multislice(path,data=name,#tail=istr,
#     mulslice=False,keV=200,#tilt=[tx*np.pi/180/1000,0],
#     NxNy=1024,repeat=[1,1,1],slice_thick=25,Nhk=5,
#     #TDS=True,T=300,n_TDS=15,
#     opt='sr',fopt='',v='nctr',ppopt='w',#nctrdDR',
#     #ssh=ssh,hostpath=hostpath
#     )


# padded simulation
# name = path+'IRELOH001_test.xyz'
# name = path+'IRELOH001.xyz'
pad,rep = 1,[1,1,1]
# mut.import_cif(cif_file,name,n=[0,0,1],rep=rep,pad=pad)#dopt='s')
# mut.show_grid(name,opts='xy',name=figpath+'001_xy.png',opt='p')
# mut.show_grid(path+'ireloh001.xyz',opts='xy',name=figpath+'001_xy.png',opt='p',equal=1)
# mut.show_grid(name,opts='zx',name=figpath+'001_zx.png',opt='ps')
# mut.show_grid(name,opts='zy',name=figpath+'001_zy.png',opt='ps')
# h,k = np.meshgrid(np.arange(6)*(rep[0]+2*int(pad)),np.arange(6)*(rep[1]+2*int(pad)))
# hk=[(h0,k0) for h0,k0 in zip(h.flatten(),k.flatten())];#print(hk)


# multi = mupy.Multislice(path,data=name,tail='pptest',
#     mulslice=False,keV=200,#tilt=[tx*np.pi/180/1000,0],
#     NxNy=2**12,repeat=[1,1,1],slice_thick=1,hk=hk,#Nhk=5,
#     #TDS=True,T=300,n_TDS=15,
#     opt='srp',fopt='f',v='nctr',ppopt='wuPBf',#nctrdDR',
#     ssh='badb',#hostpath=hostpath
#     )
#
multi = pp.load_multi_obj('../dat/ireloh/ireloh_t1_autoslic.pkl')
# multi = pp.load_multi_obj('../dat/ireloh/ireloh_t1_autoslic.pkl')
# multi = pp.load_multi_obj('../dat/ireloh/ireloh_pptest_autoslic.pkl')
# multi.postprocess(ssh_alias='badb',ppopt='uP')
# multi.pattern(Iopt='Incsr',out=False,tol=1e-3,
#     # gs=1.3,caxis=[-6.2,0],
#     caxis=[0,0.01],
#     rings=[0.5,1],lw=2,
#     imOpt='cv',axPos='V',cmap='gray',
#     xylims=[-1,1,-1,1], opt='ps',name=figpath+'001_SC.png')
#

def plot_transmission_vs_structure_factor():
    plts=[]
    #### Structure factor
    if 'F' in opts:
        (qxF,qyF,qzF),Fhkl = multi.get_structure_factor(hkl=None, hklMax=hklM, sym=1, v='')
        If = np.abs(Fhkl)**2 # perfect Bragg condition and flat ewald sphere

        sF = np.s_[:,hklM,hklM] # slice in the [0,0,1] direction
        If[hklM,hklM,hklM] = 1e-6
        If/=If[hklM+2,hklM,hklM]
        plts+=[[qyF[sF],np.log10(If[sF]),'g-o','$F$']]
        # dsp.stddisp(im=[ h[s],k[s],np.log10(I[s])],xylims=[-hklM,hklM,-hklM,hklM],caxis=[0,5],
        #     imOpt='cv',axPos='V',cmap='binary')

    #### transmission function
    if 'T' in opts:
        qx,qy,It = multi.pattern(Iopt='Ics',out=True,tol=1e-3,gs=1.3,caxis=[-6.2,0],rings=[0.5,1,2],lw=2,
            imOpt='cv',axPos='V',cmap='binary',opt='c')
        # multi.pattern(Iopt='Incslg',out=False,tol=1e-3,gs=1.3,caxis=[-6.2,0],rings=[0.5,1,2],lw=2,
        #     imOpt='cv',axPos='V',cmap='binary',opt='p')
        Nx = int(qx.shape[0]/2)
        st = np.s_[:,Nx] # slice in the [0,0,1] direction
        It[Nx,Nx] = 1e-6
        It/=It[Nx+2,Nx]
        plts += [[qy[st],np.log10(It[st]),'b-o','$T_{cell}$']]

        ##small slice thickness
        multi1 = mupy.Multislice(path,data=name,tail='t1',
            mulslice=False,keV=200,#tilt=[tx*np.pi/180/1000,0],
            NxNy=1024,repeat=[1,1,1],slice_thick=1,Nhk=5,
            #TDS=True,T=300,n_TDS=15,
            opt='sr',fopt='',v='nctr',ppopt='w',#nctrdDR',
            #ssh=ssh,hostpath=hostpath
            )
        qx,qy,It1 = multi1.pattern(Iopt='Ics',out=True,tol=1e-3,gs=1.3,caxis=[-6.2,0],rings=[0.5,1,2],lw=2,
            imOpt='cv',axPos='V',cmap='binary',opt='c')
        # multi.pattern(Iopt='Incslg',out=False,tol=1e-3,gs=1.3,caxis=[-6.2,0],rings=[0.5,1,2],lw=2,
        #     imOpt='cv',axPos='V',cmap='binary',opt='p')
        Nx = int(qx.shape[0]/2)
        st = np.s_[:,Nx] # slice in the [0,0,1] direction
        It1[Nx,Nx] = 1e-6
        It1/=It1[Nx+2,Nx]
        plts += [[qy[st],np.log10(It1[st]),'c-o','$T_1$']]



        dsp.stddisp(plts,labs=['$q$','$I$'],xylims=[-2,2,-10,1])

# dsp.stddisp(im=[ h[s],k[s],np.log10(T[s])],xylims=[-hklM,hklM,-hklM,hklM],caxis=[0,5],
#     imOpt='cv',axPos='V',cmap='binary')
