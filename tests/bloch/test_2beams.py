from utils import*
# from blochwave import bloch2D as bl          ;imp.reload(bl)
# from multislice import mupy_utils as mut     ;imp.reload(mut)
# import multislice.multi_2D as ms             ;imp.reload(ms)
plt.close('all')
opts = 'BM' #(Bloch) (Multislice)

K = 5.024937810560445   # exact Bragg condition for hk=(1,0) beam
keV = cst.lam2keV(1/K)

def setup():
    fig,axB=dsp.create_fig()
    cmap='jet'
    thick=1000  #thickness
    eps = 0.01  #strength of potential
    n = 10      #u=[1,n]
    e0=1
    # e0=0.98,0.965
    #### Bloch
    if 'B' in opts:
        file = 'dat/p1.txt'
        # p1 = mut.import_wallpp(file)
        # p1.plot_unit_cells(opts='uAa',nh=3,nk=3)
        b0  = bl.Bloch2D(file,keV=keV,u=[1,n],thick=thick,eps=e0*eps,
            Nmax=8,Smax=0.2,solve=1)
        # b0.show_beams(opts='B',fz=np.abs)
        # b0.show_ewald()
        # b0.show_beams(opts='S')
        b0.show_beams_vs_thickness(thicks=(0,thick,1000),strong=['I'],m={'I':100},
            linespec='--',marker='',cm=cmap,ax=axB,opt='')
        # b0.G('Sw')

        idx = b0.get_Istrong(Icols=['I'],out=1)
        ibs = b0.get_refl(idx)
        print('pendulosung thickness')
        print(b0.get_Xig())
setup()
