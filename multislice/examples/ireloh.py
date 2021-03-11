# from utils.handler3D import handler_3d
import importlib as imp
from utils import*                          ; imp.reload(dsp)
import multislice.multislice as mupy        ; imp.reload(mupy)
import multislice.postprocess as pp         ; imp.reload(pp)
import multislice.mupy_utils as mut         ; imp.reload(mut)
import multislice.rotating_crystal as rcc   ; imp.reload(rcc)
from crystals import Crystal
plt.close('all')
path = '../dat/ireloh/'
figpath = '../docs_fig/ireloh/'
cif_file = path+'ireloh.cif'
opts = 'FT'  # F(structure factor) E(Ewald) T(tilts) B(Beams)
hklM = 20

rotx = lambda a:np.array([[1,0,0],[0,np.cos(a),np.sin(a)],[0,-np.sin(a),np.cos(a)]])
rotz = lambda a:np.array([[np.cos(a),np.sin(a),0],[-np.sin(a),np.cos(a),0],[0,0,1]])

imV = mut.Image_viewer('/home/tarik/Documents/data/ireloh/IRELOH_ED_Dataset_1/n14_a004_0521.cbf')
fig,ax = imV.show_image(stack=5,caxis=[0,50],cmap='binary',opt='p')

# mut.show_image('/home/tarik/Documents/data/ireloh/IRELOH_ED_Dataset_1/n14_a004_0484.cbf')
# exit()
# cbf_v = mut.Viewer_cbf(exp_path='/home/tarik/Documents/data/ireloh/IRELOH_ED_Dataset_1/',figpath=figpath)

rsa = np.array([
    [  5.75197579,   1.08642157,   5.77468802],
    [ -1.09493538,  -9.45325930,   2.86911892],
    [ 12.49046736,  -4.94063233, -11.5118351 ],])
rsa_ = mut.get_reciprocal(rsa)
# rsa_ :  np.array([
#  [ 0.09721788 ,  0.05915872,  0.06860005],
#  [-0.01449509 , -0.07170057,  0.07450839],
#  [ 0.04263771 , -0.03194969, -0.0268509 ],])
lat_rec = np.linalg.norm(rsa_,axis=1)

# get the indexing matrix A-1
A = np.array([
[0.08540147986530082, -0.010093029325705745, 0.03984388435141725],
[0.016954218646467147, -0.09572669491816205, -0.015506584465560233],
[0.08490755747184861, 0.029266310687469856, -0.036979330900490597],
])
Ainv = np.linalg.inv(A)

#identified reflections
# xh0s = np.array([[-0.05,0.63],[-0.25,0.14],[-0.48,0.16,0],[-0.16,-0.37,0],[-0.49,-0.18,0]])
xh0s = np.array([[-0.49,-0.18],[-0.16,-0.37],[-0.558,0.329],[0.311,0.257],[-0.25,0.14]])
xh0s = np.hstack([xh0s,np.zeros((xh0s.shape[0],1))])*0.3
x0,y0 = 0.07*0.3,-0.07*0.3 #270,270,
xh = np.array(xh0s)-np.array([x0,y0,0])

phi = np.deg2rad(-29.014)
xh_r = rotx(-phi).dot(xh.T).T
# h = rsa_.dot(xh)/lat_rec**2;#print(h)
hs = Ainv.dot(xh_r.T).T
hstr = [ ','.join(['%.1f' %i for i in h ]) for h in hs]


#plot reflections and axi
txts = [[x[0]+x0,x[1]+y0,l,c] for x,l,c in zip(rsa_,['$a^{*}$','$b^{*}$','$c^{*}$'],['r']*3)]
arrs = [[x0,y0,x[0],x[1],c] for x,c in zip(rsa_,['r']*3)]
txts += [[xh0[0],xh0[1],'(%s)' %h_str,'k'] for h_str,xh0 in zip(hstr,xh0s)]
#resolution rings
t = np.linspace(0,2*np.pi,500)
ct,st = np.cos(t),np.sin(t)
plts = [[R*ct,R*st,'r',''] for R in 1/np.array([1.32,0.93,0.76,0.66])]
#display
imV = mut.Image_viewer('/home/tarik/Documents/data/ireloh/IRELOH_ED_Dataset_1/n14_a004_0521.cbf')
# fig,ax = imV.show_image(stack=5,caxis=[0,50],cmap='binary',opt='')
# dsp.stddisp(plts,ax=ax,fig=fig,texts=txts,arrows=arrs,labs=['$q_x$','$q_y$'],
#     opt='p',xylims=0.2)

xy_px = np.array(imV.s2px(xh0s)).T
print('xy_px:',xy_px)
print('xy_mm :',xy_px*imV.pxy*1000)


def show_setup(name,gopts=['xy','yz','xz'],**kwargs):
    xyz_name = dsp.basename(name)
    for opts in gopts:
        figname=figpath+xyz_name.replace('.xyz','')+opts+'.png'
        mut.show_grid(name,opts,name=figname,**kwargs)

    # dials.frame_orientations integrated.expt

    # n= np.array([-1249.0467359091645, 494.06323332809154, 1151.1835101484317])/200
    # n = 10*np.array([-1.51089556, -0.51374498, 0.65074564])
    # array([-0.08540403, -0.02904403,  0.0367669 ])
    angles = np.deg2rad([134.61, 106.78, 49.40])
    crys = Crystal.from_cif(cif_file)
    lat_vec = np.array(crys.lattice_vectors)
    lat_params = np.array(crys.lattice_parameters[:3])


    n= np.array([-0.08540851, -0.02904118, 0.03678561])*lat_params
    ## original structure
    # name = path+'ireloh111.xyz'
    # rcc.import_cif(cif_file,xyz=name)
    #
    # fig,ax = mut.show_grid3(name,opt='',ms=50)
    # rcc.show_unit_cell(ax,a=0.2,c='b',lw=2,lat_params=crys.lattice_parameters[:3])
    # rcc.show_trihedron(ax,uvw=np.array([n]) ,x0=[0,0,0],cs=['k--'],clabs=['beam'],lw=3,rc=0.1,h=0.2,opt='')
    # rcc.show_trihedron(ax,uvw=lat_vec,x0=[0,0,0],cs=['r--','g--','b--'],clabs=['$a$','$b$','$c$'],lw=3,rc=0.1,h=0.2,
    #     labs=['$x$','$y$','$z$'],xylims=[-10,15,-10,15,-5,20],view=[10,-20],name=figpath+'1_484_xyz.png',opt='p')
    #
    # ## rotated structure
    ez=[0,0,-1]

    uvw  = rcc.orient_crystal(lat_vec,n_u=n/np.linalg.norm(n),ez=ez)
    n_rot= rcc.orient_crystal(n      ,n_u=n/np.linalg.norm(n),ez=ez)

    np.rad2deg(np.arccos(uvw.dot([0,0,-1])/lat_params))


    a = 0
    # a = np.deg2rad(225)

    uvw = np.array([rotz(a).dot(u) for u in uvw])
    # n_rot = rotz(a).dot(n_rot)

    mut.import_cif(cif_file,xyz=name,n=n)
    #
    # fig,ax = mut.show_grid3(name,opt='',ms=50)
    # rcc.show_unit_cell(ax,n=n,a=0.2,c='b',lw=2,lat_params=crys.lattice_parameters[:3])
    # rcc.show_trihedron(ax,uvw=np.array([n_rot]) ,x0=[0,0,0],cs=['k'],clabs=['beam'],lw=3,rc=0.1,h=0.2,opt='')
    # rcc.show_trihedron(ax,uvw=uvw    ,x0=[0,0,0],cs=['r','g','b'],clabs=['$a$','$b$','$c$'],lw=3,rc=0.1,h=0.2,
    #     labs=['$x$','$y$','$z$'],xylims=[-10,20,-10,20,-10,20],view=[20,-100],name=figpath+'1_484_rot.png',opt='p')


    # fig,ax = mut.show_grid3(name,opt='',ms=50)
    # name = path+'ireloh111_rotated484.xyz'
    # rcc.show_unit_cell(ax,n=n,a=0.2,c='b',lw=2,lat_params=crys.lattice_parameters[:3])
    # rcc.show_trihedron(ax,uvw=np.array([n_rot]) ,x0=[0,0,0],cs=['k'],clabs=['beam'],lw=3,rc=0.1,h=0.2,opt='')
    # rcc.show_trihedron(ax,uvw=uvw    ,x0=[0,0,0],cs=['r','g','b'],clabs=['$a$','$b$','$c$'],lw=3,rc=0.1,h=0.2,
    #     labs=['$x$','$y$','$z$'],xylims=[-10,20,-10,20,-10,20],view=[20,10],name=figpath+'1_484_dav.png',opt='p')

    rsa = np.array([
        [  5.75197579,   1.08642157,   5.77468802],
        [ -1.09493538,  -9.45325930,   2.86911892],
        [ 12.49046736,  -4.94063233, -11.5118351 ],])
    rsa_ = mut.get_reciprocal(rsa)

    # name = path+'ireloh111_SRFU.xyz'
    a,b,c = crys.lattice_parameters[:3]
    fig,ax=dsp.stddisp(rc='3d')
    # fig,ax = mut.show_grid3(name,opt='',ms=50)
    # rcc.show_unit_cell(ax,n=n,ez=ez,a=0.2,c='b',lw=2,lat_params=[a,b,c])
    # rcc.show_unit_cell(ax,a=0.05,c='b',lw=2,lat_params=[a,b,c],ls='--')
    rcc.show_trihedron(ax,uvw=np.array([10*n,10*n_rot,[10,0,0]]) ,x0=[0,0,0],cs=['k-.','k--','k-'],clabs=['beam','$b_z$','phi'],lw=3,rc=0.1,h=0.2)
    # rcc.show_trihedron(ax,uvw=np.array() ,x0=[0,0,0],cs=['k'],clabs=['beam'],lw=3,rc=0.1,h=0.2,opt='')
    # rcc.show_trihedron(ax,uvw=rsa    ,x0=[0,0,0],cs=['r-','g-','b-']   ,clabs=['$a$','$b$','$c$'],lw=3,rc=0.1,h=0.2)
    rcc.show_trihedron(ax,uvw=rsa_   ,x0=[0,0,0],cs=['r--','g--','b--'],clabs=['$a^{*}$','$b^{*}$','$c^{*}$'],lw=3,rc=0.1,h=0.2)
    # rcc.show_trihedron(ax,uvw=uvw1   ,x0=[0,0,0],cs=['r','g','b'],clabs=['$a1$','$b1$','$c1$'],lw=3,rc=0.1,h=0.2)
    # rcc.show_trihedron(ax,uvw=uvw    ,x0=[0,0,0],cs=['r','g','b'],clabs=['$a$','$b$','$c$'],lw=3,rc=0.1,h=0.2,

    dsp.stddisp(ax=ax,fig=fig,labs=['$x$','$y$','$z$'],xylims=np.array([-10,20,-10,20,-10,20]),view=[20,10],name=figpath+'1_484_dav.png',opt='p')

    # hdl = handler_3d(fig,persp=False)

# name = path+'ireloh001.xyz'
# name = path+'ireloh484.xyz'
# name = path+'IRELOH001_test.xyz'
# name = path+'IRELOH001.xyz'
# name = path+'ireloh_rotated484.xyz'
name = path+'ireloh111_rotated484.xyz'
# name = path+'ireloh_rotated900.xyz'
# name = path+'ireloh111_rotated484_TD.xyz'

# show_setup(name,gopts=[],opt='sc')


def xyz_gif(rpath,dpath):
    rpath = '/home/tarik/Documents/git/ccp4/reports/ireloh/xyz/'
    dpath='/home/tarik/Documents/git/ccp4/src/david_scripts/ED_processing/IRELOH_ED_Dataset_1-dials/'
    for i in range(500,1351,50) :
        print(i)
        name = dpath+'ireloh_rotated%d.xyz' %i
        mut.show_grid(name,'xy',name=rpath+'xy/1_xy%s.png' %(str(i).zfill(4)),opt='sc',xylims=[100,450,100,450])
        mut.show_grid(name,'yz',name=rpath+'yz/1_yz%s.png' %(str(i).zfill(4)),opt='sc',xylims=[100,500,0,310])
        mut.show_grid(name,'xz',name=rpath+'xz/1_xz%s.png' %(str(i).zfill(4)),opt='sc',xylims=[100,450,0,310])


# mut.show_grid(name,'xz')#,xylims=[100,450,100,450],name=figpath+'ireloh_rotated484_xy.png',opt='sc')



def run_simu(name):
    pad,rep = 1,[1,1,1]
    # mut.import_cif(cif_file,name,n=[0,0,1],rep=rep,pad=pad)#dopt='s')

    # h,k = np.meshgrid(np.arange(-6,7)*(rep[0]+2*int(pad)),np.arange(-6,7)*(rep[1]+2*int(pad)))
    h,k = np.meshgrid(np.arange(-30,31),np.arange(-30,31))
    hk=[(h0,k0) for h0,k0 in zip(h.flatten(),k.flatten())];#print(hk)

    multi = mupy.Multislice(path,data=name,tail='900',
        mulslice=False,keV=200,#tilt=[tx*np.pi/180/1000,0],
        NxNy=2**12,repeat=[1,1,1],slice_thick=1,hk=hk,#Nhk=5,
        #TDS=True,T=300,n_TDS=15,
        opt='srp',fopt='f',v='nctr',ppopt='wuPBf',#nctrdDR',
        ssh='badb',#hostpath=hostpath
        )
# run_simu(name)


def postpp(multi):
    fig,ax = multi.pattern(Iopt='Incs',out=False,Nmax=256,#tol=1e-3,
        # gs=1.3,caxis=[-6.2,0],
        caxis=[0,0.001],
        rings=[0.25,0.5],lw=2,
        imOpt='cv',axPos='V',cmap='gray',
        xylims=0.4*np.array([-1,1,-1,1]),name=figpath+'900_I.png',opt='')

    qx,qy,I = multi.pattern(Iopt='Incs',out=True,Nmax=256,opt='c')

    Nmax = 256
    q = np.sqrt((qx-qx.min())**2+(qy-qy.min())**2)

    nx,ny,nd = 2,2,4
    ss = [np.s_[Nmax+int(i*Nmax/nx),:] for i in range(-nx+1,nx)]
    ss+= [np.s_[:,Nmax+int(i*Nmax/ny)] for i in range(-ny+1,ny)]
    ss+=[np.s_[range(2*Nmax-1,-1,-1),range(2*Nmax)],np.s_[range(2*Nmax-1,-1,-1),range(2*Nmax-1,-1,-1)]]
    sl = ['$q_y=%d$' %i for i in range(-nx+1,nx)]+['$q_x=%d$' %i for i in range(-ny+1,ny)]+['$q_x=-q_y$','$q_x=q_y$']

    sc=dsp.getCs('Spectral',len(ss))
    plts=[[qx[s],qy[s],c,l] for s,l,c in zip(ss,sl,sc)]
    dsp.stddisp(plts,fig=fig,ax=ax,lw=2)

    dsp.stddisp([[q[s],I[s],c,l] for s,l,c in zip(ss,sl,sc)],lw=2,xylims=['y',0,0.001])

    # multi.beam_vs_thickness(iBs=range(10),cm='Spectral',name=figpath+'900_B.svg',opt='sp')
    # multi.beam_vs_thickness(tol=1e-3,cm='Spectral',name=figpath+'900_B.svg',opt='sp',xylims=['y',0,0.005])

# multi = pp.load_multi_obj('../dat/ireloh/ireloh_t1_autoslic.pkl')
multi = pp.load_multi_obj('../dat/ireloh/ireloh_484_autoslic.pkl')
# multi = pp.load_multi_obj('../dat/ireloh/ireloh_900_autoslic.pkl')
# multi = pp.load_multi_obj('../dat/ireloh/ireloh_t1_autoslic.pkl')
# multi = pp.load_multi_obj('../dat/ireloh/ireloh_pptest_autoslic.pkl')
# multi.postprocess(ssh_alias='badb',ppopt='uB')


# postpp(multi)



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
