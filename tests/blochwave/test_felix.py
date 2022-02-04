from utils import*
from crystals import Crystal
from blochwave import bloch as bl;imp.reload(bl)
from scattering import scattering_factors as scatf  ;imp.reload(scatf)
from utils import pytest_util
import pytest,os,copy
plt.close('all')
keV=200

try:
    out,ref,dir = pytest_util.get_path(__file__)
    path=out+'/felix/'
    b = bl.Bloch(cif_file='GaAs',u=[-1,1,0],keV=keV,Nmax=10,Smax=0.2,path=out,solve=0)

    cif='GaAs_felix.cif'
    if not os.path.exists(path+'eigenvals.txt'):
        bf=copy.copy(b)
        bf._solve_Felix(cif,nbeams=200,thicks=(10,20,10))

    A   = np.loadtxt(path+'eigenvals.txt')
    hkl = np.array(A[:,:3],dtype=int)
    hkl_str = [str(tuple(h)) for h in hkl]
    b.solve(hkl=hkl,Smax=0,Nmax=10)#Smax=0.02,Nmax=7)
    # idx = b.get_beam(refl=hkl_str)
except Exception as e:
    print(e)
    pass


@pytest.mark.new
@pytest_util.cmp_ref(__file__)
def test_solve_felix():
    bf=copy.copy(b)
    bf._solve_Felix(cif,nbeams=200,thicks=(10,250,10))
    return bf.gammaj

@pytest.mark.lvl2
def test_coords():
    dfa = pd.read_csv(path+'atoms.txt',sep="  *",names=['Z','x','y','z'],engine='python')
    crys = b.crys
    xyz  = pd.DataFrame()
    xyz['Z'] = np.array([a.atomic_number for a in crys.atoms])
    xyz[['x','y','z']] = [a.coords_cartesian for a in crys.atoms] #,columns=['x','y','z'])

# @pytest.mark.lvl2
# def test_scattering_factors():
#     abcd = np.load('/home/tarik/Documents/git/ccp4/src/electron-diffraction/scattering/data/abcd.npy')
#     aGa,aAs = abcd[[31,33]]
#     a_Ga = np.array([6.255284640e-01,1.100056500e-01,2.053029010e00,2.410957860e00,2.896081200e-01,4.786857360e01,
#       2.079105940e-01,3.278072240e-01,3.450796170e-01,7.431390610e-01,6.556342980e-03,3.094113690e-02])
#     a_As = np.array([7.778752180e-01,1.507331570e-01,5.938481500e-01,1.428822090e02,1.959187510e00,1.747503390e00,
#       1.798802260e-01,3.318008520e-01,8.632672220e-01,5.854902740e00,9.590534270e-03,2.337775690e-02])
#     print(aGa-a_Ga, aAs-a_As)

@pytest.mark.lvl1
def test_structure_factor():
    df  = pd.read_csv(path+'StructureFactors.txt',sep="  *",names=['h','k','l','qx','qy','qz','Fr','Fi'],engine='python')
    df.index=[str(tuple(hkl)) for hkl in df[['h','k','l']].values]
    df=df.loc[hkl_str] #

    df['Rg'] = np.linalg.norm(df[['qx','qy','qz']].values,axis=1)
    df['q']  = df.Rg/(2*np.pi)
    df['qb'] = b.df_G.q.loc[hkl_str]
    # Za = np.array(b.pattern[:,-1],dtype=int)
    # fj = scatf.get_fe(Za,df.q)

    df['Ug'] = b.df_G.loc[df.index,'Vg']
    # df['Ug'] = np.sum(fj,axis=1)
    df['Ufelix'] = df.Fr+1J*df.Fi
    df['Uf'] = abs(df.Ufelix)
    df['Udiff'] = abs(df.Ug-df.Ufelix)#/df.Uf
    # print(b.df_G.loc[df.index,['qx','qy','qz']])
    # print(df.iloc[:20][['Rg','Ufelix','Ug']])
    print(df.iloc[:20][['Ufelix','Ug','Udiff']])
    idx,cols = (df.Udiff>0.01) &(df.Uf>1e-8),['Ufelix','Ug','Udiff','Uf']
    print(df.loc[idx,cols])
    dmax=df.Udiff.max()
    print(dmax)
    assert dmax<1e-4

def get_H():
    H = np.loadtxt(path+'matrix.txt')/100
    H = H[:,3::2]+1J*H[:,4::2]
    return H

@pytest.mark.lvl1
def test_excitation():
    H=get_H()
    dfM = pd.DataFrame(hkl_str,columns=['hkl'],index=hkl_str)
    dfM['Sw_felix'] = np.diag(H).real
    dfM['Sw_bloch'] = b.df_G.loc[hkl_str]['Sw']*(2*np.pi)
    print('mean excitation diff:',abs(dfM.Sw_felix-dfM.Sw_bloch).mean())
    dmax = abs(dfM.Sw_felix-dfM.Sw_bloch).max()
    print(dfM);print(dmax)
    assert dmax<1e-3


@pytest.mark.lvl1
def test_matrix():
    H=get_H()
    Hdiff = abs(b.H-H)
    print('mean matrix error :' ,Hdiff.mean())
    print('max matrix error :' ,Hdiff.max())
    assert Hdiff.max()<1e-3

def get_eigen():
    g = np.loadtxt(path+'eigenvals.txt')
    C = np.loadtxt(path+'eigenvec.txt')
    g = g[:,3::2]+1J*g[:,4::2];g=g[:,0]
    C = C[:,3::2]+1J*C[:,4::2]
    return g,C


@pytest.mark.lvl1
def test_eigen():
        #felix
        g,C = get_eigen()
        idx_g  = np.argsort(g)
        g = g[idx_g]
        C = C[:,idx_g]
        #bloch
        idx_b  = np.argsort(b.gammaj)
        gammaj = b.gammaj[idx_b]
        Cj     = b.CjG[:,idx_b]

        # print('Cj.gamma.Cjdag',np.linalg.norm(Cj.dot(np.diag(gammaj)).dot(np.linalg.inv(Cj))-b.H))
        # print('C.g.Cdag',np.linalg.norm(C.dot(np.diag(g)).dot(np.linalg.inv(C))-H))
        dg = (abs(g-gammaj)/abs(g))
        dfg = pd.DataFrame(np.array([g,gammaj,dg]).T,columns=['felix','bloch','diff'])
        print('mean eigen val error percentage :',dg.mean())
        print('max eigen val error percentage :' ,dg.max())
        print(dfg)
        assert dg.mean()<5e-2

@pytest.mark.lvl1
@pytest_util.add_link(__file__)
def test_intensities():
        with open(path+'felix.inp','r') as f:l=[l.strip() for l in f.readlines()]
        z0,z1,dz=[float(s.split("=")[1]) for s in l if 'Thickness' in s]
        z = np.arange(z0,z1+1,dz)
        I = np.loadtxt(path+'intensities.txt')[:,3:]
        idx = np.arange(1,10)
        I=I[idx]            #;print(I.shape)
        #### force calculations with felix eigen solutions
        # g,C=get_eigen()
        # b.gammaj=g
        # b.CjG=C
        # print(g.shape,C.shape)
        # b.invCjG=np.linalg.inv(C)

        b._set_beams_vs_thickness(thicks=z)
        Ib = b.get_beams_vs_thickness(idx=idx)
        # Rfz = np.abs(np.abs(I).sum(axis=0)-np.abs(Ib).sum(axis=0))/np.abs(Ib).sum(axis=0)
        # dsp.stddisp([z,Rfz,'b-'],labs=['$z$','Rf'])
        # print('max Rfactor error :',np.max(Rfz))
        # idI = I>1e-4
        print(Ib.shape)
        dI = np.abs(Ib-I)
        print('error in I : max=%.2f, mean=%.2E' %(dI.max(),dI.mean()))
        # print('diff with ref:',abs(Rfz-np.load(path+'Rfz.npy')).sum())

        zA = 250

        idz=abs(z-zA).argmin()#;print(idz)
        dI = np.array([(np.abs(I1-I0)[I0>1e-3]/I0[I0>1e-3]).max() for I0,I1 in zip(I[:,:idz],Ib[:,:idz])])
        print(dI)
        print('max err I over %dA : %.2E' %(zA,dI.max()))
        # assert dI.max()<0.05

        cs = dsp.getCs('jet',I.shape[0])
        plts = [[z,If,[c,'--o'],''] for If,c in zip(I,cs)]
        plts += [[z,Ib,[c,'-x'],''] for Ib,c in zip(Ib,cs)]
        legElt = {'bloch':'k-x','felix':'k--o'}
        return dsp.stddisp(plts,lw=2,legElt=legElt,
            xylims=['x',0,zA],opt='p')

# test_excitation()
# test_matrix()
# test_structure_factor()
# test_eigen()
# test_solve_felix()
# fig,ax=test_intensities();fig.show()
