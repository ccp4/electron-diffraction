import blochwave,time,sys,pytest,os,copy
from subprocess import check_output
from pytest_html import extras
from utils import glob_colors as colors
from utils import*                      ;imp.reload(dsp)
from utils import pytest_util           ;imp.reload(pytest_util)
from blochwave import bloch             ;imp.reload(bloch)
from blochwave import util as bloch_util;imp.reload(bloch_util)
from EDutils import pets as pets_imp    ;imp.reload(pets_imp)
from EDutils import utilities as ut     ;imp.reload(ut)

plt.close('all')

out,ref,dir = pytest_util.get_path(__file__)

pts_path = 'glycine/glycine.pts'
F=3 # frame at which comparison is done
thick = 330
opts = 'k'

pets = ut.load_pkl('glycine/pets.pkl')


def read_fortran_csv(file):
    csv = file.replace('.txt','.csv')
    cmd ="sed -E 's/ +/,/g' %s | sed 's/^,//' | sed 's/,$//' > %s" %(file,csv)
    # print(cmd)
    check_output(cmd,shell=True)
    df = pd.read_csv(csv,sep=",",engine='python',index_col=None,header=None)
    return df


if 'r' in opts:
    pets.make_eldyn(F=3,thickness=thick,phisteps=1,thr=4,
        bargs={'Nmax':4},v=0)
    pets.run_dyngo()
    # pets.dyngo_path=os.path.join(pets.path,'dyngo')
    b0 = pets.run_blochwave(F=F,Nmax=4,Smax=0.02)

df = pets.get_edout(F=F,opts='')
    ### compare
    # cols = ['Sw']



v = 1
if opts:
    hkl0 = [-4,-1,-3]
    ### compare orientation
    phi=1.816000
    alpha,beta,scale = pets.dyn.loc[F,['alpha','beta','scale']].values
    Rx = lambda a:(lambda ca,sa: np.array([[1,0,0],[0,ca,-sa],[0,sa,ca]]))(np.cos(np.deg2rad(a)),np.sin(np.deg2rad(a)))
    Ry = lambda a:(lambda ca,sa: np.array([[ca,0,sa],[0,1,0],[-sa,0,ca]]))(np.cos(np.deg2rad(a)),np.sin(np.deg2rad(a)))
    Rmat= Rx(alpha+0*phi/2).dot(Ry(beta))
    R = Rmat.dot(pets.A)

    #### check K in lab frame
    Klat = -pets.dyn.loc[F,['u','v','w']].values
    Kcar = pets.lat_vec.T.dot(Klat)
    Krec = np.linalg.inv(pets.lat_vec1.T).dot(Kcar)
    Klab = R.dot(Krec)
    Ruvw = R.dot(np.linalg.inv(pets.lat_vec1.T).dot(pets.lat_vec.T))
    if 'k' in opts:
        # Klab = np.linalg.inv(Rcar).dot(Kcar)
        print('K_lab : ', Klab/np.linalg.norm(Klab))

        #getting u,v,w from A
        print(Klat)
        u = np.linalg.inv(Ruvw).dot([0,0,-1])
        u /= u[1]
        print('K_lat : ',u)
    #### compare g coordinates
    glab = R.dot(hkl0)
    if 'g' in opts:
        # hkl = df[['h','k','l']].values
        print('g_lab : ', glab)


    #### compare excitation errors
    if 'w' in opts:
        print(colors.blue+'...comparing excitation errors...'+colors.black)
        b0 = pets.load_b0()
        k0 = b0.k0
        Klab*=k0/np.linalg.norm(Klab)
        if v:print('K_lab : ', Klab)
        if v:print('g_lab : ', glab)
        # sg_lab = (k0**2 - np.linalg.norm(glab+np.array([0,0,-k0]))**2)/(2*k0)
        sg_lab = (k0**2 - np.linalg.norm(glab+Klab)**2)/(2*k0)
        gcar = pets.lat_vec1.T.dot(hkl0)
        Kcar*=k0/np.linalg.norm(Kcar)
        sg_car = (k0**2 - np.linalg.norm(gcar+Kcar)**2)/(2*k0)
        if v:print('Sg_lab=%.6f'%sg_lab)
        if v:print('sg_car=%.6f'%sg_car)

        # if v:print('gcar bw')
        # if v:print(b0.df_G.loc[str(tuple(hkl0))][['qx','qy','qz']])

        cols = ['Sw','Sw_bw']
        df_allint=pets.read_allint(F)
        hkl  = df_allint.index
        b0._set_excitation_errors(Smax=0,hkl=hkl)
        # df['Sw_bw']=np.inf
        df_allint.loc[hkl,'Sw_bw']=b0.df_G.loc[hkl,'Sw']
        print(df_allint.loc[hkl,['Sw','Sw_bw']])

    if 'M' in opts or 'I' in opts:

        mat_file = 'glycine/dyngo/mat_back.txt'
        mat = read_fortran_csv(mat_file)
        # csv = mat.replace('.txt','.csv')
        # cmd ="sed -E 's/ +/,/g' %s | sed 's/^,//' | sed 's/,$//' > %s" %(mat,csv)
        # print(cmd)
        # check_output(cmd,shell=True)
        # mat = pd.read_csv(csv,sep=",",engine='python',index_col=None,header=None)
        hkls = mat[[0,1,2]].values
        # hkl =[str(tuple(h)) for h in hkls]
        mat.index = [str(tuple(h)) for h in hkls]
        mat = mat.drop([0,1,2],axis=1)
        H = mat.values[:,::2] + 1J*mat.values[:,1::2]
        Hd = H.copy().T
        np.fill_diagonal(Hd,0)
        H += Hd

        # df.loc[str((0,0,0))] = 0
        # df.loc[hkl,'Sw_mat'] = np.diag(v1)
        b0 = pets.load_b0()
        b0.solve(Nmax=2*hkls.max(),thick=thick,hkl=hkls,
            dyngo_args={'Rmat':Rmat,'scale':scale})
        H_bw = b0.H*2*b0.k0/(2*np.pi)

        def printH(H):
            print(abs(H).min(),abs(H).max(),abs(H).mean() )
            return abs(H).argmax()
        # printH(H)
        # printH(H_bw)
        #### check differences
        tol=1e-3
        idx = (abs(H_bw)>tol) & (abs(H)>tol)
        idmax=printH((H_bw[idx]-H[idx])/H_bw[idx])
        #print(H_bw[idx][idmax],H[idx][idmax]) #a bit too strong
        i = np.unravel_index(abs(H_bw-H).argmax(),H.shape)
        print(H_bw[i],H[i],abs(H_bw[i]-H[i]))

        # np.allclose(b0.H/(2*np.pi),H)
        # idx,hkl_bw = np.array([ [i,h] for i,h in enumerate(b0.df_G.index) if h in df.index]).T
        # df.loc[hkls,'Sw_bw'] = np.inf
        # hkl = mat.index
        mat['Sw_bw'] = np.diag(H_bw.real)#[np.array(idx,dtype=int)]
        mat['Sw']    = np.diag(H).real
        mat['U_bw']  = np.sum(H_bw.real)#[np.array(idx,dtype=int)]
        mat['U']     = np.sum(H).real
        print(mat[['Sw_bw','Sw','U_bw' ,'U']])
        # print(mat[['U_bw' ,'U']])

        if 'e' in opts:
            #eigen
            df_eig = read_fortran_csv(file='glycine/dyngo/eig.txt')
            # b0 = pets.load_b0()
            # b0.solve(Nmax=13,Smax=0.05,thick=thick)
            # df['I_bw']=np.inf
            eig = df_eig.values
            gj=eig[:,0]
            Cj=eig[:,2::2]+1J*eig[:,3::2]
            mat['g']    = gj
            mat['g_bw'] = b0.gammaj*2*b0.k0/(2*np.pi)
            mat['C']    = abs(Cj).sum(axis=1)
            mat['C_bw'] = abs(b0.CjG).sum(axis=1)
            print(mat[['Sw_bw','Sw','U_bw' ,'U','g','g_bw','C','C_bw']])


        if 'I' in opts:
            df_allint=pets.read_allint(F,'O')
            hkl = df_allint.index
            mat.loc[hkl,['S']]  = df_allint.loc[hkl,'S1'] + 1J*df_allint.loc[hkl,'S2']
            mat.loc[hkl,'I']    = abs(mat.loc[hkl,'S'])**2
            mat.loc[hkl,'S_bw'] = b0.df_G.loc[hkl,'S']
            mat.loc[hkl,'I_bw'] = b0.df_G.loc[hkl,'I']
            print(mat[['S_bw','S','I','I_bw']])

            hkl = [ h for h in df.index if h in b0.df_G.index]
            df.loc[hkl,'I_bw']=b0.df_G.loc[hkl,'I']
            print(df.loc[hkl,['I','I_bw']])
            # print(b0.df_G.sort_values('I')['I'][-10:])
            # print(df.sort_values('I')['I'][-10:])

###
# eldyn = make_eldyn(pets,F=3)
# run_dyngo(eldyn)
# df = get_edout(eldyn)
