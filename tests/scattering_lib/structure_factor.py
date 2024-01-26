from utils import*
import io
from subprocess import check_output
from crystals import Crystal
from scattering import structure_factor as sf;imp.reload(sf)
from scattering import scattering_factors as scatf;imp.reload(scatf)
opts='G'
ptype='electron'
# ptype='xray'

cif_file='dat/LTA.cif'

crys = Crystal.from_cif(cif_file)
pattern = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in crys.atoms] )
lat_vec = np.array(crys.reciprocal_vectors)
df_Fhkl=sf.get_structure_factor(cif_file,hklMax=2,v='',ed=ptype=='electron')
df_Fhkl['amp']=np.abs(df_Fhkl.F)
df_Fhkl['hkl']=df_Fhkl.index


if 'G' in opts:
    cmd='gemmi sfcalc --wavelength={lam} --for={ptype} --dmin=0.9 {cif_file}'.format(
        cif_file=cif_file,
        lam={'electron':0.025,'xray':1.54059}[ptype],
        ptype=ptype,)


    out=check_output(cmd,shell=True).decode()
    df_FhklG=pd.read_csv(io.StringIO(out),names=['hkl','amp','phi'],sep="\t",index_col=0)
    df_FhklG['Fc2'] = np.abs(df_FhklG.amp*np.exp(1J*np.deg2rad(df_FhklG.phi)))**2
    hkl = [eval(','.join(hkl.split(" "))[1:]) for hkl in df_FhklG.index]
    df_FhklG['h']= [h[0] for h in hkl]
    df_FhklG['k']= [h[1] for h in hkl]
    df_FhklG['l']= [h[2] for h in hkl]
    df_FhklG.index=[str(h) for h in hkl]
    df_FhklG['hkl']=df_FhklG.index
    df_FhklG['I'] = df_FhklG.Fc2

    df_gb=pd.merge(df_FhklG,df_Fhkl,on='hkl',suffixes=('_gemmi','_bloch'),how='inner')
    df_gb=df_gb.set_index('hkl')
    print(df_gb.sort_values('amp_gemmi')[['amp_gemmi','amp_bloch']])

if 'V' in opts:
    df_FhklV=pd.read_csv('dat/vesta_lta.csv',header=0)
    df_FhklV.index=[str(tuple(h)) for h in df_FhklV[['h','k','l']].values]
    df_FhklV['hkl']=df_FhklV.index
    df_FhklV.rename(columns={'|F|':'amp'},inplace=True)

    df1=pd.merge(df_FhklG,df_FhklV,on='hkl',suffixes=('_gemmi','_vesta'),how='inner')
    df1=df1.set_index('hkl')
    print(df1.sort_values('I_vesta')[['amp_gemmi','amp_vesta','I_gemmi','I_vesta']])

if 'd' in opts:
    crys = Crystal.from_database('diamond')
    pattern = np.array([np.hstack([a.coords_fractional,a.atomic_number]) for a in crys.atoms] )
    lat_vec = np.array(crys.reciprocal_vectors)
    hkl,Fhkl=sf.structure_factor3D(pattern,lat_vec,hklMax=2,sym=0,v='')
    h,k,l = hkl
    I = np.abs(Fhkl)**2

    id0 = ~((h%2==0) & (k%2==0) & (l%2==0)) & ~((h%2==1) & (k%2==1) &(l%2==1))

    print(np.abs(I[id0]))                   #should be 0
    print(np.abs(I[~id0 & ((h+k+l)%4==0)])) #should be 64
    print(np.abs(I[~id0 & ((h+k+l)%2==1)])) #should be 32
    print(np.abs(I[~id0 & ((h+k+l)%4==2)])) #should be 0
