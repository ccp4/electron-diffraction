from utils import*
from scattering import structure_factor as sf;imp.reload(sf)

def solve_Bloch(k0,hkl,Gs,Fhkl,v=0):
    ''' Solves fast electron wave equation with Bloch wave system
    k0   : incident wave vector amplitude
    idx  : Nx2 or Nx3 - miller indices of each beam
    Gs   : Nx2 or Nx3 - reciprocal lattice beams
    Fhkl : 2Nx+1 x 2Ny+1 - 3D Potential components (Volts)
    '''
    Ug = cst.Vg2Ug(Fhkl,k0)
    N,dim = Gs.shape                    # Number of beams and dimensions 2D/3D
    z0 = np.zeros((dim)); z0[-1] = 1    # beam direction
    K0 = k0*z0                          # Incident beam
    # V0 index
    Nx,Ny = Fhkl.shape                    # nb Fhkl components in x and y directions
    i_hkl0 = (np.array(Fhkl.shape)-1)/2   # Such that Fhkl[i_hkl0] = V0
    i_hkl0 = np.array(i_hkl0,dtype=int)
    #normalized excitation errors
    Sg = (k0**2-np.linalg.norm(K0-Gs,axis=1)**2)/(2*k0)         #;print(Sg)

    if v:print(green+'\t ...Assembling matrix with %d x %dD beams ...' %(N,dim) +black)
    # H = np.diag(Sg+0J);
    H = np.zeros((N,N),dtype=complex)
    for iG,G in enumerate(hkl) :
        #Vg is accessed as Fhkl[i,j] where i,j = G-hkl
        V_iG = np.array([Fhkl[tuple(hkl_p+i_hkl0)] for hkl_p in G-hkl]) #;print(V_iG.shape)
        V_iG = V_iG/(2*k0) #;print(iG,V_iG)
        H[iG,iG] = Sg[iG]   #diagonal terms as excitation errors
        H[iG,:] += V_iG     #off diagonal terms as potential
    if v:print(H.real)
    gammaj,CjG = np.linalg.eigh(H) #;print(red+'Ek',lk,black);print(wk)
    return gammaj,CjG


###############################################################################
#### misc
###############################################################################
def get_pattern(name,N=1,opt='',**kwargs):
    # crys = Crystal.from_database(name)
    # lat_vec  = crys.reciprocal_vectors
    # pattern  = np.array([list(a.coords_fractional)+[a.atomic_number] for a in crys.atoms])
    pptype,ax,by,cz,angle = 'p1',2,2,2,90
    pattern = np.array([[1,1,1,1]])
    lat_vec = 2*np.pi*np.diag([1/ax,1/by,1/cz])
    hkl,Fhkl = sf.structure_factor3D(pattern,lat_vec,hklMax=2*N,sym=1)
    if opt:sFact.plot_structure3D(hkl,Fhkl,opt=opt,**kwargs)
    a,b,c = lat_vec
    hkl = np.array([i[N:3*N+1,N:3*N+1,N:3*N+1].flatten() for i in hkl])
    #hkl = hkl[0].flatten(),hkl[1].flatten(),hkl[2].flatten()])
    h,k,l = hkl
    kx = a[0]*h + b[0]*k + c[0]*l
    ky = a[1]*h + b[1]*k + c[1]*l
    kz = a[2]*h + b[2]*k + c[2]*l
    Gs = np.array([kx,ky,kz]).T/(2*np.pi)
    Vg = Fhkl #.flatten()
    return hkl.T,Gs,Vg/(ax*by*cz)


def get_pattern2D(N,s,**kwargs):
    ''' Setup for Bloch wave simulation
        - N : number of beams such that Vg.shape= 2N+1 x 2N+1
        - s : np.s_[...] - slices of Vg beams for the G beams to perform blochwave calculation
    '''
    pptype,ax,cz,angle = 'p1',2,2,90
    pattern = np.array([[1,1,1]])
    lat_vec = 2*np.pi*np.diag([1/ax,1/cz])
    hl,Fhl = sf.structure_factor2D(pattern,lat_vec,hkMax=2*N,sym=1,eps=0.01)   #;print(hl)
    a,c = lat_vec
    h,l = hl[0][s],hl[1][s]
    kx = a[0]*h + c[0]*l
    kz = a[-1]*h + c[-1]*l
    Gs = np.array([kx,kz]).T/(2*np.pi)
    # Fhl = hl[0]+hl[1] # Was used to test assembling
    return np.array([h,l]).T,Gs,Fhl/(ax*cz)



########################################################################
#def : test
########################################################################
def test_Ug():
    #sanity check on units here
    k0 = 5.008343065817388      #A-1
    keV = cst.lam2keV(1/k0)     #keV
    sig = cst.keV2sigma(keV)    #kVA
    Vg = 1                      #Volts
    xig = np.pi/(sig*Vg/1000)   #A
    print(xig,k0/cst.Vg2Ug(Vg,k0))


if __name__ == '__main__':
    test_Ug()
    # print('bloch wave')
