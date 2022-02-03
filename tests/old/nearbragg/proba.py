from utils import*
from math import factorial as fact
path='../../figures/dynamical_diffraction/Pcoh_kin_dual_dyn'
plt.close('all')

Nmax = 1000
N   = np.arange(Nmax)+1
sig = 0.01
ncol = 6

P = np.zeros((ncol,Nmax))
P[4] = sig**5*(N-3)*(N-4)*(N-5)*(N-6)*(N-7)/(5*4*3*2)/2
P[4] = sig**4*(N-2)*(N-3)*(N-4)*(N-5)/(4*3*2)/2
P[3] = sig**3*(N-1)*(N-2)*(N-3)/(3*2)/2
P[2] = sig**2*N*(N-1)/2 - P[3]
P[1] = sig*N-P[2]
P[0] = 1-P[1]

cs = dsp.getCs('jet',ncol)
plts=[[N,Pi,cs[i],'$P_%d$' %i] for i,Pi in enumerate(P)]
# plts+= [[N,P0,'b','P0']]
# plts+= [[N,P1,'r','P1']]
# plts+= [[N,P2,'g','P2']]
dsp.stddisp(plts,lw=2,labs=[r'$N$','P(z)'])

def scat(topt=3):
    rho   = 0.106   #atoms/A^3
    sig_e = np.array([0.001,0.0025,0.005]) #A^2
    le_s  = 1/(sig_e*rho) #A

    z = np.linspace(0,10000,1000) #A

    cs=dsp.getCs('jet',topt+2)
    l_e  = le_s[1]
    y = [1/fact(m)*(z/l_e)**m*np.exp(-z/l_e) for m in range(topt+1)]
    plts = [ [z/10,y[m],cs[m],'$%d$' %m] for m in range(topt+1)]
    plts+= [[z/10,1-np.array(y).sum(axis=0),'k--','']]
    dsp.stddisp(plts,lw=2,labs=[r'$z(nm)$','P(z)'],#title=tle,
        name=path+'%d.svg' %0, opt='p',fonts={'leg':20,'title':20})


def cross_section_sweep():
    for i,le in zip(range(le_s.size),le_s):
        Pcoh = np.exp(-z/le)
        Pkin = z/le*np.exp(-z/le)
        Pdual = 0.5*(z/le)**2*np.exp(-z/le)
        if topt:
            Ptri = 1/6*(z/le)**3*np.exp(-z/le)
            Pdyn = 1-np.exp(-z/le)*(1+z/le+0.5*(z/le)**2+(z/le)**3/6)
            S = Pcoh+Pkin+Pdual+Ptri+Pdyn
        else:
            Pdyn = 1-np.exp(-z/le)*(1+z/le+0.5*(z/le)**2)
            S = Pcoh+Pkin+Pdual+Pdyn

        tle = r'$\sigma_e=%.3f\AA^2$, $l_e$=%.d nm' %(sig_e[i],int(le/10))
        cs=dsp.getCs('jet',4)
        plts = [[z/10,Pcoh,cs[0],'$coh$'],[z/10,Pkin,cs[1],'$kin$']]
        plts+= [[z/10,Pdual,cs[2],'$double$']]
        if topt: plts+= [[z/10,Ptri ,cs[3],'$tri$']]
        plts+= [[z/10,Pdyn,'k','$dyn$']]#,[z/10,S,'k--']]
        dsp.stddisp(plts,lw=2,labs=[r'$z(nm)$','P(z)'],#title=tle,
            name=path+'%d.svg' %i, opt='p',fonts={'leg':20,'title':20})
# cross_section_sweep()
scat(topt=ncol-2)
