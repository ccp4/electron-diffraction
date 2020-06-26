from utils import*
plt.close('all')

rho   = 0.106 #atoms/A^3
sig_e = np.array([0.001,0.0025,0.005]) #A^2
le_s  = 1/(sig_e*rho) #A

z = np.linspace(0,10000,1000) #A
for i,le in zip(range(le_s.size),le_s):
    Pcoh = np.exp(-z/le)
    Pkin = z/le*np.exp(-z/le)
    Pdyn = 1-np.exp(-z/le)*(1+z/le)
    S = Pcoh+Pdyn+Pkin

    tle = r'$\sigma_e=%.3f\AA^2$, $l_e$=%.d nm' %(sig_e[i],int(le/10))
    plts = [[z/10,Pcoh,'r','$coh$'],[z/10,Pkin,'g','$kin$'],[z/10,Pdyn,'b','$dyn$'],[z/10,S,'k']]
    dsp.stddisp(plts,lw=2,labs=[r'$z(nm)$','P(z)'],title=tle,
        name='figures/Pcoh_kin_dyn%d.svg' %i, opt='ps',fonts={'leg':20,'title':20})
