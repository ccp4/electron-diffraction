from utils import*
plt.close('all')

opts='d'


dfs = np.arange(-100,101,50)
Cs = np.array([1,2])*1e6 #nm
keV = np.array([100,200])
lam = cst.keV2lam(keV)*0.1 #nm

u = np.linspace(0,6,1000) #nm^-1
Xi = lambda lam,df,Cs: np.sin((df*lam*u**2 + 0.5*Cs*lam**3*u**4)*np.pi)

if 'a' in opts:
    cs00 = dsp.getCs('Reds',dfs.size)
    fig,axs = dsp.create_fig(rc=[dfs.size,1],pad=0.0)
    for i,df in enumerate(dfs):
        plts00 =[[u,Xi(lam[1],df,Cs[0]),'b',r'$\Delta f=%dnm$' %df]]
        dsp.stddisp(plts00,fig=fig,ax=axs[i],
            lw=2,opt='',pOpt='ltG',fonts={'tick':10,'leg':10},legLoc='upper left')#title='$\lambda=%dkeV, Cs=%dmm$' %(keV[1],Cs[0]*1e-6))
    axs[0].set_title('$\lambda=%dkeV, Cs=%dmm$' %(keV[1],Cs[0]*1e-6),fontsize=15)
    axs[-1].set_xlabel('$u(nm^{-1})$',fontsize=15)
    fig.show()
if 'b' in opts:
    cs10 = dsp.getCs('Blues',dfs.size)
    cs01 = dsp.getCs('Greens',dfs.size)

    plts10 =[[u,Xi(lam[1],df,Cs[0]),cs10[i],r'$\Delta f=%d$' %df] for i,df in enumerate(dfs)]
    plts01 =[[u,Xi(lam[1],df,Cs[1]),cs01[i],r'$\Delta f=%d$' %df] for i,df in enumerate(dfs)]
    dsp.stddisp(plts10,lw=2,labs=['$u(nm^{-1})$',laby],title='$\lambda=%dkeV, Cs=%dmm$' %(keV[1],Cs[0]*1e-6))
    dsp.stddisp(plts01,lw=2,labs=['$u(nm^{-1})$',laby],title='$\lambda=%dkeV, Cs=%dmm$' %(keV[1],Cs[1]*1e-6))

if 'c' in opts:
    dfs = np.linspace(-100,-50,10)
    for i,df in enumerate(dfs):
        plts00 =[[u,Xi(lam[1],df,Cs[0]),'b',r'$\Delta f=%d$' %df]]
        dsp.stddisp(plts00,lw=2,labs=['$u(nm^{-1})$','$\sin\chi(u)$'],title='$\lambda=%dkeV, Cs=%dmm$' %(keV[1],Cs[0]*1e-6))

if 'd' in opts:
    df = -65
    betas = np.array([5,10])
    us = (betas*1e-3)/lam[1]
    plts00 =[[u,Xi(lam[1],df,Cs[0]),'b',r'$\Delta f=%dnm$' %df]]
    plts00+=[[[u,u],[-1,1],'r',''] for u in us]
    dsp.stddisp(plts00,lw=2,labs=['$u(nm^{-1})$','$\sin\chi(u)$'],title='$\lambda=%dkeV, Cs=%dmm$' %(keV[1],Cs[0]*1e-6))
