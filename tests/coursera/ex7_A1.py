from utils import*
plt.close('all')

opts = '1'

if '1' in opts:
    dfs = np.arange(-100,101,50)
    Cs = np.array([1,2])*1e6 #nm
    keV = np.array([100,200])
    lam = cst.keV2lam(keV)*0.1 #nm

    u = np.linspace(0,6,1000) #nm^-1

    Xi = lambda lam,df,Cs:(df*lam*u**2 + 0.5*Cs*lam**3*u**4)*np.pi
    if 'S' in opts:Xi = lambda lam,df,Cs: np.sin((df*lam*u**2 + 0.5*Cs*lam**3*u**4)*np.pi)

    cs00 = dsp.getCs('Reds',dfs.size)
    cs10 = dsp.getCs('Blues',dfs.size)
    cs01 = dsp.getCs('Greens',dfs.size)

    plts00 =[[u,Xi(lam[0],df,Cs[0]),cs00[i],r'$\Delta f=%d$' %df] for i,df in enumerate(dfs)]
    plts10 =[[u,Xi(lam[1],df,Cs[0]),cs10[i],r'$\Delta f=%d$' %df] for i,df in enumerate(dfs)]
    plts01 =[[u,Xi(lam[1],df,Cs[1]),cs01[i],r'$\Delta f=%d$' %df] for i,df in enumerate(dfs)]
    laby = '$\chi(u)$'
    if 'S' in opts : laby = '$\sin\chi(u)$'
    dsp.stddisp(plts00,lw=2,labs=['$u(nm^{-1})$',laby],title='$\lambda=%dkeV, Cs=%dmm$' %(keV[0],Cs[0]*1e-6))
    dsp.stddisp(plts10,lw=2,labs=['$u(nm^{-1})$',laby],title='$\lambda=%dkeV, Cs=%dmm$' %(keV[1],Cs[0]*1e-6))
    dsp.stddisp(plts01,lw=2,labs=['$u(nm^{-1})$',laby],title='$\lambda=%dkeV, Cs=%dmm$' %(keV[1],Cs[1]*1e-6))

if '2' in opts:
    lam = cst.keV2lam(200)
    beta = np.linspace(5,50,100) #mrad
    beta0 = np.array([5,10,20,40])
    plts = [[beta,1.22*lam/beta*1e3 ,'b']]
    plts +=[[beta0,1.22*lam/beta0*1e3,'b*']]
    # plts += [beta,1.22*lam/beta*1e3,'b']
    dsp.stddisp(plts,labs=[r'$\beta(mrad)$',r'$d(\AA)$'],lw=2,ms=10,xylims=[1,51,0.5,6.5])
