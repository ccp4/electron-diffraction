from utils import*

def plot_Ikin_Idyn(opt='p'):
    txi = [0.1,0.5]
    t   = [100,500]                   #A
    xi  = np.array(t)/np.array(txi) #A
    K   = 100                       #A^-1
    Ug  = K/np.array(xi)            #A^-2
    print('xi(A):',xi);print('Ug(A^-2):',Ug)
    w = np.linspace(-10,10,1000)
    Ikin = lambda t,xi:np.sinc(t/xi*w)**2
    Idyn = lambda t,xi:np.sinc(t/xi*np.sqrt(1+w**2))**2
    plts = [[w,Ikin(t[0],xi[0]),'b--',''],[w,Idyn(t[0],xi[0]),'b',''],
            [w,Ikin(t[1],xi[1]),'r--',''],[w,Idyn(t[1],xi[1]),'r','']]
    legElt = {'kin':'k--','dyn':'k',
             't=%dA' %int(t[0]):'b', 't=%dA' %int(t[1]):'r'}
    dsp.stddisp(plts,labs=[r'$w (A^{-1})$' ,'$I(w)/\sigma v_g t$'],
        lw=2,opt=opt,figsize='12',legElt=legElt,name=figpath+'kin_dyn.svg')

def plot_Ikin_Idyn_S0(opt='p'):
    K  = 100        #A^-1
    Ug = 0.1        #A^-2
    tlam = K/Ug
    t    = [10,100,500,1000,2000]
    txi  = np.array(t)/tlam; #
    print('xi   = %.2f A' %(tlam))
    print('t(A) =',t)
    print('txi  =',txi)
    if 'p' in opt:
        t = np.linspace(0,1000,1000) #A
        t0 = np.array([100,500])
        Ikin = lambda t:(pi*t/tlam)**2
        Idyn = lambda t:np.sin(pi*t/tlam)**2
        plts = [[t,Ikin(t),'k--','kinematic'],[t,Idyn(t),'k','dynamic']]
        plts+= [[t0[0],Ikin(t0[0]),'bs',''],[t0[0],Idyn(t0[0]),'bs','']]
        plts+= [[t0[1],Ikin(t0[1]),'rs',''],[t0[1],Idyn(t0[1]),'rs','']]
        dsp.stddisp(plts,labs=['$Thickness$ (A)' ,'$I(S_G=0)$'],
            lw=2,opt=opt,figsize='12',name=figpath+'kin_dyn0.svg')

####### display
def plot_Idyn():
    tzeta = [0.1,0.25,0.5,0.75,1,1.5]#,2,2.5,3];
    tzeta = [0.75,1.5]
    nTs=len(tzeta)
    cs = getCs('Reds',nTs)
    w = np.linspace(-6,6,1000)
    I = lambda tz:np.sin(pi*tz*np.sqrt(1+w**2))**2/(1+w**2)
    plts=[[w,I(tzeta[i]),cs[i],'%.2f' %(tzeta[i])] for i in range(nTs)]
    dsp.stddisp(plts,labs=[r'$\omega$','$Intensity$'],title=r'$t/\xi$',lw=2)#,xylims=[-6,6,0,1])

#plot_Idyn()
plot_Ikin_Idyn_S0(opt='s')
plot_Ikin_Idyn(opt='s')
