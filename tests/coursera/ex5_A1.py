from utils import*
plt.close('all')


df=pd.read_csv('dat/ex5_A1.csv') #Aluminum
ext = 'Extinc. / nm'

hs = np.arange(2,11,2)[:1]
s = np.linspace(-0.1,0.1,1000) # nm^-1
t = 130 #nm

xi=np.zeros(hs.shape)
for i,h in enumerate(hs):
    f = df[(df[' h ']==h) & (df['k ']==0) & (df['l ']==0)]
    xi_g =f[ext].values[0];#print(xi_g)
    xi[i] = xi_g
    Ig0 = 1/(1+xi_g**2*s**2)
    Ig = np.sin(np.pi*t*np.sqrt(1/xi_g**2+s**2))**2*Ig0
    plts = [[s,Ig,'k'],[s,Ig0,'r--']]
    dsp.stddisp(plts,xylims=['y',0,1.1*Ig.max()],title='(%d 0 0)' %h,labs=['$s(nm^{-1})$','$I_g$'],lw=2)
