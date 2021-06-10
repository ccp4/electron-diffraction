from utils import*
plt.close('all')
path='../figures/'

t = np.linspace(1,1000,1000)
K = 1/0.025
Ug,Uh = 0.1,0.08
Sg,Sh = 0,0.02

Ugh = Ug-Uh
Ug_e = Ug - Uh*Ugh/(2*K*Sh)
So_e = - (abs(Uh)**2/(2*K*Sh))/(2*K)
Sg_e = (2*K*Sg - abs(Ugh)**2/(2*K*Sh))/(2*K)
Se = Sg_e-So_e
Ig = abs(Ug_e)**2/(K**2*Se**2+abs(Ug_e)**2)*np.sin(np.pi*t/K*np.sqrt(K**2*Se**2+abs(Ug_e)**2))**2

print(2*K*Sh,Ug)
print(K/Ug,K/Ug_e)
print(Sg_e)

plts = [[t,Ig,'b','']]
dsp.stddisp(plts,labs=['Thickness($A$)','$I_g$'],lw=2)
