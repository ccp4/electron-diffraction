from utils import*
import utils.FourierUtils as fu
from utils import displayStandards as dsp
import scipy.fftpack as fft

Tf,dt = 50, 0.01
T,nT,duty = 2,1,1

#square func
t = np.arange(0,Tf,dt) - Tf/2; N=t.size
#y = fu.squareBis(t,T,nT,duty)
y = 0*np.ones((N),dtype=complex) #(1+1J)*1e-5*np.ones((N),dtype=complex);
y[list(int(N/2)+10*np.array([-1,1]))]=np.array([0,1])*0.5/dt # + 1e-5*1J


# sine cardinal
# Y = fft.fftshift(fft.fft(fft.fftshift(y))
# f = fft.fftfreq(t.size,dt)
f,Y = fu.get_FFT(t,fft.fftshift(y),dt,pOpt=False)
df,fmax = fu.fftParams(dt,Tf,'')

# inverse fourier transform of the sine cardinal
t1,y1 = fu.get_iFFT(f,Y,df,pOpt=False)
dt1,T1 = fu.fftParams(df,2*fmax,'')


plt.close('all')
fig,(ax1,ax2,ax3) = plt.subplots(ncols=3,nrows=1)
plt.tight_layout(pad=2)
# fig,(ax1,ax2,ax3) = dsp.stddisp(rc=[3,1],pad=2.5)
ax1.plot(t,np.real(y),'b', t,np.imag(y),'r');
ax2.plot(f,np.real(Y),'b',f,abs(np.imag(Y)),'r',f,np.abs(Y)**2,'g');
ax3.plot(t1,np.real(y1),'b',t1,(np.imag(y1)),'r')

fonts = {'lab':15,'tick':10}
dsp.standardDisplay(ax1,['t','y' ],pOpt='tGX',fonts=fonts)#,xylims=[-1,1,-1.2,1.2])
dsp.standardDisplay(ax2,['f','F' ],pOpt='tGX',fonts=fonts)#,xylims=[-5,5,-2,2])
dsp.standardDisplay(ax3,['t','y1'],pOpt='tGX',fonts=fonts)#,xylims=[-1,1,-1.2,1.2])
fig.show()
