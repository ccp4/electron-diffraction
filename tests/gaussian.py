import numpy as np
import scipy.fftpack as fft
import scipy.signal.windows as win
import matplotlib.pyplot as plt
plt.close('all')

N = 10000
k = 100
t = np.linspace(-1,1,N)
y = np.exp(-1J*k*t**2)
y = win.tukey(N)
fig1,ax=plt.subplots()
plt.plot(t,y.real,'b')
plt.plot(t,y.imag,'r')
fig1.show()

dt = t[1]-t[0]
f = fft.fftfreq(N,dt)
F = fft.fft(y)*dt
f = fft.fftshift(f)
F = fft.fftshift(F)
Fth = np.sqrt(np.pi/k)*np.exp(1J*(np.pi**2*f**2/k-np.pi/4))

fig2,ax=plt.subplots()
plt.plot(F.real)
# plt.plot(f,F.real,'b')#,f,Fth.real,'b--')
# plt.plot(f,F.imag,'r')#,f,Fth.imag,'r--')
fig2.show()
