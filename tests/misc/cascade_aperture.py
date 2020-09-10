from scipy.fftpack import ifftn
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import utils.displayStandards as dsp
#
# N = 30
# #f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
# f, (ax1, ax2) = plt.subplots(2, 1, sharex='col', sharey='row')
# plt.tight_layout(w_pad=0.75,h_pad=0.2)
# xf = np.zeros((N,N))
# i=3
# xf[0, i] = 1
# xf[0, N-i] = 1
# Z = ifftn(xf)
# ax1.imshow(xf, cmap=cm.Reds)
# ax2.imshow(np.real(Z), cmap=cm.binary)
#
#
# plt.show()

B = 0.1
x = np.linspace(0,10,1000)
x0s = [0,0.5,0.75]
n = np.arange(10)

T = np.exp(-B*x**2)
S = np.zeros(n.shape,dtype=complex)
Sf = np.zeros(x.shape,dtype=complex)
ST = np.zeros(x.shape,dtype=complex)
for x0 in x0s : S += np.exp(2J*np.pi*x0*n)
for Si,ni in zip(S,n) : Sf += Si*np.exp(-(x-ni)**2/0.01)
# for x0 in x0s : ST += T*np.exp(2J*np.pi*x0*x)
ST=T*Sf


N = np.stack([n,n])
S = np.stack([np.zeros(n.shape),np.abs(S)**2])
plts = [[x,T,'g--','$T$'],[N,S,'r--',''],[x,np.abs(Sf)**2,'m--','$S_f$'] ,[x,np.abs(ST)**2,'b','$ST$']]
dsp.stddisp(plts,labs=['$x$','$F$'],lw=2)
