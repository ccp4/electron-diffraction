from utils import*
import scipy.signal as sig
plt.close('all')


dfI=pd.read_csv('dat/ex5_2A_I.csv') #Silicon 2-beam intensities
dfS=pd.read_csv('dat/ex5_2A_Si_200kVdat.csv') #Silicon Simu
h,k,l,ext = ' h ','k ','l ','Extinc. / nm'
cols=[h,k,l,ext]
ds = 'Distance (nm)'
bf = 'Bright-field intensity (norm)'
df = 'Dark-field intensity (norm)'

t,Ib,Id = dfI[ds],dfI[bf],dfI[df]
plts = [[t,Ib,'r','$I_b$'],[t,Id,'b','$I_d$']]


idx = sig.find_peaks(-Ib,width=10)[0]
# idx = sig.find_peaks(-Id,width=10)[0]
plts += [[t[idx],Ib[idx],'gs']]
# plts += [[t[idx],Id[idx],'gs']]
dsp.stddisp(plts,labs=['$t(nm)$','$I$'],lw=2)

xi = np.diff(t[idx]);print(xi)

print(dfS[cols])
