from utils import*


tab_path='/home/ronan/Downloads/'
elts = ['H','C','N','O','S','P']
cH,cC=0.8,0.6
cs=[(cH,cH,cH),(cC,cC,cC),(0,0,1),(1,0,0),(0,1,1),(1,0,1)]


df=pd.read_html(tab_path+'Atomic_form_factors.html')[-1]
cols=df.iloc[0];df=df[1:];df.columns=cols
df=df.set_index('Element')
df=df.loc[elts].astype(float)
df['color']=cs

q = np.linspace(0,30,100)
q2 = ((q/(4*pi))**2)[:,None] #A

Ai = df[['a%d' %i for i in range(1,5)]].values
Bi = df[['b%d' %i for i in range(1,5)]].values
C  = df['c'].values


df['fq'] = [(ai*np.exp(-bi*q2)).sum(axis=1) + c  for ai,bi,c in zip(Ai,Bi,C)]

plts=[[q,elt.fq,elt.color,'$%s$' %elt.name ] for index,elt in df.iterrows()]
stddisp(plts,labs=['$q(A)$','$f_j(e)$'],legLoc='upper right')