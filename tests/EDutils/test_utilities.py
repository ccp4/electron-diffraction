from utils import*
from EDutils import utilities as ut;imp.reload(ut)
plt.close('all')

def test_sweep():
    class C:
        def __init__(self,name,path,a,b=1,c=''):
            if not isinstance(c,str):
                raise TypeError('c must be of type str')
            self.a=a
            self.b=b
            self.c=c
            self.name=name
            self.path=path
            self.save()

        def _get_pkl(self):
            return ut.get_pkl(path=self.path,name=self.name)
        def save(self):
            ut.save(self,file=self.get_pkl())
        def __str__(self):
            return 'a=%d, b=%d, c=%s' %(self.a,self.b,self.c)

    df1 = ut.sweep_var(C,path='dat/sweep',params='a',vals=np.arange(2),tag='single')
    df2 = ut.sweep_var(C,path='dat/sweep',params=['a','b'],vals=np.random.rand(3,2),tag='double')
    df3 = ut.sweep_var(C,path='dat/sweep',params=['a','c'],vals=[[1,'test1'],[2,'test2']],tag='mixed')

    print('df1[0 ] : ',ut.load_pkl(df1.iloc[0].pkl))
    print('df2[1 ] : ',ut.load_pkl(df2.iloc[1].pkl))
    print('df3[-1] : ',ut.load_pkl(df3.iloc[-1].pkl))

def test_uvw():
    u0 = [1,0,1]
    args = {'show':True,'h3d':False,'xylims':1}
    ut.get_uvw_cont(u0,u1=[1,1,1] ,nframes=100   ,title='cont'      ,**args)
    ut.get_uvw_rock(e0=u0,e1=[1,0],npts=20,deg=1 ,title='rock_theta',**args)
    ut.get_uvw_rock(e0=u0,e1=[0,1],npts=20,deg=1 ,title='rock_phi'  ,**args)
    ut.get_uvw_CBED(u0=[1,1,1],cart=False ,npts=[10,20],deg=20  ,title='CBED'      ,**args)


# test_sweep()
# test_uvw()
