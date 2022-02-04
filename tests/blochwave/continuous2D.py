from utils import*
from blochwave import bloch2D as bl          ;imp.reload(bl)
from blochwave import bloch as bl0           ;imp.reload(bl0)# load_Bloch
from wallpp import config as cg      ;imp.reload(cg)
plt.close('all')

class bloch2D_Cont:
    def __init__(self,u,tag='test',path='dat',run=False,bargs={}):
        self.u=u
        self.path=path
        self.tag=tag
        if run:self._run(bargs)
        self._set_df_G()

    def _run(self,bargs):
        for i,u0 in enumerate(self.u):
            name = '%s_%s' %(self.tag,str(i).zfill(3))
            b0  = bl.Bloch2D(file=self.tag,u=u0,path=self.path,name=name,**bargs)
            hk = [str((h,k)) for h,k in zip(*b0.df_G[['h','k']].values.T)]
            b0.df_G.index = hk
            b0.save()

    def _set_df_G(self):
        b = self.load_b()

        hk = np.unique(np.hstack([b0.df_G.index.values for b0 in b]))
        df_G=pd.DataFrame(np.zeros((hk.size,3)),columns=['qx','qy','I'],index=hk)
        for i,b0 in enumerate(b):
            df_G.loc[b0.df_G.index,['qx','qy']] = b0.df_G[['qx','qy']].values
            df_G.loc[b0.df_G.index,'I'] += b0.df_G.I.values
        df_G.loc[str((0,0)),'I'] = 0
        self.df_G=df_G


    #####################################################################
    #### available
    #####################################################################
    def load_b(self):
        return [bl0.load_Bloch(self.path,
            '%s_%s' %(self.tag,str(i).zfill(3)),v=0) for i,u0 in enumerate(self.u)]

    def show_I(self,thick=None,fz=abs,ms=80,pOpt='im',caxis=[],
        **kwargs):
        if thick:self.set_I(thick)
        if isinstance(fz,str):
            if fz=='log':fz = lambda x:np.log10(np.maximum(x,1e-10))
            else:fz=abs
        qx,qy,I = self.df_G[['qx','qy','I']].values.T
        if not len(caxis)==2:caxis = [fz(I).min(),fz(I).max()]
        dsp.stddisp(scat=[qx,qy,ms,fz(I)],labs=['$q_x$','$q_y$'],
            title='$z=%d\AA$' %self.thick,caxis=caxis,
            pOpt=pOpt,**kwargs)

    def set_I(self,thick):
        self.thick=thick
        b = self.load_b()
        self.df_G.I=0
        for b0 in b:
            b0.set_thickness(thick)
            self.df_G.loc[b0.df_G.index,'I'] += b0.df_G.I.values
        self.df_G.loc[str((0,0)),'I'] = 0

t = np.linspace(0,np.pi,91)
u = np.array([np.cos(t),np.sin(t)]).T
pattern = np.array([[0.2,0.1,1],[0.1,0.3,2]])
pp = {'a':5,'b':5,'alpha':90,'pattern':pattern,'interp':False,'fract':True}
bargs = {'crys':pp,'keV':200,'thick':500,'Nmax':20,'Smax':0.025,'solve':1}

for pp_type in cg.pp_types[2:]:
    pp['pp_type']=pp_type
    bargs['crys']=pp
    b2d = bloch2D_Cont(u,tag=pp_type,path='dat',run=True,bargs=bargs)

    thicks = np.arange(5,101,5)
    # thicks = np.hstack([np.arange(5,26,5),np.arange(30,101,10), np.arange(120,201,20),np.arange(250,501,50)])
    for thick in thicks:#[:2]:
        b2d.show_I(thick=thick,fz='log',xylims=4,
            name='figures/cont2D_%s_z%d.png' %(pp_type,thick),opt='cs')


# fig,ax = dsp.stddisp(opt='')
# for b0 in b: b0.show_ewald(ax=ax,legOpt=0,opt='')
# fig.show()
