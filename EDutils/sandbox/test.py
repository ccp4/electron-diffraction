# from EDutils import __version__
from . import __version__

print(__version__)

##too slow
class Shortcut_finder:
    def __init__(self,df_keys):
        self.df_keys=df_keys
        self.fig,self.ax = dsp.stddisp(figsize=(0.3,0.2),pOpt='')#,texts=['type shortcut'])
        cid = self.fig.canvas.mpl_connect('key_press_event', self)
    def __call__(self,event):
        dfk = self.is_key(event.key)#;print(dfk)
        self.ax.cla()
        self.ax.axis('off')
        self.ax.text(0.5,0.5,'Shortcut : %s\n' %event.key     ,color='k',fontsize=30,horizontalalignment='center')
        if any(dfk.index) :
            txt = dfk.index[0].replace('_',' ')
            self.ax.text(0.5,0.25,txt,color='g',fontsize=30,horizontalalignment='center')
        self.fig.canvas.draw()
    def is_key(self,key):return self.df_keys.loc[self.df_keys.key==key]
    def is_cmd(self,cmd):return self.df_keys.iloc[[cmd in k for k in  self.df_keys.index]]
