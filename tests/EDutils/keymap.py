from utils import*                  #;imp.reload(dsp)
from EDutils import viewer_config as vw_cfg   ;imp.reload(vw_cfg)
plt.close('all')

key_map = vw_cfg.KeyBinder()
txt = 'type a key to see the mapping on the console'
fig,ax = dsp.stddisp(pOpt='',texts=[0.5,0.5,txt,'g'],figsize='21')#,[0.5,0.25],)
cid = fig.canvas.mpl_connect('key_press_event', key_map)
fig.show()
