import pickle,easygui,os
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
from utils import glob_colors as colors

class KeyBinder:
    def __init__(self,name=''):
        if not name:
            name=os.path.join(os.path.dirname(__file__),'data','viewer')
        d = pd.read_pickle(name+'_dict.pkl')
        self.df_keys = pd.read_pickle(name+'_df.pkl')

        self.df_generic = self.df_keys.loc[d['generic']]
        self.df_modes   = self.df_keys.loc[d['modes']]
        self.df_orient  = self.df_keys.loc[d['orient']+d['dorient']]
        self.df_frames  = self.df_keys.loc[d['frame'] +d['dframe']]
        self.df_thicks  = self.df_keys.loc[d['thick'] +d['dthick']]
        self.df_bright  = self.df_keys.loc[d['brightness']+d['dbrightness']+d['rbrightness']]
        self.df_rot     = self.df_keys.loc[d['rot']]
        self.df_vals    = self.df_keys.loc[d['vals']]
        # self.df_bloch   = self.df_keys.loc[d['bloch']]
        self.df_pets    = self.df_keys.loc[d['dsp']+d['opt']]
        self.df_multi   = self.df_keys.loc[d['multi']]
        # self.df_pets['alias'].iloc[:4] = Pnum_keys
        self.dict_k = self.df_keys.key.to_dict()
        self.dict_a = self.df_keys.alias.to_dict()

        for type_k,list_k in d.items():
            self.__dict__[type_k+'_keys'] = list(self.df_keys.loc[list_k].key.values)

        # print(self.dsp_keys)
        self.pets_keys = self.dsp_keys+self.opt_keys
        self.opt_chars = dict(zip(self.opt_keys,'hkrl'))
        self.dsp_chars = dict(zip(self.dsp_keys,'PMBKVSL'))

        #### Keys that trigger the show
        self.show_keys = self.modes_keys.copy()
        self.show_keys += self.pets_keys
        self.show_keys += self.rot_keys+self.thick_keys+self.brightness_keys+self.rbrightness_keys
        self.show_keys +=[self.dict_k[k] for k in ['settings','hold','solve_bloch']]

    def is_key(self,key):return self.df_keys.loc[(self.df_keys.key==key) | (self.df_keys.alias==key)]
    def is_cmd(self,cmd):return self.df_keys.iloc[[cmd in k for k in  self.df_keys.index]]

    def __call__(self,event):
        self.shortcut_call(event.key)

    def shortcut_call(self,key):
        combs = ['control','shift','alt',
         'ctrl+shift','ctrl+alt','ctrl+alt+shift',
         'shift+control','shift+super','shift+ctrl+super','shift+alt+control',
         'alt+control','alt+shift']
        if not key in combs :
            dfk = self.is_key(key)
            if any(dfk.index):print(colors.green+'found key'+colors.black);print(dfk)
            else:print(colors.red+'no key binding for %s' %key+colors.black)


    ###################################
    ##### Help
    ###################################
    def show_help_cli(self):
        print(colors.green+'Shortcuts : '+colors.black)
        print(colors.blue+'Generic :                    '+colors.black);print(self.df_generic)
        print(colors.blue+'Modes :                      '+colors.black);print(self.df_modes)
        print(colors.blue+'Frames (Frames mode) :       '+colors.black);print(self.df_frames)
        print(colors.blue+'Orientation (Rotate mode) :  '+colors.black);print(self.df_orient)
        print(colors.blue+'Thickness :                  '+colors.black);print(self.df_thicks)
        print(colors.blue+'Brightness :                 '+colors.black);print(self.df_bright)
        print(colors.blue+'Show values :                '+colors.black);print(self.df_vals)
        print(colors.blue+'Display:                     '+colors.black);print(self.df_pets)
        # print(colors.blue+'Display (Rotate mode): '+colors.black);print(self.df_bloch)

    def show_help_gui(self):
        # msg = str(self.df_keys)
        msg =''
        msg+= '\t\t\t SHORTCUTS : \n'
        msg+= 'Generic :     '+'\n'
        msg+=str(self.df_generic)+'\n'
        msg+= 'Modes :       '+'\n'
        msg+=str(self.df_modes)+'\n'
        msg+= 'Frames (Frames mode) :      '+'\n'
        msg+=str(self.df_frames)+'\n'
        msg+= 'Orientation (Rotate mode) : '+'\n'
        msg+=str(self.df_orient)+'\n'
        msg+= 'Thickness :   '+'\n'
        msg+=str(self.df_thicks)+'\n'
        msg+= 'Brightness :  '+'\n'
        msg+=str(self.df_bright)+'\n'
        msg+= 'Show values : '+'\n'
        msg+=str(self.df_vals)+'\n'
        # msg+= 'Display (Rotate mode): '+'\n'
        # msg+=str(self.df_bloch)+'\n'
        msg+= 'Display: '+'\n'
        msg+=str(self.df_pets)+'\n'
        easygui.msgbox(msg, title='help')

    def show_help_simple(self):
        print(self.df_keys)

################################################################################
#### keymapping
################################################################################
def save_keys(name):
    dict_k,dict_a = {},{}

    generic  = ['new_fig'   ,'settings' ,'help_cli','help_gui','help_simple' ]
    gen_keys = ['ctrl+enter','enter'    ,'ctrl+2'  ,'ctrl+1'  ,'ctrl+ '      ]
    generic += ['save_config','save_config_as','open_config','reload_config',]
    gen_keys+= ['ctrl+s'     ,'ctrl+S'        ,'ctrl+o'     ,'backspace'    ,]
    generic += ['save_fig','save_bloch','solve_bloch']
    gen_keys+= ['s'       ,'ctrl+alt+S','ctrl+B'     ]
    generic += ['hold','find_shortcut','find_command']
    gen_keys+= [' '   ,'ctrl+P'       ,'ctrl+p'      ]
    dict_k.update(dict(zip(generic,gen_keys)))

    multi      = ['open_multi','solve_multi','set_multi','run_multi','show_multi','show_tags' ,'log_info']
    multi_keys = ['ctrl+O'    ,'m'          ,'ctrl+m'   ,'ctrl+M'   ,'ctrl+alt+m','ctrl+alt+M','ctrl+l'  ]
    dict_k.update(dict(zip(multi,multi_keys)))

    #modes
    modes = ['frames_mode','rotate_mode']
    # modes_keys = ['ctrl+alt+'+c for c in '12']
    # modes_keys = ['ctrl+shift+'+c for c in '12']
    modes_keys = ['!','"']
    dict_k.update(dict(zip(modes,modes_keys)))

    #orientation
    orient = ['elev_up','elev_down','azim_up','azim_down']
    orient_keys  = ['up','down','right','left']
    dict_k.update(dict(zip(orient,orient_keys)))
    dorient = ['d'+s for s in orient]
    dorient_keys = ['ctrl+'+s for s in orient_keys]
    dict_k.update(dict(zip(dorient,dorient_keys)))

    #rotate in frame
    rots  = ['rot_up','rot_down','rotS_up','rotS_down']
    rot_keys  = ['r','R','ctrl+r','ctrl+R']
    dict_k.update(dict(zip(rots,rot_keys)))
    # drots = ['d'+s for s in rots]
    # dorient_keys = ['ctrl+'+s for s in orient_keys]
    # dict_k.update(dict(zip(drots,dorient_keys)))

    #frames
    frames = ['frame_'+s for s in ['up','down']] #'right','down','left']]
    dframes = ['inc_'+s for s in frames]
    frame_keys    = ['up','down']
    dframe_keys   = ['ctrl+'+s for s in frame_keys]
    frame_keys_a  = ['right','left']
    # dframe_keys_a = ['ctrl+'+s for s in frame_keys_a]
    dict_k.update(dict(zip(frames,frame_keys)))
    dict_k.update(dict(zip(dframes,dframe_keys)))
    dict_a.update(dict(zip(frames,frame_keys_a)))
    frame_keys = frame_keys.copy()+frame_keys_a

    #thickness
    thicks  = ['thick_up','thick_down']
    dthicks = ['d'+s for s in thicks]
    thick_keys  = ['shift+ctrl+up','shift+ctrl+down']
    dthick_keys = ['ctrl+t','ctrl+T']
    dict_k.update(dict(zip(thicks,thick_keys)))
    dict_k.update(dict(zip(dthicks,dthick_keys)))

    #brightness
    bright  = ['brightness '+s for s in ['down','up']]
    dbright = ['d'+s for s in bright]
    rbright = [s+'brightness reset' for s in ['','d']]
    brightness_keys  = ['pagedown','pageup']
    dbrightness_keys = ['ctrl+'+k for k in brightness_keys]
    rbrightness_keys = ['shift+ctrl+'+k for k in brightness_keys]
    dict_k.update(dict(zip(bright,brightness_keys)))
    dict_k.update(dict(zip(dbright,dbrightness_keys)))
    dict_k.update(dict(zip(rbright,rbrightness_keys)))

    #vals
    vals = ['hkl','I','z','u','v','beams_vs_thickness','Ig','Vg','Sw']
    vals = ['show_'+s for s in vals]
    vals_keys = ['ctrl+'+c for c in 'hizuvbIVE']
    dict_k.update(dict(zip(vals,vals_keys)))

    # #display Bloch
    # Fs    = ['L','Sw','Vg','S','I','Ig']#,'E']
    # fopt  = ['m','L' ,'l' ,'m','m','m' ]#,'m']
    # bloch = ['Lattice','Excitation error','Potential','Scattering','Intensities','Kinematic']
    # bloch = ['show '+s for s in bloch]
    # F_keys    = [''+c for c in 'LGVSIK']
    # Fnum_keys = [''+c for c in '123456']
    # df_F = pd.DataFrame.from_dict(dict(zip(['names','key','alias','F','fopt'],[bloch,F_keys,Fnum_keys,Fs,fopt])))
    # df_F = df_F.set_index('names')

    #display Frames
    dict_k['show_Exp']='0'
    petsD = ['Processed','Multislice','Bloch','Kinematic','Potential','Excitation error','Lattice']
    petsD = ['show_'+s for s in petsD]
    dsp_keys   = ['P','M','B','K','V','S','L']
    dsp_keys_a = [''+c for c in '1234567']
    dsp_chars = dict(zip(dsp_keys,'PMBKVSL'))
    dict_k.update(dict(zip(petsD,dsp_keys)))
    dict_a.update(dict(zip(petsD,dsp_keys_a)))
    #options
    petsO = ['hkl_exp','hkl_kin','rings','legend']
    petsO = ['show_'+s for s in petsO]
    opt_keys = ['j','k','g','l']
    opt_chars = dict(zip(opt_keys,'hkrl'))
    dict_k.update(dict(zip(petsO,opt_keys)))


    #### Data frame for the help command
    df_keys = pd.DataFrame.from_dict(dict_k,orient='index',columns=['key'])
    # df_keys = df_keys.append(df_F[['key','alias']])

    df_keys['alias'] = ['']*df_keys.shape[0]
    for k,v in dict_a.items():df_keys.loc[k,'alias'] = v
    df_keys['mode']  = ['any']*df_keys.shape[0]
    df_keys.loc[orient+dorient,'mode'] = 'rotate'
    df_keys.loc[frames+dframes,'mode'] = 'frames'

    d = dict(zip(['generic','modes','frame','dframe','orient','dorient','rot','thick','dthick','brightness','dbrightness','rbrightness','vals','dsp'  ,'opt'  ,'multi'],
                 [ generic , modes , frames, dframes, orient , dorient , rots, thicks, dthicks, bright     , dbright     , rbright     , vals , petsD , petsO , multi ]
                 ))

    df_keys.to_pickle(name+'_df.pkl')
    with open(name+'_dict.pkl','wb') as f:pickle.dump(d,f,pickle.HIGHEST_PROTOCOL)
    print(colors.yellow+name+colors.green+' saved'+colors.black)
    # return df_keys,d

if __name__=='__main__':
    save_keys(name='data/viewer')
