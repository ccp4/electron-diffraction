import importlib as imp
import pandas as pd,numpy as np,easygui
from utils import displayStandards as dsp #;imp.reload(dsp)


def show_frame(opts='Pqr',mag=10,rot=0,hkl_idx=[],
    df_pets=None,im_multi=None,df_bloch=None,
    **kwargs):
    ''' Show a frame with information specified by opts
    - opts : P(proc), M(Multislice), B(Bloch), K(kin), V(Vg), S(Sg), L(lat)
        h(hkl_exp),r(rings), k(hkl_kin) l(legend)
    '''

    labs,qmax,wm = ['$q_x$','$q_y$'],0,0.05
    plts,scat,txts,im,legElt = [[0,0,'b+','']],(),[],None,None

    if 'P' in opts and isinstance(df_pets,pd.DataFrame):
        qx_p,qy_p,I = df_pets[['qx','qy','I']].values.T
        if qx_p.size:break
        scat += ([qx_p,qy_p,I,cp,mp],)
        qmax = np.ceil(max(abs(qx_p).max(),abs(qy_p).max()))

        if 'h' in opts:
            hx,kx,lx = df_pets[['hx','kx','lx']].values.T
            sx = lambda x:['%d' %round(x),'%.1f' %x][abs(x-round(x))>0.06]
            txts += [[x+0*wm,y+wm,'(%s,%s,%s)' %(sx(h0),sx(k0),sx(l0)),(0.5,)*3] for x,y,h0,k0,l0 in zip(qx_p,qy_p,hx,kx,lx)]

    if isinstance(df_bloch,pd.DataFrame):
        qx_b,qy_b = df_bloch[['px','py']].values.T
        if rot:
            ct,st = np.cos(np.deg2rad(rot)),np.sin(np.deg2rad(rot))
            qx_b,qy_b = ct*qx_b-st*qy_b,st*qx_b+ct*qy_b
        for c,C in dfb.iterrows() :
            if c in opts:
                I = C.fz(df_bloch[C.F])
                scat += ( [qx_b,qy_b,I/I.max()*mag,C.color,C.marker],)
        bs = [c for c in dfb.index if c in opts]
        if 'k' in opts and bs and any(hkl_idx):
            h,k,l,qx0,qy0 = df_bloch.iloc[hkl_idx][['h','k','l','qx','qy']].values.T
            if rot:qx0,qy0 = ct*qx0-st*qy0,st*qx0+ct*qy0
            ctxt = dfb.loc[bs[0]].color
            txts += [[x,y+wm,'(%d,%d,%d)' %(h0,k0,l0),ctxt] for x,y,h0,k0,l0 in zip(qx0,qy0,h,k,l)] #,I_b) if I>Itol]
        qmax = np.ceil(max(qmax,qx_b.max(),qy_b.max()))

    if 'M' in opts and im_multi:
        qx_m,qy_m,I=im_multi
        qmax = np.ceil(max(qmax,qx_m.max(),qy_m.max()))
        im=[qx_m,-qy_m,I]

    if 'r' in opts:
        t = np.linspace(0,np.pi*2,100)
        ct,st = np.cos(t),np.sin(t)
        plts+=[[i*ct,i*st,'m','',0.5] for i in np.arange(0.25,qmax,0.25)]
        plts+=[[i*ct,i*st,'m','',2] for i in np.arange(1,qmax)]


    if 'l' in opts:
        legElt = {k:v for c,(k,v) in zip('PMBKVSL',legElts.items()) if c in opts}

    if not 'fonts' in kwargs.keys():kwargs['fonts'] = {'text':15}
    if not 'xylims' in kwargs.keys():kwargs['xylims'] = qmax
    dsp.stddisp(plts,ms=20,scat=scat,im=im,texts=txts,bgcol=None,gridOn=0,
        labs=labs,sargs={'alpha':0.75},legElt=legElt,**kwargs)

########################################################################
#### misc
########################################################################
def get_fz(opts):
    '''mapping for function to apply
    - opts : 'r'(real) 'i'(imag) 'a'(angle) 'm'(mag) 'l'(log10(|Fhkl|+1)) '2'(mag^2) 'L'(logM)
    '''
    if not opts:opts='m'
    keys   =['r','i','a','m','l','2','L']
    fs     = [np.real,np.imag,np.angle,np.abs,logF,abs2,logM]
    fs_str = ['real part','imag part','phase','magnitude','$\log_{10}(|F|+1)$','$|F|^2$','$-\log_{10}$']
    fz     = dict(zip(keys,fs))[opts]
    fz_str = dict(zip(keys,fs_str))[opts]
    return fz,fz_str

abs2 = lambda F:np.abs(F)**2
logF = lambda F:np.log10(np.abs(F)+1)
logM = lambda F:-np.log10(np.maximum(F,1e-5))

dfb = { 'B':['I' ,get_fz('m')[0],'b','o','Intensity Bloch$I_{g}$'],
        'K':['Ig',get_fz('m')[0],'g','o','Intensity kinematic $I_{g,kin}$'],
        'V':['Vg',get_fz('l')[0],'c','d','Potential $V_{g}$'],
        'S':['Sw',get_fz('L')[0],'m','o','Excitation error $\zeta_{g}$'],
        'L':['L' ,get_fz('m')[0],'k','s','Lattice']}
dfb = pd.DataFrame.from_dict(dfb,orient='index',columns=['F','fz','color','marker','legend'])
cp,mp = (0.5,)*3,'o'
legElts={'Processed':[cp,mp],'Multislice':'r*'}
legElts.update({C.legend:[C.color,C.marker] for c,C in dfb.iterrows()})


def multenterbox(msg,title,fieldValues,fieldNames):
    fieldValues = easygui.multenterbox(msg, title, fieldNames,fieldValues)
    while True:
        if fieldValues is None:
            break
        errs = list()
        for n, v in zip(fieldNames, fieldValues):
            if v.strip() == "":errs.append('"{}" is a required field.'.format(n))
        if not len(errs):
            break
        fieldValues = easygui.multenterbox("\n".join(errs), title, fieldNames, fieldValues)
    if fieldValues:
        return dict(zip(fieldNames,fieldValues))
