import importlib as imp
import tifffile,os,glob,pickle5,subprocess
import numpy as np,pandas as pd
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from utils import displayStandards as dsp   #;imp.reload(dsp)
from utils import glob_colors as colors     #;imp.reload(colors)
from utils import handler3D as h3D          #;imp.reload(h3D)
from . import utilities as ut               ;imp.reload(ut)

class Rocking:
    def __init__(self,Simu,uvw,tag,path,Sargs):
        """ simulate rocking curve

        Parameters
        ----------
        Simu
            Simulator object
        uvw
            orientation vectors
        tag
            a tag identifying the sweep
        Sargs
            Simulator constructor arguments
        """
        self.path = path
        self.tag  = tag
        self.uvw  = uvw
        self.Sargs = Sargs
        self.df = ut.sweep_var(Simu,params='u',vals=uvw,tag=tag,path=path,**Sargs)
        ts = np.arange(uvw.shape[0])
        self.ts       = ts
        self.df['ts'] = ts
        self.Iz_dyn = {}
        self.Iz_kin = {}
        self.save(v=1)

    ###########################################################################
    #### compute
    ###########################################################################
    def _integrate_rocking(self,refl=[],new=0):
        if new:self.Iz_dyn,self.Iz_kin = {},{}

        # refl,nbs = self.get_beams(cond=cond,refl=refl)
        nbs=len(refl)
        z,nzs = self._get_z()
        hkl  = [h for h in refl if not h in self.Iz_dyn.keys()] #;print(hkl)
        refl = [eval(h) for h in hkl]

        nbs,nts = len(hkl),self.ts.size
        if nbs:
            Iz_dyn  = dict(zip(hkl, np.zeros((nbs,nzs)) ))
            Iz_kin  = dict(zip(hkl, np.zeros((nbs,nzs)) ))
            for i in range(nts):
                sim_obj = self.load(i)
                idx = sim_obj.get_beam(refl=refl,cond='',index=True)
                if idx:
                    hkl0 = sim_obj.df_G.iloc[idx].index #[str(tuple(h)) for h in sim_obj.get_hkl()[idx]]
                    for idB,hkl_0 in zip(idx,hkl0):
                        Iz_dyn[hkl_0] += sim_obj.Iz[idB,:]
                        Iz_kin[hkl_0] += sim_obj.Iz_kin[idB,:]

            self.Iz_dyn.update(Iz_dyn)
            self.Iz_kin.update(Iz_kin)
            self.save()
            print(colors.green+'rock.Iz updated'+colors.black)

    def get_rocking(self,iZs:Optional[slice]=-1,
        zs:Optional[Iterable[float]]=None,
        refl:Sequence[tuple]=[],cond:str=''):
        """Get intensities at required beams and thicknesses

        Parameters
        ----------
        zs
            selected thicknesses
        iZs
            slice into thicknesses
        refl
            reflections to get
        cond
            condition to apply to selecte reflections

        Returns
        -------
            z,dict
                z,{beam:I(z)}
        """
        iZs,nzs  = self._get_iZs(iZs,zs)        #;print(iZs)
        z0 = self.load(0).z.copy()[iZs]
        print('setting thickness to %dA' %z0[-1])
        self.do('set_thickness',thick=z0,v=0)

        refl,nbs = self.get_beams(cond=cond,refl=refl)  #;print(refl)
        nts = self.ts.size
        I = {}
        for h in refl : I[str(h)]=np.nan*np.ones((nts,nzs))

        #gather the intensities
        for i in range(nts):
            # print(colors.red,i,colors.black)
            sim_obj = self.load(i)
            hkl0  = sim_obj.get_beam(refl=refl,index=False)
            if hkl0:
                idx  = sim_obj.get_beam(refl=refl,index=True)
                # hkl0 = [str(tuple(h)) for h in sim_obj.get_hkl()[idx]]
                for idB,hkl_0 in zip(idx,hkl0):
                    I[hkl_0][i,:] = np.array(sim_obj.Iz[idB,iZs])
        z = self.load(0).z.copy()[iZs]
        return z,I


    ###########################################################################
    #### Display
    ###########################################################################
    # def QQplot(self,zs=None,iZs=10,refl=[],cond='',
    #     int_opt=True,cmap='Spectral',**kwargs):
    #     if int_opt:self.integrate_rocking(refl,cond)
    #     iZs,nzs  = self._get_iZs(iZs,zs)    #;print(iZs)
    #     z  = self.load(0).z.copy()[iZs]
    #
    #     refl = self.get_beams(refl=refl,cond=cond)
    #     refl = [str(tuple(h)) for h in refl]
    #     Iz_dyn = np.array([self.Iz_dyn[h].copy()[iZs] for h in refl])
    #     Iz_kin = np.array([self.Iz_kin[h].copy()[iZs] for h in refl])
    #     # print(Iz_dyn.shape)
    #     iB = np.argsort(np.sum(Iz_dyn,axis=1))[-1]
    #     Iz_dyn/= Iz_dyn[iB,:]
    #     Iz_kin/= Iz_kin[iB,:]
    #     cs = dsp.getCs(cmap,nzs) #; print(len(cs),Iz_dyn.shape,Iz_kin.shape)
    #
    #     plts=[[I_kin,I_dyn,[cs[i],'o'],r'$z=%d \AA$' %z0] for i,(z0,I_dyn,I_kin) in enumerate(zip(z,Iz_dyn.T,Iz_kin.T))]
    #     plts+=[ [[0,1],[0,1],[(0.5,)*3,'--'],''] ]
    #     dsp.stddisp(plts,labs=['$I_{kin}$','$I_{dyn}$'],sargs={'alpha':0.5},**kwargs)

    def plot_integrated(self,refl,new:bool=False,cm='Spectral',**kwargs):
        """plot the integrated rocking curves for selected beams as function of z

        Parameters
        ----------
        refl
            beams to consider
        new
            update integration
        """
        self._integrate_rocking(refl=refl,new=new)
        nbs = len(refl)
        z = self.load(0).z

        cs = dsp.getCs(cm,nbs)
        plts = [[z,self.Iz_dyn[h],cs[i],'%s' %h] for i,h in enumerate(refl)]
        dsp.stddisp(plts,labs=[r'$z(\AA)$','$I_{int}$'],**kwargs)

    def plot_rocking(self,cmap='viridis',x:str='Sw',
        cond='',refl=[],iZs=-1,zs=None,
        **kwargs):
        """plot rocking curve for set of selected beams at thickness zs

        Parameters
        ----------
        cond,refl,iZs,zs
            select beams and thicknesses (see :meth:`~Rocking.get_rocking`)
        x
            to display on x axis
        kwargs
            args to pass to dsp.stddisp
        """
        z,I = self.get_rocking(cond=cond,refl=refl,zs=zs,iZs=iZs)
        refl,plts = list(I.keys()),[]           #;print(refl)
        xlab = {'frame':'frame','Sw':r'$S_w(\AA^{-1})$','theta':r'$\theta(deg)$'}[x]

        nbs,nzs = len(refl),z.size
        if nbs>nzs:
            cs,ms = dsp.getCs(cmap,nbs), dsp.markers
            legElt = { '%s' %refl0:[cs[i],'-'] for i,refl0 in enumerate(refl)}
            for iz,zi in enumerate(z):
                legElt.update({'$z=%d A$' %(zi):['k',ms[iz]+'-']})
                for i,refl0 in enumerate(refl):
                    df=self.get_frames(refl0,cols=['Sw'])
                    # plts += [[df.Sw,I[refl0][df.frame,iz],[cs[i],ms[iz]+'-'],'']]
                    plts += [[df[x],I[refl0][df.frame,iz],[cs[i],ms[iz]+'-'],'']]
        else:
            # rocking for different thicknesses
            cs,ms = dsp.getCs(cmap,nzs),  dsp.markers
            legElt = { '%s' %refl0:['k','-'+ms[i]] for i,refl0 in enumerate(refl)}
            # self.get_frames(hkl,iTs=slice(0,None))
            for i,refl0 in enumerate(refl):
                df=self.get_frames(refl0,cols=['Sw'])
                for iz,zi in enumerate(z):
                    # plts += [[ts,I[refl0][:,iz],[cs[iz],ms[i]+'-'],''] for i,refl0 in enumerate(refl)]
                    plts += [[df[x],I[refl0][df.frame,iz],[cs[iz],ms[i]+'-'],'']]
            legElt.update({'$z=%d A$' %(zi):[cs[iz],'-'] for iz,zi in enumerate(z) })
        return dsp.stddisp(plts,labs=[xlab,'$I$'],legElt=legElt,**kwargs)


    def Sw_vs_theta(self,refl=[[0,0,0]],cond='',thick=None,fz=abs,opts='',
        iTs=slice(0,None),ts=None,
        cm='Spectral',figname='',**kwargs):
        """Displays Sw and I for a range of angle simulations at given thickness

        Parameters
        ----------
        refl,cond
            selected reflections
        thick
            thickness
        fz
            functor for Sw
        Iopt
            plot I
        """
        Iopt = 'I' in opts
        iTs,nts = self._get_iTs(iTs,ts)
        xlab,ts = r'$\theta$',self.ts.copy()[iTs]
        if 'f' in opts:
            xlab,ts = 'frame',np.arange(1,self.ts.size+1)[iTs]       #;print(iTs,ts)

        if Iopt:
            if thick:
                if not self.load(0).thick==thick:
                    # self.do('_set_excitation_errors',Smax=0.025)
                    # self.do('solve',thick=thick,Smax=0.025)
                    self.do('set_thickness',thick=thick)
            else:
                thick = self.load(0).thick
        refl,nbs = self.get_beams(cond,refl)                        #;print(refl)

        Sw = pd.DataFrame(np.ones((nts,nbs)),columns=[str(h) for h in refl])
        if Iopt:I  = pd.DataFrame(np.zeros((nts,nbs)),columns=[str(h) for h in refl])
        for i,name in enumerate(self.df.index[iTs]):
            b = self.load(i) #;print(i)
            hkl0 = b.get_beam(refl=refl,cond=cond,opt=0)
            # hkl0 = [str(tuple(h)) for h in b.get_hkl()[idx]]
            Sw.loc[i,hkl0] = b.df_G.loc[hkl0,'Sw'].values
            if Iopt:I.loc[i,hkl0] = b.df_G.loc[hkl0,'I'].values

        #locate minimum excitation errors
        iSmin = np.argmin(abs(Sw).values.T,axis=1) #locate minimums
        Sw[Sw==1]=np.nan
        SwE = fz(Sw.values.T)

        if 'i' in opts:
            dsp.stddisp(im=[SwE],pOpt='im',labs=[xlab,'$beam$'],caxis=[fz(1e-2),fz(1e-6)],
                cmap='Reds',title='Excitation error',name=figname+'_Sw.svg',**kwargs)
            # print(refl)

        else:
            cs,txts = dsp.getCs(cm,nbs),[]
            plts = [[ts,Sw0,[cs[i],'-o'],''] for i,Sw0 in enumerate(SwE)]
            if 't' in opts:txts = [[ts[idx],SwE[i,idx],'%s' %str(h),cs[i]] for i,(h,idx) in enumerate(zip(refl,iSmin))]
            dsp.stddisp(plts,texts=txts,labs=[xlab,'$S_w$'],name=figname+'_Sw.svg',**kwargs)

        if Iopt:
            IE = I.values.T
            if 'i' in opts:
                dsp.stddisp(im=[IE],pOpt='im',labs=[xlab,'$beam$'],caxis=[0,0.2],
                    cmap='YlGnBu',title='Intensity',name=figname+'_Sw.svg',**kwargs)
            else:
                if 'w' in opts:
                    xlab='$S_w$'
                    plts = [[Sw0,I0,[cs[i],'-o'],''] for i,(Sw0,I0) in enumerate(zip(SwE,IE))]
                    txts = [[0,IE[i,idx],'%s' %str(h),cs[i]] for i,(h,idx) in enumerate(zip(refl,iSmin))]
                else:
                    plts = [[ts,I0,[cs[i],'-o'],''] for i,I0 in enumerate(IE)]
                    txts = [[ts[idx],IE[i,idx],'%s' %str(h),cs[i]] for i,(h,idx) in enumerate(zip(refl,iSmin))]
                dsp.stddisp(plts,texts=txts,labs=[xlab,'$I$'],
                    title=r'thickness=$%d A$' %thick,name=figname+'_I.svg',**kwargs)
            # return IE
        # return Sw

    ###########################################################################
    #### gets
    ###########################################################################
    def get_frames(self,hkl:str,cols:Sequence[str]=['Sw','I']):
        """get frames in which beam hkl appears

        Parameters
        ---------
        hlk
            single reflection to gather information about
        cols
            info to fetch

        Returns
        -------
        DataFrame with info for that reflection
        """
        frames = [ i for i,name in enumerate(self.df.index)
            if hkl in self.load(i).df_G.index]

        df = pd.DataFrame(frames,columns=['frame'])
        df[cols]=0

        for i,f in enumerate(frames):
            vals = self.load(f).df_G.loc[hkl,cols].values #;print(f,vals)#df.iloc[i,cols])
            df.loc[i,cols] = vals.real
        return df

    def get_beams(self,cond='',refl=[]):
        if cond:
            refl = []
            for i,name in enumerate(self.df.index):
                b = self.load(i)
                refl += b.get_beam(cond=cond,index=False)
        refl = np.unique(refl)           #;print(refl)
        if not isinstance(refl[0],str):
            refl=[str(tuple(h)) for h in refl]
        nbs = len(refl)#.size;print(nbs,refl)
        print('total number of beams:%d' %nbs)
        return refl,nbs

    def _get_iTs(self,iTs,ts):
        t = self.ts
        if isinstance(ts,float) or isinstance(ts,int):ts = [ts]
        if isinstance(ts,list) or isinstance(ts,np.ndarray):
            iTs = [np.argmin(abs(t-t0)) for t0 in ts]
        if isinstance(iTs,int):iTs=[iTs]
        if not type(iTs) in [list,np.ndarray,slice]:iTs = slice(0,None,1)
        nts = t[iTs].size
        return iTs,nts

    def _get_iZs(self,iZs,zs):
        z = self.load(0).z
        if isinstance(zs,float) or isinstance(zs,int):zs = [zs]
        if isinstance(zs,list) or isinstance(zs,np.ndarray):
            iZs = [np.argmin(abs(z-z0)) for z0 in zs]
        if isinstance(iZs,int):iZs=[iZs]
        nzs = z[iZs].size
        return iZs,nzs

    def _get_z(self):
        z = self.load(0).z
        nzs = z.size
        return z,nzs

    def _get_ts(self,i,ts):
        if type(ts) in [float,int]:i=np.argmin(abs(self.ts-ts))
        return i,self.ts[i]

    ###########################################################################
    #### misc
    ###########################################################################
    def set_tag(self,tag):
        nts   = self.ts.size
        pad   = int(np.ceil(np.log10(nts)))
        names = [ '%s_%s%s' %(tag,'u',str(i).zfill(pad)) for i in range(nts)]
        self.df.index = names
        cmd ="cd %s; rename 's/%s/%s/' %s*.pkl df_%s.pkl rock_%s.pkl" %(self.path,self.tag,tag,self.tag,self.tag,self.tag)
        p = subprocess.Popen(cmd,shell=True);p.wait()
        self.tag = tag
        for i,name in enumerate(names):
            sim_obj = load_pkl(os.path.join(self.path,name+'.pkl'))
            sim_obj.set_name(name=name,path=self.path)
            sim_obj.save()
            self.df.loc[name,'pkl'] = sim_obj.get_pkl()
        df_file = os.path.join(self.path,'df_%s.pkl' %tag)
        self.df.to_pickle(df_file)
        print(colors.green+'Dataframe saved : '+colors.yellow+df_file+colors.black)
        self.save()

    def save(self,v=1):
        """save this object"""
        file = rock_name(self.path,self.tag)
        with open(file,'wb') as out :
            pickle5.dump(self, out, pickle5.HIGHEST_PROTOCOL)
        if v:print(colors.green+"object saved\n"+colors.yellow+file+colors.black)

    def load(self,i:int=0,ts:float=None):
        """load one of the simu

        Parameters
        ----------
        i
            index of the simu
        ts
            value
        """
        i,ts = self._get_ts(i,ts)
        file = self.df.iloc[i].pkl
        sim_obj = ut.load_pkl(file)
        return sim_obj

    def do(self,f,v=True,**args):
        """ apply function to all simus
        """
        for i in range(self.df.shape[0]):
            obj = self.load(i)
            obj.__getattribute__(f)(**args)
            obj.save(v=v)


def rock_name(path,tag):
    return os.path.join(path,'rock_%s.pkl' %(tag))

def load_rock(path,tag):
    return load_pkl(file=exp.rock_name(path,tag))
