import importlib as imp
import tifffile,os,glob,pickle5,subprocess
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np,pandas as pd
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Sequence, Union
from utils import displayStandards as dsp   #;imp.reload(dsp)
from utils import glob_colors as colors     #;imp.reload(colors)
from utils import handler3D as h3D          #;imp.reload(h3D)
from scipy.integrate import trapz
from scipy.stats import linregress
from . import utilities as ut               ;imp.reload(ut)

class Rocking:
    def __init__(self,Simu:object,
            uvw:list,tag:str,path:str,
            Sargs:dict,frames:list=None,
            params=[],vals=[],verbose=True):
        """ simulate rocking curve

        Parameters
        ----------
        Simu
            Simulator object
        uvw
            list of orientation vectors
        tag
            a tag identifying the sweep
        Sargs
            Simulator constructor arguments common to all simulations
        """
        self.path = path
        self.tag  = tag
        self.uvw  = uvw
        self.Sargs = Sargs
        self.n_simus  = uvw.shape[0]

        params+=['u']
        vals+=[uvw]
        if any(frames):
            params+=['frame']
            vals+=[frames]

        #reshape vals into [val_param_1,..,val_param_N] x n_simus
        _vals = []#[[v] for v in vals[0]]
        for i in  range(self.n_simus):
            _vals.append([row[i] for row in vals ])


        self.df = ut.sweep_var(Simu,params=params,vals=_vals,tag=tag,path=path,**Sargs)
        ts            = np.arange(self.n_simus)
        self.ts       = ts
        self.df['ts'] = ts
        self.Iz_dyn = {}
        self.Iz_kin = {}
        self._build_index(v=verbose)
        self.save(v=verbose)

    ###########################################################################
    #### compute
    ###########################################################################
    def _build_index(self,v=True):
        frames = [ i for i,name in enumerate(self.df.index)]
        beams = np.unique(np.hstack([self.load(i).df_G.index.values for i,name in enumerate(self.df.index)]))
        self.nbeams = beams.size
        if v:print('...building index...')
        self.beams = pd.DataFrame(index=beams,columns=['Frame','Sw'])
        for h in beams:
            self.beams.loc[h,'Frame']=[]
            self.beams.loc[h,'Sw']=[]
        for i,f in enumerate(self.df.index):
            for h,c in self.load(i).df_G.iterrows():
                self.beams.loc[h,'Frame'].append(i)
                self.beams.loc[h,'Sw'   ].append(c.Sw.real)
        self.df['nbeams']=[self.load(i).df_G.index.shape[0] for i,f in enumerate(self.df.index)]
        self.beams['nframes' ] = [len(r) for r in self.beams.Frame]
        self.beams['f_min' ] = [np.min(r.Frame) for h,r in self.beams.iterrows()]
        self.beams['f_cen' ] = [r.Frame[np.abs(r.Sw).argmin()] for h,r in self.beams.iterrows()]
        self.beams['f_max' ] = [np.max(r.Frame) for h,r in self.beams.iterrows()]
        self.beams['f_range'] = [(r.f_min,r.f_cen,r.f_max) for h,r in self.beams.iterrows()]
        self.beams['Sw_min'] = [np.min(r.Sw) for h,r in self.beams.iterrows()]
        self.beams['Sw_cen'] = [np.min(np.abs(r.Sw)) for h,r in self.beams.iterrows()]
        self.beams['Sw_max'] = [np.max(r.Sw) for h,r in self.beams.iterrows()]

    def get_full_refl(self,Swm=0.01):
        hkl = self.beams.loc[[ (r.Sw_min<-Swm) & (r.Sw_max>Swm)
            for h,r in self.beams.iterrows()]].index
        return hkl

    def show_beams(self,hkl=[],cols=['f_range','Sw_min','Sw_cen','Sw_max'],n=-1):
        e='{:>7.2e}'
        formats = {'Sw_min':e,'Sw_cen':e,'Sw_max':e}
            # 'f_range':'{:>6.2f}'}
        df = self.beams.sort_values('Sw_cen')
        if any(hkl):
            df=df.loc[hkl]
        if n:
            df=df[:n]
        print(df[cols].to_string(
                formatters={k: v.format for k, v in formats.items()}))

    def reset_int(self):
        self.Iz_dyn,self.Iz_kin = {},{}

    def _integrate_rocking(self,refl=[],new=0):
        if not new:
            # refl,nbs = self.get_beams(cond=cond,refl=refl)
            nbs=len(refl)
            refl  = [h for h in refl if not h in self.Iz_dyn.keys()] #;print(hkl)

        hkl = [eval(h) for h in refl]
        z,nzs = self._get_z()
        nbs,nts = len(hkl),self.ts.size
        if nbs:
            Iz_dyns  = dict(zip(refl, [np.zeros((nts,nzs)) for h in refl] ))
            Iz_kins  = dict(zip(refl, [np.zeros((nts,nzs)) for h in refl] ))
            for i in range(nts):
                sim_obj = self.load(i)
                idx = sim_obj.get_beam(refl=refl,cond='',index=True)
                if idx:
                    hkl0 = sim_obj.df_G.iloc[idx].index #[str(tuple(h)) for h in sim_obj.get_hkl()[idx]]
                    for idB,hkl_0 in zip(idx,hkl0):
                        Iz_dyns[hkl_0][i,:] = sim_obj.Iz[idB,:]
                        Iz_kins[hkl_0][i,:] = sim_obj.Iz_kin[idB,:]

            Iz_dyn  = dict()
            Iz_kin  = dict()
            for h in refl:
                df_b = self.beams.loc[h]
                if len(df_b.Sw)>1:
                    s = np.sign(df_b.Sw[1]-df_b.Sw[0])
                    Iz_dyn[h]=[s*trapz(Iz_dyns[h][df_b.Frame,iz],df_b.Sw) for iz in range(nzs)]
                    Iz_kin[h]=[s*trapz(Iz_kins[h][df_b.Frame,iz],df_b.Sw) for iz in range(nzs)]
                else:
                    Iz_dyn[h] = [Iz_dyns[h][df_b.Frame,iz] for iz in range(nzs)]
                    Iz_kin[h] = [Iz_kins[h][df_b.Frame,iz] for iz in range(nzs)]
            self.Iz_dyn.update(Iz_dyn)
            self.Iz_kin.update(Iz_kin)
            self.save()
            print(colors.green+'rock.Iz updated'+colors.black)


    def get_rocking(self,iZs:Optional[slice]=-1,
        zs:Optional[Iterable[float]]=None,
        refl:Sequence[tuple]=[],cond:str='',opts:str='',n=0,
        kin:bool=False,
        v=False):
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
        opts
            F(include friedel), O(include central beam)
        kin
            Also include kinematic values if True(default False)
        Returns
        -------
            z,dict
                z,{beam:I(z)}
        """
        iZs,nzs  = self._get_iZs(iZs,zs)        #;print(iZs)
        z0 = self.load(0).z.copy()[iZs][-1]

        refl,nbs = self.get_beams(cond=cond,refl=refl,opts=opts)  #;print(refl)
        nts = self.ts.size
        I = {}
        # for h in refl : I[str(h)]=np.nan*np.ones((nts,nzs))
        for h in refl :
            if kin:
                I[str(h)]=[np.zeros((nts,nzs)),np.zeros((nts,nzs))]
            else:
                I[str(h)]=np.zeros((nts,nzs))

        if v : print("gathering the intensities")
        for i in range(nts):
            # print(colors.red,i,colors.black)
            sim_obj = self.load(i)
            hkl0  = sim_obj.get_beam(refl=refl,index=False)
            if hkl0:
                idx  = sim_obj.get_beam(refl=refl,index=True)
                # hkl0 = [str(tuple(h)) for h in sim_obj.get_hkl()[idx]]
                for idB,hkl_0 in zip(idx,hkl0):
                    if kin:
                        I[hkl_0][0][i,:] = np.array(sim_obj.Iz[idB,iZs])
                        I[hkl_0][1][i,:] = np.array(sim_obj.Iz_kin[idB,iZs])
                    else:
                        I[hkl_0][i,:] = np.array(sim_obj.Iz[idB,iZs])
        if n and nbs>n:
            print('keeping only %d strongest beams' %n)
            df_Imax = pd.DataFrame.from_dict({h:Ib[:,0].max() for h,Ib in I.items()},orient='index',columns=['I'])
            hkls = df_Imax.sort_values('I',ascending=False)[:n].index
            I = {h:I[h] for h in hkls}

        z = self.load(0).z.copy()[iZs]
        return z,I

    
    def get_rocking_old(self,iZs:Optional[slice]=-1,
        zs:Optional[Iterable[float]]=None,
        refl:Sequence[tuple]=[],cond:str='',opts:str='',n=0,
        kin:bool=False,
        v=False):
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
        opts
            F(include friedel), O(include central beam)
        kin
            Also include kinematic values if True(default False)
        Returns
        -------
            z,dict
                z,{beam:I(z)}
        """
        iZs,nzs  = self._get_iZs(iZs,zs)        #;print(iZs)
        z0 = self.load(0).z.copy()[iZs][-1]
        # if not self.load(0).thick==z0:
        #     print('setting thickness to %dA' %z0)
        #     self.do('_set_I',verbose=0,iZ=iZs[-1])

        refl,nbs = self.get_beams(cond=cond,refl=refl,opts=opts)  #;print(refl)
        nts = self.ts.size
        I = {}
        # for h in refl : I[str(h)]=np.nan*np.ones((nts,nzs))
        for h in refl :
            if kin:
                I[str(h)]=[np.zeros((nts,nzs)),np.zeros((nts,nzs))]
            else:
                I[str(h)]=np.zeros((nts,nzs))

        if v : print("gathering the intensities")
        for i in range(nts):
            # print(colors.red,i,colors.black)
            sim_obj = self.load(i)
            hkl0  = sim_obj.get_beam(refl=refl,index=False)
            if hkl0:
                idx  = sim_obj.get_beam(refl=refl,index=True)
                # hkl0 = [str(tuple(h)) for h in sim_obj.get_hkl()[idx]]
                for idB,hkl_0 in zip(idx,hkl0):
                    if kin:
                        I[hkl_0][0][i,:] = np.array(sim_obj.Iz[idB,iZs])
                        I[hkl_0][1][i,:] = np.array(sim_obj.Iz_kin[idB,iZs])
                    else:
                        I[hkl_0][i,:] = np.array(sim_obj.Iz[idB,iZs])
        if n and nbs>n:
            print('keeping only %d strongest beams' %n)
            df_Imax = pd.DataFrame.from_dict({h:Ib[:,0].max() for h,Ib in I.items()},orient='index',columns=['I'])
            hkls = df_Imax.sort_values('I',ascending=False)[:n].index
            I = {h:I[h] for h in hkls}

        z = self.load(0).z.copy()[iZs]
        return z,I

    def integrate(self,thick=None):
        if not thick:
            return self._integrate_all()
        else:
            return self._integrate_thick(thick)

    def _integrate_thick(self,thick):
        b0 = self.load(0)
        if not b0.thick==thick:
            self.do('set_thickness',thick=thick)
        df_int = pd.DataFrame(0,index=self.beams.index,columns=['I'])
        for i in range(self.n_simus):
            b0 = self.load(i)
            df_int.loc[b0.df_G.index,'I'] += b0.df_G.I
        # self.save()
        return df_int

    def _integrate_all(self):
        b0 = self.load(0)
        thicks=['%.1fA' %z for z in b0.z ]
        self.z=b0.z
        self.df_int = pd.DataFrame(0,index=self.beams.index,columns=thicks)

        for i in range(self.n_simus):
            b0 = self.load(i)
            self.df_int.loc[b0.df_G.index] += b0.Iz
        self.save()
        # print(self.z.shape,self.df_int.values.shape)
        return self.df_int

    def Rfactor(self,df_exp):
        refl = self.df_int.loc[self.df_int.index.isin(df_exp.index)].index
        I_exp = df_exp.loc[refl,'I'].values             #;print(I_exp)
        prescale=I_exp.mean()

        z=self.df_int.columns
        self.R = pd.DataFrame(index=z,columns=['scale','r_value','Rfac'])
        Isum = I_exp.sum()
        for z0 in z :
            I_sim = self.df_int.loc[refl,z0].values*prescale     #;print(I_sim)
            scale, intercept, r_value, p_value, std_err = linregress(I_sim, I_exp)
            Rfac = abs(I_sim*scale - I_exp).sum()/Isum
            self.R.loc[z0] = [scale,r_value,Rfac]

        self.save()
        return self.R
    ###########################################################################
    #### Display
    ###########################################################################
    def show_excitation_map(self,cmaps=('Reds','Reds_r'),hkls=[],
            sw_color=lambda Sw:Sw,vm=0.05,vmin=-0.005,vmax=0.005,
            figs=(20,5),nb_max=80,sw_min=1e-3,sort='nframes',**kwargs,
        ):
        # cmaps=('PuRd','OrRd_r')
        # swc,vm = lambda Sw:-np.sign(Sw)*np.log10(np.maximum(np.abs(Sw),1e-10)),4
        cm = mcolors.LinearSegmentedColormap.from_list('my_cmap', np.vstack((
            plt.get_cmap(cmaps[0])(np.linspace(0, 1, 128)),
            plt.get_cmap(cmaps[1])(np.linspace(0, 1, 128)),
        )))
        if vm:
            vmin=-vm
            vmax=vm


        if not any(hkls):
            hkls  = self.beams.index
        if len(hkls)>nb_max:
            hkls = self.beams.loc[self.beams.Sw_cen< sw_min].index
            print('%d beam within sw_min=%.2e' %(len(hkls),sw_min))
            if len(hkls)>nb_max:
                print('too many beams %d/%d : keeping only top %d' %(len(hkls),nb_max,nb_max))
                hkls = self.beams.sort_values('Sw_cen')[:nb_max].index

        # sort by how long they sitck around
        df = self.beams.loc[hkls]
        if sort :
            df = self.beams.loc[hkls].sort_values('nframes',ascending=False)
        hkls    = df.index
        df['i'] = np.arange(len(hkls))

        fig,ax = dsp.create_fig(figsize=figs)
        ##highlight beams with full rocking curves
        for h,r in df.iterrows():
            ax.scatter([r.i]*len(r.Frame),r.Frame,40,sw_color(r.Sw),marker='s',cmap=cm,vmin=vmin,vmax=vmax)
        hkl_f = df.index[df.index.isin(self.get_full_refl(Swm=vmax))]
        idx,f_cen = df.loc[hkl_f,['i','f_cen']].values.T
        ax.plot(idx,f_cen,'gs',markersize=8,mew=2,mfc='none',mec='g')
        # ax.scatter(idx,f_cen,50,'g',marker='o',linewidths=3)

        ax.tick_params(axis='x',direction='in',labelrotation=90);
        ax.set_xticks(list(range(len(hkls))));
        ax.set_xticklabels(list(hkls));
        args=dict(labs=['','Frames'],axPos=[0.05,0.3,0.93,0.65],
            xylims=['x',-1,len(hkls)])
        args.update(kwargs)
        dsp.standardDisplay(ax=ax,**args);
        return fig,ax

    def plot_integrated(self,refl,new:bool=False,cm='Spectral',kin=False,
        **kwargs):
        """plot the integrated rocking curves for selected beams as function of z

        Parameters
        ----------
        refl
            beams to consider
        kin
            True to overlay kinematic integration
        new
            update integration
        """
        self._integrate_rocking(refl=refl,new=new)
        nbs = len(refl)
        z = self.load(0).z
        Iz = np.array([self.Iz_dyn[h] for i,h in enumerate(refl)])

        cs,legElt = dsp.getCs(cm,nbs),{}
        plts = [[z,Iz[i],cs[i],'%s' %h] for i,h in enumerate(refl)]
        if kin:
            Iz = np.array([self.Iz_kin[h] for i,h in enumerate(refl)])
            plts += [[z,Iz[i],[cs[i],'--'],''] for i,h in enumerate(refl)]
            legElt={'dyn':'k-','kin':'k--'}
        fig,ax=dsp.stddisp(plts,labs=[r'$z(\AA)$','$I_{int}$'],
            legElt=legElt**kwargs)
        vw=Rock_viewer(self,fig,ax,z,Iz,refl)
        return fig,ax


    def plot_rocking(self,cmap='viridis',x:str='Sw',
        cond='',refl=[],opts:str='',iZs=-1,zs=None,n:int=0,
        kin=False,v=False,ms='',ls='',
        **kwargs):
        """plot rocking curve for set of selected beams at thickness zs

        Parameters
        ----------
        cond,refl,iZs,zs,opts
            select beams and thicknesses (see :meth:`~Rocking.get_rocking`)
        x
            to display on x axis (Frame,Sw,theta)
        kwargs
            args to pass to dsp.stddisp
        """

        z,I = self.get_rocking(cond=cond,refl=refl,opts=opts,zs=zs,iZs=iZs,n=n,kin=kin)

        refl,plts = list(I.keys()),[]           #;print(refl)
        xlab = {'Frame':'frame','Sw':r'$S_w(\AA^{-1})$','theta':r'$\theta(deg)$'}[x]

        if v:print('gathering plots')
        nbs,nzs = len(refl),z.size
        if not ls:
            ls={'kin':'--','dyn':'-'}        
        if nbs>=nzs:
            cs = dsp.getCs(cmap,nbs)
            if not ms:ms=dsp.markers
            legElt = { '%s' %refl0:[cs[i],'-'] for i,refl0 in enumerate(refl)}
            for iz,zi in enumerate(z):
                legElt.update({'$z=%d A$' %(zi):['k',ms[iz]+'-']})
                for i,refl0 in enumerate(refl):
                    df_b=self.beams.loc[refl0]
                    if kin :
                        plts += [[df_b[x],I[refl0][0][df_b.Frame,iz],[cs[i],ms[i]+ls['dyn']],'']]
                        plts += [[df_b[x],I[refl0][1][df_b.Frame,iz],[cs[i],ms[i]+ls['kin']],'']]
                    else:
                        plts += [[df_b[x],I[refl0][df_b.Frame,iz],[cs[i],ms[iz]+ls['dyn']],'']]

        else:
            # rocking for different thicknesses
            cs = dsp.getCs(cmap,nzs)
            if not ms:ms=dsp.markers
            legElt = { '%s' %refl0:['k','-'+ms[i]] for i,refl0 in enumerate(refl)}
            # self.get_frames(hkl,iTs=slice(0,None))
            for i,refl0 in enumerate(refl):
                for iz,zi in enumerate(z):
                    df_b=self.beams.loc[refl0]
                    if kin :
                        plts += [[df_b[x],I[refl0][0][df_b.Frame,iz],[cs[iz],ms[i]+ls['dyn']],'']]
                        plts += [[df_b[x],I[refl0][1][df_b.Frame,iz],[cs[iz],ms[i]+ls['kin']],'']]
                    else:
                        plts += [[df_b[x],I[refl0][df_b.Frame,iz],[cs[iz],ms[i]+ls['dyn']],'']]
            legElt.update({'$z=%d A$' %(zi):[cs[iz],'-'] for iz,zi in enumerate(z) })

        if kin :
            legElt['kin'] = 'k--'
        # print('displaying')
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

        # dsp.stddisp([[r.Frame,r.Sw,[c,'-o'],h] for c,(h,r) in zip(dsp.getCs('jet',len(refl)),rock.beams.loc[refl].iterrows()) ] )

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

    def get_beams(self,cond='',refl=[],opts='',n=None,v=True):
        if cond:
            refl = []
            for i,name in enumerate(self.df.index):
                b = self.load(i)
                refl += b.get_beam(cond=cond,index=False)
        refl = np.unique(refl)           #;print(refl)
        if not isinstance(refl[0],str):
            refl=[str(tuple(h)) for h in refl]
        if not 'F' in opts:
            refl=ut.remove_friedel_pairs(refl)
        if 'O' in opts:
            refl=np.hstack([refl,str((0,0,0))])

        if n and len(refl)>n:refl=refl[:n]
        nbs = len(refl)#.size;print(nbs,refl)
        if v :print('total number of beams:%d' %nbs)

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
        # file = os.path.join(self.path,os.path.basename(self.df.iloc[i].pkl))
        sim_obj = ut.load_pkl(file)
        return sim_obj

    def change_path(self,path):
        if not self.path==path:
            self.path=path
            self.df.pkl=[s.replace(os.path.dirname(s),self.path) for s in self.df.pkl]
            self.save()
            print(self.df.pkl)
            for i in range(self.df.shape[0]):
                obj = self.load(i)
                obj.path=path
                obj.save(v=0)

    def do(self,f,verbose=True,**args):
        """ apply function to all simus
        """
        for i in range(self.df.shape[0]):
            obj = self.load(i)
            if verbose:print('rock %d' %i)
            obj.__getattribute__(f)(**args)
            obj.save(v=verbose)


def rock_name(path,tag):
    return os.path.join(path,'rock_%s.pkl' %(tag))

def load_rock(path,tag):
    return load_pkl(file=exp.rock_name(path,tag))

class Rock_viewer:
    def __init__(self,rock,fig,ax,z,Iz,refl):
        self.rock=rock
        self.fig = fig
        self.ax = ax
        self.z  = z
        self.Iz = Iz
        self.refl=refl
        cid = self.fig.canvas.mpl_connect('key_press_event', self)

    def __call__(self, event):
        # print(event.key)
        keys={
            'single':'enter' ,
            'thicks':'z'     ,
        }
        if event.key == keys['single']:
            z,I=dsp.plt.ginput(n=1)[0]
            iz = abs(self.z-z).argmin()
            ih = abs(I-self.Iz[:,iz]).argmin()
            print(iz,ih)
            refl,iZs = [self.refl[ih]],[iz]

        elif event.key == keys['thicks']:
            dats=dsp.plt.ginput(n=-1)
            iZs = [abs(self.z-z).argmin() for z,I in dats]
            I,z = dats[0]
            ih  = abs(I-self.Iz[:,iZs[0]]).argmin()
            refl,iZs = [self.refl[ih]],iZs

        if event.key in keys.values():
            self.rock.plot_rocking(refl=refl,iZs=iZs)
