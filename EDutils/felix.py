import importlib as imp
import numpy as np
import pandas as pd
import os,glob
from utils import displayStandards as dsp
from utils import plt
from . import dyn_utils as dyn
from . import utilities as ut;imp.reload(ut)
from . import pets as pt;imp.reload(pt)

def load_felix(path):
    return ut.load_pkl(os.path.join(path,'felix.pkl'))

class Felix():
    '''Felix class for generating full simulations

    Parameters :
    ------------
        path
            simulation folder (must contain reflprofiles_strong.dat and *dyn.cif_pets )
        xtal
            crystal name
        ww
            width of the extracted rocking curve +/-
        keV
            wavelength
        frames_per_sim
            frames per simulation
        wid
            radius of the LACBED pattern in pixels (the patterns are 2*wid x 2*wid pixels)
    '''
    def __init__(self,path,xtal,
            ww=40,keV=200,frames_per_sim=20,wid=200):
        self.path = os.path.realpath(path)
        self.xtal = xtal
        self.ww   = ww
        self.wid  = wid
        self.nwx  = 2*self.wid
        self.file_profiles = os.path.join(self.path,'reflprofiles_strong.dat')
        self.file_dyn = glob.glob(os.path.join(self.path,'*dyn.cif_pets'))[0]
        self.paths = {
            'rock_exp' : os.path.join(self.path,'rock_curves'),
            'optim'    : os.path.join(self.path,'optimization'),
            'ref'      : os.path.join(self.path,'DWF_refinement'),
            'rock'     : os.path.join(self.path,'plots'),
            }
        for f in self.paths:
            if not os.path.exists(f):os.mkdir(f)
        self.import_rocks()
        print('... load pets info ... ')
        self.load_dyn()
        print('... integrating exp reflections ... ')
        self.integrate_exp_intensity()
        self.save()

    def save(self):
        ut.save_pkl(self,file=os.path.join(self.path,'felix.pkl'))

    def load_dyn(self):
        #  Now use dyn.cif_pets to set up the felix simulation input files
        dyn = pt.load_dyn(self.file_dyn)
        # print(dyn)
        self.__dict__.update(dyn)
        # self.u,self.v,self.w,self.n_frames,self.frame_alpha,self.cell =

        # doc = cif.read_file(self.file_dyn)
        # block = doc.sole_block()
        #
        # # Extract data:
        # # first the frames and their orientation for the felix.inp file
        # # frame number fn
        # fn = np.array([int(i) for i in list(block.find_loop('_diffrn_zone_axis_id'))])
        # self.n_frames = len(fn)  # no. of frames
        # # incident beam orientation u,v,w
        # self.u = np.array([float(i) for i in list(block.find_loop('_diffrn_zone_axis_u'))])
        # self.v = np.array([float(i) for i in list(block.find_loop('_diffrn_zone_axis_v'))])
        # self.w = np.array([float(i) for i in list(block.find_loop('_diffrn_zone_axis_w'))])
        #
        # # Frame angles
        # self.frame_alpha = np.array([float(i) for i in
        #                         list(block.find_loop('_diffrn_zone_axis_alpha'))])
        # # Unit cell
        # cell_a = float(block.find_value('_cell_length_a'))
        # cell_b = float(block.find_value('_cell_length_b'))
        # cell_c = float(block.find_value('_cell_length_c'))
        # self.cell=[cell_a,cell_b,cell_c]

    def load_refl_df(self):
        cols = list(self.df_refl.columns.values)
        with open (self.file_profiles,'r') as f:lines = f.readlines()
        #fix touching columns
        for i in range(len(lines)) : lines[i] = lines[i][:96]+' '+lines[i][96:]
        # get reflection block lines
        l_refl = [i for i,l in enumerate(lines) if '#' in l]+[len(lines)]
        for l1,l2 in zip(l_refl[:-1],l_refl[1:]):
            with open ('tmp.csv','w') as f:f.write(''.join(lines[l1+1:l2-2]))
            df=pd.read_csv('tmp.csv',sep=' +',engine='python',
                names=['h','k','l']+cols,
                index_col=False)
            # hkl = str(tuple(np.array(df.loc[0,['h','k','l']],dtype=int)))
            hkl = ','.join(df.loc[0,['h','k','l']].values)
            self.df_refl.loc[hkl,cols[:-2]] = [list(df[c]) for c in cols[:-2]]
            self.df_refl.loc[hkl,'frame'] = list(np.array(df['frame'],dtype=int))
            self.df_refl.loc[hkl,'omega'] = list(df['omega'])
            # print(hkl)

    def load_refl(self):
        f = open(self.file_profiles, "r")
        # number of reflection
        n = ([])
        h_list = ([])  # h
        k_list = ([])  # k
        l_list = ([])  # l
        dstar_list = ([])  # reciprocal d-spacing
        s_list = ([])  # deviation parameter, A^-1
        rsg_list = ([])  # the ratio between the excitation error of the reflection
        # and the maximum excitation error spanned by the oscilation of the frame
        # Thus, if |Rsg|=0, the reflection is in the diffraction condition exactly
        # in the middle of the frame. If |Rsg|=1, the reflection is in exact Bragg
        # condition at the edge of the wedge covered by the frame.
        Iobs_list = ([])  # measured intensity
        sigma_list = ([])  # SD, as in I/sigma
        Ifit_list = ([])  # fitted intensity
        frame_list = ([])  # frame number
        frame_max = 0  # max frame number
        omega_list = ([])  # rotation angle
        while(True):
            # read next line
            line = f.readline()
            # empty line is EOF
            if not line:
                break
            if (line.find("#") >= 0):  # we have a new reflection
                frame = ([])  # list of frames for this reflection
                omega = ([])  # list of rotations for this reflection
                s = ([])  # list of deviation parameters for this reflection
                Iobs = ([])  # observed intensities for this reflection
                Ifit = ([])  # fit intensities for this reflection
                rsg = ([])  # rsg, see above
                sigma = ([])  # dunno
                n.append(int(line[line.find("#")+1:]))
                # first line of data for this reflection- never empty, we hope
                line = f.readline()
                h_list.append(int(line[0:4]))
                k_list.append(int(line[4:8]))
                l_list.append(int(line[8:12]))
                dstar_list.append(float(line[12:26]))
                s.append(float(line[26:40]))
                rsg.append(float(line[40:54]))
                Iobs.append(float(line[54:70]))
                sigma.append(float(line[70:84]))
                Ifit.append(float(line[84:96]))
                frame.append(int(line[96:100]))
                omega.append(float(line[100:109]))
                if (int(line[96:100]) > frame_max):
                    frame_max = int(line[96:100])
                # second & subsequent lines of data for this reflection
                line = f.readline()
                while (line[0] != "\n"):
                    s.append(float(line[26:40]))
                    rsg.append(float(line[40:54]))
                    Iobs.append(float(line[54:70]))
                    sigma.append(float(line[70:84]))
                    Ifit.append(float(line[84:96]))
                    frame.append(int(line[96:100]))
                    omega.append(float(line[100:109]))
                    if (int(line[96:100]) > frame_max):
                        frame_max = int(line[96:100])
                    line = f.readline()
                # add each parameter to their arrays
                # frame_list.append([frame])
                # omega_list.append([omega])
                # s_list.append([s])
                # Iobs_list.append([Iobs])
                # Ifit_list.append([Ifit])
                # rsg_list.append([rsg])
                # sigma_list.append([sigma])
                # print(tuple(h_list[-1],k_list[-1],l_list[-1]))
                refl = str((h_list[-1],k_list[-1],l_list[-1]))
                self.df_refl.loc[refl,'dstar']  = dstar_list[-1]
                self.df_refl.loc[refl,'s']      = s
                self.df_refl.loc[refl,'rsg']    = rsg
                self.df_refl.loc[refl,'Iobs']   = Iobs
                self.df_refl.loc[refl,'sigma']  = sigma
                self.df_refl.loc[refl,'Ifit']   = Ifit
                self.df_refl.loc[refl,'frame']  = frame
                self.df_refl.loc[refl,'omega']  = omega

    def load_rock(self):
        # Iobs is ordered according to deviation parameter s
        # reorder according to frame ID, and extract the peak over a frame range 2*ww+1
        # peak_f = ([])  # frame numbers for the extracted range
        # peak_s = ([])  # deviation parameters s for the extracted range
        # peak_I = ([])  # intensities for the extracted range
        # width of the extracted rocking curve +/-
        ww = self.ww
        # n_refl = len(self.h)
        # for i in range(n_refl):
        for i,r in self.df_refl.iterrows():
            # frame list for this reflection
            # c = np.array(sum(r.frame, []))
            c = np.array(r.frame,dtype=int)
            # Iobs for this reflection
            # a = np.array(sum(r.Iobs, []))
            a = np.array(r.Iobs)
            # s for this reflection
            # b = np.array(sum(r.s, []))
            b = np.array(r.s)
            # re-order according to frame number
            # print(a.shape,argsort(c))
            f1 = c[np.argsort(c)]
            a1 = a[np.argsort(c)]
            b1 = b[np.argsort(c)]
            # extract the peak in a 2*ww-frame window around the max
            # the different ways here just limit the window to the actual frame list
            pp = np.argmax(a1)  # pp=peak position in this sub-list
            if(pp-ww < 0):  # peak is close to the start
                ff = f1[:pp+ww]  # frame list
                ss = b1[:pp+ww]  # s list
                II = a1[:pp+ww]  # intensity list
            elif(pp+ww > len(a1)):  # peak is close to the end
                ff = f1[pp-ww:]
                ss = b1[pp-ww:]
                II = a1[pp-ww:]
            else:
                ff = f1[pp-ww:pp+ww]
                ss = b1[pp-ww:pp+ww]
                II = a1[pp-ww:pp+ww]
            # peak_f.append(ff)
            # peak_s.append(ss)
            # peak_I.append(II)
            refl = r.name
            # print(refl,ff,ss,II)
            # self.df_refl.loc[refl,['peak_f','peak_s','peak_I']] = [ff,ss,II]
            self.df_peak.loc[refl,'f'] = ff
            self.df_peak.loc[refl,'s'] = ss
            self.df_peak.loc[refl,'I'] = II


    def import_rocks(self):
        cols = ['dstar','s','rsg','Iobs','sigma','Ifit','frame','omega']
        self.df_refl=pd.DataFrame(columns=cols,dtype=object)
        self.load_refl()
        print('... ordering ...')
        self.df_peak = pd.DataFrame(columns=['f','s','I'],dtype=object,)
        self.load_rock()

    def plot_exp_rock(self,opt='s',show=False):
        outpath = self.paths['rock_exp']
        fi = np.argsort(self.df_int['fmax'])
        # print(fi.values)
        # if not show:opt+='c'
        n = int(np.ceil(np.log10(self.df_int.shape[0])))
        for i,(refl,r) in enumerate(self.df_peak.iloc[fi].iterrows()):
            name = '%s.%s.png' %(str(i).zfill(n),
                refl[1:-1].replace(', ','_'))
            x = np.array(r.f)
            y = np.array(r.I)
            dsp.stddisp([x,y,'b-o'],labs=['Frame','Intensity'],title=refl,
                xylims=[x.min(),x.max()+1,1e2*np.floor(y.min()/1e2),1e2*np.ceil(y.max()/1e2)],
                xyTicks=[5,y.max()/10],
                name = os.path.join(outpath,name) ,opt=opt)
            if show:plt.show()
            plt.close()

    def integrate_exp_intensity(self,thresh = 600,plot_rock=False,save_rock=False):
        f_max = ([])  # frame number for the peak
        I_int_expt_raw = ([])  # raw integrated intensities
        I_int_expt_Lor = ([])  # Lorentz scaled integrated intensities
        ds = ([])  # delta s for the peak
        outpath = self.paths['rock']

        for i,r in self.df_peak.iterrows():
            x = np.array(r.f)
            # default is that background is the minimum value
            back = np.where(r.I >= thresh, min(r.I), r.I)
            # alternative default: background is threshold value
            # back = np.where(peak_I[i]>=thresh, 1500, peak_I[i])
            ItB = back*1.0
            # we assume all data points below thresh are background
            if (min(r.I) < thresh):
                # fit a straight line to the background
                ItB = np.where(r.I >= thresh, np.nan, r.I)
                # squeeze & remove the points above background
                x1 = x[~np.isnan(ItB)]
                b1 = ItB[~np.isnan(ItB)]
                if (len(x1) > 1):
                    # least squares linear fit
                    fit = np.polyfit(x1, b1, 1)
                    # background interpolated under the peak
                    back = x*fit[0]+fit[1]
            # Intensities for the plot
            It0 = np.where(r.I < thresh, np.nan, r.I-back)
            # Intensities for integration (can't have NaNs)
            ItI = np.where(r.I < thresh, 0, r.I-back)
            # Integrated intensity without Lorentz scaling
            I_int_expt_raw.append(sum(ItI))
            # delta s per frame for this reflection
            deltas = abs((r.s[1]-r.s[0]))
            ds.append(deltas)
            # Apply Lorentz scaling, multiply by delta s
            I_int_expt_Lor.append(sum(ItI)*deltas)
            # centroid of thresholded rocking curve
            centr = np.dot(r.f, ItI)/np.sum(ItI)
            # frame number for the peak
            f_max.append(centr)
            # f_max.append(ff[np.argmax(II)])  # simple maximum, deprecated
            if (plot_rock):
                name = r.name
                fig,ax = dsp.create_fig(figsize=(6, 4))
                # ax = fig.add_subplot(111)
                ax.plot(x, r.I, color='grey')
                ax.bar(x, It0, color='green')
                ax.plot(x, back, color='red')
                ax.scatter(x, ItB, s=3, color='red')
                ax.set_xlabel('Frame')
                ax.set_ylabel('Intensity')
                plt.suptitle(name)
                plt.show
                if (save_rock):
                    plotname = os.path.join(outpath, name + '.png')
                    plt.savefig(plotname)
                    plt.close()

        print('integration done')
        self.df_int = pd.DataFrame(
            np.array([f_max,ds,I_int_expt_raw,I_int_expt_Lor]).T,
            columns=['fmax','ds','Iint_raw','Iint_lor'],
            index=self.df_refl.index)

    def correct_rock(self):
       return

    ###########################################################################
    ###########################################################################
    ### Generate and launch simulation
    ###########################################################################
    ###########################################################################
    def write_felix_inp(self,thicks,sig_fig=1e5):
       return


    ###########################################################################
    ###########################################################################
    ### post process : LACBED patterns,simulated integration
    ###########################################################################
    ###########################################################################
    def list_LACBED_hkl(self):
       return

    def import_sims(self,nwx=None):
        cols = ['sim','thick','nwx','refl','file']
        self.df_sims = pd.DataFrame(columns=cols)
        # print()
        sims_dir = glob.glob(os.path.join(self.path,self.xtal,'*/'))
        sims = np.sort([int(s.split('/')[-2]) for s in sims_dir])
        for sim,sim_d in zip(sims,sims_dir):
            folders = glob.glob(os.path.join(sim_d,'*/'))
            for d in folders:
                # print(d.split('/'))#[-1].split('_')[2:])
                thick,nwxs = d.split('/')[-2].split('_')[2:]
                # nwx0 = tuple([int(i) for i in nwxs[:-1].split('x')])
                nwx0 = int(nwxs[:-1].split('x')[0])
                if not type(nwx)==type(None):
                    if not int(nwx0)==nwx:
                        continue
                files = glob.glob(os.path.join(d,'*'))
                refls = [f.split('_')[-1][:-4] for f in files]
                hkl = [(',-'.join(','.join(r.split('+')).split('-'))[1:]) for r in refls]
                n = len(files)
                self.df_sims = pd.concat([
                    self.df_sims,
                    pd.DataFrame(np.array([[sim]*n,[thick]*n,[nwx0]*n,hkl,files]).T,columns=cols),
                    ])
        # print('done')

    def show_lacbed(self,show=0,**kwargs):
        for i,r in self.df_sims.iterrows():
            I = np.reshape(np.fromfile(r.file,dtype=np.float64),(int(r.nwx),)*2)
            # tle = os.path.basename(r.file)
            tle='refl=%s, thick=%s, sim %s,' %(r.sim,'(%s)' %r.refl,r.thick)
            dsp.stddisp(im=[I],pOpt='im',title=tle,
                name=r.file.replace('.bin','.png'),**kwargs)
            if show:plt.show()


    def find_best_row_sim(self):
        ''' Find the best row for each simulation'''
        return
    def extract_sim_integrated_intensity(self):
        return
    def calculate_dynamical_Rfactors(self):
        return
