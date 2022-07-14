# -*- coding: utf-8 -*-
# %% 0: Subroutines
"""
Created on Mon Jul 12 17:16:01 2021

@author: Richard

Import data from a reflprofiles PETS file
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from gemmi import cif


def strucfac(hkl, B):
    # atom coordinates
    aSi = 5.432
    r = np.array([0.0, 0.0, 0.0,
                  0.5, 0.5, 0.0,
                  0.0, 0.5, 0.5,
                  0.5, 0.0, 0.5,
                  0.25, 0.25, 0.25,
                  0.75, 0.75, 0.25,
                  0.25, 0.75, 0.75,
                  0.75, 0.25, 0.75])
    r = r.reshape(8, 3)
    gmag = (1/aSi)*np.dot(hkl, hkl)**0.5
    # calculate scattering factor
    aSi = ([1.065438920E+00, 1.201436910E-01, 1.809152630E-01])
    bSi = ([1.041184550E+00, 6.871133680E+01, 8.875339260E-02])
    cSi = ([1.120656200E+00, 3.054528160E-02, 1.599635020E+00])
    dSi = ([3.700626190E+00, 2.140978970E-01, 9.990966380E+00])
    fSi = 0
    for i in range(3):
        fSi = fSi + aSi[i]/((gmag**2)+bSi[i])
        + cSi[i]*np.exp(-(dSi[i]*gmag**2))

    # structure factor
    F = np.sum(fSi*np.exp(-B*gmag*gmag/4)*np.exp(2*np.pi*1j*np.dot(r, hkl)))
    return(F)


def writefelixinp(u, v, w, x1, x2, x3, npix, lpha, th_start, th_end, th_step):
    # currently this just writes a basic Felix file
    # but expect more input arguments to be added

    fel = open('felix.inp', 'w')

    fel.write('# Input file for Felix version :master: v1.03 :14Mar2019:\n')
    fel.write('# ------------------------------------\n\n')
    fel.write('# ------------------------------------\n\n')
    fel.write('# control flags\n')
    fel.write('IWriteFLAG                = 1\n')
    fel.write('IImageFLAG                = 1\n')
    fel.write('IScatterFactorMethodFLAG  = 0\n')
    fel.write('IBlochMethodFLAG          = 0\n')
    fel.write('IMaskFLAG                 = 1\n')
    fel.write('IHolzFLAG                 = 0\n')
    fel.write('IAbsorbFLAG               = 2\n')
    fel.write('IAnisoDebyeWallerFlag     = 0\n')
    fel.write('IByteSize                 = 8\n\n')
    fel.write('# radius of the beam in pixels\n')
    fel.write('IPixelCount               = '+npix+'\n\n')
    fel.write('# beam selection criteria\n')
    fel.write('IMinReflectionPool        = 200\n')
    fel.write('IMinStrongBeams           = 100\n')
    fel.write('IMinWeakBeams             = 0\n\n')
    fel.write('# crystal settings\n')
    fel.write('RDebyeWallerConstant      = 0.0\n')
    fel.write('RAbsorptionPer            = 10.0\n\n')
    fel.write('# microscope settings\n')
    fel.write('ROuterConvergenceAngle    = '+lpha+'\n')
    # Incident beam - from the input
    fel.write('IIncidentBeamDirection    = ['+u+','+v+','+w+']\n')
    # x-direction as calculated from cross product
    fel.write('IXDirection               = ['+x1+','+x2+','+x3+']\n')
    # not very sure about this - assume it is the beam direction
    # fel.write('INormalDirection          = ['+u+','+v+','+w+']\n')
    fel.write('INormalDirection          = [-1,0,1]\n')
    fel.write('RAcceleratingVoltage (kV) = 200.0\n')
    fel.write('RAcceptanceAngle          = 0.0\n\n')
    fel.write('# Image Output Options\n')
    fel.write('RInitialThickness        = '+str(th_start)+'\n')
    fel.write('RFinalThickness          = '+str(th_end)+'\n')
    fel.write('RDeltaThickness          = '+str(th_step)+'\n')
    fel.write('RPrecision               = 0.00002\n\n')
    fel.write('#Refinement Specific Flags\n')
    fel.write('IRefineModeFLAG          = S\n')  # simulation only
    # all these lines are irrelevant for sim only
    fel.write('IWeightingFLAG           = 0\n')
    fel.write('IRefineMethodFLAG        = 3\n')
    fel.write('ICorrelationFLAG         = 2\n')
    fel.write('IImageProcessingFLAG     = 4\n')
    fel.write('RBlurRadius              = 0.0\n')
    fel.write('INoofUgs                 = 40\n')
    fel.write('IAtomicSites             = (1)\n')
    fel.write('IPrint                   = 0\n')
    fel.write('RSimplexLengthScale      = 0.1\n')
    fel.write('RExitCriteria            = 0.000001\n')

    fel.close()


def hklfile(xtal, hkl):
    # make a filename corresponding to a felix output e.g. Si_+0-2+7.bin
    if (hkl[0] < 0):
        h_str = str(hkl[0])
    else:
        h_str = '+'+str(hkl[0])
    if (hkl[1] < 0):
        k_str = str(hkl[1])
    else:
        k_str = '+'+str(hkl[1])
    if (hkl[2] < 0):
        l_str = str(hkl[2])
    else:
        l_str = '+'+str(hkl[2])
    fname = xtal + '_' + h_str+k_str+l_str + '.bin'

    return fname


def gaussian(r, sigma):
    return np.exp(-(r*r/(2*sigma)))


def conv2D(imagein, kernel):
    wx = len(imagein)
    wc = len(kernel)
    cc = int(wc/2)
    imgout = imagein*0
    for i in range(wc):
        for j in range(wc):
            imgout[i:wx-wc+i, j:wx-wc+j] += imagein[cc:wx+cc-wc,
                                                    cc:wx+cc-wc]*kernel[i, j]
    return imgout


def optimise_orientation(xtal,sim, local_hkl_list, local_fmax_list, inpath, outpath,
                         fram, kernel, plot, save):
    # Get the peak positions (frame number) for reflections in a given LACBED
    # pattern and vary the x-coordinate of the origin so that the peaks
    # lie along a line parallel to x
    # number of DF LACBED patterns
    if(plot):
        print('Simulation:%d' %sim)
    # os.chdir(inpath)
    start = -9+inpath[-10:].find('x')
    wx = int(inpath[start:])  # the patterns are wx x wx pixels
    sim_lacbed = ([])
    # load the LACBED patterns for this simulation
    for hkl in local_hkl_list:
        fnom = hklfile(xtal, hkl)
        lacbed_name = os.path.join(inpath,fnom)
        print(lacbed_name)
        if(os.path.exists(lacbed_name)):
            if '+0+0+0' not in lacbed_name:
                sim_lacbed.append(conv2D(
                    np.fromfile(lacbed_name, dtype=np.float64)
                    .reshape(wx, wx), kernel))
        else:
            print(fnom+' not found')

    # now find out if the offset in the frame indexing and the best line
    # that matches the relative positions of the maxima in different frames.
    # We assume that the x-direction is correct.
    # The calculation is quick so just go through a range and find the
    # one which gives the closest to a straight line along x
    fom_list = ([])
    best_fom = 1.0e+10  # a big number, we want a minimum
    frame_range = 10  # how far we look for the origin +/-
    # NB first frame starts at x=wx/4, each frame covers wx/(2*fram) pixels
    # number of first frame is 2*sim*fram +1 (e.g. fram=10 gives 1,21,41...)
    frf = frame_range*fram
    for offset in range(-frf, frf):
        row_list = ([])
        # if(plot):
        #     plt.figure()
        #     plt.title(str(offset))
        # get the figure of merit - i.e. the sum of squares of
        # deviations of peak positions from a horizontal line
        # looping over the reflections in this simulation
        for j in range(len(sim_lacbed)):
            # pixel coords of the maximum for this reflection
            # we have sub-frame accuracy since fmax is a centroid
            col = int(wx/4 + fram*(local_fmax_list[j] - 2*sim*fram)+offset+0.5)
            # print(str(offset)+', '+str(j)+'; max '+str(local_fmax_list[j])
            # +' =column '+str(col)+'\n')
            # intensity for the specified column
            I_sim = sim_lacbed[j][:, col]
            # if(plot):
            #     plt.plot(I_sim)
            # the row corresponding to the centroid intensity of this column
            # if PETS gave a perfect answer this would be the central row of
            # the simulated image for all reflections
            row_list.append(np.dot(list(range(wx)), I_sim)/np.sum(I_sim))
        # if(plot) & (save):
        #     plotname = outpath + '\\' + str(sim) + '_' +str(offset)+ '.png'
        #     plt.savefig(plotname)
        #     plt.close()
        # find the sum of squares of the deviation from the average row
        # this is our figure of merit fom_row
        mean_row = np.mean(row_list)
        fom_row = np.dot((row_list-mean_row), (row_list-mean_row))
        if (fom_row < best_fom):
            best_fom = fom_row
            best_row = int(np.rint(mean_row))
            best_offset = offset/fram
        fom_list.append(fom_row)
    # plot if required
    if(plot):
        plt.figure()
        plt.plot(np.arange(-frame_range, frame_range, 1/fram), fom_list)
        plt.title(str(sim))
        plt.xlabel('Frame offset')
        plt.ylabel('Figure of merit')
        # plt.show
        print('best row = '+str(best_row)+', offset = '
              + str(best_offset) + ' frames')
    # save plot if required
    if(save):
        plotname = os.path.join(outpath,str(sim) + '.png')
        plt.savefig(plotname)
        plt.close()

    return best_row, best_offset
