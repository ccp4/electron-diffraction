''' Computes path length for all pixel in a detector with respect 
to :
    Incident beam with potentially spatial and spectral broadening
    Atomic distribution with potentially non zero B factor
features :
    - Fraunhofer diffraction framework
    - Loss of coherency option
'''
#nearBragg program
#main                                                   #58
#parse arguments                                        #150
#display help                                           #400
#allocate                                               #450
#handle parameters                                      #500
#load lattice points                                    #750
#initializing                                           #815

#The mega loop
for yp,xp in zip(ypixels,xpixels):                      #950
    Fa,Fb = 0,0
    for suby,subz in zip(oversample,oversample):        #976
        pixel=[distance,Ydet-Ybeam,Xdet-Xbeam]
        steps = divsteps*dispsteps*depthsteps*oversample*oversample
        for lam in wavelengths:                         #1010
            for hdiv,vdiv in solid_angle:               #1015
                for source_tic in source_pts:           #1024
                    for i in atoms:                     #1042
                        phase=2*pi*(source_to_atom_path+atom_to_pixel_path)/lam + phsft
                        S0,S = atom-source, pixel-atom
                        stol=norm(S-S0)/lam*1e10
                        DWF = occ[i]*exp(-Bfac[i]*stol**2)
                        Fa+=DWF*cos(phase)/source_to_atom_path/atom_to_pixel_path
                        Fb+=DWF*sin(phase)/source_to_atom_path/atom_to_pixel_path
                if not coherent : I += Fa**2 + Fb**2; Fa,Fb=0,0
            if not coherent : I += Fa**2 + Fb**2; Fa,Fb=0,0
    if coherent :                                       #1120
        Fa += cosimage[j]*steps
        Fb += sinimage[j]*steps
        cosimage[j] = Fa/steps
        sinimage[j] = Fb/steps
        Fa /= sqrt(steps);
        Fb /= sqrt(steps);
    I += Fa**2 + Fb**2
    floatimage[j] += I/steps*omega_Rsqr_pixel*fluence*r_e_sqr
    j+=1

# write floatimage into 'floatimage.bin'                #1193
# convert and dump into file 'intimage.img'             #1211
intimage = float2gray(floatimage)
np.savetxt('intimage.img',intimage)
