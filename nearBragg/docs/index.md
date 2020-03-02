# Near bragg

An algorithm computing the interference patterns of an assembly of scatterers based
on their path difference to a detector.

Source file format : x,y,z, O,B,P

## Parameters :

Type                | parameters
----- --------------|-----------
**detector**        |distance
                    |detsize,detpixels,pixel
**beam spectral**   |lambda,dispersion,dispteps
**beam spatial**    |divergence,divsteps,hdivrange,vdivrange,hdivsteps,vdivsteps
                    |source_distance,source_depth,depthstep
**misc**            |oversample,curvedet,point_pixel, Xbeam,Ybeam
**features**        |coherent, nopolar + occ, bfactor, phase-shift through file parsing
**io**              |file,cosfile,sinfile,floatfile,intifle
                    |scale, roi
**print**           |noprogress, printout, accumulate

## Code structure


```
double fluence = 1.25932015286227e+29;//incident flux in particles/m^2 */
double r_e_sqr = 7.94079248018965e-30;//classical radius of the electron
double steps = divsteps*dispsteps*depthsteps*oversample*oversample;//total number of sub-steps to normalize over 
double airpath = sqrt(pixel_X*pixel_X+pixel_Y*pixel_Y+pixel_Z*pixel_Z);//solid angle subtended by a pixel from the origin: (pix/airpath)^2*cos(2theta)
double omega_Rsqr_pixel = pixel*pixel*distance/airpath;//m^2
    DWF = occ[i]*exp(-Bfac[i]*stol**2)
    Fa+=DWF*cos(phase)/source_to_atom_path/atom_to_pixel_path
    I += Fa**2 + Fb**2;// normalized intensity
floatimage[j] += I/steps*omega_Rsqr_pixel*fluence*r_e_sqr;
```

## Monte-Carlo style 

Mutliple and inelastic scattering :

- Stochasticly using the theoretical scattering rate.
- The differential scattering cross section can be used for more accurately account for the scattering factor of individual atoms.
- Inelastic scattering loss of coherency and spectral broadening can be represented through rigorous bookeeping.
