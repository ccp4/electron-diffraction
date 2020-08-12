/* general scattering simulator						-James Holton		3-5-16

example:

gcc -O -O -o nearBragg nearBragg.c -lm
./nearBragg -file lattice_points.txt -lambda 1 -dispersion 0.1 -dispstep 3 -distance 100  -detsize 100 pixel 0.1 \
  -hdiv 0.28 -hdivstep 0.02 -vdiv 0.28 -vdivstep 0.02  -source_distance 10000 

lattice positions and wavelength (lambda) should be provided in Angstrom, three numbers per line
detector distance, detsize and pixel size in mm
divergence in mrad
dispersion in percent

 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.141592653589793
#define twoPI 6.283185307179586
double RTD = 180.0/PI;

/* rotate a 3-vector in space applied in order phix,phiy,phiz*/
double *rotate(double *v, double *new, double phix, double phiy, double phiz);
/* rotate a 3-vector about a unit vector axis */
double *rotate_axis(double *v, double *new, double *axis, double phi);

/* generate unit vector in random direction */
float uniform3Ddev(float *dx, float *dy, float *dz, long *idum);
/* random deviate with Poisson distribution */
float poidev(float xm, long *idum);
/* random deviate with Gaussian distribution */
float gaussdev(long *idum);
/* random deviate with Lorentzian distribution */
float lorentzdev(long *idum);
/* random deviate with triangle-shaped distribution */
float triangledev(long *idum);
/* random deviate with exponential distribution (>0) */
float expdev(long *idum);
/* random deviate with uniform distribution */
float ran1(long *idum);

char *infilename;
FILE *infile = NULL;
char line[1024];
char *token;
const char delimiters[] = " \t\r,;:!";
const char numberstuf[] = "0123456789-+.EGeg";
char *floatfilename = "floatimage.bin\0";
char *sinfilename = "sinimage.bin\0";
char *cosfilename = "cosimage.bin\0";
char *intfilename = "intimage.img\0";
char *Ifilename = "I.txt\0";

FILE *outfile = NULL;

int main(int argc, char** argv)
{
    int printout = 0;
    int printout_ypixel,printout_xpixel=-1;
    int accumulate = 0;
    double r;

    /* Thomson cross section */
    double r_e_sqr = 7.94079248018965e-30;
    /* incident x-ray fluence in photons/m^2 */
    double fluence = 1.25932015286227e+29;
    /* polarization factor */
    double polar,costwotheta;

    double test;


    double distance = 100.0e-3;
    double detsize_x = 102.4e-3;
    double detsize_y = 102.4e-3;
    double fdet_vector[4]  = {0,0,0,1};
    double sdet_vector[4]  = {0,0,-1,0};
    double odet_vector[4]  = {0,1,0,0};
    double pix0_vector[4]  = {0,0,0,0};
    double beam_vector[4]  = {0,1,0,0};
    double polar_vector[4] = {0,0,0,1};
    double spindle_vector[4] = {0,0,0,1};
    int curved_detector = 0;
    double pixel = 0.1e-3;
    double subpixel;
    double Xdet,Ydet,Xbeam=-1e99,Ybeam=-1e99,Rdet;
    int xpixel,ypixel,xpixels=0,ypixels=0,pixels;
    long progress_pixel,progress_pixels;
    int progress_meter=1;
    int debug=0;
    int oversample = -1,recommended_oversample,suby,subz;
    double xtalsize_max,xtalsize_a,xtalsize_b,xtalsize_c,reciprocal_pixel_size;
    int roi_xmin=-1,roi_xmax=-1,roi_ymin=-1,roi_ymax=-1;
    double airpath,omega_pixel,omega_Rsqr_pixel,omega_sum;
    int point_pixel = 0;
    int nopolar = 0;

    double stol,twotheta,theta;
    double detector_twotheta = -1e99;

    double lambda,dispersion=0.0,dispstep=-1,lambda0 = 1.0e-10;
    double phi0=0.0;
    double hdiv,hdivstep=-1.0,hdivrange= -1.0;
    double vdiv,vdivstep=-1.0,vdivrange= -1.0;
    int    round_div = 1;
    double ddepth,depthstep=-1.0,depthrange= 0;
    double source_distance = 10.0,source_depth=0.0,depth_step=-1;
    int far_source = 0;
        
    int n,atoms=1000000;
    int hdiv_tic,vdiv_tic,disp_tic,source_tic;
    int steps,divsteps=-1,hdivsteps=-1,vdivsteps=-1,dispsteps=-1,depthsteps=-1;
    int i,j,k;
    double *atomX,*atomY,*atomZ,*occ,*Bfac,*phsft;
    float *floatimage;
    float *sinimage;
    float *cosimage;
    unsigned short int *intimage;
    double X,Y,Z,O,B,P,DWF;
    double maxX=-1e99,maxY=-1e99,maxZ=-1e99,minX=1e99,minY=1e99,minZ=1e99;
    double ShannX,ShannY;
    double S_X,S_Y,S_Z,S0_X,S0_Y,S0_Z,s_x,s_y,s_z;

        
    /* x-ray endpoints */
    double vector[4];
    double newvector[4];
    double pixel_X,pixel_Y,pixel_Z;
    double source_X,source_Y,source_Z;

    double source_to_atom_path,atom_to_pixel_path,phase,Fa,Fb,I;
    int coherent = 0;
    double max_I = 0.0;
    double max_I_x = 0.0,max_I_y = 0.0;
    double intfile_scale = 0.0;
    double pgm_scale = 0.0;
    double sum,sumsqr,avg,rms,rmsd;
    int sumn=0;
    int overloads = 0;        
        
    long seed;
        
//    seed = -time((time_t *)0);
//    printf("GOTHERE seed = %u\n",seed);
      
        
    /* check argument list */
    for(i=1; i<argc; ++i)
    {
        if(argv[i][0] == '-')
        {
            /* option specified */
            if(strstr(argv[i], "-atoms") && (argc > (i+1)))
            {
                atoms = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-Xbeam") && (argc > (i+1)))
            {
                Xbeam = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-Ybeam") && (argc > (i+1)))
            {
                Ybeam = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-fdet_vector") && (argc > (i+3)))
            {
                fdet_vector[1] = atof(argv[i+1]);
                fdet_vector[2] = atof(argv[i+2]);
                fdet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-sdet_vector") && (argc > (i+3)))
            {
                sdet_vector[1] = atof(argv[i+1]);
                sdet_vector[2] = atof(argv[i+2]);
                sdet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-odet_vector") && (argc > (i+3)))
            {
                odet_vector[1] = atof(argv[i+1]);
                odet_vector[2] = atof(argv[i+2]);
                odet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-beam_vector") && (argc > (i+3)))
            {
                beam_vector[1] = atof(argv[i+1]);
                beam_vector[2] = atof(argv[i+2]);
                beam_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-polar_vector") && (argc > (i+3)))
            {
                polar_vector[1] = atof(argv[i+1]);
                polar_vector[2] = atof(argv[i+2]);
                polar_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-spindle_axis") && (argc > (i+3)))
            {
                spindle_vector[1] = atof(argv[i+1]);
                spindle_vector[2] = atof(argv[i+2]);
                spindle_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-pix0_vector") && (argc > (i+3)))
            {
                pix0_vector[0] = 1.0;
                pix0_vector[1] = atof(argv[i+1]);
                pix0_vector[2] = atof(argv[i+2]);
                pix0_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-divergence") && (argc > (i+1)))
            {
                hdivrange = vdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivrange") && (argc > (i+1)))
            {
                hdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-vdivrange") && (argc > (i+1)))
            {
                vdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivstep") && (strlen(argv[i]) == 9) && (argc > (i+1)))
            {
                hdivstep = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivsteps") && (argc > (i+1)))
            {
                hdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-vdivstep") && (strlen(argv[i]) == 9) && (argc > (i+1)))
            {
                vdivstep = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-vdivsteps") && (argc > (i+1)))
            {
                vdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-divsteps") && (argc > (i+1)))
            {
                hdivsteps = vdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-distance") && (argc > (i+1)))
            {
                distance = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-twotheta") && (argc > (i+1)))
            {
                detector_twotheta = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-source_dist") && (argc > (i+1)))
            {
		if(strstr(argv[i+1], "far")) {
		    far_source = 1;
		    source_distance = 0.0;
		    source_depth = 0.0;
		}else{
                    source_distance = atof(argv[i+1])/1000.0;
		}
            }
            if(strstr(argv[i], "-source_depth") && (argc > (i+1)))
            {
                 source_depth = atof(argv[i+1])*1e-10;
            }
            if(strstr(argv[i], "-depthsteps") && (argc > (i+1)))
            {
                depthsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detsize") && (strlen(argv[i]) == 8) && (argc > (i+1)))
            {
                detsize_x = atof(argv[i+1])/1000.0;
                detsize_y = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detsize_x") && (argc > (i+1)))
            {
                detsize_x = atof(argv[i+1])/1000.0;
            }
             if(strstr(argv[i], "-detsize_y") && (argc > (i+1)))
            {
                detsize_y = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detpixels") && (strlen(argv[i]) == 10) && (argc > (i+1)))
            {
                xpixels = ypixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detpixels_x") && (argc > (i+1)))
            {
                xpixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detpixels_y") && (argc > (i+1)))
            {
                ypixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-curved_det") && (argc > (i+1)))
            {
                curved_detector = 1;
            }
            if(strstr(argv[i], "-pixel") && (argc > (i+1)))
            {
                pixel = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-point_pixel") )
            {
                point_pixel = 1;
            }
            if(strstr(argv[i], "-nopolar") )
            {
                nopolar = 1;
            }
            if(strstr(argv[i], "-oversample") && (argc > (i+1)))
            {
                oversample = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-roi") && (argc > (i+4)))
            {
                roi_xmin = atoi(argv[i+1]);
                roi_xmax = atoi(argv[i+2]);
                roi_ymin = atoi(argv[i+3]);
                roi_ymax = atoi(argv[i+4]);
            }
            if((strstr(argv[i], "-lambda") || strstr(argv[i], "-wave")) && (argc > (i+1)))
            {
                lambda0 = atof(argv[i+1])/1.0e10;
            }
           if(strstr(argv[i], "-energy") && (argc > (i+1)))
            {
                lambda0 = (12398.42/atof(argv[i+1]))/1.0e10;
            }
            if(strstr(argv[i], "-fluence") && (argc > (i+1)))
            {
                fluence = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-dispersion") && (argc > (i+1)))
            {
                dispersion = atof(argv[i+1])/100.0;
            }
            if(strstr(argv[i], "-dispsteps") && (argc > (i+1)))
            {
                dispsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-phi") && strlen(argv[i])==4 && (argc > (i+1)))
            {
                phi0 = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-file") && (argc > (i+1)))
            {
                infilename = argv[i+1];
		infile = fopen(infilename,"r");
            }
            if(strstr(argv[i], "-floatfile") && (argc > (i+1)))
            {
                floatfilename = argv[i+1];
            }
            if(strstr(argv[i], "-sinfile") && (argc > (i+1)))
            {
                sinfilename = argv[i+1];
            }
            if(strstr(argv[i], "-cosfile") && (argc > (i+1)))
            {
                cosfilename = argv[i+1];
            }
            if(strstr(argv[i], "-intfile") && (argc > (i+1)))
            {
                intfilename = argv[i+1];
            }
            if(strstr(argv[i], "-Ifilename") && (argc > (i+1)))
            {
                Ifilename = argv[i+1];
            }
            if(strstr(argv[i], "-scale") && (argc > (i+1)))
            {
                intfile_scale = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-coherent") )
            {
		/* turn off incoherent addition */
                coherent = 1;
            }
            if(strstr(argv[i], "-printout") )
            {
		/* turn on console printing */
                printout = 1;
            }
            if(strstr(argv[i], "-noprogress") )
            {
		/* turn off progress meter */
                progress_meter = 0;
            }
            if(strstr(argv[i], "-progress") )
            {
		/* turn on progress meter */
                progress_meter = 1;
            }
            if(strstr(argv[i], "-debug") )
            {
                /* turn on debug messages */
                debug = 1;
            }
            if(strstr(argv[i], "-printout_pixel") && (argc > (i+2)))
            {
                printout_xpixel = atoi(argv[i+1]);
                printout_ypixel = atoi(argv[i+2]);
            }
            if(strstr(argv[i], "-accumulate") )
            {
		/* accumulate from old float-file */
                accumulate = 1;
            }
        }
    }

    printf("nearBragg assumption-free diffraction simulator - James Holton 3-19-15\n");

    if(infile == NULL){
	printf("usage: nearBragg -file atom_list.txt\n");
	printf("options:\n");\
	printf("\t-file filename.txt  \ttext file containing point scatterer coordinates in Angstrom relative to the origin.  The x axis is the x-ray beam and Y and Z are parallel to the detector Y and X coordinates, respectively\n");
        printf("\t-distance        \tdistance from origin to detector center in mm\n");
        printf("\t-detsize         \tdetector size in mm.  may also use -detsize_x -detsize_y\n");
        printf("\t-detpixels       \tdetector size in pixels.  may also use -detpixels_x -detpixels_y\n");
        printf("\t-pixel           \tdetector pixel size in mm.\n");
        printf("\t-Xbeam           \timage X coordinate of direct-beam spot (mm). (default: center)\n");
        printf("\t-Ybeam           \timage Y coordinate of direct-beam spot (mm). )default: center)\n");
        printf("\t-twotheta       \trotation of detector about spindle axis (deg). )default: 0)\n");
        printf("\t-lambda          \tincident x-ray wavelength in Angstrom. may also use -energy in eV\n");
        printf("\t-dispersion      \tspectral dispersion: delta-lambda/lambda in percent\n");
        printf("\t-dispsteps       \tnumber of wavelengths in above range\n");
        printf("\t-divergence      \thorizontal and vertical angular spread of source points in mrad\n");
        printf("\t-divsteps        \tnumber of points across source\n");
        printf("\t-square_div      \tfull divergence grid (default: round off corners)\n");
        printf("\t-hdivrange       \thorizontal angular spread of source points in mrad\n");
        printf("\t-vdivrange       \tvertical angular spread of source points in mrad\n");
        printf("\t-hdivsteps       \tnumber of source points in the horizontal\n");
        printf("\t-vdivsteps       \tnumber of source points in the vertical\n");
        printf("\t-detsize_x          \tdetector size in x direction (mm)\n");
        printf("\t-detsize_y          \tdetector size in y direction (mm)\n");
	printf("\t-detpixels_x        \tdetector size in x direction (pixels)\n");
	printf("\t-detpixels_y        \tdetector size in y direction (pixels)\n");
        printf("\t-point_pixel        \tturn off solid-angle correction for square flat pixels\n");
	printf("\t-curved_det         \tall detector pixels same distance from sample (origin)\n");
        printf("\t-nopolar            \tturn off the polarization correction\n");
        printf("\t-oversample         \tnumber of sub-pixels per pixel. Default: 1\n");

	printf("\t-roi xmin xmax ymin ymax\tonly render pixels within a set range. Default: all detector\n");
        printf("\t-floatfile          \tname of binary pixel intensity output file (4-byte floats)\n");
	printf("\t-sinfile            \tname of binary imaginary structure factor output file (4-byte floats)\n");
	printf("\t-cosfile            \tname of binary real structure factor output file (4-byte floats)\n");
        printf("\t-intfile            \tname of smv-formatted output file.\n");
	printf("\t-scale              \tscale factor for intfile. Default: fill dynamic range\n");
        printf("\t-coherent           \tcoherently add everything, even different wavelengths. Not the default\n");
	printf("\t-accumulate         \timport contents of floatfile or sinfile/cosfile and add to them. Not the default\n");

        printf("\t-noprogress         \tturn off the progress meter\n");
        printf("\t-printout           \tprint pixel values out to the screen\n");

        printf("\t-source_distance    \tdistance from origin to source (mm) (can be \"far\" to save time)\n");
	printf("\t-source_depth       \tdistance from front of source to back (mm)\n");
	printf("\t-depthstep          \tnumber of source points along path to sample. Default: 1\n");
exit(9);
    }

    /* how many atoms do we need? */
    atoms = 0;
    while ( fgets ( line, sizeof line, infile ) != NULL ) {
	token = line;
	token += strspn(token,delimiters);
	if(strcmp(token,"\n")==0) {
	    if(debug) printf("blank\n");
	    continue;
	}
	++atoms;
    }
    printf("counted %d atoms in %s\n",atoms,infilename);
    rewind(infile);

    /* allocate memory */
    atomX = calloc(atoms+10,sizeof(double));
    atomY = calloc(atoms+10,sizeof(double));
    atomZ = calloc(atoms+10,sizeof(double));
    occ   = calloc(atoms+10,sizeof(double));
    Bfac  = calloc(atoms+10,sizeof(double));
    phsft = calloc(atoms+10,sizeof(double));

    if(xpixels) {
	detsize_x = pixel*xpixels;
    }
    if(ypixels) {
	detsize_y = pixel*ypixels;
    }
    xpixels = ceil(detsize_x/pixel);
    ypixels = ceil(detsize_y/pixel);
    pixels = xpixels*ypixels;
    floatimage = calloc(pixels+10,sizeof(float));
    sinimage = calloc(pixels+10,2*sizeof(float));
    cosimage = calloc(pixels+10,2*sizeof(float));
    intimage = calloc(pixels+10,sizeof(unsigned short int));

    /* defaults? */
    if(Xbeam <= -1e99) Xbeam = (detsize_x - pixel)/2.0;
    if(Ybeam <= -1e99) Ybeam = (detsize_y + pixel)/2.0;
    if(detector_twotheta > -1e99) {
        rotate_axis(fdet_vector,newvector,spindle_vector,detector_twotheta);
        fdet_vector[1] = newvector[1];
        fdet_vector[2] = newvector[2];
        fdet_vector[3] = newvector[3];
        rotate_axis(sdet_vector,newvector,spindle_vector,detector_twotheta);
        sdet_vector[1] = newvector[1];
        sdet_vector[2] = newvector[2];
        sdet_vector[3] = newvector[3];
        rotate_axis(odet_vector,newvector,spindle_vector,detector_twotheta);
        odet_vector[1] = newvector[1];
        odet_vector[2] = newvector[2];
        odet_vector[3] = newvector[3];
    }
    if(pix0_vector[0] == 0.0) {
        pix0_vector[1] = -Xbeam*fdet_vector[1]-Ybeam*sdet_vector[1]+distance*odet_vector[1];
        pix0_vector[2] = -Xbeam*fdet_vector[2]-Ybeam*sdet_vector[2]+distance*odet_vector[2];
        pix0_vector[3] = -Xbeam*fdet_vector[3]-Ybeam*sdet_vector[3]+distance*odet_vector[3];
    }
    if(detector_twotheta <= -1e99) detector_twotheta = 0;
    if(roi_xmin < 0) roi_xmin = 0;
    if(roi_xmax < 0) roi_xmax = xpixels;
    if(roi_ymin < 0) roi_ymin = 0;
    if(roi_ymax < 0) roi_ymax = ypixels;
    progress_pixels = (roi_xmax-roi_xmin+1)*(roi_ymax-roi_ymin+1);


    if(hdivsteps <= 0){
        /* auto-select number of steps */
	if(hdivrange < 0.0) {
            /* auto-select range */
            if(hdivstep <= 0.0) {
	        /* user doesn't care about anything */
	        hdivsteps = 1;
		hdivrange = 0.0;
		hdivstep = 0.0;
	    } else {
		/* user specified stepsize and nothing else */
		hdivrange = hdivstep;
		hdivsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(hdivstep <= 0.0) {
	        /* range specified, but nothing else */
                hdivstep = hdivrange;
                hdivsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                hdivsteps = ceil(hdivrange/hdivstep);
	    }
	}
    } else {
	/* user-specified number of steps */
	if(hdivrange < 0.0) {
            /* auto-select range */
            if(hdivstep <= 0.0) {
	        /* user cares only about number of steps */
		hdivrange = 1.0;
		hdivstep = hdivrange/hdivsteps;
	    } else {
		/* user doesn't care about range */
		hdivrange = hdivstep;
		hdivsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(hdivstep <= 0.0) {
	        /* range and steps specified */
		if(hdivsteps <=1 ) hdivsteps = 2;
                hdivstep = hdivrange/(hdivsteps-1);
            } else {
                /* everything specified */
	    }
	}
    }

    if(vdivsteps <= 0){
        /* auto-select number of steps */
	if(vdivrange < 0.0) {
            /* auto-select range */
            if(vdivstep <= 0.0) {
	        /* user doesn't care about anything */
	        vdivsteps = 1;
		vdivrange = 0.0;
		vdivstep = 0.0;
	    } else {
		/* user specified stepsize and nothing else */
		vdivrange = vdivstep;
		vdivsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(vdivstep <= 0.0) {
	        /* range specified, but nothing else */
                vdivstep = vdivrange;
                vdivsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                vdivsteps = ceil(vdivrange/vdivstep);
	    }
	}
    } else {
	/* user-specified number of steps */
	if(vdivrange < 0.0) {
            /* auto-select range */
            if(vdivstep <= 0.0) {
	        /* user cares only about number of steps */
		vdivrange = 1.0;
		vdivstep = vdivrange/vdivsteps;
	    } else {
		/* user doesn't care about range */
		vdivrange = vdivstep;
		vdivsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(vdivstep <= 0.0) {
	        /* range and steps specified */
		if(vdivsteps <=1 ) vdivsteps = 2;
                vdivstep = vdivrange/(vdivsteps-1);
            } else {
                /* everything specified */
	    }
	}
    }
    

    if(dispsteps <= 0){
        /* auto-select number of steps */
	if(dispersion < 0.0) {
            /* auto-select range */
            if(dispstep <= 0.0) {
	        /* user doesn't care about anything */
	        dispsteps = 1;
		dispersion = 0.0;
		dispstep = 0.0;
	    } else {
		/* user specified stepsize and nothing else */
		dispersion = dispstep;
		dispsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(dispstep <= 0.0) {
	        /* range specified, but nothing else */
                dispstep = dispersion;
                dispsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                dispsteps = ceil(dispersion/dispstep);
	    }
	}
    } else {
	/* user-specified number of steps */
	if(dispersion < 0.0) {
            /* auto-select range */
            if(dispstep <= 0.0) {
	        /* user cares only about number of steps */
		dispersion = 1.0;
		dispstep = dispersion/dispsteps;
	    } else {
		/* user doesn't care about range */
		dispersion = dispstep;
		dispsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(dispstep <= 0.0) {
	        /* range and steps specified */
		if(dispsteps <=1 ) dispsteps = 2;
                dispstep = dispersion/(dispsteps-1);
            } else {
                /* everything specified */
	    }
	}
    }
        

    if(depthsteps <= 0){
        /* auto-select number of steps */
	if(source_depth < 0.0) {
            /* auto-select range */
            if(depthstep <= 0.0) {
	        /* user doesn't care about anything */
	        depthsteps = 1;
		source_depth = 0.0;
		depthstep = 0.0;
	    } else {
		/* user specified stepsize and nothing else */
		source_depth = depthstep;
		depthsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(depthstep <= 0.0) {
	        /* range specified, but nothing else */
                depthstep = source_depth;
                depthsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                depthsteps = ceil(source_depth/depthstep);
	    }
	}
    } else {
	/* user-specified number of steps */
	if(source_depth < 0.0) {
            /* auto-select range */
            if(depthstep <= 0.0) {
	        /* user cares only about number of steps */
		source_depth = 1.0;
		depthstep = source_depth/depthsteps;
	    } else {
		/* user doesn't care about range */
		source_depth = depthstep;
		depthsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(depthstep <= 0.0) {
	        /* range and steps specified */
		if(depthsteps <=1 ) depthsteps = 2;
                depthstep = source_depth/(depthsteps-1);
            } else {
                /* everything specified */
	    }
	}
    }
    

    /* sanity checks */
    if(hdivrange <= 0.0 || hdivstep <= 0.0 || hdivsteps <= 0) {
        hdivsteps = 1;
        hdivrange = 0.0;
        hdivstep = 0.0;
    }
    if(vdivrange <= 0.0 || vdivstep <= 0.0 || vdivsteps <= 0) {
        vdivsteps = 1;
        vdivrange = 0.0;
        vdivstep = 0.0;
    }
    if(dispersion <= 0.0 || dispstep <= 0.0 || dispsteps <= 0) {
        dispsteps = 1;
        dispersion = 0.0;
        dispstep = 0.0;
    }
    if(source_depth <= 0.0 || depthstep <= 0.0 || depthsteps <= 0) {
        depthsteps = 1;
        source_depth = 0.0;
        depthstep = 0.0;
    }

  
    /* load the lattice points */
    n = 0;
    printf("reading %s\n",infilename);
    while ( fgets ( line, sizeof line, infile ) != NULL ) { /* read a line */

       // initialize

       X=Y=Z=0.0;
       O=1.0;
       B=0.0;
       P=-90;

	token = line;
	token += strspn(token,delimiters);
	if(strcmp(token,"\n")==0) {
	    if(debug) printf("blank\n");
	    continue;
	}
	++atoms;

	while (strcmp(token,"\n")!=0) 
        {
	    X=atof(token);

	    token += strspn(token,numberstuf);
	    if (strcmp(token,"\n")==0) continue;
	    token += strcspn(token,delimiters);
	    token += strspn(token,delimiters);
	    if (strcmp(token,"\n")==0) continue;

	    Y=atof(token);

	    token += strspn(token,numberstuf);
	    if (strcmp(token,"\n")==0) continue;
	    token += strcspn(token,delimiters);
	    token += strspn(token,delimiters);
	    if (strcmp(token,"\n")==0) continue;

	    Z=atof(token);

	    token += strspn(token,numberstuf);
	    if (strcmp(token,"\n")==0) continue;
	    token += strcspn(token,delimiters);
	    token += strspn(token,delimiters);
	    if (strcmp(token,"\n")==0) continue;

	    O=atof(token);

	    token += strspn(token,numberstuf);
	    if (strcmp(token,"\n")==0) continue;
	    token += strcspn(token,delimiters);
	    token += strspn(token,delimiters);
	    if (strcmp(token,"\n")==0) continue;

	    B=atof(token);

	    token += strspn(token,numberstuf);
	    if (strcmp(token,"\n")==0) continue;
	    token += strcspn(token,delimiters);
	    token += strspn(token,delimiters);
	    if (strcmp(token,"\n")==0) continue;

	    P=atof(token);

	    break;
	}
		
	if(debug) printf("initializing: %d %lg %lg %lg with Z=%lg B=%lg\n",n,X,Y,Z,O,B);
	atomX[n] = X/1e10; atomY[n] = Y/1e10; atomZ[n] = Z/1e10;
	occ[n]=O ; Bfac[n]=B; phsft[n]=P/RTD;
	//printf("initialized: %d at %lg %lg %lg with Z=%lg B=%lg\n",n,X,Y,Z,O,B);

	/* rough estimate of object size? */
	if(X>maxX)maxX=X;
	if(Y>maxY)maxY=Y;
	if(Z>maxZ)maxZ=Z;
	if(X<minX)minX=X;
	if(Y<minY)minY=Y;
	if(Z<minZ)minZ=Z;

        ++n;
	if(n > atoms) {
	    printf("WARNING: more than %d atoms in %s\n",atoms,infilename);
	    break;
	}
    }
    fclose(infile);

    /* check size of object vs pixel size */
    xtalsize_a = ((maxX-minX)*1e-10);
    xtalsize_b = ((maxY-minY)*1e-10);
    xtalsize_c = ((maxZ-minZ)*1e-10);
    xtalsize_max = xtalsize_a;
    if(xtalsize_max < xtalsize_b) xtalsize_max = xtalsize_b;
    if(xtalsize_max < xtalsize_c) xtalsize_max = xtalsize_c;
    reciprocal_pixel_size = lambda0*distance/pixel;
    recommended_oversample = ceil(3.0 * xtalsize_max/reciprocal_pixel_size);
    if(recommended_oversample <= 0) recommended_oversample = 1;
    if(oversample <= 0) {
	oversample = recommended_oversample;
	printf("auto-selected %d-fold oversampling\n",oversample);
    }
    if(oversample < recommended_oversample)
    {
	printf("WARNING: maximum dimension of sample is %g A\n",xtalsize_max*1e10);
	printf("         but reciprocal pixel size is %g A\n", reciprocal_pixel_size*1e10 );
        printf("         intensity may vary significantly across a pixel!\n");
	printf("         recommend -oversample %d to work around this\n",recommended_oversample);
    }


    /* load old float-image? */
    if(accumulate && ! coherent)
    {
	infile = fopen(floatfilename,"r");
 	if(infile == NULL)
	{
	    printf("WARNING: no %s float-file to read!\n",floatfilename);
	    printf("setting initial image to all zeroes...\n");
	}
	else
	{
            printf("importing intensity: %s\n",floatfilename);
	    if(pixels != fread(floatimage,sizeof(float),pixels,infile) ) {perror("ERROR: fread() did not read all pixels");};
	    fclose(infile);
	}
    }
    if(accumulate && coherent)
    {
	infile = fopen(sinfilename,"r");
 	if(infile == NULL)
	{
	    printf("WARNING: no %s sin-file to read!\n",sinfilename);
	    printf("setting initial Im(F) to all zeroes...\n");
	}
	else
	{
            printf("importing sin: %s\n",sinfilename);
	    if(pixels != fread(sinimage,sizeof(float),pixels,infile) ) {perror("ERROR: fread() did not read all pixels");};
	    fclose(infile);
	}
	infile = fopen(cosfilename,"r");
 	if(infile == NULL)
	{
	    printf("WARNING: no %s cos-file to read!\n",cosfilename);
	    printf("setting initial Re(F) to all zeroes...\n");
	}
	else
	{
            printf("importing cos: %s\n",cosfilename);
	    if(pixels != fread(cosimage,sizeof(float),pixels,infile) ) {perror("ERROR: fread() did not read all pixels");};
	    fclose(infile);
	}
    }

    /* count divsteps sweep over solid angle of beam divergence */
    divsteps = 0;
    for(hdiv_tic=0;hdiv_tic<hdivsteps;++hdiv_tic){
        for(vdiv_tic=0;vdiv_tic<vdivsteps;++vdiv_tic){                
	    hdiv = hdivstep * hdiv_tic - hdivrange/2.0 ;
	    vdiv = vdivstep * vdiv_tic - vdivrange/2.0 ;
	    /* force an elliptical divergence */
	    test = (hdiv*hdiv-hdivstep*hdivstep/4.0*(1-hdivsteps%2))/hdivrange/hdivrange ;
	    test += (vdiv*vdiv-vdivstep*vdivstep/4.0*(1-vdivsteps%2))/vdivrange/vdivrange ;
	    if( round_div && test*4.0 > 1.1) continue;
                
	    ++divsteps;
            printf("divergence deviation: %g %g\n",hdiv,vdiv);
        }
    }

    /* print out wavelength steps with sweep over spectral dispersion */
    for(disp_tic=0;disp_tic<dispsteps;++disp_tic){
	lambda = lambda0 * ( 1.0 + dispstep * disp_tic - dispersion/2.0 ) ;
	printf("lambda%d = %.15g\n",disp_tic,lambda);
    }

    /* print out source steps with sweep over all source planes */
    for(source_tic=0;source_tic<depthsteps;++source_tic){
	ddepth = depthstep * source_tic - source_depth/2.0 ;
	printf("depth%d = %g A\n",source_tic,ddepth*1e10);
    }


    /* total number of sub-steps to normalize over */
    steps = divsteps*dispsteps*depthsteps*oversample*oversample;
    subpixel = pixel/oversample;

    atoms = n;
    printf("  %d atoms\n",atoms);
    printf("  wave=%g meters +/- %g%% in %d steps\n",lambda0,dispersion*100,dispsteps);
    printf("  distance=%g detsize=%gx%g  pixel=%g meters (%dx%d pixels)\n",distance,detsize_x,detsize_y,pixel,xpixels,ypixels);
    printf("  Xbeam=%g Ybeam=%g\n",Xbeam,Ybeam);
    printf("  roi: %d < z < %d && %d < y < %d\n",roi_xmin,roi_xmax,roi_ymin,roi_ymax);
    printf("  hdivrange=%g hdivstep=%g  radians\n",hdivrange,hdivstep);
    printf("  vdivrange=%g vdivstep=%g  radians\n",vdivrange,vdivstep);
    printf("  %d divergence steps\n",divsteps);
    if(far_source){
	printf("  source is VERY far away\n");
    }else{
        printf("  source at %g meters and %g deep in %d depth steps\n",source_distance,source_depth,depthsteps);
    }
    printf("  %d pixel oversample steps\n",oversample);
    printf("  coherent source: %d\n",coherent);

 


    /* sweep over detector */   
    sum = sumsqr = 0.0;
    j = 0;
    progress_pixel = 0;
    omega_sum = 0.0;
    for(ypixel=0;ypixel<ypixels;++ypixel){
        for(xpixel=0;xpixel<xpixels;++xpixel){

            if(xpixel < roi_xmin || xpixel > roi_xmax || ypixel < roi_ymin || ypixel > roi_ymax) {
	   	++j; continue;
	    }

	    /* reset photon count */
	    I = 0;

	    /* reset coherent addition to input value (or zero if no input value) */
	    /* (this will be re-asserted below for each incoherent beam) */
	    Fa=Fb=0.0;

	    for(suby=0;suby<oversample;++suby){
		for(subz=0;subz<oversample;++subz){

		    Xdet = subpixel*(xpixel*oversample + subz) + subpixel/2.0;
		    Ydet = subpixel*(ypixel*oversample + suby) + subpixel/2.0;
//		    Xdet = pixel*xpixel;
//		    Ydet = pixel*ypixel;

		    /* construct detector pixel position */
//		    pixel_X = distance;
//		    pixel_Y = Ydet-Ybeam;
//		    pixel_Z = Xdet-Xbeam;
                    pixel_X = Xdet*fdet_vector[1]+Ydet*sdet_vector[1]+pix0_vector[1];
                    pixel_Y = Xdet*fdet_vector[2]+Ydet*sdet_vector[2]+pix0_vector[2];
                    pixel_Z = Xdet*fdet_vector[3]+Ydet*sdet_vector[3]+pix0_vector[3];
		    if(curved_detector) {
			/* construct detector pixel that is always "distance" from the origin */
			vector[1] = distance; vector[2]=0 ; vector[3]=0;
			rotate(vector,newvector,0,pixel_Z/distance,pixel_Y/distance);
			pixel_X = newvector[1];
			pixel_Y = newvector[2];
			pixel_Z = newvector[3];
		    }


		    /* solid angle subtended by a pixel from the origin: (pix/airpath)^2*cos(2theta) */
		    airpath = sqrt(pixel_X*pixel_X+pixel_Y*pixel_Y+pixel_Z*pixel_Z);
//		    omega_pixel = pixel*pixel/airpath/airpath*distance/airpath;
		    omega_Rsqr_pixel = pixel*pixel*distance/airpath;
	
		    if(point_pixel) omega_Rsqr_pixel = 1.0;

		    omega_sum += omega_Rsqr_pixel;

		    /* sweep over wavelengths */
    		    for(disp_tic=0;disp_tic<dispsteps;++disp_tic){
    		       lambda = lambda0 * ( 1.0 + dispstep * disp_tic - dispersion/2.0 ) ;

		       /* sweep over solid angle of beam divergence */
	               for(hdiv_tic=0;hdiv_tic<hdivsteps;++hdiv_tic){
                         for(vdiv_tic=0;vdiv_tic<vdivsteps;++vdiv_tic){                
	                    hdiv = hdivstep * hdiv_tic - hdivrange/2.0 ;
	                    vdiv = vdivstep * vdiv_tic - vdivrange/2.0 ;
	                    /* force an elliptical divergence */
	                    test = (hdiv*hdiv-hdivstep*hdivstep/4.0*(1-hdivsteps%2))/hdivrange/hdivrange ;
	                    test += (vdiv*vdiv-vdivstep*vdivstep/4.0*(1-vdivsteps%2))/vdivrange/vdivrange ;
	                    if( round_div && test*4.0 > 1.1) continue;

		                for(source_tic=0;source_tic<depthsteps;++source_tic){
	                            ddepth = depthstep * source_tic - source_depth/2.0 ;

				    /* construct source position (flat, coherently-emitting plane) */
				    if(! far_source){
					source_X = -source_distance+ddepth;
					source_Y = atan(hdiv)*source_distance;
					source_Z = atan(vdiv)*source_distance;
				    }

				    /* construct source position (curved, coherently-emitting sphere) */
				    //vector[1] = -source_distance; vector[2]=0 ; vector[3]=0;
				    //rotate(vector,newvector,0,0,-phi/2-(theta0+dtheta));
				    //source_X = newvector[1];
				    //source_Y = newvector[2];
				    //source_Z = newvector[3];

				    /* consider each atom in the crystal */
				    for(i=0;i<atoms;++i){

					/* distance from source to atom to pixel */
					if(far_source) {
					    source_to_atom_path = atomX[i];
					}else{
					    source_to_atom_path = sqrt((source_X-atomX[i])*(source_X-atomX[i])+(source_Y-atomY[i])*(source_Y-atomY[i])+(source_Z-atomZ[i])*(source_Z-atomZ[i]));
					}
					atom_to_pixel_path  = sqrt((pixel_X-atomX[i])*(pixel_X-atomX[i])+0*(pixel_Y-atomY[i])*(pixel_Y-atomY[i])+(pixel_Z-atomZ[i])*(pixel_Z-atomZ[i]));

					/* calculate how many radians (2*pi*cycles) from source to detector point */
					phase = twoPI*(source_to_atom_path+atom_to_pixel_path)/lambda + phsft[i];
					//printf("%d %g\n",i,phase);
					if(far_source){
					    source_to_atom_path=1.0;
					}

					DWF=occ[i];
					if(Bfac[i] != 0.0) {
					    /* need to know the scattering vector length for this */
					    if(far_source){
						S0_X=1.0;S0_Y=S0_Z=0.0;
					    }else{
						S0_X=(atomX[i]-source_X)/source_to_atom_path;
						S0_Y=(atomY[i]-source_Y)/source_to_atom_path;
						S0_Z=(atomZ[i]-source_Z)/source_to_atom_path;
					    }
					    S_X =(pixel_X-atomX[i])/atom_to_pixel_path;
					    S_Y =(pixel_Y-atomY[i])/atom_to_pixel_path;
					    S_Z =(pixel_Z-atomZ[i])/atom_to_pixel_path;
					    /* B factors are in Angstrom units */
					    s_x=(S_X-S0_X)/(lambda*1e10);
					    s_y=(S_Y-S0_Y)/(lambda*1e10);
					    s_z=(S_Z-S0_Z)/(lambda*1e10);
					    stol = sqrt(s_x*s_x+s_y*s_y+s_z*s_z)/2.0;
					    DWF *= exp(-Bfac[i]*stol*stol);
					}

					/* add up sine waves to form structure factor of lattice */
					/* also account for inverse square law */
					Fa += DWF*cos(phase)/source_to_atom_path/atom_to_pixel_path;
					Fb += DWF*sin(phase)/source_to_atom_path/atom_to_pixel_path;
					if( debug && printout && ((xpixel==printout_xpixel && ypixel==printout_ypixel) || printout_xpixel < 0) )
					{
					  printf("%d  %g %g %g -> %g %g %g %g\n",i,atomX[i],atomY[i],atomZ[i],phase,cos(phase),Fa,Fb);
					}
				    }
				    /* end of atom loop */

				}
				/* end of source loop */

				/* convert amplitudes into intensity (photons per steradian) */
				if( ! coherent ) {
				    I += (Fa*Fa + Fb*Fb);
				    /* store Fa,Fb in coherent image? */
				    /* reset back to zero */
				    Fa=Fb=0.0;
				}
			    }
			    /* end of vdiv loop */
			}
			/* end of hdiv loop */

			if( ! coherent ) {
			    I += (Fa*Fa + Fb*Fb);
			    /* store Fa,Fb in coherent image? */
			    /* reset back to zero */
			    Fa=Fb=0.0;
			}
		    }
		    /* end of lambda loop */
		}
		/* end of subz loop */
	    }
	    /* end of suby loop */

	    /* accumulate the coherent image */
	    if( coherent )
	    {
		/* Fa and Fb will only be non-zero if -coherent was specified */
		/* cosimage and sinimage will only be non-zero if -accumulate was specified */

		/* pretend cosimage and sinimage were accumulated with the same oversampling */
	    	Fa += cosimage[j]*steps;
	    	Fb += sinimage[j]*steps;
		/* save the average complex structure factor */
	    	cosimage[j] = Fa/steps;
	    	sinimage[j] = Fb/steps;

		/* now enforce conservation of energy: 
		   final intensity must be independent of pixel/source steps */
		Fa /= sqrt(steps);
		Fb /= sqrt(steps);
	    }

	    /* convert amplitudes into intensity (photons per steradian) */
	    I += (Fa*Fa + Fb*Fb);

	    /* accumulate the intensity image */
	    /* in units of F^2*omega if fluence was not specified */
	    /* or in photons/pixel if fluence was speficied */
	    floatimage[j] += I ;///steps*omega_Rsqr_pixel*fluence*r_e_sqr;

	    if(floatimage[j] > max_I) {
	        max_I = floatimage[j];
	        max_I_x = Xdet;
	        max_I_y = Ydet;
	    }
	    sum += floatimage[j];
            sumsqr += floatimage[j]*floatimage[j];
            ++sumn;

	    if( printout )
	    {
		if((xpixel==printout_xpixel && ypixel==printout_ypixel) || printout_xpixel < 0)
		{
		    twotheta = atan2(sqrt(pixel_Y*pixel_Y+pixel_Z*pixel_Z),pixel_X);
		    stol = sin(twotheta/2.0)/(lambda0*1e10);
		    printf("%4d %4d : stol = %g\n", xpixel,ypixel,stol);
		    printf(" Fa=%g  Fb=%g   I = %g\n", Fa,Fb,I);
		    printf("cos=%g sin=%g\n", cosimage[j],sinimage[j]);
		    printf("I/sr   %15.10g\n", floatimage[j]/(omega_Rsqr_pixel/airpath/airpath));
		    printf("omega  %15.10g\n", omega_Rsqr_pixel/airpath/airpath);
		    printf("pixel  %15.10g\n", floatimage[j]);
		}
	    }
	    else
	    {
		if(progress_meter && progress_pixels/100 > 0){
		    if(progress_pixel % ( progress_pixels/20 ) == 0 || ((10*progress_pixel<progress_pixels || 10*progress_pixel>9*progress_pixels) && (progress_pixel % (progress_pixels/100) == 0))){
			printf("%ld%% done\n",progress_pixel*100/progress_pixels);
		    }
		}
	    }
	    ++j;
	    ++progress_pixel;
	}
	/* end of X pixel loop */
    }
    /* end of Y pixel loop */
    printf("\n");

    printf("solid angle subtended by detector = %g steradian ( %g%% sphere)\n",omega_sum,100*omega_sum/4/PI);

    /* do some stats? */
    if(sumn<=0) sumn=1;
    avg = sum/sumn;
    if(sumn<=1) sumn=2;
    rms = sqrt(sumsqr/(sumn-1));

    printf("writing %s as %d %lu-byte floats\n",floatfilename,pixels,sizeof(float));
    outfile = fopen(floatfilename,"wb");
    fwrite(floatimage,sizeof(float),pixels,outfile);
    fclose(outfile);

    FILE* fu = fopen(Ifilename,"wb");
    j = 0;
    for(ypixel=0;ypixel<ypixels;++ypixel){
      for(xpixel=0;xpixel<xpixels;++xpixel){
        if(xpixel < roi_xmin || xpixel > roi_xmax || ypixel < roi_ymin || ypixel > roi_ymax) {
	  ++j; continue;
	}
	fprintf(fu,"%g ", floatimage[j]);
	++j;
      }
      fprintf(fu,"\n");
    }    
    fclose(fu);
    
    /* no point in writing these if they are all zero? */
    if( coherent )
    {
	printf("writing %s as %d %lu-byte floats\n",sinfilename,pixels,sizeof(float));
	outfile = fopen(sinfilename,"wb");
	fwrite(sinimage,sizeof(float),pixels,outfile);
	fclose(outfile);
	printf("writing %s as %d %lu-byte floats\n",cosfilename,pixels,sizeof(float));
	outfile = fopen(cosfilename,"wb");
	fwrite(cosimage,sizeof(float),pixels,outfile);
	fclose(outfile);
    }

    /* output as ints */   
    j = 0;
    printf("max_I = %g  at %g %g\n",max_I,max_I_x,max_I_y);
    printf("mean= %g rms= %g\n",avg,rms);
    if(intfile_scale <= 0.0){
	intfile_scale = 1.0;
	if(max_I > 0.0) intfile_scale = 55000.0/max_I;
    }
    printf("intfile_scale = %g\n",intfile_scale);
    for(ypixel=0;ypixel<ypixels;++ypixel){
      for(xpixel=0;xpixel<xpixels;++xpixel){
        if(xpixel < roi_xmin || xpixel > roi_xmax || ypixel < roi_ymin || ypixel > roi_ymax) {
	   ++j; continue;
        }
	test = floatimage[j] *intfile_scale+40.0;
	if(test > 65535.0) test = 65535.0;
	intimage[j] = (unsigned short int) ( test );
//	printf("%d %d = %d\n",xpixel,ypixel,intimage[j]);
	++j;
      }
    }

    printf("writing %s as %lu-byte integers\n",intfilename,sizeof(unsigned short int));
    outfile = fopen(intfilename,"wb");
    fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=little_endian;\nTYPE=unsigned_short;\n");
    fprintf(outfile,"SIZE1=%d;\nSIZE2=%d;\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",xpixels,ypixels,pixel*1000.0,distance*1000.0);
    fprintf(outfile,"WAVELENGTH=%g;\nBEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",lambda0*1e10,Xbeam*1000.0,(detsize_y-Ybeam)*1000);
    fprintf(outfile,"PHI=%g;\nOSC_START=%g;\nOSC_RANGE=0;\n",phi0*RTD,phi0*RTD);
    fprintf(outfile,"DETECTOR_SN=000;\n");
    fprintf(outfile,"BEAMLINE=fake;\n");
    fprintf(outfile,"}\f");
    while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    fwrite(intimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);


    return 0;
}


double *rotate(double *v, double *new, double phix, double phiy, double phiz) {
    
    double rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz;
    double new_x,new_y,new_z,rotated_x,rotated_y,rotated_z;
    
    new_x=v[1];
    new_y=v[2];
    new_z=v[3];
    
    if(phix != 0){
        /* rotate around x axis */
        //rxx= 1;         rxy= 0;         rxz= 0;
        ryx= 0;         ryy= cos(phix); ryz=-sin(phix);
        rzx= 0;         rzy= sin(phix); rzz= cos(phix);
        
        rotated_x = new_x;
        rotated_y = new_y*ryy + new_z*ryz;
        rotated_z = new_y*rzy + new_z*rzz;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }
    
    if(phiy != 0) {
        /* rotate around y axis */
        rxx= cos(phiy); rxy= 0;         rxz= sin(phiy);
        //ryx= 0;         ryy= 1;         ryz= 0;
        rzx=-sin(phiy); rzy= 0;         rzz= cos(phiy);
        
        rotated_x = new_x*rxx + new_y*rxy + new_z*rxz;
        rotated_y = new_y;
        rotated_z = new_x*rzx + new_y*rzy + new_z*rzz;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }
    
    if(phiz != 0){
        /* rotate around z axis */
        rxx= cos(phiz); rxy=-sin(phiz); rxz= 0;
        ryx= sin(phiz); ryy= cos(phiz); ryz= 0;
        //rzx= 0;         rzy= 0;         rzz= 1;
        
        rotated_x = new_x*rxx + new_y*rxy ;
        rotated_y = new_x*ryx + new_y*ryy;
        rotated_z = new_z;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }
    
    new[1]=new_x;
    new[2]=new_y;
    new[3]=new_z;
    
    return new;
}



/* rotate a point about a unit vector axis */
double *rotate_axis(double *v, double *new, double *axis, double phi) {

    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double dot = (axis[1]*v[1]+axis[2]*v[2]+axis[3]*v[3])*(1.0-cosphi);
    double temp[4];

    temp[1] = axis[1]*dot+v[1]*cosphi+(-axis[3]*v[2]+axis[2]*v[3])*sinphi;
    temp[2] = axis[2]*dot+v[2]*cosphi+(+axis[3]*v[1]-axis[1]*v[3])*sinphi;
    temp[3] = axis[3]*dot+v[3]*cosphi+(-axis[2]*v[1]+axis[1]*v[2])*sinphi;
    new[1]=temp[1]; new[2]=temp[2]; new[3]=temp[3];

    return new;
}


