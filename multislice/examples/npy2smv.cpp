//g++ -o npy2smv npy2smv.cpp  -lcnpy -lz --std=c++11
//usage :
// ./npy2smv test.npy test.smv
#include "cnpy.h"
#include <stdio.h>
#include <stdlib.h>
#define PI 3.141592653589793
double RTD = 180.0/PI;

int main(int argc, char *argv[])
{
  char npyfilename[100];// = "test.npy\0";
  char intfilename[100];// = "test.smv\0";
  FILE *outfile = NULL;
  short unsigned  int *intimage;
  double* floatimage;
  double max_I = 0.0;
  double max_I_x = 0.0,max_I_y = 0.0;
  double intfile_scale = 0.0;
  int xpixels,ypixels,pixels;

  double lambda0=0.025e-10;
  double distance = 100.0e-3;
  double detsize_x = 102.4e-3;
  double detsize_y = 102.4e-3;
  double pixel = 0.1e-3;
  double Xbeam=-1e99,Ybeam=-1e99;
  double phi0=0.0,phi=0.0;

  if (argc !=3){
    printf("wrong usage : 2 args needed \n"
           "./npy2smv input.npy output.smv\n");
    return -1;
  }
  else {
    sprintf(npyfilename,argv[1]);printf("input:%s\n",npyfilename);
    sprintf(intfilename,argv[2]);printf("output:%s\n",intfilename);
  }

  /*Load npy array*/
  cnpy::NpyArray arr = cnpy::npy_load(npyfilename);
  floatimage = arr.data<double>();
  xpixels = arr.shape[0];
  ypixels = arr.shape[1];
  pixels = xpixels*ypixels;
  printf("shape : %d x %d \n", xpixels,ypixels);
  // for (int i=0;i<pixels;i++)
  //     printf("data[i] = %.2f\n", floatimage[i]);


  /* output as ints james holton code snippet  */
  if(intfile_scale <= 0.0){
  	intfile_scale = 1.0;
  	if(max_I > 0.0) intfile_scale = 55000.0/max_I;
  }
  printf("intfile_scale = %g\n",intfile_scale);
  double test=0;
  intimage = (unsigned short int*)calloc(pixels+10,sizeof(unsigned short int));
  for (int j=0;j<pixels;j++){
  	test = floatimage[j] *intfile_scale+40.0;
    if(test > 65535.0) test = 65535.0;
  	intimage[j] = (unsigned short int) ( test );
    // printf("data[j] = %.2f\n", floatimage[j]);
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
