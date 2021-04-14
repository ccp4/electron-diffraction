#include "slicelib.hpp"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int Z=1;                       //Atomic number
double r0=0.0,rmax=2.0,dr=-1.0;//radius in Angstrom
int nrs=100;                   //number of points
char outf[100]="vatom.txt";    //outfile name
char func='a';                 //function called (a(vatom),z(vzatom),L(vzLUT))
double (*vfunc)(int,double);
int compute=1;
/*####################################################################
                            utils
######################################################################*/
void parse(int argc, char **argv){
  char c;
  while ((c = getopt (argc, argv, "Z::N::r::R::d::f::o::")) != -1) {
     switch (c) {
       case 'Z': Z    = atoi(optarg);break;
       case 'r': r0   = atof(optarg);break;
       case 'R': rmax = atof(optarg);break;
       case 'd': dr   = atof(optarg);break;
       case 'N': nrs  = atoi(optarg);break;
       case 'o': if(optarg) sprintf(outf,"%s",optarg);break;
       case 'f': func=optarg[0];break;
     }
   }
   if (dr<0.0) dr=(rmax-r0)/(double)nrs;
   switch (func) {
     case 'a': vfunc=&vatom     ;break;
     case 'z': vfunc=&vzatom    ;break;
     case 's': {vfunc=&vzatomSpline;compute=0;}break;
     case 'L': vfunc=&vzatomLUT ;break;
   }

   printf("Z=%d, N=%d,r0=%.2f,dr=%.2f,rmax=%.2f,func=%c,outf=%s\n",
         Z,nrs,r0,dr,rmax,func,outf);
}
void write_v(double *r, double *v){
  FILE *fp=fopen(outf,"w");
  for (int i=0; i<nrs; i++){
    fprintf(fp,"%g %g \n",r[i],v[i]);
  }
  fclose(fp);
  printf("successfully written to %s\n",outf);
}

/*####################################################################
                            Main
######################################################################*/
void compute_vfunc(){
  double *r = new double[nrs];
  double *v = new double[nrs];
  double ri=r0,r_rsq=0.0;
  for (int i=0; i<nrs; i++){
    ri  += dr;
    r[i] = ri;
    r_rsq = (func=='L')?ri*ri:ri;//printf("%.2f ", r_rsq);
    v[i] = vfunc(Z, r_rsq);
  }

  write_v(r,v);
  delete[] r,v;
}
int main(int argc, char **argv){
  parse(argc,argv);
  if(compute)
    compute_vfunc();
  else
    vfunc(Z,0.0);
  return 0;
}
