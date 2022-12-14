// Simple C example program that explores the CLT  
// by Kevin E Bassler, originally written 7/2021, last modified by Yi-Kai May 25, 2022
//
// Increments are +/- 1. 
// Program uses Numerical Recipies (NR) arrays, random numbers, 
// and illustrates linear histogram binning.
//
// Version 2 adds calculation and output of error bars (showing uncertainity range).
// 
// Parameters:
// M - size of ensemble, i.e. number of walks in ensemble
// N - number of steps in each walk (should be an even integer)
// seed - random number seed (should be a [large] negative integer)
// W - width of each bin (should be an integer that is an odd multiple of 2)
// nSIG - number of standard deviations of error (sigmas) in outputed uncertainity bars (should be a positive real number)
//

#include <stdio.h>	  // standard io library
#include <stdlib.h>	  // standard library with lots of functions
#include <math.h>	  // standard math library
#define NRANSI		  // needed for NR
#include "my_nrutil.h"    // include NR header files


float ran2(long *idum);	  // typecast ran2

int main(int argc, char *argv[]){   // argc and argv used for command line input
  long M, N;
  int W;
  long seed;
  double x, delta, iii;
  long bin, nbins, offsetbin;
  long i,j;
  double value, errors, nSIG;
  long minbin, maxbin;

  float *histo;		// dimension pointer to int (needed for NR ivector subroutine)

  FILE *outfile;	// pointer to filename

  if (argc == 7){		// require argc be equal to number of command line entries
    M = atol(argv[1]);		// read long variable from command line spot 1
    N = atol(argv[2]);
    seed = atol(argv[3]);
    W = atoi(argv[4]);		// read int variable from command line spot 4
    nSIG = atof(argv[5]);		// read double variable from command line spot 5
  } else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: CLTex.x M N seed W nSIG outfile \n");  // correct input syntax
    return 1;
  }
  
  nbins = (N+W/2)/W;					// max possible number of bins required
  histo = vector(-nbins,nbins);			// use NR subroutine to allocate memory for the array histo[bin], where bin in range [-nbins, nbins]
  for(bin=-nbins; bin<=nbins; ++bin) histo[bin]=0.0;	// initialize histo[] array

  for (i=0; i<M; ++i){					// loop over walks
    x = 0;						// initialize position of walker
    for (j=0; j<N; ++j){				// loop over steps in walk
      // Next two lines create a random increment delta between -1 to +1  
      iii = 2*ran2(&seed);				// iii is randomly number between 0 to 2
      delta = iii-1;					// convert iii to number between -1 to 1
      x += delta;					// increment x by delta
    }							
    offsetbin = (x+W/2+nbins*W)/W;			// shifted value of bin for x
    bin = offsetbin - nbins;				// correct value of bin for x
    ++histo[bin];					// record bin of x in histo
  }

// next six lines finds min and max extent of non-zero bins in histo
  minbin = 0;
  maxbin = 0;
  for (bin = 1; bin <= nbins; ++bin){
    if (histo[-bin] > 0) minbin = -bin;
    if (histo[bin] > 0) maxbin = bin;
  }

  outfile = fopen(argv[6],"w");				// get output filename from command line + open output file
  for (bin=minbin; bin<=maxbin; ++bin){		// loop over non-zero bins for output
    value = ((double) histo[bin])/(W*((double) M));	// normalized real-valued value of histo entry
    errors = nSIG*sqrt((double) histo[bin])/(W*((double) M));	// normalized real-valued value of histo entry
    fprintf(outfile,"%ld, %le, %le \n",W*bin,value,errors);		// write normalized value to output file
  }
  fclose(outfile);					// close output file

  free_vector(histo,-nbins,nbins);			// NR subroutine to free allocated memory 

  return 0;
}

// below is a NR random number generator. It generated float numbers evenly over range [0,1)
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software *1(.|a. */
