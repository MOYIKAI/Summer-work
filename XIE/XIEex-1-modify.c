// Simple C example program that explores the XIE model with "heat-bath" dynamics
// by Kevin E Bassler, originally written 9/2021, last modified Sept. 5, 2021
//
// Program uses Numerical Recipies (NR) arrays, random numbers, 
// and uses simple graphics to show the evolution of the adjacency matrix.
//
// Parameters:
// NI - number of introverts
// NE - number of extroverts
// seed - random number seed (should be a [large] negative integer)
//
//

#include <stdio.h>	      // standard io library
#include <stdlib.h>	      // standard library with lots of functions
#include <math.h>	        // standard math library
#define NRANSI		        // needed for NR
#include "my_nrutil.h"    // include NR header files

long N;     // Size of adjacency matrix A
long seed;  // Random number seed
int **A;		// dimension pointer to pointer to int (needed for NR imatrix subroutine)

void initA(); // initialized matrix A
void plotA(); // plot A on Terminal 
float ran2(long *idum);	  // typecast ran2

int main(int argc, char *argv[]){   // argc and argv used for command line input
  long NI, NE;                      // numbers of introvert and extrovert
  long firstnode, secondnode, rdn, col, row; // Index of matrix
  long t, tmax, plot;                     // how many time steps the system evlove
  long ii, ee, clkie, clkei, tclk;          // different kinds of link
  long bin, nbins, maxbin, minbin;  // varaiables relate to histogram

  int *histoe; // dimension pointer to int (needed for NR ivector subroutine)
  int *histoi; // dimension pointer to int (needed for NR ivector subroutine)
  int *histot;  // dimension pointer to int (needed for NR ivector subroutine)

  FILE *outfile;	// pointer to filename

  if (argc == 7){		// require argc be equal to number of command line entries
    NI = atol(argv[1]);		// read long variable from command line spot 1
    NE = atol(argv[2]);
    seed = atol(argv[3]);
    tmax = atol(argv[4]);
    plot = atol(argv[5]);
  } else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: XIEex.x NI NE seed tmax plot(0/1) outfile \n");  // correct input syntax
    
    return 1;
  }
  
  N = NI+NE;					// number of total nodes
  A = imatrix(0,N-1,0,N-1);			// use NR subroutine to allocate memory for the array A[i][j], where i and j in range [0, N-1]
  
  initA();					// initialize adjacency matrix A
  if (plot == 1) {  // show the initial state
    fprintf(stderr,"Initial configuration: \n");
    plotA();}       

  for (t=1; t<=tmax; ++t){
    firstnode = N*ran2(&seed);
    rdn = (N-1)*ran2(&seed);
    secondnode = (rdn + firstnode + 1)%N;
    if (firstnode < NI) {
      A[firstnode][secondnode] = 0; // 0 means cutting link
      A[secondnode][firstnode] = 0; // 0 means cutting link
    } else {
      A[firstnode][secondnode] = 1; // 1 means linking
      A[secondnode][firstnode] = 1; // 1 means linking
    }
  }

  if (plot == 1){ // show the final state
    fprintf(stderr,"On MC Step number %ld: \n", t);
    fprintf(stderr,"Nodes chosen for update: %ld and %ld; New configuration: \n",firstnode,secondnode);
    plotA();}
  
  // simple check if it is steady state
  ii = 0;
  ee = 0;
  for (row=0; row<NI; ++row){
    for (col=0; col<NI; ++col){
      if (A[row][col] == 0){ ++ii;}
    }
  }
  for (row=NI; row<N; ++row){
    for (col=NI; col<N; ++col){
      if (A[row][col] == 1){ ++ee;}
    }
  }

  nbins = NI*NE;   // maximum possible number of bins required(max cross links)
  histot = ivector(0,nbins); //use NR subroutine to allocate memory for the array histot[bin], where bin in range [0, nbins]
  for (bin=0; bin<=nbins; ++bin) { histot[bin]=0; } //initialize histo[] array
  histoi = ivector(0,nbins); //use NR subroutine to allocate memory for the array histoi[bin], where bin in range [0, nbins]
  for (bin=0; bin<=nbins; ++bin) { histoi[bin]=0; } //initialize histo[] array
  histoe = ivector(0,nbins); //use NR subroutine to allocate memory for the array histo[bin], where bin in range [0, nbins]
  for (bin=0; bin<=nbins; ++bin) { histoe[bin]=0; } //initialize histo[] array
  
  //Below's for loops are finding the histogram of cross links
  // First get the cross links of introvert to extrovert
  
  tclk =  0; // initialized number of total cross links
  clkie = 0; // initialized number of I-E cross links

  for (row=0; row<NI; ++row){
    for (col=NI; col<N; ++col){
      if (A[row][col] == 1){ ++tclk; ++clkie;}
    }
    ++histoi[clkie]; ++histot[tclk]; tclk = 0; clkie =0;
  }
  // Second get the cross links of extrovert to introvert
  tclk = 0; // initialized number of total cross links
  clkei =0; // initialized number of E-I cross links

  for (row=NI; row<N; ++row){
    for (col=0; col<NI; ++col){
      if (A[row][col] == 1){ ++tclk; ++clkei;}
    }
    ++histoe[clkei]; ++histot[tclk]; tclk = 0; clkei =0;
  }

  // next six lines finds max extent of non-zero bins in histo
  maxbin = 0;
  minbin = 1;
  for (bin = 1; bin <= nbins; ++bin){
    if (histot[bin] > 0){ maxbin = bin;}
    if (histoe[bin] > 0){ maxbin = bin;}
    if (histoi[bin] > 0){ maxbin = bin;}
  }

  //files for steady state
  if (ii == NI*NI && ee == (NE-1)*NE){
    // next five lines output the histogram file
    outfile = fopen(argv[6],"w");
    for (bin=minbin; bin<=maxbin; ++bin){
      fprintf(outfile, "%ld, %ld, %ld, %ld \n", bin, histot[bin], histoi[bin], histoe[bin]);}
    fclose(outfile);
    fprintf(stderr, "Steady state \n");
  }

  else {
    outfile = fopen(argv[6],"w");
    for (bin=minbin; bin<=maxbin; ++bin){
      fprintf(outfile, "%ld, %ld, %ld, %ld \n", bin, histot[bin], histoi[bin], histoe[bin]);}
    fclose(outfile); 
    fprintf(stderr, "Not steady state \n");
    fprintf(stderr, "%ld, %ld \n", ii,ee);
    }

  free_ivector(histot,0,nbins);      // Free allocated memory for histot
  free_ivector(histoi,0,nbins);      // Free allocated memory for histoi
  free_ivector(histoe,0,nbins);      // Free allocated memory for histoe
  free_imatrix(A,0,N-1,0,N-1);			// NR subroutine to free allocated memory 

  return 0;
}

void initA(){
  long i,j;
  float roll;
  for (i=0; i<N; ++i){
    A[i][i] = 0;
    for (j=i+1; j<N; ++j){
      roll = ran2(&seed);
      if (roll < 0.5){
        A[i][j] = 1;
	      A[j][i] = 1;
      } 
      else {
        A[i][j] = 0;
	      A[j][i] = 0;
      }
    }
  }
  return;
}

void plotA(){
  long i,j;
  for (i=0; i<N; ++i){
    for (j=0; j<N; ++j){
      if (A[i][j] == 0) fprintf(stderr,"o ");
      else fprintf(stderr,"+ ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n\n");

  return;
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
