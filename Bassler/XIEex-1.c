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

#include <stdio.h>	  // standard io library
#include <stdlib.h>	  // standard library with lots of functions
#include <math.h>	  // standard math library
#define NRANSI		  // needed for NR
#include "my_nrutil.h"    // include NR header files

long N;
long seed;
int **A;		// dimension pointer to pointer to int (needed for NR imatrix subroutine)

void initA();
void plotA();
float ran2(long *idum);	  // typecast ran2

int main(int argc, char *argv[]){   // argc and argv used for command line input
  long NI, NE;
  long firstnode, secondnode, rdn;
  long t, tmax;

  FILE *outfile;	// pointer to filename

  if (argc == 5){		// require argc be equal to number of command line entries
    NI = atol(argv[1]);		// read long variable from command line spot 1
    NE = atol(argv[2]);
    seed = atol(argv[3]);
    tmax = atol(argv[4]);
  } else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: XIEex.x NI NE seed tmax \n");  // correct input syntax
    
    return 1;
  }
  
  N = NI+NE;					// number of total nodes
  A = imatrix(0,N-1,0,N-1);			// use NR subroutine to allocate memory for the array A[i][j], where i and j in range [0, N-1]
  initA();					// initialize adjacency matrix A

  fprintf(stderr,"Initial configuration: \n");
  plotA();
  getchar();

  for (t=1; t<=tmax; ++t){
    firstnode = N*ran2(&seed);
    rdn = (N-1)*ran2(&seed);
    secondnode = (rdn + firstnode + 1)%N;
    if (firstnode < NI) {
      A[firstnode][secondnode] = 0;
      A[secondnode][firstnode] = 0;
    } else {
      A[firstnode][secondnode] = 1;
      A[secondnode][firstnode] = 1;
    }
    fprintf(stderr,"On MC Step number %ld: \n", t);
    fprintf(stderr,"Nodes chosen for update: %ld and %ld; New configuration: \n",firstnode,secondnode);
    plotA();
    getchar();
  }

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
      } else {
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
