// Simple C example program that performs a Metropolis Monte Carlo Simulation with "Glauer" dynamics of 2 dimensional Ising model on a square lattice with periodic boundary conditions.
// There is no magnetic field in this simulation.
// by Kevin E Bassler, originally written 11/2021, last modified Nov. 14, 2021
//
// Program uses Numerical Recipies (NR) arrays, random numbers, 
// and uses simple graphics to show the evolution of the lattice.
//
// Parameters:
// L - linear extent of the lattice. (N = L*L total spins)
// T - temperature (real number)
// seed - random number seed (should be a [large] negative integer)
//
//

#include <stdio.h>	  // standard io library
#include <stdlib.h>	  // standard library with lots of functions
#include <math.h>	  // standard math library
#define NRANSI		  // needed for NR
#include "my_nrutil.h"    // include NR header files

long L, N;
double T;
long seed;
int **A;		// dimension pointer to pointer to int (needed for NR imatrix subroutine)
double *Boltz;

void initA();
void plotA();
int MCmove(int i, int j);
float ran2(long *idum);	  // typecast ran2

int main(int argc, char *argv[]){   // argc and argv used for command line input
  int rint, i, j;
  int accepted;
  long t, tmax;

  FILE *outfile;	// pointer to filename

  if (argc == 5){		// require argc be equal to number of command line entries
    L = atol(argv[1]);		// read long variable from command line spot 1
    T = atof(argv[2]);		// read double variable
    seed = atol(argv[3]);
    tmax = atol(argv[4]);
  } else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: Ising2Dex.x L T seed tmax \n");  // correct input syntax
    
    return 1;
  }
  
  N = L*L;					// number of total nodes
  A = imatrix(0,L-1,0,L-1);			// use NR subroutine to allocate memory for the array A[i][j], where i and j in range [0, N-1]
  initA();					// initialize lattice A
  Boltz = dvector(-2,2);
  for (i=-2; i<=2; ++i) Boltz[i] = exp(-4.0*i/T);

  fprintf(stderr,"Initial configuration: \n");
  plotA();
  getchar();

  for (t=1; t<=tmax; ++t){
    rint = N*ran2(&seed);
    i = rint/L;
    j = rint%L;

    accepted = MCmove(i,j);

    fprintf(stderr,"On MC Move number %ld: \n", t);
    if (accepted == 1) fprintf(stderr,"Lattice site chosen for move: (%d, %d); Move accepted? Yes \n",i,j);
    else fprintf(stderr,"Lattice site chosen for move: (%d, %d); Move accepted? No \n",i,j);
    plotA();
    getchar();
  }

  free_imatrix(A,0,L-1,0,L-1);			// NR subroutine to free allocated memory 
  free_dvector(Boltz,-2,2);

  return 0;
}

int MCmove(int i, int j){
  double roll;
  int ip,im,jp,jm;
  int nsum;

  ip = (i+1)%L;
  im = (i-1+L)%L;
  jp = (j+1)%L;
  jm = (j-1+L)%L;

  nsum = (A[ip][j] + A[i][jp] + A[im][j] + A[i][jm])*A[i][j];

  roll = ran2(&seed);
  if (roll < Boltz[nsum/2]) {
    A[i][j] *= -1;
    return 1;
  } else return 0;
}

void initA(){
  long i,j;
  float roll;
  for (i=0; i<L; ++i){
    for (j=0; j<L; ++j){
      roll = ran2(&seed);
      if (roll < 0.5){
        A[i][j] = 1;
      } else {
        A[i][j] = -1;
      }
    }
  }
  return;
}

void plotA(){
  long i,j;
  for (i=0; i<L; ++i){
    for (j=0; j<L; ++j){
      if (A[i][j] == -1) fprintf(stderr,"- ");
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
