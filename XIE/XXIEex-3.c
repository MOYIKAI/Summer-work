// Example program that explores the XIE model with "extreme" dynamics
// by Kevin E Bassler, originally written 10/2021, last modified Oct. 18, 2021
//
// Program uses Numerical Recipies (NR) arrays, random numbers, 
// and uses simple graphics to show the evolution of the adjacency matrix.
//
// It also shows the use of lists.
//
// Parameters:
// NI - number of introverts
// NE - number of extroverts
// seed - random number seed (should be a [large] negative integer)
//
//

#include <string.h>   // standard string lib
#include <stdio.h>	  // standard io library
#include <stdlib.h>	  // standard library with lots of functions
#include <math.h>	  // standard math library
#define NRANSI		  // needed for NR
#include "my_nrutil.h"    // include NR header files

long NI, NE;
long N;
long seed;
int **A;		// dimension pointer to pointer to int (needed for NR imatrix subroutine) for Adjacency matrix
int **ConList;		// dimension pointer to pointer to int for Connection Lists (one list for every node == matrix)
int **ListLoc;		// p2p to int for List Locations
int *ListLength;	// pointer to int for List Lengths
int *histot;      // pointer to histogram
char *buffer;

void initA();
void plotA();
void plotLists();
float ran2(long *idum);	  // typecast ran2
void initLists();
int MCupdate(int firstnode);

int main(int argc, char *argv[]){   // argc and argv used for command line input
  long firstnode, secondnode;
  long t, tmax;
  long bin, nbins, maxbin, minbin, row, col, ii, ee;  // varaiables relate to histogram
  char sn[15];

  FILE *outfile;	// pointer to filename

  if (argc == 6){		// require argc be equal to number of command line entries
    NI = atol(argv[1]);		// read long variable from command line spot 1
    NE = atol(argv[2]);
    seed = atol(argv[3]);
    tmax = atol(argv[4]);
  } else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: XXIEex.x NI NE seed tmax outfile\n");  // correct input syntax
    
    return 1;
  }
  
  N = NI+NE;					// number of total nodes
  A = imatrix(0,N-1,0,N-1);			// use NR subroutine to allocate memory for the array A[i][j], where i and j in range [0, N-1]
  ConList = imatrix(0,N-1,0,N-2);
  ListLoc = imatrix(0,N-1,0,N-1);
  ListLength = ivector(0,N-1);
  nbins = NI*NE;   // maximum possible number of bins required(max cross links)
  histot = ivector(0,nbins); //use NR subroutine to allocate memory for the array histot[bin], where bin in range [0, nbins]
  initA();					// initialize adjacency matrix A
  initLists();


  for (t=1; t<=tmax; ++t){
    firstnode = N*ran2(&seed);
    secondnode = MCupdate(firstnode);
    if (t>200000000 && t%20000 == 0) {

      for (row=0; row<nbins; ++row){ histot[row] = 0;} // Initializted the histogram

      // Getting hisogram
      for (row=0; row<N; ++row){
      if (row<NI){++histot[ListLength[row]];}
        else { ++histot[NE - ListLength[row]];}
      }

      // finds max extent of non-zero bins in histo
      maxbin = 0;
      minbin = 1;
      for (bin = 1; bin <= nbins; ++bin){if (histot[bin] > 0){ maxbin = bin;}}
      
      // create a string for outfile
      buffer = malloc(strlen(argv[5]) + 15);
      strcpy(buffer, argv[5]);
      sprintf(sn, "%d", t);
      strcat(buffer, sn);
        //fprintf(stderr, "%s\n", buffer);
      // saving files
      outfile = fopen(buffer,"w");
      for (bin=minbin; bin<=maxbin; ++bin){
        fprintf(outfile, "%ld, %ld \n", bin, histot[bin]);}
      fclose(outfile);
      }  
  }

  free(buffer);
  free_ivector(histot,0,nbins);
  free_imatrix(A,0,N-1,0,N-1);			// NR subroutine to free allocated memory 
  free_imatrix(ConList,0,N-1,0,N-2);
  free_imatrix(ListLoc,0,N-1,0,N-1);
  free_ivector(ListLength,0,N-1);

  return 0;
}

int MCupdate(int firstnode){
  int i,j,l,secondnode;
  // For Introvert listlength = 0 means no cross links. For Extrovert 0 means they link to all other nodes
  if (ListLength[firstnode] == 0){ secondnode = -1;}

  else {
    j = ListLength[firstnode]*ran2(&seed); // j is related to cross link length "more cross link <=> bigger j" 0 =< j <= max listlength-1
    secondnode = ConList[firstnode][j];    // second node is not randomly choose, choose from the connection one

    if (j < (ListLength[firstnode] - 1)){ // if j is not at the maxlength
      l = ConList[firstnode][ListLength[firstnode]-1];
      ConList[firstnode][j] = l;
      ListLoc[firstnode][l] = j;
    }
    
    --ListLength[firstnode]; // First node has been chosen. For Extrovert -- is adding link, For Introvert -- is cutting link.

    // If fistnode is Introvert cut link in A
    if (firstnode < NI) { 
      A[firstnode][secondnode] = 0;
      A[secondnode][firstnode] = 0;

      // Below statement is considering the second node is I or E
      // For I2I is -1 Listlegnth to second node. 
      if (secondnode < NI) {
        i = ListLoc[secondnode][firstnode];
	      if (i < (ListLength[secondnode] -1)){
	        l = ConList[secondnode][ListLength[secondnode]-1];
	        ConList[secondnode][i] = l;
	        ListLoc[secondnode][l] = i;
          }
        --ListLength[secondnode];
      }
      // For I2E is +1 Listlength to second node
      else {
        ConList[secondnode][ListLength[secondnode]] = firstnode;
	      ListLoc[secondnode][firstnode] = ListLength[secondnode];
        ++ListLength[secondnode];
      }
    }
    // Else is Extrovert add link in A
    else {
      A[firstnode][secondnode] = 1;
      A[secondnode][firstnode] = 1;
      if (secondnode >= NI) {
        i = ListLoc[secondnode][firstnode];
	      if (i < (ListLength[secondnode] -1)){
	        l = ConList[secondnode][ListLength[secondnode]-1];
	        ConList[secondnode][i] = l;
	        ListLoc[secondnode][l] = i;
	      }
        --ListLength[secondnode];
      } 
      else {
        ConList[secondnode][ListLength[secondnode]] = firstnode;
	      ListLoc[secondnode][firstnode] = ListLength[secondnode];
        ++ListLength[secondnode];
      }
    }
  }
  return secondnode;
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

void initLists(){
  long i,j,jj;
  // If A[i][jj] List length
  for (i=0; i<NI; ++i){
    ListLength[i] = 0;
    for (j=1; j<N; ++j){
      jj = (i+j)%N;
      if (A[i][jj] == 1){
        ConList[i][ListLength[i]] = jj;
	      ListLoc[i][jj] = ListLength[i];
        ++ListLength[i]; // In Introvert add 1 means adding link in Listlength
      }
    }
  }

  for (i=NI; i<N; ++i){
    ListLength[i] = 0;
    for (j=1; j<N; ++j){
      jj = (i+j)%N;
      if (A[i][jj] == 0){
        ConList[i][ListLength[i]] = jj;
	      ListLoc[i][jj] = ListLength[i];
        ++ListLength[i]; // In Extrovert add 1 means cutting link in Listlength
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
      else if (A[i][j] == 1) fprintf(stderr,"+ ");
      else fprintf(stderr,"E ");
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n\n");

  return;
}

void plotLists(){
  long i,j;

  fprintf(stderr,"ConList: \n");
  for (i=0; i<N; ++i){
    fprintf(stderr,"For node %ld, Length = %d: ",i,ListLength[i]);
    for (j=0; j<ListLength[i]; ++j){
      fprintf(stderr,"%d ",ConList[i][j]);
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n\n");

  fprintf(stderr,"ListLoc: \n");
  for (i=0; i<N; ++i){
    for (j=0; j<N; ++j){
      fprintf(stderr,"%d ",ListLoc[i][j]);
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
