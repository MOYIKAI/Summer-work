#include <stdio.h>	  // standard io library
#include <string.h>   // standard string library
#include <stdlib.h>	  // standard library with lots of functions
#include <math.h>	    // standard math library
#define NRANSI		    // needed for NR
#include "my_nrutil.h"    // include NR header files

#define MAX_LINES 500000

int main(int argc, char *argv[]){ 
  int nums[MAX_LINES]; // create an integer array with max line 500000
  long i,j,lnum;
  float errors, nSIG;
  long max, min;
  long offsetbin, bin, nbins, W;

  int *histo;

  FILE *infile;   //pointer to inputfile
  FILE *outfile;	// pointer to outputfile

  if (argc == 5){
    W = atol(argv[1]);      // read long variable from command line spot 1
    nSIG = atof(argv[2]);		// read double variable from command line spot 2
  }
  else {			// error in value of argc 
    fprintf(stderr,"\n Initialization error:\n");
    fprintf(stderr,"Usage: histo.x W nSIG infile outfile\n");  // correct input syntax 
    return 1;
  }
  

  if (infile = fopen(argv[3],"r"))
  {
    min = 0;  // min num in file
    max = 0;  // max num in file
    i = 0;
    j = 0;
    lnum = 0; // numbers of element in num array

    // Read all lines in the file and save into num array
    while (fscanf(infile, "%d", &nums[i]) != EOF ){
      if (nums[i] < min) min = nums[i];
      if (nums[i] > max) max = nums[i];
      i++;
      lnum++;
    }
    fclose(infile);

    nbins = (max-min+1)/2*W;		                      // max possible number of bins required
    histo = ivector(-nbins,nbins);                    // use NR subroutine to allocate memory for the array histo[bin], where bin in range [-nbins, nbins]
    for(bin=-nbins; bin<=nbins; ++bin){histo[bin]=0;} // initialize histo[] array

    // Get historgram
    for (j=0; j<=lnum; j++){bin = nums[j];++histo[bin];}

    // Print the histogram value out
    outfile = fopen(argv[4],"w");
    for (bin=-nbins; bin<=nbins; ++bin){
      errors = nSIG*sqrt((float) histo[bin]);	// normalized real-valued value of histo entry
      fprintf(outfile, "%d, %d, %f\n", bin, histo[bin], errors);
    }
    fclose(outfile);

    free_ivector(histo,-nbins,nbins);
  }
  else{ fprintf(stderr,"No such input file \n"); return 1;}

  return 0;
}