#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "mystat.h"
#include "myio.h"
#include "probterms.h"
#include "simfunctions.h"


//reads allele freqs from files as specified
void readAlleleFreqs(char *path, char **STRs, int *nalleles, int npops, double ***AFs, int **samplesizes) {
  int i, j, k;
  FILE *fin;
  char filename[100], tmp[100];

  //run through STR files
  for(i=0; i<13; i++) {
    //get filename
    sprintf(filename, "%s/%s.dat", path, STRs[i]);
    //printf("Opening input file: %s\n", filename);

    //open file
    fin = fopen(filename, "r"); 
    if(fin == NULL) {
      printf("Unable to locate input file: %s\n", filename);
      exit(1);
    }
    
    //read STR name
    fscanf(fin, "%s", tmp);

    //read population names
    for(j=0; j<npops; j++) {
      fscanf(fin, "%s", tmp);
    }      
    
    //run through alleles
    for(j=0; j<nalleles[i]; j++) {
      //read allele name
      fscanf(fin, "%s", tmp);
      //printf("allele is %s \n", tmp);

      //run through pops
      for(k=0; k<npops; k++) {
	fscanf(fin, "%lf", &AFs[i][j][k]);
      }
      //printDouArr(AFs[i][j], npops);
    }

    fclose(fin);
  }

  //ALERT: THIS BELOW SHOULD ONLY BE ON WHEN AF FILES ARE OUT OF 1, NOT 100  (/100)
  //so that means use it for BudowleSheaNiezgodaChakraborty
  /*
  for(i=0; i<13; i++) {
    for(j=0; j<nalleles[i]; j++) {
      for(k=0; k<npops; k++) {
	AFs[i][j][k] = AFs[i][j][k]*100.0;
      }
    }
  }
  */

  //get sample sizes file name
  sprintf(filename, "%s/samplesize.dat", path);
  
  //open sample sizes file
  fin = fopen(filename, "r"); 
  if(fin == NULL) {
    printf("Unable to locate input file: %s\n", filename);
    exit(1);
  }
  
  //read header with 'samplesize' and population names
  for(i=0; i<(npops+1); i++) {
    fscanf(fin, "%s", tmp);
  }

  //read sample sizes for each pop
  for(j=0; j<13; j++) {
    fscanf(fin, "%s", tmp);
    for(i=0; i<npops; i++) {
      fscanf(fin, "%d", &samplesizes[j][i]);
    }
  }
}


//reads reference Y haplotypes from files as specified
void readRefHaps(char *path, int npops, int ***refHaps, int *samplesizes, int nSTRs, char **popnames) {
  int i, j, k;
  FILE *fin;
  char pathfile[1000];

  //run through populations
  for(i=0; i<npops; i++) { 
    //open pop file
    sprintf(pathfile, "%s/%s.dat", path, popnames[i]);
    fin = fopen(pathfile, "r"); 
    if(fin == NULL) {
      printf("Unable to locate input file: %s\n", pathfile);  exit(1);
    }

    //run through individuals
    for(j=0; j<samplesizes[i]; j++) {
      //run through loci
      for(k=0; k<nSTRs; k++) {
	fscanf(fin, "%d", &refHaps[i][j][k]);
      }
    }
    fclose(fin);
  }
}

//reads reference Y haplotypes from files as specified
void readYSampleSizes(char *path, int npops, int *samplesizes) {
  int i;
  FILE *fin;
  char tmp[100], pathfile[1000];

  //get sample sizes file name
  sprintf(pathfile, "%s/samplesize.dat", path);
  
  //open sample sizes file
  fin = fopen(pathfile, "r");
  if(fin == NULL) {
    printf("Unable to locate input file: %s\n", pathfile);  exit(1);
  }
  
  //read header with 'samplesize' and population names
  for(i=0; i<(npops); i++)
    fscanf(fin, "%s", tmp);

  //read sample sizes for each pop
  for(i=0; i<npops; i++)
    fscanf(fin, "%d", &samplesizes[i]);

  fclose(fin);
}





///////////////////////////////////////////////////////////////////////////////
// main
int main(int argc, char** argv) {
  int i, j, k, l, m, n;
  char path[1000] = "", Ypath[1000]="", filename[1500], relationship[50];
  double ***alleleFreqs, *tmparr1, *tmparr2;
  int **samplesizes, *Ysamplesizes, *modYsamplesizes;
  int **indiv1, **indiv2, *Yindiv1, *Yindiv2;
  double *ks;
  double *LRs;
  int *nalleles;
  int ***refHaps, totYindivs;
  double **estYThetas, **YhapFreqs, **modYhapFreqs;
  int *popmap;
  double **FPRs;
  int *nSameYchr;

  //constants
  int npops = 5;
  int nSTRs = 13;
  double ndatabase = 1824085.0;
  int seed = time(NULL);
  double simtheta = 0.010;
  double calctheta = 0.010;
  int niter = 10000;
  char *STRs[] = {"CSF1PO", "FGA", "TH01", "TPOX", "vWA", "D3S1358", "D5S818", 
		  "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11"};
  int nYpops = 13;
  int nYSTRs = 16;  
  char *Ypopnames[] = {"African_American", "Asian", "Asian_Indian", "Caucasian", "Chinese", 
		       "Filipino", "Hispanic", "Japanese", "Malay", "Native_American", 
		       "Sub-sharan_African", "Thai", "Vietnamese"};

  srandom(seed);
  



  //parsing command line arguments  
  if(argc < 3) {
    printf("Missing input file(s). \n\tUsage: ./autYsimUnrelPopPair -f PATH -g YPATH -n nLOCI -s simtheta -c calctheta -i nITER\n");
    exit(0);
  }

  for(i=1; i<argc; i++) {
    if(argv[i][0] == '-'){
      switch(argv[i][1]) {
      case 'f' :
	(void) strcpy(path, &argv[i][3]); break;
      case 'g' :
	(void) strcpy(Ypath, &argv[i][3]); break;
      case 'n' :
	nSTRs = atoi(&argv[i][3]); break;
      case 's' :
	simtheta = atof(&argv[i][3]); break;
      case 'c' :
	calctheta = atof(&argv[i][3]); break;
      case 'i' :
	niter = atoi(&argv[i][3]); break;
      default :
	printf("Unknown option %s\n",argv[i]);
      }
    }
  }

  //malloc and fill-in nalleles
  nalleles = (int*) malloc(13*sizeof(int));
  //ALERT: use the first version for BudowleMoretti and the second for BudowleSheaNiezgodaChakraborty
  nalleles[0]=11;  nalleles[1]=23; nalleles[2]=9; nalleles[3]=7; nalleles[4]=10;
  nalleles[5]=9; nalleles[6]=10; nalleles[7]=9; nalleles[8]=11; nalleles[9]=9;
  nalleles[10]=9; nalleles[11]=16; nalleles[12]=23;
  /*
  nalleles[0]=14;  nalleles[1]=31; nalleles[2]=10; nalleles[3]=10; nalleles[4]=12; 
  nalleles[5]=12; nalleles[6]=12; nalleles[7]=18; nalleles[8]=11; nalleles[9]=12; 
  nalleles[10]=11; nalleles[11]=20; nalleles[12]=34;
  */


  //malloc samplesizes
  samplesizes = (int**) malloc(13*sizeof(int*));
  for(i=0; i<13; i++)
    samplesizes[i] = (int*) malloc(npops*sizeof(int));
  Ysamplesizes = (int*) malloc(nYpops*sizeof(int));
  modYsamplesizes = (int*) malloc(nYpops*sizeof(int));

  //malloc indiv1 and indiv2 [STR][allele]
  indiv1 = (int**) malloc(nSTRs*sizeof(int*));
  indiv2 = (int**) malloc(nSTRs*sizeof(int*));
  for(j=0; j<nSTRs; j++) {
    indiv1[j] = (int*) malloc(2*sizeof(int));
    indiv2[j] = (int*) malloc(2*sizeof(int));
  }

  //malloc ks
  ks = (double*) malloc(3*sizeof(double*));

  //malloc tmparrs
  tmparr1 = (double*) malloc(nYpops*sizeof(double));  
  tmparr2 = (double*) malloc(nYpops*sizeof(double));

  //malloc LRs[assumed pop number]
  LRs = (double*) malloc(npops*sizeof(double));

  //malloc FPRs[sim pop1][sim pop 2
  FPRs = (double**) malloc(npops*sizeof(double*));
  for(i=0; i<npops; i++) 
    FPRs[i] = (double*) malloc(npops*sizeof(double));

  //malloc alleleFreqs[STR number][allele][population]
  alleleFreqs = (double***) malloc(13*sizeof(double**));
  for(i=0; i<13; i++) {
    alleleFreqs[i] = (double**) malloc(nalleles[i]*sizeof(double*));
    for(j=0; j<nalleles[i]; j++) {
      alleleFreqs[i][j] = (double*) malloc(npops*sizeof(double));
    }
  }


  //read in sample sizes for Y haps
  readYSampleSizes(Ypath, nYpops, Ysamplesizes);
  //printf("ysamplesizes are: \n");
  //printIntArr(Ysamplesizes, nYpops);

  //init totYindivs
  totYindivs = 0;
  for(i=0; i<nYpops; i++) 
    totYindivs += Ysamplesizes[i];

  //malloc refHaps[population][indiv][STR]
  refHaps = (int***) malloc(nYpops*sizeof(int**));
  for(i=0; i<nYpops; i++) {
    refHaps[i] = (int**) malloc(Ysamplesizes[i]*sizeof(int*));
    for(j=0; j<Ysamplesizes[i]; j++) {
      refHaps[i][j] = (int*) malloc(nYSTRs*sizeof(int));
    }
  }

  //malloc estYThetas[pop1][pop2]
  estYThetas = (double**) malloc(nYpops*sizeof(double*));
  for(i=0; i<nYpops; i++) 
    estYThetas[i] = (double*) malloc(nYpops*sizeof(double));

  //malloc YhapFreqs[population][haplotype]
  YhapFreqs = (double**) malloc(nYpops*sizeof(double*));
  for(i=0; i<nYpops; i++) 
    YhapFreqs[i] = (double*) malloc(totYindivs*sizeof(double));

  //malloc modYhapFreqs[population][haplotype]
  modYhapFreqs = (double**) malloc(nYpops*sizeof(double*));
  for(i=0; i<nYpops; i++) 
    modYhapFreqs[i] = (double*) malloc(totYindivs*sizeof(double));

  //malloc and set popmap to be popmap[auto pop i] = Y pop i
  popmap = (int*) malloc(npops*sizeof(int));
  popmap[0] = 12;
  popmap[1] = 0;
  popmap[2] = 3;
  popmap[3] = 6;
  popmap[4] = 9;

  //malloc Yindivs
  Yindiv1 = (int*) malloc(1*sizeof(int));
  Yindiv2 = (int*) malloc(1*sizeof(int));

  //read in allele frequencies
  readAlleleFreqs(path, STRs, nalleles, npops, alleleFreqs, samplesizes);

  //read in reference haplotypes
  readRefHaps(Ypath, nYpops, refHaps, Ysamplesizes, nYSTRs, Ypopnames);

  //estimate theta for Y chr simulations
  estimateTheta(refHaps, nYpops, Ysamplesizes, nYSTRs, estYThetas);
  //print2dDouArr(estYThetas, nYpops, nYpops);

  //get Y hap frequencies
  calcHapFreqs(refHaps, nYpops, Ysamplesizes, nYSTRs, YhapFreqs);
  sprintf(filename, "results/autYsimUnrelPopPair/%s/YhapFreqs1.res", path);
  write2dDouArr(YhapFreqs, nYpops, totYindivs, filename); 

  //loop through relationships
  for(n=0; n<2; n++) {
    //set relationship as p-o and sib
    switch(n) {
    case 0: //set parent-offspring
      ks[0] = 0.0; ks[1] = 1.0; ks[2] = 0.0; 
      sprintf(relationship, "p-o");
      break;
    case 1: //set sib
      ks[0] = .25; ks[1] = .5; ks[2] = .25; 
      sprintf(relationship, "sib");
      break;
    }
 
    //init FPRs 
    for(i=0; i<npops; i++)
      for(j=0; j<npops; j++) 
	FPRs[i][j] = 0.0;
    
    /*
    nSameYchr = (int*) malloc(npops*sizeof(int));
    for(i=0; i<npops; i++) 
      nSameYchr[i] = 0;
    */

    for(k=0; k<niter; k++) {
      //run through each population pair, making rand indivs and calcing LR for a match  
      for(i=0; i<npops; i++) {
	for(j=i; j<npops; j++) {
	  //make two unrel individuals assuming memberships in each population
	  if(i!=j) {
	    thetaRandIndiv(alleleFreqs, 13, nalleles, i, simtheta, indiv1);
	    thetaRandIndiv(alleleFreqs, 13, nalleles, j, simtheta, indiv2);
	    thetaRandIndivYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, Yindiv1, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	    thetaRandIndivYchr(YhapFreqs[popmap[j]], modYhapFreqs[popmap[j]], totYindivs, Yindiv2, Ysamplesizes, modYsamplesizes, popmap[j], nYpops);
	  } else {
	    thetaRandUnrelIndivPair(alleleFreqs, 13, nalleles, i, simtheta, indiv1, indiv2);
	    //thetaRandUnrelIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, estYThetas[popmap[i]][popmap[i]], Yindiv1, Yindiv2, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	    randUnrelIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, Yindiv1, Yindiv2, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	    /*if(Yindiv1[0] == Yindiv2[0]) { 
	      nSameYchr[i]++;
	      //printf("got same Ychr in pop %d\n", i);
	      }*/
	  }

	  LR(indiv1, indiv2, alleleFreqs, npops, nalleles, samplesizes, tmparr1, calctheta, ks);
	  for(l=0; l<npops; l++) 
	    LRs[l] = tmparr1[l];

	  YLR(Yindiv1[0], Yindiv2[0], modYhapFreqs, nYpops, modYsamplesizes, tmparr1);
	  for(l=0; l<npops; l++)
	    LRs[l] *= tmparr1[popmap[l]];

	  //check if would be a positive
	  if(isfinite(LRs[1]) && isfinite(LRs[2]) && isfinite(LRs[3])
	     && !isnan(LRs[1]) && !isnan(LRs[2]) && !isnan(LRs[3])
	     && LRs[1]/ndatabase>.1 && LRs[2]/ndatabase>.1 && LRs[3]/ndatabase>.1
	     && (LRs[1]/ndatabase>1 || LRs[2]/ndatabase>1 || LRs[3]/ndatabase>.1)) {
	    //printf("!!!got positive with pop %d and %d, with Yhaps %d and %d, aut genos:\n", i, j, Yindiv1[0], Yindiv2[0]);
	    //print2dIntArr(indiv1, 13, 2);
	    //print2dIntArr(indiv2, 13, 2);
	    //printf("LRs:\n");
	    //printDouArr(LRs, npops);
	    FPRs[i][j]++;
	  }

	} //end pop j loop
      } //end pop i loop
    }//end niter loop

    //print nSameYchrs
    //printIntArr(nSameYchr, npops);



    //normalize FPRs 
    for(i=0; i<npops; i++)
      for(j=i; j<npops; j++) 
	FPRs[i][j] /= niter;

    //write results to files    
    sprintf(filename, "results/autYsimUnrelPopPair/%s/fpr/%s_%diter_simtheta%1.3f_calctheta%1.3f_FPRtab_rep2.res", path, relationship, niter, simtheta, calctheta);
    write2dDouArr(FPRs, npops, npops, filename);
  }//end relationship loop

  //free memory
  for(i=0; i<nSTRs; i++) {
    for(j=0; j<nalleles[i]; j++) 
      free(alleleFreqs[i][j]);
    free(alleleFreqs[i]);
  }
  free(alleleFreqs);
  
  free(nalleles);

  for(i=0; i<nSTRs; i++)
    free(samplesizes[i]);
  free(samplesizes);

  for(j=0; j<nSTRs; j++) {
      free(indiv1[j]);
      free(indiv2[j]);
  }
  free(indiv1);
  free(indiv2);

  free(ks);
  free(tmparr1);
  free(tmparr2);

  free(LRs);

  free(Yindiv1);
  free(Yindiv2);

  return 0;
}
