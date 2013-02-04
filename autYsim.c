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

//subsamples so that all AFs are from same nindivs
void sampleAlleleFreqs(int *nalleles, int npops, double ***AFs, int **samplesizes) {
  int i, j, k, l, m;
  int nA, curA, *obsAs;  
  int subsampsize = 100;
  double *curAFs;

  //run through STR files
  for(i=0; i<13; i++) {
    nA = nalleles[i];

    //malloc obsAs and curAFs
    obsAs = (int*) malloc(nA*sizeof(int));
    curAFs = (double*) malloc(nA * sizeof(double));

    for(k=0; k<npops; k++) {
      //init obsAs and curAFs
      for(j=0; j<nA; j++) {
	obsAs[j] = 0;
	curAFs[j] = AFs[i][j][k];
      }

      //get subsample
      for(l=0; l<subsampsize; l++) {
	curA = randValFromDist(curAFs, nA);
	obsAs[curA]++;
	for(m=0; m<nA; m++) {
	  if(m==curA) {
	    curAFs[m] = ((curAFs[m]*samplesizes[i][k]/100.0)-1.0)/((samplesizes[i][k]-1)/100.0);
	  } else {
	    curAFs[m] = ((curAFs[m]*samplesizes[i][k]/100.0))/((samplesizes[i][k]-1)/100.0);
	  }
	}
	samplesizes[i][k]--;
      }

      //save in AFs
      for(j=0; j<nA; j++) {
	AFs[i][j][k] = 100.0*obsAs[j]/subsampsize;
      }
    }
    free(obsAs);
    free(curAFs);
  }

  //set samplesizes
  for(j=0; j<13; j++) {
    for(i=0; i<npops; i++) {
      samplesizes[j][i] = subsampsize;
    }
  }
}

//subsamples so that all Yhaps are from same nindivs
void sampleYHapFreqs(int ***fullRefHaps, int nYpops, int *Ysamplesizes, int nYSTRs, int ***sampledRefHaps, int newsamplesize) {
  int i, j, k;
  int *fullindices;

  //run through pop groups
  for(i=0; i<nYpops; i++) {
    //make array of indices for this pop
    fullindices = malloc(Ysamplesizes[i]*sizeof(int));
    for(j=0; j<Ysamplesizes[i]; j++) 
      fullindices[j] = j;

    //shuffle indices to randomly choose sample
    shuffle(fullindices, Ysamplesizes[i]);

    //copy sampled Yhaps to sampleRefHaps
    if(i!=10) {
      for(j=0; j<newsamplesize; j++) {
	for(k=0; k<nYSTRs; k++) {
	  sampledRefHaps[i][j][k] = fullRefHaps[i][fullindices[j]][k];
	}      
      }
    }

    free(fullindices);
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


//reads reference Y haplotypes sample sizes
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
  double ***alleleFreqs, *pmatches, *tmparr1, *tmparr2;
  int **samplesizes, *Ysamplesizes, *modYsamplesizes;
  int **indiv1, **indiv2, ***indivs, *Yindivs;
  double *ks, *ks2;
  double **LRs;
  double **powLRs;
  FILE *fout;
  int *nalleles;
  int ***refHaps, ***sampledRefHaps, totYindivs;
  double **estYThetas, **YhapFreqs, **modYhapFreqs;
  int *popmap;


  //constants
  int npops = 5;
  int nSTRs = 13;
  int seed = 0; //time(NULL);
  double simtheta = 0.010;
  double calctheta = 0.010;
  double z = 1.96; 
  int niter = 10000;
  //int nallele = 10;
  //int samplesize = 200;
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
    printf("Missing input file(s). \n\tUsage: ./autYsim -f PATH -g YPATH -n nLOCI -s simtheta -c calctheta -i nITER\n");
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
  nalleles[0]=11;  nalleles[1]=23; nalleles[2]=9; nalleles[3]=7; nalleles[4]=10;
  nalleles[5]=9; nalleles[6]=10; nalleles[7]=9; nalleles[8]=11; nalleles[9]=9;
  nalleles[10]=9; nalleles[11]=16; nalleles[12]=23;

  //malloc samplesizes
  samplesizes = (int**) malloc(13*sizeof(int*));
  for(i=0; i<13; i++)
    samplesizes[i] = (int*) malloc(npops*sizeof(int));
  Ysamplesizes = (int*) malloc(nYpops*sizeof(int));
  modYsamplesizes = (int*) malloc(nYpops*sizeof(int));

  //malloc indivs [indiv][STR][allele]
  indivs = (int***) malloc(2*sizeof(int**));
  for(i=0; i<2; i++) {
    indivs[i] = (int**) malloc(nSTRs*sizeof(int*));
    for(j=0; j<nSTRs; j++) 
      indivs[i][j] = (int*) malloc(2*sizeof(int));
  }

  //malloc ks
  ks = (double*) malloc(3*sizeof(double*));
  ks2 = (double*) malloc(3*sizeof(double*));

  pmatches = (double*) malloc(npops*sizeof(double));

  //malloc tmparrs
  tmparr1 = (double*) malloc(nYpops*sizeof(double));  
  tmparr2 = (double*) malloc(nYpops*sizeof(double));

  //malloc LRs[assumed pop number][sim pop number]
  LRs = (double**) malloc(npops*sizeof(double*));
  for(j=0; j<npops; j++) {
    LRs[j] = (double*) malloc(npops*sizeof(double));
  }

  //malloc powLRs[assumed pop number][sim pop number]
  powLRs = (double**) malloc(npops*sizeof(double*));
  for(j=0; j<npops; j++) {
    powLRs[j] = (double*) malloc(npops*sizeof(double));
  }

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

  //malloc sampledRefHaps[population][indiv][STR]
  sampledRefHaps = (int***) malloc(nYpops*sizeof(int**));
  for(i=0; i<nYpops; i++) {
    sampledRefHaps[i] = (int**) malloc(103*sizeof(int*));
    for(j=0; j<103; j++) {
      sampledRefHaps[i][j] = (int*) malloc(nYSTRs*sizeof(int));
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
  Yindivs = (int*) malloc(2*sizeof(int));

  //read in allele frequencies
  readAlleleFreqs(path, STRs, nalleles, npops, alleleFreqs, samplesizes);

  //ALERT: subsample AFs for a try
  //print2dDouArr(alleleFreqs[1], nalleles[1], npops);
  //sampleAlleleFreqs(nalleles, npops, alleleFreqs, samplesizes);
  //print2dDouArr(alleleFreqs[1], nalleles[1], npops);

  //read in reference haplotypes
  readRefHaps(Ypath, nYpops, refHaps, Ysamplesizes, nYSTRs, Ypopnames);

  //ALERT: doing subsample of Yhaps, also added 'sub' to all file names to distinguish them
  //print2dIntArr(refHaps[1], Ysamplesizes[1], nYSTRs);
  sampleYHapFreqs(refHaps, nYpops, Ysamplesizes, nYSTRs, sampledRefHaps, 103);
  for(i=0; i<nYpops; i++) {
    Ysamplesizes[i] = 103;
  }
  refHaps = sampledRefHaps;
  //print2dIntArr(refHaps[0], 103, nYSTRs);


  //estimate theta for Y chr simulations
  estimateTheta(refHaps, nYpops, Ysamplesizes, nYSTRs, estYThetas);
  //print2dDouArr(estYThetas, nYpops, nYpops);

  //get Y hap frequencies
  calcHapFreqs(refHaps, nYpops, Ysamplesizes, nYSTRs, YhapFreqs);
  //print2dDouArr(YhapFreqs, nYpops, Ysamplesizes);

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

    //run fpr simulations
    for(k=0; k<niter; k++) {
      //run through each population, making rand indivs and calcing LR for a match  
      for(i=0; i<npops; i++) {
	//make two random individuals assuming memberships in this population
	ks2[0]=1.0; ks2[1]=0.0; ks2[2]=0.0;
	thetaRandIndivPair(alleleFreqs, 13, nalleles, i, simtheta, indivs, ks2);
	indiv1 = indivs[0];
	indiv2 = indivs[1];

	LR(indiv1, indiv2, alleleFreqs, npops, nalleles, samplesizes, tmparr1, calctheta, ks);
	for(l=0; l<npops; l++) 
	  LRs[l][i] = tmparr1[l];

	//make two unrel individuals Y haps assuming memberships in this population
	//thetaRandIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, estYThetas[popmap[i]][popmap[i]], Yindivs, 1.0, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	randIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, Yindivs, 1.0, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);

	YLR(Yindivs[0], Yindivs[1], modYhapFreqs, nYpops, modYsamplesizes, tmparr1);
	for(l=0; l<npops; l++)
	  LRs[l][i] *= tmparr1[popmap[l]];
      }

      //write results to files    
      sprintf(filename, "results/autYsim/%s/fpr/%s_%diter_simtheta%1.3f_calctheta%1.3f/subLRs_iter%d.res", path, relationship, niter, simtheta, calctheta, k);
      write2dDouArr(LRs, npops, npops, filename);
    }
    
    //run p-o simulations
    for(k=0; k<niter; k++) {
      //run through each population, making relatives and calcing the LR for a match  
      for(i=0; i<npops; i++) {
	ks2[0]=0.0; ks2[1]=1.0; ks2[2]=0.0;
	thetaRandIndivPair(alleleFreqs, nSTRs, nalleles, i, simtheta, indivs, ks2);
	indiv1 = indivs[0];
	indiv2 = indivs[1];
	LR(indiv1, indiv2, alleleFreqs, npops, nalleles, samplesizes, tmparr1, calctheta, ks);
	for(l=0; l<npops; l++)
	  powLRs[l][i] = tmparr1[l];
	
	//make two unrel individuals Y haps assuming memberships in this population
	thetaRandIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, estYThetas[popmap[i]][popmap[i]], Yindivs, 0.0, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	YLR(Yindivs[0], Yindivs[1], modYhapFreqs, nYpops, modYsamplesizes, tmparr1);
	for(l=0; l<npops; l++)
	  powLRs[l][i] *= tmparr1[popmap[l]];
      }

      //write results to files    
      sprintf(filename, "results/autYsim/%s/p-o/%s_%diter_simtheta%1.3f_calctheta%1.3f/subLRs_iter%d.res", path, relationship, niter, simtheta, calctheta, k);
      write2dDouArr(powLRs, npops, npops, filename);
    }

    //run sibling simulations
    for(k=0; k<niter; k++) {
      //run through each population, making relatives and calcing the LR for a match  
      for(i=0; i<npops; i++) {
	ks2[0]=0.25; ks2[1]=0.5; ks2[2]=0.25;
	thetaRandIndivPair(alleleFreqs, nSTRs, nalleles, i, simtheta, indivs, ks2);
	indiv1 = indivs[0];
	indiv2 = indivs[1];
	LR(indiv1, indiv2, alleleFreqs, npops, nalleles, samplesizes, tmparr1, calctheta, ks);
	for(l=0; l<npops; l++)
	  powLRs[l][i] = tmparr1[l];
	
	//make two unrel individuals Y haps assuming memberships in this population
	thetaRandIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, estYThetas[popmap[i]][popmap[i]], Yindivs, 0.0, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	YLR(Yindivs[0], Yindivs[1], modYhapFreqs, nYpops, modYsamplesizes, tmparr1);
	for(l=0; l<npops; l++)
	  powLRs[l][i] *= tmparr1[popmap[l]];
      }

      //write results to files    
      sprintf(filename, "results/autYsim/%s/sib/%s_%diter_simtheta%1.3f_calctheta%1.3f/subLRs_iter%d.res", path, relationship, niter, simtheta, calctheta, k);
      write2dDouArr(powLRs, npops, npops, filename);
    }


    //run cousin simulations
    for(k=0; k<niter; k++) {
      //run through each population, making relatives and calcing the LR for a match  
      for(i=0; i<npops; i++) {
	ks2[0]=.75; ks2[1]=.25; ks2[2]=0.0;
	thetaRandIndivPair(alleleFreqs, nSTRs, nalleles, i, simtheta, indivs, ks2);
	indiv1 = indivs[0];
	indiv2 = indivs[1];
	LR(indiv1, indiv2, alleleFreqs, npops, nalleles, samplesizes, tmparr1, calctheta, ks);
	for(l=0; l<npops; l++)
	  powLRs[l][i] = tmparr1[l];
	
	//make two unrel individuals Y haps assuming memberships in this population
	thetaRandIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, estYThetas[popmap[i]][popmap[i]], Yindivs, 0.0, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	YLR(Yindivs[0], Yindivs[1], modYhapFreqs, nYpops, modYsamplesizes, tmparr1);
	for(l=0; l<npops; l++)
	  powLRs[l][i] *= tmparr1[popmap[l]];
      }

      //write results to files    
      sprintf(filename, "results/autYsim/%s/cousin/%s_%diter_simtheta%1.3f_calctheta%1.3f/subLRs_iter%d.res", path, relationship, niter, simtheta, calctheta, k);
      write2dDouArr(powLRs, npops, npops, filename);
    }

    //run half-cousin simulations
    for(k=0; k<niter; k++) {
      //run through each population, making relatives and calcing the LR for a match  
      for(i=0; i<npops; i++) {
	ks2[0]=.875; ks2[1]=.125; ks2[2]=0.0;
	thetaRandIndivPair(alleleFreqs, nSTRs, nalleles, i, simtheta, indivs, ks2);
	indiv1 = indivs[0];
	indiv2 = indivs[1];
	LR(indiv1, indiv2, alleleFreqs, npops, nalleles, samplesizes, tmparr1, calctheta, ks);
	for(l=0; l<npops; l++)
	  powLRs[l][i] = tmparr1[l];
	
	//make two unrel individuals Y haps assuming memberships in this population
	thetaRandIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, estYThetas[popmap[i]][popmap[i]], Yindivs, 0.0, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	YLR(Yindivs[0], Yindivs[1], modYhapFreqs, nYpops, modYsamplesizes, tmparr1);
	for(l=0; l<npops; l++)
	  powLRs[l][i] *= tmparr1[popmap[l]];
      }

      //write results to files    
      sprintf(filename, "results/autYsim/%s/halfcousin/%s_%diter_simtheta%1.3f_calctheta%1.3f/subLRs_iter%d.res", path, relationship, niter, simtheta, calctheta, k);
      write2dDouArr(powLRs, npops, npops, filename);
    }

    //run second cousin simulations
    for(k=0; k<niter; k++) {
      //run through each population, making relatives and calcing the LR for a match  
      for(i=0; i<npops; i++) {
	ks2[0]=.9375; ks2[1]=.0625; ks2[2]=0.0;
	thetaRandIndivPair(alleleFreqs, nSTRs, nalleles, i, simtheta, indivs, ks2);
	indiv1 = indivs[0];
	indiv2 = indivs[1];
	LR(indiv1, indiv2, alleleFreqs, npops, nalleles, samplesizes, tmparr1, calctheta, ks);
	for(l=0; l<npops; l++)
	  powLRs[l][i] = tmparr1[l];
	
	//make two unrel individuals Y haps assuming memberships in this population
	thetaRandIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, estYThetas[popmap[i]][popmap[i]], Yindivs, 0.0, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	YLR(Yindivs[0], Yindivs[1], modYhapFreqs, nYpops, modYsamplesizes, tmparr1);
	for(l=0; l<npops; l++) 
	  powLRs[l][i] *= tmparr1[popmap[l]];
      }

      //write results to files    
      sprintf(filename, "results/autYsim/%s/2cousin/%s_%diter_simtheta%1.3f_calctheta%1.3f/subLRs_iter%d.res", path, relationship, niter, simtheta, calctheta, k);
      write2dDouArr(powLRs, npops, npops, filename);
    }

    //run half-sib simulations
    for(k=0; k<niter; k++) {
      //run through each population, making relatives and calcing the LR for a match  
      for(i=0; i<npops; i++) {
	ks2[0]=.5; ks2[1]=.5; ks2[2]=0.0;
	thetaRandIndivPair(alleleFreqs, nSTRs, nalleles, i, simtheta, indivs, ks2);
	indiv1 = indivs[0];
	indiv2 = indivs[1];
	LR(indiv1, indiv2, alleleFreqs, npops, nalleles, samplesizes, tmparr1, calctheta, ks);
	for(l=0; l<npops; l++)
	  powLRs[l][i] = tmparr1[l];

	//make two unrel individuals Y haps assuming memberships in this population
	thetaRandIndivPairYchr(YhapFreqs[popmap[i]], modYhapFreqs[popmap[i]], totYindivs, estYThetas[popmap[i]][popmap[i]], Yindivs, 0.0, Ysamplesizes, modYsamplesizes, popmap[i], nYpops);
	YLR(Yindivs[0], Yindivs[1], modYhapFreqs, nYpops, modYsamplesizes, tmparr1);
	for(l=0; l<npops; l++)
	  powLRs[l][i] *= tmparr1[popmap[l]];
      }

      //write results to files    
      sprintf(filename, "results/autYsim/%s/halfsib/%s_%diter_simtheta%1.3f_calctheta%1.3f/subLRs_iter%d.res", path, relationship, niter, simtheta, calctheta, k);
      write2dDouArr(powLRs, npops, npops, filename);
    }
    
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

  for(i=0; i<2; i++) {
    for(j=0; j<nSTRs; j++) 
      free(indivs[i][j]);
    free(indivs[i]);
  }
  free(indivs);

  free(ks);
  free(pmatches);
  free(tmparr1);
  free(tmparr2);

  for(j=0; j<npops; j++) {
    free(LRs[j]);
  }
  free(LRs);

  for(j=0; j<npops; j++) {
    free(powLRs[j]);
  }
  free(powLRs);

  return 0;
}
