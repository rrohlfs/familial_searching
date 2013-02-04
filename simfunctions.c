// some functions to simulate individuals 

#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include "myio.h"
#include "mystat.h"
#include "probterms.h"

//returns an index for a random value from a the given distribution, note, dist is out of 100
int randValFromDist(double *dist, int n) {
  int k;
  double curVal, tot;

  curVal = 100.0 * (random()/((double)RAND_MAX+1.0));

  //same case if dist sums to 1.0
  tot = 0.0;
  for(k=0; k<n; k++) 
    tot += dist[k];
  if(tot == 1.0) {
    tot=0.0; k=0;
    while(k<n && tot<=curVal) {
      tot += dist[k]*100.0;
      k++;
    }
    return (k-1);
  }

  // if dist sums to 100
  tot = 0.0;
  k=0;
  while(k<n && tot<=curVal) {
    tot += dist[k];
    k++;
  }
  return (k-1);
}

//returns an index for a random value from a the given distribution, but dist is proportions, not percents
int randValFromDistProp(double *dist, int n) {
  int k;
  double curVal, tot;

  curVal = (random()/((double)RAND_MAX+1.0));

  tot=0.0;
  k=0.0;
  while(k<n && tot<=curVal) {
    tot += dist[k];
    k++;
  }
  return (k-1);
}


//returns nvals distinct index for random values from a the given distribution, note, dist is out of 100
void randValsFromDist(double *dist, int n, int *vals, int nvals) {
  switch(nvals) {
  case 1:
    vals[0] = randValFromDist(dist, n);
    break;
  case 2:
    do {
      vals[0] = randValFromDist(dist, n);
      vals[1] = randValFromDist(dist, n);
    } while(vals[0]==vals[1]);
    shuffle(vals, nvals);
    break;
  case 3:
    do {
      vals[0] = randValFromDist(dist, n);
      vals[1] = randValFromDist(dist, n);
      vals[2] = randValFromDist(dist, n);
    } while(vals[0]==vals[1] || vals[0]==vals[2] || vals[1]==vals[2]);
    shuffle(vals, nvals);
    break;
  case 4:
    do {
      vals[0] = randValFromDist(dist, n);
      vals[1] = randValFromDist(dist, n);
      vals[2] = randValFromDist(dist, n);
      vals[3] = randValFromDist(dist, n);
    } while(vals[0]==vals[1] || vals[0]==vals[2] || vals[0]==vals[3] || vals[1]==vals[2] || vals[1]==vals[3] || vals[2]==vals[3]);
    shuffle(vals, nvals);
    break;
  }
}

//use Dirichlet model to modify AFs after having drawn some
void modifyAFdist(double *AFs, int nalleles, int *obsAs, int nobsAs, double *newAFs, double theta) {
  int i, j, ni;
  
  for(i=0; i<nalleles; i++) {
    //get ni (number of this allele drawn already)
    ni = 0; 
    for(j=0; j<nobsAs; j++) {
      if(obsAs[j] == i) 
	ni++;
    }

    //calculate new AF for allele i
    newAFs[i] = 100*(ni*theta + (1.0-theta)*AFs[i]/100)/(1.0 + (nobsAs-1)*theta);
   }
}




//returns one individual generated randomly from population allele freq distribution assuming theta
// note: indiv is the allele indicies, not the allele names
void thetaRandIndiv(double ***AFs, int nSTRs, int *nalleles, int npop, double theta, int **indiv) {
  int i, j;
  int nA;
  int *obsAs;
  double *modAFs, *curAFs;

  //malloc obsAs
  obsAs = (int*) malloc(1*sizeof(int));

  //run through each STR independently
  for(i=0; i<nSTRs; i++) {
    nA = nalleles[i];

    //malloc curAFs and modAFs
    modAFs = (double*) malloc(nA * sizeof(double));
    curAFs = (double*) malloc(nA * sizeof(double));
    for(j=0; j<nA; j++) {
      curAFs[j] = AFs[i][j][npop];
    }
    
    //decide first allele
    indiv[i][0] = randValFromDist(curAFs, nA);

    //modify AF distribution
    obsAs[0] = indiv[i][0];
    modifyAFdist(curAFs, nA, obsAs, 1, modAFs, theta);

    //decide second allele
    indiv[i][1] = randValFromDist(modAFs, nA);

    //free what needs freed
    free(modAFs);
    free(curAFs);
  }

  //free what needs freed
  free(obsAs);

}

//returns two individuals generated randomly from population allele freq distribution assuming theta, and some given kinship coeffs
// note: indiv is the allele indicies, not the allele names
void thetaRandIndivPair(double ***AFs, int nSTRs, int *nalleles, int npop, double theta, int ***indivs, double *ks) {
  int i, j;
  int nA, nshare;
  int *obsAs;
  double *modAFs, *curAFs;

  //malloc obsAs
  obsAs = (int*) malloc(4*sizeof(int));

  //run through each STR independently
  for(i=0; i<nSTRs; i++) {
    nA = nalleles[i];

    //malloc curAFs and modAFs
    modAFs = (double*) malloc(nA * sizeof(double));
    curAFs = (double*) malloc(nA * sizeof(double));
    for(j=0; j<nA; j++) {
      curAFs[j] = AFs[i][j][npop];
    }
    
    //decide degree of allele sharing
    nshare = randValFromDist(ks, 3);

    //decide first allele
    indivs[0][i][0] = randValFromDist(curAFs, nA);

    //modify AF distribution
    obsAs[0] = indivs[0][i][0];
    modifyAFdist(curAFs, nA, obsAs, 1, modAFs, theta);

    //decide second allele
    indivs[0][i][1] = randValFromDist(modAFs, nA);

    if(nshare == 2) {
      indivs[1][i][0] = indivs[0][i][0];
      indivs[1][i][1] = indivs[0][i][1];
    } else {
      //modify AF distribution
      obsAs[1] = indivs[0][i][1];
      modifyAFdist(curAFs, nA, obsAs, 2, modAFs, theta);
      
      //decide third allele
      indivs[1][i][0] = randValFromDist(modAFs, nA);

      if(nshare == 1) {
	if((random()/((double)RAND_MAX+1.0)) < .5)
	  indivs[1][i][1] = indivs[0][i][0];
	else
	  indivs[1][i][1] = indivs[0][i][1];
      } else { 
	//modify AF distribution    
	obsAs[2] = indivs[1][i][0];
	modifyAFdist(curAFs, nA, obsAs, 3, modAFs, theta);
	
	//decide fourth allele
	indivs[1][i][1] = randValFromDist(modAFs, nA);
      }
    }

    //free what needs freed
    free(modAFs);
    free(curAFs);
  }

  //free what needs freed
  free(obsAs);
}

//returns two unrelated individuals generated randomly from population allele freq distribution assuming theta
// note: indiv is the allele indicies, not the allele names
void thetaRandUnrelIndivPair(double ***AFs, int nSTRs, int *nalleles, int npop, double theta, int **indiv1, int **indiv2) {
  int i, j;
  int nA;
  int *obsAs;
  double *modAFs, *curAFs;

  //malloc obsAs
  obsAs = (int*) malloc(4*sizeof(int));

  //run through each STR independently
  for(i=0; i<nSTRs; i++) {
    nA = nalleles[i];

    //malloc curAFs and modAFs
    modAFs = (double*) malloc(nA * sizeof(double));
    curAFs = (double*) malloc(nA * sizeof(double));
    for(j=0; j<nA; j++) {
      curAFs[j] = AFs[i][j][npop];
    }

    //decide first allele
    indiv1[i][0] = randValFromDist(curAFs, nA);

    //modify AF distribution
    obsAs[0] = indiv1[i][0];
    modifyAFdist(curAFs, nA, obsAs, 1, modAFs, theta);

    //decide second allele
    indiv1[i][1] = randValFromDist(modAFs, nA);

    //modify AF distribution
    obsAs[1] = indiv1[i][1];
    modifyAFdist(curAFs, nA, obsAs, 2, modAFs, theta);
      
    //decide third allele
    indiv2[i][0] = randValFromDist(modAFs, nA);

    //modify AF distribution    
    obsAs[2] = indiv2[i][0];
    modifyAFdist(curAFs, nA, obsAs, 3, modAFs, theta);
    
    //decide fourth allele
    indiv2[i][1] = randValFromDist(modAFs, nA);

    //free what needs freed
    free(modAFs);
    free(curAFs);
  }

  //free what needs freed
  free(obsAs);
}

//returns two unrelated individuals generated randomly from population allele freq distribution assuming theta, given first indiv
// note: indiv is the allele indicies, not the allele names
void thetaRandUnrelSecondIndiv(double ***AFs, int nSTRs, int *nalleles, int npop, double theta, int **indiv1, int **indiv2) {
  int i, j;
  int nA;
  int *obsAs;
  double *modAFs, *curAFs;

  //malloc obsAs
  obsAs = (int*) malloc(4*sizeof(int));

  //run through each STR independently
  for(i=0; i<nSTRs; i++) {
    nA = nalleles[i];

    //malloc curAFs and modAFs
    modAFs = (double*) malloc(nA * sizeof(double));
    curAFs = (double*) malloc(nA * sizeof(double));
    for(j=0; j<nA; j++) {
      curAFs[j] = AFs[i][j][npop];
    }

    //modify AF distribution given first indiv
    obsAs[0] = indiv1[i][0];
    obsAs[1] = indiv1[i][1];
    modifyAFdist(curAFs, nA, obsAs, 2, modAFs, theta);
      
    //decide third allele
    indiv2[i][0] = randValFromDist(modAFs, nA);

    //modify AF distribution    
    obsAs[2] = indiv2[i][0];
    modifyAFdist(curAFs, nA, obsAs, 3, modAFs, theta);
    
    //decide fourth allele
    indiv2[i][1] = randValFromDist(modAFs, nA);

    //free what needs freed
    free(modAFs);
    free(curAFs);
  }

  //free what needs freed
  free(obsAs);
}

//use Dirichlet model to modify AFs after having drawn some where AFs are given as proportion, not percent
void modifyAFdistProp(double *AFs, int nalleles, int *obsAs, int nobsAs, double *newAFs, double theta) {
  int i, j, ni;
  
  for(i=0; i<nalleles; i++) {
    //get ni (number of this allele drawn already)
    ni = 0; 
    for(j=0; j<nobsAs; j++) {
      if(obsAs[j] == i) 
	ni++;
    }

    //calculate new AF for allele i
    newAFs[i] = (ni*theta + (1.0-theta)*AFs[i])/(1.0 + (nobsAs-1)*theta);
   }
}

//remove individual from hap freq distribution
void removeIndivFromHF(double *hapFreqs, double *newHapFreqs, int nalleles, int indivallele, int *nindivs, int *newnindivs, int popi, int npops) {
  int i;
  double ni;
  
  for(i=0; i<nalleles; i++) {
    ni = 1.0*hapFreqs[i]*nindivs[popi];
    if(i == indivallele) 
      newHapFreqs[i] = (ni-1.0)/(nindivs[popi]-1.0);
    else
      newHapFreqs[i] = ni/(nindivs[popi]-1.0);
  }
  
  for(i=0; i<npops; i++) {
    if(i == popi)
      newnindivs[i] = nindivs[i]-1;
    else
      newnindivs[i] = nindivs[i];
  }
}

//returns two individuals generated randomly from Y-chr ref haplotype list  assuming theta, and some given kinship coeffs
//note: indiv is the allele indicies, not the allele names
void thetaRandIndivPairYchr(double *hapFreqs, double *newHapFreqs, int nalleles, double theta, int *indivs, double k0, int* nindivs, int* modnindivs, int popi, int npops) {
  //decide first allele
  indivs[0] = randValFromDistProp(hapFreqs, nalleles);
  removeIndivFromHF(hapFreqs, newHapFreqs, nalleles, indivs[0], nindivs, modnindivs, popi, npops); 
    
  //decide second allele
  if(k0 == 0.0 || (random()/((double)RAND_MAX+1.0)) < theta) {
    indivs[1] = indivs[0];
  } else {
    indivs[1] = randValFromDistProp(newHapFreqs, nalleles);
    removeIndivFromHF(newHapFreqs, newHapFreqs, nalleles, indivs[1], nindivs, modnindivs, popi, npops); 
  }
}

void randIndivPairYchr(double *hapFreqs, double *newHapFreqs, int nalleles, int *indivs, double k0, int* nindivs, int* modnindivs, int popi, int npops) {
  //decide first allele
  indivs[0] = randValFromDistProp(hapFreqs, nalleles);
  removeIndivFromHF(hapFreqs, newHapFreqs, nalleles, indivs[0], nindivs, modnindivs, popi, npops); 
    
  //decide second allele
  indivs[1] = randValFromDistProp(newHapFreqs, nalleles);
  removeIndivFromHF(newHapFreqs, newHapFreqs, nalleles, indivs[1], nindivs, modnindivs, popi, npops); 
}

//returns two unrelated individuals generated randomly from Y-chr ref haplotype list  assuming theta
//note: indiv is the allele indicies, not the allele names
void thetaRandUnrelIndivPairYchr(double *hapFreqs, double *newHapFreqs, int nalleles, double theta, int *indiv1, int *indiv2, int* nindivs, int* modnindivs, int popi, int npops) {
  //decide first allele
  indiv1[0] = randValFromDistProp(hapFreqs, nalleles);
  removeIndivFromHF(hapFreqs, newHapFreqs, nalleles, indiv1[0], nindivs, modnindivs, popi, npops);
    
  //decide second allele
  if((random()/((double)RAND_MAX+1.0)) < theta) {
    indiv2[0] = indiv1[0];
  } else {
    indiv2[0] = randValFromDistProp(newHapFreqs, nalleles);
    removeIndivFromHF(newHapFreqs, newHapFreqs, nalleles, indiv2[0], nindivs, modnindivs, popi, npops);
  }
}

//returns two unrelated individuals generated randomly from Y-chr ref haplotype list
//note: indiv is the allele indicies, not the allele names
void randUnrelIndivPairYchr(double *hapFreqs, double *newHapFreqs, int nalleles, int *indiv1, int *indiv2, int* nindivs, int* modnindivs, int popi, int npops) {
  //decide first allele
  indiv1[0] = randValFromDistProp(hapFreqs, nalleles);
  removeIndivFromHF(hapFreqs, newHapFreqs, nalleles, indiv1[0], nindivs, modnindivs, popi, npops);

  //decide second allele
  indiv2[0] = randValFromDistProp(newHapFreqs, nalleles);
  removeIndivFromHF(newHapFreqs, newHapFreqs, nalleles, indiv2[0], nindivs, modnindivs, popi, npops);
}



//returns two unrelated individuals generated randomly from Y-chr ref haplotype list  assuming theta, given first indiv
//note: indiv is the allele indicies, not the allele names
void thetaRandUnrelSecondIndivYchr(double *hapFreqs, double *newHapFreqs, int nalleles, double theta, int *indiv1, int *indiv2, int* nindivs, int* modnindivs, int popi, int npops) {

  //decide second allele
  if((random()/((double)RAND_MAX+1.0)) < theta) {
    indiv2[0] = indiv1[0];
  } else {
    indiv2[1] = randValFromDistProp(newHapFreqs, nalleles);
    removeIndivFromHF(hapFreqs, newHapFreqs, nalleles, indiv2[0], nindivs, modnindivs, popi, npops);
  }
}




//returns one individual generated randomly from Y-chr ref haplotype list  assuming theta, and some given kinship coeffs
//note: indiv is the allele indicies, not the allele names
void thetaRandIndivYchr(double *hapFreqs, double *newHapFreqs, int nalleles, int *indiv, int* nindivs, int* modnindivs, int popi, int npops) {
  indiv[0] = randValFromDistProp(hapFreqs, nalleles);

  removeIndivFromHF(hapFreqs, newHapFreqs, nalleles, indiv[0], nindivs, modnindivs, popi, npops); 
}

