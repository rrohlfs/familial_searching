// some functions for probabilities of genotypes given relatedness and derivatives of those same probabilities

#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include "myio.h"

//add the two genotypes to the AFs in new AFs
void addIndivsToAFs(int **geno1, int **geno2, double ***AFs, double ***newAFs, int npops, int *nalleles, int **npeople) {
  int i, j, k;

  //copy AFs values to newAFs and multiply to get num alleles (instead of freqs)
  for(i=0; i<13; i++)
    for(j=0; j<nalleles[i]; j++) 
      for(k=0; k<npops; k++) 
	newAFs[i][j][k] = AFs[i][j][k]*2.0*npeople[i][k];
  
  //add genos
  for(i=0; i<13; i++) {
    for(k=0; k<npops; k++) {
      newAFs[i][geno1[i][0]][k] += 1.0;
      newAFs[i][geno1[i][1]][k] += 1.0;
      newAFs[i][geno2[i][0]][k] += 1.0;
      newAFs[i][geno2[i][1]][k] += 1.0;
    }
  }
  
  //divide out newAFs to get freqs again
  for(i=0; i<13; i++)
    for(k=0; k<npops; k++)    
      for(j=0; j<nalleles[i]; j++) 
	newAFs[i][j][k] = newAFs[i][j][k] / (2.0*npeople[i][k] + 4.0);
}

//calcs the probability of these two genotypes given k0,k1,k2 and assuming membership in the same (and each) population group
void pGenosGivenK(int **geno1, int **geno2, double ***AFs, int npops, double *pmatches, double theta, double *ks) {
  int i, k;
  int A11, A12, A21, A22;
  double pA11, pA12, pA21, pA22;

  //init pmatches
  for(i=0; i<npops; i++)
    pmatches[i] = 1.0;

  //run through all loci and pop groups, producting probabilities
  for(k=0; k<npops; k++) {
    for(i=0; i<13; i++) {

      //printDouArr(pmatches, npops);

     //get alleles
      A11 = geno1[i][0];
      A12 = geno1[i][1];
      A21 = geno2[i][0];
      A22 = geno2[i][1];
      
      //get allele frequencies
      pA11 = AFs[i][A11][k]/100.0;
      pA12 = AFs[i][A12][k]/100.0;
      pA21 = AFs[i][A21][k]/100.0;
      pA22 = AFs[i][A22][k]/100.0;
      
      //calc probability of these genotypes given relationship for different genotype cases
      //if they're homo homo with both shared
      if((A11 == A12) && (A12 == A21) && (A21 == A22)) { 
	pmatches[k] *= (ks[2]*(pA11*(theta+(1.0-theta)*pA11)) 
			+ ks[1]*(pA11*(theta+(1.0-theta)*pA11))*(2.0*theta+(1.0-theta)*pA11)/(1.0+theta)
			+ ks[0]*(pA11*(theta+(1.0-theta)*pA11))*(2.0*theta+(1.0-theta)*pA11)*(3.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(1.0+2.0*theta)));
      }
      //if they're homo homo with none shared
      else if((A11 == A12) && (A21 == A22)) { 
	pmatches[k] *= ks[0]*pA11*pA21*(1.0-theta)*(theta+(1.0+theta)*pA11)*(theta+(1.0+theta)*pA21)/((theta+1.0)*(2.0*theta+1.0));
      }
      //if they're homo het with one shared
      else if((A11 == A12) && (A12 == A21)) { 
	pmatches[k] *= (ks[1]*pA11*pA22*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
			+ ks[0]*2.0*pA11*pA22*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if((A11 == A12) && (A12 == A22)) { 
	pmatches[k] *= (ks[1]*pA11*pA21*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
			+ ks[0]*2.0*pA11*pA21*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if((A21 == A22) && (A21 == A11)) { 
	pmatches[k] *= (ks[1]*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
			+ ks[0]*2.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if((A21 == A22) && (A21 == A12)) { 
	pmatches[k] *= (ks[1]*pA12*pA11*(1.0-theta)*(theta+(1.0-theta)*pA12)/(1.0+theta)
			+ ks[0]*2.0*pA12*pA11*(1.0-theta)*(theta+(1.0-theta)*pA12)*(2.0*theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      //if they're het het with both shared
      else if((A11 == A21) && (A12 == A22)) { 
	pmatches[k] *= (ks[2]*2.0*pA11*pA12*(1.0-theta)
			+ ks[1]*pA11*pA12*(1.0-theta)*(2.0*theta+(1.0-theta)*(pA11+pA12))/(1.0+theta)
			+ ks[0]*4.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if((A11 == A22) && (A12 == A21)) { 
	pmatches[k] *= (ks[2]*2.0*pA11*pA12*(1.0-theta)
			+ ks[1]*pA11*pA12*(1.0-theta)*(2.0*theta+(1.0-theta)*(pA11+pA12))/(1.0+theta)
			+ ks[0]*4.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      //if they're het het with one shared
      else if(A11 == A21) { 
	pmatches[k] *= (ks[1]*pA11*pA12*pA22*pow((1.0-theta),2.0)/(1.0+theta)
			+ ks[0]*4.0*pA11*pA12*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if(A11 == A22) { 
	pmatches[k] *= (ks[1]*pA11*pA12*pA21*pow((1.0-theta),2.0)/(1.0+theta)
			+ ks[0]*4.0*pA11*pA12*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if(A12 == A21) { 
	pmatches[k] *= (ks[1]*pA12*pA11*pA22*pow((1.0-theta),2.0)/(1.0+theta)
			+ ks[0]*4.0*pA12*pA11*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if(A12 == A22) { 
	pmatches[k] *= (ks[1]*pA12*pA11*pA21*pow((1.0-theta),2.0)/(1.0+theta)
			+ ks[0]*4.0*pA12*pA11*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      //if they're homo het with none shared
      else if(A11 == A12) { 
 	pmatches[k] *= ks[0]*2.0*pA11*pA21*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0));
      }
      else if(A21 == A22) { 
 	pmatches[k] *= ks[0]*2.0*pA11*pA12*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA21)/((1.0+theta)*(2.0*theta+1.0));
      }
      //if they're het het with none shared
      else {
	pmatches[k] *= ks[0]*pA11*pA12*pA21*pA22*pow((1.0-theta),3.0)/((1.0+theta)*(2.0*theta+1.0));
      }      
    }
  }
}



//calcs the log probability of these two genotypes given k0,k1,k2 and assuming membership in the same (and each) population group
void lnpGenosGivenK(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *pmatches, double theta, double *ks) {
  int i, k;
  int A11, A12, A21, A22;
  double pA11, pA12, pA21, pA22;

  //init pmatches
  for(i=0; i<npops; i++)
    pmatches[i] = 0.0;

  //run through all loci and pop groups, producting probabilities
  for(k=0; k<npops; k++) {
    for(i=0; i<13; i++) {

     //get alleles
      A11 = geno1[i][0];
      A12 = geno1[i][1];
      A21 = geno2[i][0];
      A22 = geno2[i][1];
      
      //get allele frequencies
      pA11 = AFs[i][A11][k]/100.0;
      pA12 = AFs[i][A12][k]/100.0;
      pA21 = AFs[i][A21][k]/100.0;
      pA22 = AFs[i][A22][k]/100.0;

      //calc probability of these genotypes given relationship for different genotype cases
      //if they're homo homo with both shared
      if((A11 == A12) && (A12 == A21) && (A21 == A22)) { 
	pmatches[k] += log(ks[2]*(pA11*(theta+(1.0-theta)*pA11)) 
			+ ks[1]*(pA11*(theta+(1.0-theta)*pA11))*(2.0*theta+(1.0-theta)*pA11)/(1.0+theta)
			+ ks[0]*(pA11*(theta+(1.0-theta)*pA11))*(2.0*theta+(1.0-theta)*pA11)*(3.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(1.0+2.0*theta)));
      }
      //if they're homo homo with none shared
      else if((A11 == A12) && (A21 == A22)) { 
	pmatches[k] += log(ks[0]*pA11*pA21*(1.0-theta)*(theta+(1.0+theta)*pA11)*(theta+(1.0+theta)*pA21)/((theta+1.0)*(2.0*theta+1.0)));
      }
      //if they're homo het with one shared
      else if((A11 == A12) && (A12 == A21)) { 
	pmatches[k] += log(ks[1]*pA11*pA22*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
			+ ks[0]*2.0*pA11*pA22*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if((A11 == A12) && (A12 == A22)) { 
	pmatches[k] += log(ks[1]*pA11*pA21*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
			+ ks[0]*2.0*pA11*pA21*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if((A21 == A22) && (A21 == A11)) { 
	pmatches[k] += log(ks[1]*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
			+ ks[0]*2.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if((A21 == A22) && (A21 == A12)) { 
	pmatches[k] += log(ks[1]*pA12*pA11*(1.0-theta)*(theta+(1.0-theta)*pA12)/(1.0+theta)
			+ ks[0]*2.0*pA12*pA11*(1.0-theta)*(theta+(1.0-theta)*pA12)*(2.0*theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      //if they're het het with both shared
      else if((A11 == A21) && (A12 == A22)) { 
	pmatches[k] += log(ks[2]*2.0*pA11*pA12*(1.0-theta)
			+ ks[1]*pA11*pA12*(1.0-theta)*(2.0*theta+(1.0-theta)*(pA11+pA12))/(1.0+theta)
			+ ks[0]*4.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if((A11 == A22) && (A12 == A21)) { 
	pmatches[k] += log(ks[2]*2.0*pA11*pA12*(1.0-theta)
			+ ks[1]*pA11*pA12*(1.0-theta)*(2.0*theta+(1.0-theta)*(pA11+pA12))/(1.0+theta)
			+ ks[0]*4.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      //if they're het het with one shared
      else if(A11 == A21) { 
	pmatches[k] += log(ks[1]*pA11*pA12*pA22*pow((1.0-theta),2.0)/(1.0+theta)
			+ ks[0]*4.0*pA11*pA12*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if(A11 == A22) { 
	pmatches[k] += log(ks[1]*pA11*pA12*pA21*pow((1.0-theta),2.0)/(1.0+theta)
			+ ks[0]*4.0*pA11*pA12*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if(A12 == A21) { 
	pmatches[k] += log(ks[1]*pA12*pA11*pA22*pow((1.0-theta),2.0)/(1.0+theta)
			+ ks[0]*4.0*pA12*pA11*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if(A12 == A22) { 
	pmatches[k] += log(ks[1]*pA12*pA11*pA21*pow((1.0-theta),2.0)/(1.0+theta)
			+ ks[0]*4.0*pA12*pA11*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
      }
      //if they're homo het with none shared
      else if(A11 == A12) { 
 	pmatches[k] += log(ks[0]*2.0*pA11*pA21*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
      }
      else if(A21 == A22) { 
 	pmatches[k] += log(ks[0]*2.0*pA11*pA12*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA21)/((1.0+theta)*(2.0*theta+1.0)));
      }
      //if they're het het with none shared
      else {
	pmatches[k] += log(ks[0]*pA11*pA12*pA21*pA22*pow((1.0-theta),3.0)/((1.0+theta)*(2.0*theta+1.0)));
      }
    }
  }
}






//calcs the probability of these two single-locus genotypes given k0,k1,k2 and assuming membership in the same (and each) population group
void pGenoGivenK(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *pmatches, double theta, double *ks) {
  int k;
  int A11, A12, A21, A22;
  double pA11, pA12, pA21, pA22;

  //init pmatches
  for(k=0; k<npops; k++)
    pmatches[k] = 1.0;

  //run through all loci and pop groups, producting probabilities
  for(k=0; k<npops; k++) { 
    //get alleles
    A11 = geno1[locusnum][0];
    A12 = geno1[locusnum][1];
    A21 = geno2[locusnum][0];
    A22 = geno2[locusnum][1];
    
    //get allele frequencies
    pA11 = AFs[locusnum][A11][k]/100.0;
    pA12 = AFs[locusnum][A12][k]/100.0;
    pA21 = AFs[locusnum][A21][k]/100.0;
    pA22 = AFs[locusnum][A22][k]/100.0;

    //calc probability of these genotypes given relationship for different genotype cases
    //if they're homo homo with both shared
    if((A11 == A12) && (A12 == A21) && (A21 == A22)) { 
      pmatches[k] *= (ks[2]*(pA11*(theta+(1.0-theta)*pA11)) 
		      + ks[1]*(pA11*(theta+(1.0-theta)*pA11))*(2.0*theta+(1.0-theta)*pA11)/(1.0+theta)
		      + ks[0]*(pA11*(theta+(1.0-theta)*pA11))*(2.0*theta+(1.0-theta)*pA11)*(3.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(1.0+2.0*theta)));
    }
    //if they're homo homo with none shared
    else if((A11 == A12) && (A21 == A22)) { 
      pmatches[k] *= ks[0]*pA11*pA21*(1.0-theta)*(theta+(1.0+theta)*pA11)*(theta+(1.0+theta)*pA21)/((theta+1.0)*(2.0*theta+1.0));
    }
    //if they're homo het with one shared
    else if((A11 == A12) && (A12 == A21)) { 
      pmatches[k] *= (ks[1]*pA11*pA22*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
		      + ks[0]*2.0*pA11*pA22*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
    }
    else if((A11 == A12) && (A12 == A22)) { 
      pmatches[k] *= (ks[1]*pA11*pA21*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
		      + ks[0]*2.0*pA11*pA21*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
    }
    else if((A21 == A22) && (A21 == A11)) { 
      pmatches[k] *= (ks[1]*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)/(1.0+theta)
		      + ks[0]*2.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
    }
    else if((A21 == A22) && (A21 == A12)) { 
      pmatches[k] *= (ks[1]*pA12*pA11*(1.0-theta)*(theta+(1.0-theta)*pA12)/(1.0+theta)
		      + ks[0]*2.0*pA12*pA11*(1.0-theta)*(theta+(1.0-theta)*pA12)*(2.0*theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
    }
    //if they're het het with both shared
    else if((A11 == A21) && (A12 == A22)) { 
      pmatches[k] *= (ks[2]*2.0*pA11*pA12*(1.0-theta)
		      + ks[1]*pA11*pA12*(1.0-theta)*(2.0*theta+(1.0-theta)*(pA11+pA12))/(1.0+theta)
		      + ks[0]*4.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
    }
    else if((A11 == A22) && (A12 == A21)) { 
      pmatches[k] *= (ks[2]*2.0*pA11*pA12*(1.0-theta)
		      + ks[1]*pA11*pA12*(1.0-theta)*(2.0*theta+(1.0-theta)*(pA11+pA12))/(1.0+theta)
		      + ks[0]*4.0*pA11*pA12*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
    }
    //if they're het het with one shared
    else if(A11 == A21) { 
      pmatches[k] *= (ks[1]*pA11*pA12*pA22*pow((1.0-theta),2.0)/(1.0+theta)
		      + ks[0]*4.0*pA11*pA12*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
    }
    else if(A11 == A22) { 
      pmatches[k] *= (ks[1]*pA11*pA12*pA21*pow((1.0-theta),2.0)/(1.0+theta)
		      + ks[0]*4.0*pA11*pA12*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0)));
    }
    else if(A12 == A21) { 
      pmatches[k] *= (ks[1]*pA12*pA11*pA22*pow((1.0-theta),2.0)/(1.0+theta)
		      + ks[0]*4.0*pA12*pA11*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
    }
    else if(A12 == A22) { 
      pmatches[k] *= (ks[1]*pA12*pA11*pA21*pow((1.0-theta),2.0)/(1.0+theta)
		      + ks[0]*4.0*pA12*pA11*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((1.0+theta)*(2.0*theta+1.0)));
    }
    //if they're homo het with none shared
    else if(A11 == A12) { 
      pmatches[k] *= ks[0]*2.0*pA11*pA21*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((1.0+theta)*(2.0*theta+1.0));
    }
    else if(A21 == A22) { 
      pmatches[k] *= ks[0]*2.0*pA11*pA12*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA21)/((1.0+theta)*(2.0*theta+1.0));
    }
    //if they're het het with none shared
    else {
      pmatches[k] *= ks[0]*pA11*pA12*pA21*pA22*pow((1.0-theta),3.0)/((1.0+theta)*(2.0*theta+1.0));
    }
  }

}







//calcs the probability of these two single-locus genotypes given k0,k1,k2 and assuming membership in the same (and each) population group
void deltapGenoGivenK(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *deltaps, double theta, double *ks, int allelenum) {
  int k;
  int A11, A12, A21, A22;
  double pA11, pA12, pA21, pA22;

  //init pmatches
  for(k=0; k<npops; k++)
    deltaps[k] = 1.0;

  //run through all loci and pop groups, producting probabilities
  for(k=0; k<npops; k++) {
    //get alleles
    A11 = geno1[locusnum][0];
    A12 = geno1[locusnum][1];
    A21 = geno2[locusnum][0];
    A22 = geno2[locusnum][1];
    
    //get allele frequencies
    pA11 = AFs[locusnum][A11][k]/100.0;
    pA12 = AFs[locusnum][A12][k]/100.0;
    pA21 = AFs[locusnum][A21][k]/100.0;
    pA22 = AFs[locusnum][A22][k]/100.0;

    //calc probability of these genotypes given relationship for different genotype cases
    //if they're homo homo with both shared
    if((A11 == A12) && (A12 == A21) && (A21 == A22)) { 
      deltaps[k] = ks[2]*(theta + 2.0*pA11*(1.0-theta)) 
	+ ks[1]*(pA11*(theta + pA11*(1.0-theta))*(1.0-theta)/(1+theta)+(2.0*theta+(1.0-theta)*pA11)*(theta+2.0*pA11*(1.0-theta))/(1.0+theta))
	+ ks[0]*((pA11*(theta+(1.0-theta)*pA11)*(2.0*theta+(1.0-theta)*pA11)*(1.0-theta)/((1.0+theta)+(1.0+2.0*theta)))
		 + (pA11*(theta + pA11*(1.0-theta))*(1.0-theta)/(1+theta)+(2.0*theta+(1.0-theta)*pA11)*(theta+2.0*pA11*(1.0-theta))/(1.0+theta))*(3.0*theta+(1.0-theta)*pA11)/(1.0+2.0*theta));
    }
    //if they're homo homo with none shared
    else if((A11 == A12) && (A21 == A22)) { 
      if(allelenum==1 || allelenum==2) 
	deltaps[k] = ks[0]*(pA21*(1.0-theta)*(theta+(1.0-theta)*pA21)*(theta+2.0*pA11*(1.0-theta)))/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[0]*(pA11*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+2.0*pA21*(1.0-theta)))/((theta+1.0)*(2.0*theta+1.0));
    }
    //if they're homo het with one shared
    else if((A11 == A12) && (A12 == A21)) { 
      if(allelenum==1 || allelenum==2 || allelenum==3)
	deltaps[k] = ks[1]*(pA22*(1.0-theta)*(theta+2.0*pA11*(1.0-theta)))/(theta+1.0)
	  + ks[0]*(2.0*pA22*(1.0-theta)*(2.0*pow(theta,2.0)+6.0*theta*pA11*(1.0-theta)+3.0*pow(pA11,2.0)*pow((1.0-theta),2.0)))/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[1]*(pA11*(1.0-theta)*(theta+pA11*(1.0-theta)))/(theta+1.0)
	  + ks[0]*(2.0*pA11*(1.0-theta)*(theta+pA11*(1.0-theta))*(2.0*theta+(1.0-theta)*pA11))/((theta+1.0)*(2.0*theta+1.0));
    }
    else if((A11 == A12) && (A12 == A22)) { 
      if(allelenum==1 || allelenum==2 || allelenum==4)
	deltaps[k] = ks[1]*(pA21*(1.0-theta)*(theta+2.0*pA11*(1.0-theta)))/(theta+1.0)
	  + ks[0]*(2.0*pA21*(1.0-theta)*(2.0*pow(theta,2.0)+6.0*theta*pA11*(1.0-theta)+3.0*pow(pA11,2.0)*pow((1.0-theta),2.0)))/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[1]*(pA11*(1.0-theta)*(theta+pA11*(1.0-theta)))/(theta+1.0)
	  + ks[0]*(2.0*pA11*(1.0-theta)*(theta+pA11*(1.0-theta))*(2.0*theta+(1.0-theta)*pA11))/((theta+1.0)*(2.0*theta+1.0));
    }
    else if((A21 == A22) && (A21 == A11)) { 
      if(allelenum==1 || allelenum==3 || allelenum==4)
	deltaps[k] = ks[1]*(pA12*(1.0-theta)*(theta+2.0*pA11*(1.0-theta)))/(theta+1.0)
	  + ks[0]*(2.0*pA12*(1.0-theta)*(2.0*pow(theta,2.0)+6.0*theta*pA11*(1.0-theta)+3.0*pow(pA11,2.0)*pow((1.0-theta),2.0)))/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[1]*(pA11*(1.0-theta)*(theta+pA11*(1.0-theta)))/(theta+1.0)
	  + ks[0]*(2.0*pA11*(1.0-theta)*(theta+pA11*(1.0-theta))*(2.0*theta+(1.0-theta)*pA11))/((theta+1.0)*(2.0*theta+1.0));
    }
    else if((A21 == A22) && (A21 == A12)) { 
      if(allelenum==2 || allelenum==3 || allelenum==4)
	deltaps[k] = ks[1]*(pA11*(1.0-theta)*(theta+2.0*pA21*(1.0-theta)))/(theta+1.0)
	  + ks[0]*(2.0*pA11*(1.0-theta)*(2.0*pow(theta,2.0)+6.0*theta*pA21*(1.0-theta)+3.0*pow(pA21,2.0)*pow((1.0-theta),2.0)))/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[1]*(pA21*(1.0-theta)*(theta+pA21*(1.0-theta)))/(theta+1.0)
	  + ks[0]*(2.0*pA21*(1.0-theta)*(theta+pA21*(1.0-theta))*(2.0*theta+(1.0-theta)*pA21))/((theta+1.0)*(2.0*theta+1.0));

    }
    //if they're het het with both shared
    else if((A11 == A21) && (A12 == A22)) { 
      if(allelenum==1 || allelenum==3)
	deltaps[k] = ks[2]*2.0*pA12*(1.0-theta)
	  + ks[1]*pA12*(1.0-theta)*(2.0*theta+(1.0-theta)*(2.0*pA11+pA12))/(theta+1.0)
	  + ks[0]*4.0*pA12*(1.0-theta)*(theta+(1.0-theta)*pA12)*(theta+2.0*pA11*(1.0-theta))/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[2]*2.0*pA11*(1.0-theta)
	  + ks[1]*pA11*(1.0-theta)*(2.0*theta+(1.0-theta)*(2.0*pA12+pA11))/(theta+1.0)
	  + ks[0]*4.0*pA11*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+2.0*pA12*(1.0-theta))/((theta+1.0)*(2.0*theta+1.0));
    }
    else if((A11 == A22) && (A12 == A21)) { 
        if(allelenum==1 || allelenum==4)
	  deltaps[k] = ks[2]*2.0*pA12*(1.0-theta)
	    + ks[1]*pA12*(1.0-theta)*(2.0*theta+(1.0-theta)*(2.0*pA11+pA12))/(theta+1.0)
	    + ks[0]*4.0*pA12*(1.0-theta)*(theta+(1.0-theta)*pA12)*(theta+2.0*pA11*(1.0-theta))/((theta+1.0)*(2.0*theta+1.0));
	else
	  deltaps[k] = ks[2]*2.0*pA11*(1.0-theta)
	    + ks[1]*pA11*(1.0-theta)*(2.0*theta+(1.0-theta)*(2.0*pA12+pA11))/(theta+1.0)
	    + ks[0]*4.0*pA11*(1.0-theta)*(theta+(1.0-theta)*pA11)*(theta+2.0*pA12*(1.0-theta))/((theta+1.0)*(2.0*theta+1.0));
    }
    //if they're het het with one shared
    else if(A11 == A21) { 
      if(allelenum==1 || allelenum==3)
	deltaps[k] = ks[1]*pA12*pA22*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA12*pA22*pow((1.0-theta),2.0)*(theta+2.0*(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
      else if(allelenum==2)
	deltaps[k] = ks[1]*pA11*pA22*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA11*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[1]*pA11*pA12*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA11*pA12*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
    }
    else if(A11 == A22) { 
      if(allelenum==1 || allelenum==4)
	deltaps[k] = ks[1]*pA12*pA21*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA12*pA21*pow((1.0-theta),2.0)*(theta+2.0*(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
      else if(allelenum==2)
	deltaps[k] = ks[1]*pA11*pA21*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA11*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[1]*pA11*pA12*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA11*pA12*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
    }
    else if(A12 == A21) { 
      if(allelenum==2 || allelenum==3)
	deltaps[k] = ks[1]*pA11*pA22*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA11*pA22*pow((1.0-theta),2.0)*(theta+2.0*(1.0-theta)*pA12)/((theta+1.0)*(2.0*theta+1.0));
      else if(allelenum==2)
	deltaps[k] = ks[1]*pA12*pA22*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA12*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[1]*pA11*pA12*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA11*pA12*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((theta+1.0)*(2.0*theta+1.0));
    }
    else if(A12 == A22) { 
      if(allelenum==2 || allelenum==4)
	deltaps[k] = ks[1]*pA11*pA21*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA11*pA21*pow((1.0-theta),2.0)*(theta+2.0*(1.0-theta)*pA12)/((theta+1.0)*(2.0*theta+1.0));
      else if(allelenum==1)
	deltaps[k] = ks[1]*pA12*pA21*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA12*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((theta+1.0)*(2.0*theta+1.0));
      else   
	deltaps[k] = ks[1]*pA11*pA12*pow((1.0-theta),2.0)/(theta+1.0)
	  + ks[0]*4.0*pA11*pA12*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA12)/((theta+1.0)*(2.0*theta+1.0));
    }
    //if they're homo het with none shared
    else if(A11 == A12) { 
      if(allelenum==1 || allelenum==2)
	deltaps[k] = ks[0]*2.0*pA21*pA22*pow((1.0-theta),2.0)*(theta+2.0*(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
      else if(allelenum==3)
	deltaps[k] = ks[0]*2.0*pA11*pA22*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[0]*2.0*pA11*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA11)/((theta+1.0)*(2.0*theta+1.0));
    }
    else if(A21 == A22) { 
      if(allelenum==3 || allelenum==4)
	deltaps[k] = ks[0]*2.0*pA11*pA12*pow((1.0-theta),2.0)*(theta+2.0*(1.0-theta)*pA21)/((theta+1.0)*(2.0*theta+1.0));
      else if(allelenum==1)
	deltaps[k] = ks[0]*2.0*pA12*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA21)/((theta+1.0)*(2.0*theta+1.0));
      else
	deltaps[k] = ks[0]*2.0*pA11*pA21*pow((1.0-theta),2.0)*(theta+(1.0-theta)*pA21)/((theta+1.0)*(2.0*theta+1.0));
    }
    //if they're het het with none shared
    else {
      if(allelenum==1)
	deltaps[k] = ks[0]*pA12*pA21*pA22*pow((1.0-theta),3.0)/((theta+1.0)*(2.0*theta+1));
      if(allelenum==2)
	deltaps[k] = ks[0]*pA11*pA21*pA22*pow((1.0-theta),3.0)/((theta+1.0)*(2.0*theta+1));
      if(allelenum==3)
	deltaps[k] = ks[0]*pA11*pA12*pA22*pow((1.0-theta),3.0)/((theta+1.0)*(2.0*theta+1));
      if(allelenum==4)
	deltaps[k] = ks[0]*pA11*pA12*pA21*pow((1.0-theta),3.0)/((theta+1.0)*(2.0*theta+1));
    }
  }
}




//calc the derivative of the log liklihood ratio at one locus by one allele
void delta_ln_LR_given_l_i(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *dlnLRli, double theta, double *ksP, int allelenum) {
  int k;
  double *ksD, *pD, *pP, *dpD, *dpP;
  
  ksD = (double*) malloc(3*sizeof(double));
  ksD[0]=1.0;   ksD[1]=0.0;   ksD[2]=0.0;
  
  pD = (double*) malloc(npops*sizeof(double));
  pP = (double*) malloc(npops*sizeof(double));
  dpD = (double*) malloc(npops*sizeof(double));
  dpP = (double*) malloc(npops*sizeof(double));
  
  pGenoGivenK(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, pP, theta, ksP); 
  pGenoGivenK(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, pD, theta, ksD);
  deltapGenoGivenK(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dpP, theta, ksP, allelenum);
  deltapGenoGivenK(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dpD, theta, ksD, allelenum);
  
  for(k=0; k<npops; k++) {
    dlnLRli[k] = dpP[k]/pP[k] - dpD[k]/pD[k];
  }
    
  free(ksD);
  free(pD);
  free(pP);
  free(dpD);
  free(dpP);
}




//find the sum of squares of delta ln LR times var(p_i) for every possible i
void sum_square_dlnLRl(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *sumsqs, double theta, double *ksP) {
  int k;
  int A11, A12, A21, A22;
  double pA11, pA12, pA21, pA22;
  double *dlnLRli;

  dlnLRli = (double*) malloc(npops*sizeof(double));

  //init sumsqs
  for(k=0; k<npops; k++)
    sumsqs[k] = 0.0;

  //get alleles
  A11 = geno1[locusnum][0];
  A12 = geno1[locusnum][1];
  A21 = geno2[locusnum][0];
  A22 = geno2[locusnum][1];
  
  //if they're homo homo with both shared
  if((A11 == A12) && (A12 == A21) && (A21 == A22)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're homo homo with none shared
  else if((A11 == A12) && (A21 == A22)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA21 = AFs[locusnum][A21][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA21*(1.0-pA21)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're homo het with one shared
  else if((A11 == A12) && (A12 == A21)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA22 = AFs[locusnum][A22][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA22*(1.0-pA22)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if((A11 == A12) && (A12 == A22)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA21 = AFs[locusnum][A21][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA21*(1.0-pA21)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if((A21 == A22) && (A21 == A11)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if((A21 == A22) && (A21 == A12)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're het het with both shared
  else if((A11 == A21) && (A12 == A22)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    pA22=dlnLRli[0];
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }

  }
  else if((A11 == A22) && (A12 == A21)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're het het with one shared
  else if(A11 == A21) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA22 = AFs[locusnum][A22][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA22*(1.0-pA22)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if(A11 == A22) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA21 = AFs[locusnum][A21][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA21*(1.0-pA21)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if(A12 == A21) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA22 = AFs[locusnum][A22][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA22*(1.0-pA22)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if(A12 == A22) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA21 = AFs[locusnum][A21][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA21*(1.0-pA21)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're homo het with none shared
  else if(A11 == A12) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA21 = AFs[locusnum][A21][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA21*(1.0-pA21)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA22 = AFs[locusnum][A22][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA22*(1.0-pA22)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if(A21 == A22) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA21 = AFs[locusnum][A21][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA21*(1.0-pA21)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're het het with none shared
  else {
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 1);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      sumsqs[k] = pow(dlnLRli[k],2.0)*(pA11*(1.0-pA11)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA12 = AFs[locusnum][A12][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA12*(1.0-pA12)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA21 = AFs[locusnum][A21][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA21*(1.0-pA21)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA22 = AFs[locusnum][A22][k]/100.0;
      sumsqs[k] += pow(dlnLRli[k],2.0)*(pA22*(1.0-pA22)*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }

  free(dlnLRli);
}



//find the sum of products of delta ln LR times cov(p_i,p_j) for every possible i j pair
void sum_prod_dlnLRl(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *sums, double theta, double *ksP) {
  int k;
  int A11, A12, A21, A22;
  double pA11, pA12, pA21, pA22;
  double *dlnLRli1, *dlnLRli2, *dlnLRli3, *dlnLRli4;

  //malloc delta ln LR given l,i  1,2,3,4
  dlnLRli1 = (double*) malloc(npops*sizeof(double));
  dlnLRli2 = (double*) malloc(npops*sizeof(double));
  dlnLRli3 = (double*) malloc(npops*sizeof(double));
  dlnLRli4 = (double*) malloc(npops*sizeof(double));

  //get alleles
  A11 = geno1[locusnum][0];
  A12 = geno1[locusnum][1];
  A21 = geno2[locusnum][0];
  A22 = geno2[locusnum][1];
    
  //if they're homo homo with both shared
  if((A11 == A12) && (A12 == A21) && (A21 == A22)) { 
    for(k=0; k<npops; k++) {
      sums[k] = 0.0;
    }
  }
  //if they're homo homo with none shared
  else if((A11 == A12) && (A21 == A22)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA21 = AFs[locusnum][A21][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're homo het with one shared
  else if((A11 == A12) && (A12 == A21)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA22 = AFs[locusnum][A22][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if((A11 == A12) && (A12 == A22)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 3);   
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA21 = AFs[locusnum][A21][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if((A21 == A22) && (A21 == A11)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }	
  }
  else if((A21 == A22) && (A21 == A12)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're het het with both shared
  else if((A11 == A21) && (A12 == A22)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if((A11 == A22) && (A12 == A21)) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're het het with one shared
  else if(A11 == A21) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli3, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      pA22 = AFs[locusnum][A22][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli1[k]*dlnLRli3[k]*(-1.0*pA11*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli2[k]*dlnLRli3[k]*(-1.0*pA12*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if(A11 == A22) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli3, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      pA21 = AFs[locusnum][A21][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli1[k]*dlnLRli3[k]*(-1.0*pA11*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli2[k]*dlnLRli3[k]*(-1.0*pA12*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if(A12 == A21) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli4, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      pA22 = AFs[locusnum][A22][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli1[k]*dlnLRli3[k]*(-1.0*pA11*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli2[k]*dlnLRli3[k]*(-1.0*pA12*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if(A12 == A22) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli3, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      pA21 = AFs[locusnum][A21][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli1[k]*dlnLRli3[k]*(-1.0*pA11*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli2[k]*dlnLRli3[k]*(-1.0*pA12*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're homo het with none shared
  else if(A11 == A12) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli3, theta, ksP, 3);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli4, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA21 = AFs[locusnum][A21][k]/100.0;
      pA22 = AFs[locusnum][A22][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli1[k]*dlnLRli3[k]*(-1.0*pA11*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli2[k]*dlnLRli3[k]*(-1.0*pA21*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  else if(A21 == A22) { 
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli3, theta, ksP, 3);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      pA21 = AFs[locusnum][A21][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli1[k]*dlnLRli3[k]*(-1.0*pA11*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli2[k]*dlnLRli3[k]*(-1.0*pA12*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }
  //if they're het het with none shared
  else {
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli1, theta, ksP, 1);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli2, theta, ksP, 2);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli3, theta, ksP, 3);
    delta_ln_LR_given_l_i(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, dlnLRli4, theta, ksP, 4);
    for(k=0; k<npops; k++) {
      pA11 = AFs[locusnum][A11][k]/100.0;
      pA12 = AFs[locusnum][A12][k]/100.0;
      pA21 = AFs[locusnum][A21][k]/100.0;
      pA22 = AFs[locusnum][A22][k]/100.0;
      sums[k] = dlnLRli1[k]*dlnLRli2[k]*(-1.0*pA11*pA12*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli1[k]*dlnLRli3[k]*(-1.0*pA11*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli1[k]*dlnLRli4[k]*(-1.0*pA11*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli2[k]*dlnLRli3[k]*(-1.0*pA12*pA21*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli2[k]*dlnLRli4[k]*(-1.0*pA12*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]))
	+ dlnLRli3[k]*dlnLRli4[k]*(-1.0*pA21*pA22*((2.0*npeople[locusnum][k]-1.0)*theta+1.0)/(2.0*npeople[locusnum][k]));
    }
  }

  //free stuff
  free(dlnLRli1);
  free(dlnLRli2);
  free(dlnLRli3);
  free(dlnLRli4);
}



//calc the variance of the log LR at one locus
void var_ln_LR_given_l(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, int locusnum, double *vars, double theta, double *ksP) {
  int k;
  double *sum1, *sum2; 
  
  //malloc sums
  sum1 = (double*) malloc(npops*sizeof(double));
  sum2 = (double*) malloc(npops*sizeof(double));

  //calc sums
  sum_square_dlnLRl(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, sum1, theta, ksP);
  sum_prod_dlnLRl(geno1, geno2, AFs, npops, nalleles, npeople, locusnum, sum2, theta, ksP);

  for(k=0; k<npops; k++) {
    vars[k] = sum1[k] + sum2[k];
  }

  //free
  free(sum1);
  free(sum2);
}

//calculates the variance of the log LR computed using all loci
void var_ln_LR(int **geno1, int **geno2, double ***AFs, int npops, int* nalleles, int **npeople, double *vars, double theta, double *ksP) {
  int i, k, nloci=13;
  double *varl;
  int n;

  varl = (double*) malloc(npops*sizeof(double));
  
  for(k=0; k<npops; k++)
    vars[k] = 0.0;

  for(i=0; i<nloci; i++) {
    var_ln_LR_given_l(geno1, geno2, AFs, npops, nalleles, npeople, i, varl, theta, ksP);
    for(k=0; k<npops; k++) {
      vars[k] += varl[k];


      
      if(varl[k]<0) {
	printf("weird varl is %f  \n", varl[k]);
	printf("genos are %d %d, %d %d\n", geno1[i][0], geno1[i][1], geno2[i][0], geno2[i][1]);
	for(n=0; n<nalleles[i]; n++) {
	  //printf("%f ", AFs[i][n][k]);
	}
	//printf("\n");
	

      }
    
      
    }
  }

  free(varl);
}

//calculates the log LR computed using all loci
void ln_LR(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *lnLRs, double theta, double *ksp) {
  int i;
  double *ksd, *lnp_hp, *lnp_hd;

  ksd = (double*) malloc(3*sizeof(double));
  lnp_hp = (double*) malloc(npops*sizeof(double));
  lnp_hd = (double*) malloc(npops*sizeof(double));

  ksd[0]=1.0; ksd[1]=0.0; ksd[2]=0.0;

  lnpGenosGivenK(geno1, geno2, AFs, npops, nalleles, npeople, lnp_hp, theta, ksp);
  lnpGenosGivenK(geno1, geno2, AFs, npops, nalleles, npeople, lnp_hd, theta, ksd);

  for(i=0; i<npops; i++)
    lnLRs[i] = lnp_hp[i] - lnp_hd[i];

  free(ksd);
  free(lnp_hp);
  free(lnp_hd);
}
	
//calculates the LR computed using all loci
void LR(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *LRs, double theta, double *ksp) {
  int i, j, k;
  double*** newAFs;

  //malloc newAFs and init
  newAFs = (double***) malloc(13*sizeof(double**));
  for(i=0; i<13; i++) {
    newAFs[i] = (double**) malloc(nalleles[i]*sizeof(double*));
    if(newAFs[i] == NULL)
      printf("ran out of memory!\n");
    for(j=0; j<nalleles[i]; j++) {
      newAFs[i][j] = (double*) malloc(npops*sizeof(double*));
      if(newAFs[i][j] == NULL)
	printf("ran out of memory!\n");
      for(k=0; k<npops; k++) {
	newAFs[i][j][k] = -1.0;
      }
    }
  }

  addIndivsToAFs(geno1, geno2, AFs, newAFs, npops, nalleles, npeople);

  ln_LR(geno1, geno2, newAFs, npops, nalleles, npeople, LRs, theta, ksp);

  for(i=0; i<npops; i++)
    LRs[i] = exp(LRs[i]);

  //free newAFs
  for(i=0; i<13; i++) {
    for(j=0; j<nalleles[i]; j++) {
      free(newAFs[i][j]);
    }
    free(newAFs[i]);
  }
  free(newAFs);
}



//calculates confidence intervals for the LR
void CIs(int **geno1, int **geno2, double ***AFs, int npops, int *nalleles, int **npeople, double *lci, double *uci, double theta, double *ksP, double z) {
  int i, j, k;
  double *LRs, *vars;
  double ***newAFs;
  //int m, n;

  //malloc LRs and vars
  LRs = (double*) malloc(npops*sizeof(double));
  vars = (double*) malloc(npops*sizeof(double));
  
  //malloc newAFs and init
  newAFs = (double***) malloc(13*sizeof(double**));
  for(i=0; i<13; i++) {
    newAFs[i] = (double**) malloc(nalleles[i]*sizeof(double*));
    if(newAFs[i] == NULL)
      printf("ran out of memory!\n");
    for(j=0; j<nalleles[i]; j++) {
      newAFs[i][j] = (double*) malloc(npops*sizeof(double*));
      if(newAFs[i][j] == NULL)
	printf("ran out of memory!\n");
      for(k=0; k<npops; k++) {
	newAFs[i][j][k] = -1.0;
      }
    }
  }


  LR(geno1, geno2, AFs, npops, nalleles, npeople, LRs, theta, ksP);

  addIndivsToAFs(geno1, geno2, AFs, newAFs, npops, nalleles, npeople);

  var_ln_LR(geno1, geno2, newAFs, npops, nalleles, npeople, vars, theta, ksP);


  for(i=0; i<npops; i++) {
    uci[i] = ( LRs[i] * exp(z*pow(vars[i], 0.5)));
    lci[i] = ( LRs[i] / exp(z*pow(vars[i], 0.5)));

    if(isnan(lci[i])) {
      printf("lCI is nan!, var is %f\n", vars[i]);
    }
  }

  if(uci[0]<lci[0])
    printf("upper-lower switch here\n");

  //printf("LRs:  ");
  //printDouArr(LRs, npops);
  //printf("vars: ");
  //printDouArr(vars, npops);

  free(LRs);
  free(vars);
  //free newAFs
  for(i=0; i<13; i++) {
    for(j=0; j<nalleles[i]; j++) {
      free(newAFs[i][j]);
    }
    free(newAFs[i]);
  }
  free(newAFs);
}

//returns 1 if they completely match, 0 otherwise
int checkHapMatch(int *hap1, int *hap2, int nloci) {
  int i, match=1;

  for(i=0; i<nloci; i++) {
    if(hap1[i] != hap2[i]) 
      match = 0;
  }

  return(match);
}

//returns number of hap1 haps in refHaps in popi
int countHaps(int *hap, int ***refHaps, int popi, int nindivs, int nloci) {
  int i, j;
  int match, nhaps = 0;
  
  for(i=0; i<nindivs; i++) {
    match = 1;
    for(j=0; j<nloci; j++) 
      if(hap[j] != refHaps[popi][i][j]) {
	match = 0;
	break;
      }

    if(match == 1)
      nhaps++;
  }

  return(nhaps);
}

//calculates the Y haplotype confidence limit as suggested by SWGDAM
void YCLhap(int *hap1, int *hap2, int ***refHaps, int npops, int *nindivs, int nloci, double *YCLs) {
  int i;
  int curNHaps;
  double phap;
  
  if(checkHapMatch(hap1, hap2, nloci)==0) {
    for(i=0; i<npops; i++) 
      YCLs[i] = 0;
    return;
  }

  for(i=0; i<npops; i++) {
    curNHaps = countHaps(hap1, refHaps, i, nindivs[i], nloci);

    if(curNHaps == 0)
      YCLs[i] = 1.0 - pow(0.05,(1.0/nindivs[i]));
    else {
      phap = 1.0*curNHaps/nindivs[i];
      YCLs[i] = phap + 1.96*pow(((phap)*(1.0-phap))/(nindivs[i]), 0.5);
    }
  }
}



//calculates the Y haplotype confidence limit as suggested by SWGDAM, given indexed haps instead of full haps
void YCL(int hap1, int hap2, double **hapFreqs, int npops, int *nindivs, double *YCLs) {
  int i;
  int curNHaps;
  double phap;
  
  if(hap1 != hap2) {
    for(i=0; i<npops; i++) 
      YCLs[i] = 0;
    return;
  }

  for(i=0; i<npops; i++) {
    curNHaps = 1.0*hapFreqs[i][hap1]*nindivs[i];

    if(curNHaps == 0)
      YCLs[i] = 1.0 - pow(0.05,(1.0/nindivs[i]));
    else {
      phap = 1.0*curNHaps/nindivs[i];
      YCLs[i] = phap + 1.96*pow(((phap)*(1.0-phap))/(nindivs[i]), 0.5);
    }
  }
}

//calculates the Y haplotype confidence limit as suggested by SWGDAM, given indexed haps instead of full haps, but uses cumulative hap freqs
void YCLall(int hap1, int hap2, double **hapFreqs, int npops, int *nindivs, double *YCLs) {
  int i;
  int nHaps, totnindivs;
  double phap;
  
  if(hap1 != hap2) {
    for(i=0; i<npops; i++) 
      YCLs[i] = 0;
    return;
  }

  nHaps = 0;
  totnindivs = 0;
  for(i=0; i<npops; i++) {
    nHaps += 1.0*hapFreqs[i][hap1]*nindivs[i];
    totnindivs += nindivs[i];
  }

  if(nHaps == 0)
    for(i=0; i<npops; i++) 
      YCLs[i] = 1.0 - pow(0.05,(1.0/totnindivs));
  else {
    for(i=0; i<npops; i++) {
      phap = 1.0*nHaps/totnindivs;
      YCLs[i] = phap + 1.96*pow(((phap)*(1.0-phap))/(totnindivs), 0.5);
    }
  }
}

//returns the LR for this Y hap pair
void YLR(int hap1, int hap2, double **hapFreqs, int npops, int *nindivs, double *YLRs) {
  int i;
  double *YCLs;

  YCLs = (double*) malloc(npops*sizeof(double));

  YCLall(hap1, hap2, hapFreqs, npops, nindivs, YCLs);

  for(i=0; i<npops; i++) {
    YLRs[i] = 1.0/YCLs[i];
  }

  free(YCLs);
}

//count haps and calc freqs for each population
void calcHapFreqs(int ***refHaps, int npops, int *nindivs, int nloci, double **hapFreqs) {
  int i, j, k, l, hapFreqsj;
  int numHap, popHaps, curmatch, totindivs;
  int **doneHaps;

  //malloc/init doneHaps
  doneHaps = (int**) malloc(npops*sizeof(int*));
  for(i=0; i<npops; i++) {
    doneHaps[i] = (int*) malloc(nindivs[i]*sizeof(int));
    if(doneHaps[i] == NULL)  printf("ran out of memory!\n");
    for(j=0; j<nindivs[i]; j++)
      doneHaps[i][j] = 0;
  }
  
  //init totindivs
  totindivs = 0;
  for(i=0; i<npops; i++) 
    totindivs += nindivs[i];
  
  //init hapFreqs  
  for(i=0; i<npops; i++) 
    for(j=0; j<totindivs; j++) 
      hapFreqs[i][j] = 0.0;

  //run through all haplotypes counting observations
  hapFreqsj = 0;
  for(i=0; i<npops; i++) {
    for(j=0; j<nindivs[i]; j++) {
      if(doneHaps[i][j] == 0) {
	//get total number of haplotypes in all pop groups to check singletons
	numHap = 0;
	for(k=0; k<npops; k++) 
	  numHap += countHaps(refHaps[i][j], refHaps, k, nindivs[k], nloci);
	
	//if singleton, set hapFreqs and doneHaps
	if(numHap == 1) {
	  hapFreqs[i][hapFreqsj] = 1.0/nindivs[i];
	  doneHaps[i][j] = 1;
	} else {	//handle replicates
	  for(k=i; k<npops; k++) {
	    popHaps = 0;
	    for(l=0; l<nindivs[k]; l++) {
	      curmatch = checkHapMatch(refHaps[i][j], refHaps[k][l], nloci);
	      if(curmatch == 1) {
		popHaps += 1;
		doneHaps[k][l] = 1;
	      }
	    }
	    hapFreqs[k][hapFreqsj] = 1.0*popHaps/nindivs[k];
	  }
	}
	hapFreqsj++;
      }
    }
  }


  for(i=0; i<npops; i++)
    free(doneHaps[i]);
  free(doneHaps);
}

//estimate theta between each pop group as outlined in weir and cockerham
void estimateTheta(int ***refHaps, int npops, int *nindivs, int nloci, double **thetas) {
  int i, j, k;
  double topsum, botsum;
  double nlc, p1, p2, pbar, MSP, MSG;
  double **hapFreqs;
  int totindivs = 0;

  //init totindivs
  for(i=0; i<npops; i++) 
    totindivs += nindivs[i];

  //malloc hapFreqs[population][haplotype]
  hapFreqs = (double**) malloc(npops*sizeof(double*));
  for(i=0; i<npops; i++) {
    hapFreqs[i] = (double*) malloc(totindivs*sizeof(double));
    if(hapFreqs[i] == NULL)  printf("ran out of memory!\n");
  }
  
  //calc hapFreqs
  calcHapFreqs(refHaps, npops, nindivs, nloci, hapFreqs);
  //print2dDouArr(hapFreqs, npops, totindivs);  

  //run through all pop group pairs
  for(i=0; i<npops; i++) {
    for(j=0; j<npops; j++) {
      thetas[i][j] = 0.0;
      topsum = 0.0;
      botsum = 0.0;

      //run through all haplotypes (aka alleles)
      for(k=0; k<totindivs; k++) {
	nlc = (2.0*nindivs[i]*nindivs[j])/(nindivs[i]+nindivs[j]);
	p1 = hapFreqs[i][k];
        p2 = hapFreqs[j][k];
	pbar = (nindivs[i]*p1 + nindivs[j]*p2)/(nindivs[i]+nindivs[j]);

        MSP = nindivs[i]*pow((p1-pbar),2.0) + nindivs[j]*pow((p2-pbar),2.0);
        MSG = (nindivs[i]*p1*(1.0-p1) + nindivs[j]*p2*(1.0-p2))/(nindivs[i]+nindivs[j]-2.0);

	topsum += (MSP-MSG);
	botsum += (MSP+(nlc-1.0)*MSG);
      }
      thetas[i][j] = 1.0*topsum/botsum;
    }
  }

  //free
  for(i=0; i<npops; i++) 
    free(hapFreqs[i]);
  free(hapFreqs);
}
 
