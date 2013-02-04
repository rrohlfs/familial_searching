#include <stdio.h>

/* i/o functions i've written, mostly printing different arrays */

/*
//prints char* array
void printStrArr(char **arr, int length) {
  int i;
  for(i=0; i<length; i++) {
    printf("%s ", &arr[i]);
  }
  printf("\n");
}
*/

//prints char** array
void print2dCharArr(char **arr, int nrow) {
  int i;
  for(i=0; i<nrow; i++) {
      printf("%s\n", &arr[i][0]);
  }
}


//prints int array
void printIntArr(int *arr, int length) {
  int i;
  for(i=0; i<length; i++) {
    printf("%d ", arr[i]);
  }
  printf("\n");
}

//prints float array
void printFloArr(float *arr, int length) {
  int i;
  for(i=0; i<length; i++) {
    printf("%f ", arr[i]);
  }
  printf("\n");
}

//prints double array
void printDouArr(double *arr, int length) {
  int i;
  for(i=0; i<length; i++) {
    printf("%e ", arr[i]);
  }
  printf("\n");
}


//prints long double array
void printLDouArr(long double *arr, int length) {
  int i;
  for(i=0; i<length; i++) {
    printf("%Le ", arr[i]);
  }
  printf("\n");
}

//prints 2-d int array
void print2dIntArr(int **arr, int nrows, int ncols) {
  int i,j;
  for(i=0; i<nrows; i++) {
    for(j=0; j<ncols; j++) {
      printf("%d ", arr[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

//prints 2-d double array
void print2dDouArr(double **arr, int nrows, int ncols) {
  int i,j;
  for(i=0; i<nrows; i++) {
    for(j=0; j<ncols; j++) {
      printf("%f ", arr[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

//prints 2-d long double array
void print2dLDouArr(long double **arr, int nrows, int ncols) {
  int i,j;
  for(i=0; i<nrows; i++) {
    for(j=0; j<ncols; j++) {
      printf("%Lf ", arr[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

//prints 2-d float array
void print2dFloArr(float **arr, int nrows, int ncols) {
  int i,j;
  for(i=0; i<nrows; i++) {
    for(j=0; j<ncols; j++) {
      printf("%f ", arr[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}
/*
//parses a .ped file
void parsePED(char *filename) {
  FILE *fin;
  int i;

  fin = fopen(filename, "r");
  if(fin == NULL) {
    printf("Unable to locate input file: %s\n", filename);
    exit(1);
  }
  
  fscanf(fin, "WTCCC%d WTCCC%d %d %d %d %d", &indivID, &indivID, &par1, &par2, &sex, &pheno);
  while(curSNP
  fscanf(fin, "%d%d", &npeople, &nSNPs);

}
*/

//prints a 2d Double array to file
void write2dDouArr(double** arr, int nrows, int ncols, char* filename) {
  int i, j;
  FILE *fout;

  fout = fopen(filename, "w");

  for(i=0; i<nrows; i++) {
    for(j=0; j<ncols; j++) 
      fprintf(fout, "%e ", arr[i][j]);
    fprintf(fout, "\n");
  }

  fclose(fout);
}

//prints a 3d Double array to file
void write3dDouArr(double*** arr, int imax, int jmax, int kmax, char* filename) {
  int i, j, k;
  FILE *fout;

  fout = fopen(filename, "w");

  for(i=0; i<imax; i++) {
    for(j=0; j<jmax; j++) {
      for(k=0; k<kmax; k++) 
	fprintf(fout, "%e ", arr[i][j][k]);
      fprintf(fout, "\n");
    }
    fprintf(fout, "\n");
  }

  fclose(fout);
}
