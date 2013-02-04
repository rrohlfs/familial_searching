#include <stdio.h>


//print 1d arrays
//void printStrArr(char **arr, int length);
void printIntArr(int *arr, int length);
void printFloArr(float *arr, int length);
void printDouArr(double *arr, int length);
void printLDouArr(long double *arr, int length);

//print 2d arrays
void print2dIntArr(int **arr, int nrows, int ncols);
void print2dDouArr(double **arr, int nrows, int ncols);
void print2dFloArr(float **arr, int nrows, int ncols);
void print2dCharArr(char **arr, int nrow);//, int ncol);
void print2dLDouArr(long double **arr, int nrows, int ncols);

//write arrays to file
void write2dDouArr(double** arr, int nrows, int ncols, char* filename);
void write3dDouArr(double*** arr, int imax, int jmax, int kmax, char* filename);
