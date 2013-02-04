/* header file for mystat */

//numerical recipes chisq stuff
float gammln(float xx);
void gser(float *gamser, float a, float x, float *gln);
void gcf(float *gammcf, float a, float x, float *gln);
float gammp(float a, float x);

//my own simple math functions
//float fmin(float a, float b);
double dmin(double a, double b);
int imin(int a, int b);
long double factorial(int n);
long double logfactorial(int n);
int trinomial(int n, int k1, int k2, int k3);
long double logtrinomial(int n, int k1, int k2, int k3);
int choose(int n, int k);
long double logchoose(int n, int k);
int mylonground(long double a);
int myround(double a);

//quick sorts
void ldsort(long double arr[], int arr2[], int beg, int end);
void isort(int arr[], long double arr2[], long double arr3[], int beg, int end);
void ldldsort(long double arr[], long double arr2[], int beg, int end);
void ldldldsort(long double arr[], long double arr2[], long double arr3[], int beg, int end);
void ldldldldsort(long double arr[], long double arr2[], long double arr3[], long double arr4[], int beg, int end);

//shuffles
void shuffle(int *arr, int imax);

//zero count
int zeroCount2dDouArr(double **mat, int nrow, int ncol);
