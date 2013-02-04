/* stat functions i've written (or taken from nr) */

#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include "myio.h"

#define ITMAX 1000 //Maximum allowed number of iterations.
#define EPS 3.0e-7 //Relative accuracy.
#define FPMIN 1.0e-30 //Number near the smallest representable floating-point number.


/* Returns the value ln[%GÃ´ÂÂÂ%@(xx)] for xx > 0.
Internal arithmetic will be done in double precision, a nicety that you can omit if ve-gure
accuracy is good enough.  */
float gammln(float xx) {
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

/* Returns the incomplete gamma function P(a; x) evaluated by its series representation as gamser.
Also returns ln %GÃ´ÂÂÂ%@(a) as gln. */
void gser(float *gamser, float a, float x, float *gln) {
  float gammln(float xx);
  int n;
  float sum,del,ap;
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) printf("x less than 0 in routine gser\n");
    *gamser=0.0;
    return;
  } 
  else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    printf("a too large, ITMAX too small in routine gser\n");
    return;
  }
}


/* Returns the incomplete gamma function Q(a; x) evaluated by its continued fraction representation
   as gammcf. Also returns ln %GÃ´ÂÂÂ%@(a) as gln. */
void gcf(float *gammcf, float a, float x, float *gln) {
  float gammln(float xx);
  int i;
  float an,b,c,d,del,h;
  *gln=gammln(a);
  b=x+1.0-a; //Set up for evaluating continued fraction by modied Lentz's method (x5.2) with b0 = 0.
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) { //Iterate to convergence.
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) printf("a too large, ITMAX too small in gcf\n");
  *gammcf=exp(-x+a*log(x)-(*gln))*h; //Put factors in front.
}



//Returns the incomplete gamma function P(a; x).
float gammp(float a, float x) {
  void gcf(float *gammcf, float a, float x, float *gln);
  void gser(float *gamser, float a, float x, float *gln);
  float gamser,gammcf,gln;
  if (x < 0.0 || a <= 0.0) printf("Invalid arguments in routine gammp\n");
  if (x < (a+1.0)) { //Use the series representation.
    gser(&gamser,a,x,&gln);
    return gamser;
  } 
  else { //Use the continued fraction representation and take its complement.
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf; 
  }
}

/*
//my own min function
float fmin(float a, float b) {
  if(a<b)
    return a; 
  else
    return b;
}
*/

//my own min function
double dmin(double a, double b) {
  if(a<b)
    return a; 
  else
    return b;
}

//my own min function
int imin(int a, int b) {
  if(a<b)
    return a; 
  else
    return b;
}

//my own factorial function
long double factorial(int n) {
  int zero=0;
  if(n==zero) return 1;
  long double k = n;
  while(--n > zero) k*=n;
  return k;
}

//my own log factorial function
long double logfactorial(int n) {
  int zero=0;
  if(n==zero) return 0.0;
  long double k = log(n);
  while(--n > zero) k+=log(n);
  return k;
}

//returns multinomial of n choose k1, k2, k3
int trinomial(int n, int k1, int k2, int k3) {
  int ans = factorial(n)/(factorial(k1)*factorial(k2)*factorial(k3));
  return ans;  
}

//returns multinomial of n choose k1, k2, k3
long double logtrinomial(int n, int k1, int k2, int k3) {
  long double ans;

  ans = logfactorial(n) - (logfactorial(k1)+logfactorial(k2)+logfactorial(k3));
  
  return ans;
}

//returns binomial n choose k
int choose(int n, int k) {
  int ans = factorial(n)/(factorial(k)*factorial(n-k));
  return ans;
}

//returns log of binomial n choose k
long double logchoose(int n, int k) {
  long double ans = logfactorial(n) - (logfactorial(k)+logfactorial(n-k));
  return ans;
}

//rounds to the nearest integer
long long int mylonground(long double a) {
  long long int ans;

  if(a - floorl(a) < 0.5)
    ans = floorl(a);
  else 
    ans = ceill(a);

  return ans;
}

//rounds to the nearest integer
int myround(double a) {
  int ans;

  if(a - floor(a) < 0.5)
    ans = (int) floor(a);
  else 
    ans = (int) ceil(a);

  return ans;
}

void ldswap(long double *a, long double *b) {
  long double t=*a;
  *a=*b;
  *b=t;
}

void iswap(int *a, int *b) {
  int t=*a;
  *a=*b;
  *b=t;
}

void ldsort(long double arr[], int arr2[], int beg, int end) {
  if (end > beg + 1) {
    long double piv = arr[beg];
    int l = beg + 1, r = end;

    while (l < r) {
      if (arr[l] <= piv)
	l++;
      else {
	ldswap(&arr[l], &arr[--r]);
	iswap(&arr2[l], &arr2[r]);
      }
    }
    ldswap(&arr[--l], &arr[beg]);
    iswap(&arr2[l], &arr2[beg]);
    ldsort(arr, arr2, beg, l);
    ldsort(arr, arr2, r, end);
  }
}

void ldldsort(long double arr[], long double arr2[], int beg, int end) {
  if (end > beg + 1) {
    long double piv = arr[beg];
    int l = beg + 1, r = end;
    
    while (l < r) {
      if (arr[l] <= piv)
	l++;
      else {
	ldswap(&arr[l], &arr[--r]);
	ldswap(&arr2[l], &arr2[r]);
      }
    }
    ldswap(&arr[--l], &arr[beg]);
    ldswap(&arr2[l], &arr2[beg]);
    ldldsort(arr, arr2, beg, l);
    ldldsort(arr, arr2, r, end);
  }
}


void ldldldsort(long double arr[], long double arr2[], long double arr3[], int beg, int end) {
  if (end > beg + 1) {
    //printLDouArr(arr, end-beg);

    long double piv = arr[beg];
    int l = beg + 1, r = end;
    
    while (l < r) {
      if (arr[l] <= piv)
	l++;
      else {
	ldswap(&arr[l], &arr[--r]);
	ldswap(&arr2[l], &arr2[r]);
	ldswap(&arr3[l], &arr3[r]);
      }
    }
    ldswap(&arr[--l], &arr[beg]);
    ldswap(&arr2[l], &arr2[beg]);
    ldswap(&arr3[l], &arr3[beg]);
    ldldldsort(arr, arr2, arr3, beg, l);
    ldldldsort(arr, arr2, arr3, r, end);
  }
}

void ldldldldsort(long double arr[], long double arr2[], long double arr3[], long double arr4[], int beg, int end) {
  if (end > beg + 1) {
    //printLDouArr(arr, end-beg);

    long double piv = arr[beg];
    int l = beg + 1, r = end;
    
    while (l < r) {
      if (arr[l] <= piv)
	l++;
      else {
	ldswap(&arr[l], &arr[--r]);
	ldswap(&arr2[l], &arr2[r]);
	ldswap(&arr3[l], &arr3[r]);
  	ldswap(&arr4[l], &arr4[r]);
    }
    }
    ldswap(&arr[--l], &arr[beg]);
    ldswap(&arr2[l], &arr2[beg]);
    ldswap(&arr3[l], &arr3[beg]);
    ldswap(&arr4[l], &arr4[beg]);
    ldldldldsort(arr, arr2, arr3, arr4, beg, l);
    ldldldldsort(arr, arr2, arr3, arr4, r, end);
  }
}

void isort(int arr[], long double arr2[], long double arr3[], int beg, int end) {
  if (end > beg + 1) {
    int piv = arr[beg];
    int l = beg + 1, r = end;
    
    while (l < r) {
      if (arr[l] <= piv)
	l++;
      else {
	iswap(&arr[l], &arr[--r]);
	ldswap(&arr2[l], &arr2[r]);
	ldswap(&arr3[l], &arr3[r]);
      }
    }
    iswap(&arr[--l], &arr[beg]);
    ldswap(&arr2[l], &arr2[beg]);
    ldswap(&arr3[l], &arr3[beg]);
    isort(arr, arr2, arr3, beg, l);
    isort(arr, arr2, arr3, r, end);
  }
}

//shuffles elements of an int array
void shuffle(int *arr, int imax) {
  int i, swapi;
  int tmp;

  for(i=0; i<imax; i++) {
    swapi = i + (int)((imax-i)*(random()/((double)RAND_MAX+1.0)));

    tmp = arr[i];
    arr[i] = arr[swapi];
    arr[swapi] = tmp;
  }
}

//returns the number of zero entries in the 2d double matrix
int zeroCount2dDouArr(double **mat, int nrow, int ncol) {
  int i, j;
  int nzeros = 0;

  for(i=0; i<nrow; i++) {
    for(j=0; j<ncol; j++) {
      if(mat[i][j] == 0.0) {
	nzeros++;
      }
    }
  }

  return nzeros;
}
