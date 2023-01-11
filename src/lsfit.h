#ifndef __LSFIT__
#define __LSFIT__
/*these are taken from TPCCLIB http://www.turkupetcentre.net/tpcclib-doc/v2/ */
int nnls_c(double* a, const int* mda, const int* m, const int* n, double* b,
         double* x, double* rnorm, double* w, double* zz, int* index,
         int* mode); 
/*****************************************************************************/
int nnls(
  double **a,
  int m, 
  int n, 
  double *b,
  double *x, 
  double *rnorm,
  double *wp,
  double *zzp,
  int *indexp
); 
    
int nnlsWght(
  int N,
  int M,
  double **A,
  double *b,
  double *weight
);

int nnlsWghtSquared(
  int N,
  int M,
  double **A,
  double *b,
  double *sweight
); 


int bvls(
  int key,
  const /*unsigned*/ int m,
  const /*unsigned*/ int n,
  double *a,
  double *b,
  double *bl,
  double *bu,
  double *x,
  double *w,
  double *act,
  double *zz,
  int *istate,
  int *iter,
  int verbose
);

int llsqWght(
  int N,
  int M,
  double **A,
  double *a,
  double *b,
  double *weight
);

int llsqWghtSquared(
  int N,
  int M,
  double **A,
  double *a,
  double *b,
  double *sweight
); 

#endif
