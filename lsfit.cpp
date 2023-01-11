/*****************************************************************************/
//#include "tpcclibConfig.h"
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*****************************************************************************/
#include "lsfit.h"
/*****************************************************************************/
/* The routines in this file have been translated from Fortran to C by
   f2c. Additional modifications have been made to remove the
   dependencies on the f2c header file and library. The original
   Fortran 77 code accompanies the SIAM Publications printing of
   "Solving Least Squares Problems," by C. Lawson and R. Hanson and is
   freely available at www.netlib.org/lawson-hanson/all. */

/* nnls.F -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* The next line was removed after the f2c translation */
/* #include "f2c.h" */

/* The next lines were added after the f2c translation. Also swapped
   abs for nnls_abs and max for nnls_max to avoid confusion with some
   compilers. */
#include <stdio.h>
#include <math.h>
//#include "nnls.h"
#define nnls_max(a,b) ((a) >= (b) ? (a) : (b))
#define nnls_abs(x) ((x) >= 0 ? (x) : -(x))
typedef int integer;
typedef double doublereal;

/* The following subroutine was added after the f2c translation */
double d_sign(double *a, double *b)
{
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;

/*#define c__1 1
#define c__0 0
#define c__2 2*/

// minor edits by Nicolas Bonneel
int h12_(integer *mode, integer *lpivot, integer *l1, const integer *m, doublereal *u, integer *iue, doublereal *up, doublereal *c__, integer *ice, const integer *icv, integer *ncv);
doublereal diff_(doublereal* x, doublereal* y);
int g1_(doublereal *a, doublereal *b, doublereal *cterm, doublereal *sterm, doublereal *sig);


/*     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) */

/*  Algorithm NNLS: NONNEGATIVE LEAST SQUARES */

/*  The original version of this code was developed by */
/*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory */
/*  1973 JUN 15, and published in the book */
/*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
/*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */

/*     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN */
/*     N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM */

/*                      A * X = B  SUBJECT TO X .GE. 0 */
/*     ------------------------------------------------------------------ */
/*                     Subroutine Arguments */

/*     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE */
/*                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N */
/*                     MATRIX, A.           ON EXIT A() CONTAINS */
/*                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN */
/*                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY */
/*                     THIS SUBROUTINE. */
/*     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- */
/*             TAINS Q*B. */
/*     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL */
/*             CONTAIN THE SOLUTION VECTOR. */
/*     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE */
/*             RESIDUAL VECTOR. */
/*     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN */
/*             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0. */
/*             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z */
/*     ZZ()     AN M-ARRAY OF WORKING SPACE. */
/*     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N. */
/*                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS */
/*                 P AND Z AS FOLLOWS.. */

/*                 INDEX(1)   THRU INDEX(NSETP) = SET P. */
/*                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z. */
/*                 IZ1 = NSETP + 1 = NPP1 */
/*                 IZ2 = N */
/*     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING */
/*             MEANINGS. */
/*             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY. */
/*             2     THE DIMENSIONS OF THE PROBLEM ARE BAD. */
/*                   EITHER M .LE. 0 OR N .LE. 0. */
/*             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. */

/*     ------------------------------------------------------------------ */
/* Subroutine */ int nnls_(doublereal *a, const integer *mda, const integer *m, const integer *n, doublereal* b, doublereal* x, doublereal* rnorm, doublereal* w, doublereal* zz, integer* index, integer* mode)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* The following lines were commented out after the f2c translation */
    /* double sqrt(); */
    /* icnteger s_wsfe(), do_fio(), e_wsfe(); */

    /* Local variables */
    //extern doublereal diff_();
    /*static*/ integer iter;
    /*static*/ doublereal temp, wmax;
    /*static*/ integer i__, j, l;
    /*static*/ doublereal t, alpha, asave;
    /*static*/ integer itmax, izmax, nsetp;
    //extern /* Subroutine */ int g1_();
    /*static*/ doublereal dummy, unorm, ztest, cc;
    //extern /* Subroutine */ int h12_();
    /*static*/ integer ii, jj, ip;
    /*static*/ doublereal sm;
    /*static*/ integer iz, jz;
    /*static*/ doublereal up, ss;
    /*static*/ integer rtnkey, iz1, iz2, npp1;

    /* Fortran I/O blocks */
    /* The following line was commented out after the f2c translation */
    /* static cilist io___22 = { 0, 6, 0, "(/a)", 0 }; */


/*     ------------------------------------------------------------------ 
*/
/*     integer INDEX(N) */
/*     double precision A(MDA,N), B(M), W(N), X(N), ZZ(M) */
/*     ------------------------------------------------------------------ 
*/
    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --b;
    --x;
    --w;
    --zz;
    --index;

    /* Function Body */
    *mode = 1;
    if (*m <= 0 || *n <= 0) {
	*mode = 2;
	return 0;
    }
    iter = 0;
    itmax = *n * 1; //3

/*                    INITIALIZE THE ARRAYS INDEX() AND X(). */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = 0.;
/* L20: */
	index[i__] = i__;
    }

    iz2 = *n;
    iz1 = 1;
    nsetp = 0;
    npp1 = 1;
/*                             ******  MAIN LOOP BEGINS HERE  ****** */
L30:
/*                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION. 
*/
/*                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED. */

    if (iz1 > iz2 || nsetp >= *m) {
	goto L350;
    }

/*         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W(). 
*/

    i__1 = iz2;
    for (iz = iz1; iz <= i__1; ++iz) {
	j = index[iz];
	sm = 0.;
	i__2 = *m;
	for (l = npp1; l <= i__2; ++l) {
/* L40: */
	    sm += a[l + j * a_dim1] * b[l];
	}
	w[j] = sm;
/* L50: */
    }
/*                                   FIND LARGEST POSITIVE W(J). */
L60:
    wmax = 0.;
    i__1 = iz2;
    for (iz = iz1; iz <= i__1; ++iz) {
	j = index[iz];
	if (w[j] > wmax) {
	    wmax = w[j];
	    izmax = iz;
	}
/* L70: */
    }

/*             IF WMAX .LE. 0. GO TO TERMINATION. */
/*             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS. 
*/

    if (wmax <= 0.) {
	goto L350;
    }
    iz = izmax;
    j = index[iz];

/*     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P. */
/*     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID */
/*     NEAR LINEAR DEPENDENCE. */

    asave = a[npp1 + j * a_dim1];
    i__1 = npp1 + 1;
    h12_(&c__1, &npp1, &i__1, m, &a[j * a_dim1 + 1], &c__1, &up, &dummy, &
	    c__1, &c__1, &c__0);
    unorm = 0.;
    if (nsetp != 0) {
	i__1 = nsetp;
	for (l = 1; l <= i__1; ++l) {
/* L90: */
/* Computing 2nd power */
	    d__1 = a[l + j * a_dim1];
	    unorm += d__1 * d__1;
	}
    }
    unorm = sqrt(unorm);
    d__2 = unorm + (d__1 = a[npp1 + j * a_dim1], nnls_abs(d__1)) * .01;
    if (diff_(&d__2, &unorm) > 0.) {

/*        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE Z
Z */
/*        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ). */

	i__1 = *m;
	for (l = 1; l <= i__1; ++l) {
/* L120: */
	    zz[l] = b[l];
	}
	i__1 = npp1 + 1;
	h12_(&c__2, &npp1, &i__1, m, &a[j * a_dim1 + 1], &c__1, &up, &zz[1], &
		c__1, &c__1, &c__1);
	ztest = zz[npp1] / a[npp1 + j * a_dim1];

/*                                     SEE IF ZTEST IS POSITIVE */

	if (ztest > 0.) {
	    goto L140;
	}
    }

/*     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P. */
/*     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL */
/*     COEFFS AGAIN. */

    a[npp1 + j * a_dim1] = asave;
    w[j] = 0.;
    goto L60;

/*     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM */
/*     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER */
/*     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN */
/*     COL J,  SET W(J)=0. */

L140:
    i__1 = *m;
    for (l = 1; l <= i__1; ++l) {
/* L150: */
	b[l] = zz[l];
    }

    index[iz] = index[iz1];
    index[iz1] = j;
    ++iz1;
    nsetp = npp1;
    ++npp1;

    if (iz1 <= iz2) {
	i__1 = iz2;
	for (jz = iz1; jz <= i__1; ++jz) {
	    jj = index[jz];
	    h12_(&c__2, &nsetp, &npp1, m, &a[j * a_dim1 + 1], &c__1, &up, &a[
		    jj * a_dim1 + 1], &c__1, mda, &c__1);
/* L160: */
	}
    }

    if (nsetp != *m) {
	i__1 = *m;
	for (l = npp1; l <= i__1; ++l) {
/* L180: */
	    a[l + j * a_dim1] = 0.;
	}
    }

    w[j] = 0.;
/*                                SOLVE THE TRIANGULAR SYSTEM. */
/*                                STORE THE SOLUTION TEMPORARILY IN ZZ(). 
*/
    rtnkey = 1;
    goto L400;
L200:

/*                       ******  SECONDARY LOOP BEGINS HERE ****** */

/*                          ITERATION COUNTER. */

L210:
    ++iter;
    if (iter > itmax) {
	*mode = 3;
	/* The following lines were replaced after the f2c translation */
	/* s_wsfe(&io___22); */
	/* do_fio(&c__1, " NNLS quitting on iteration count.", 34L); */
	/* e_wsfe(); */
	fprintf(stdout, "\n NNLS quitting on iteration count.\n");
	fflush(stdout);
	goto L350;
    }

/*                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE. */
/*                                  IF NOT COMPUTE ALPHA. */

    alpha = 2.;
    i__1 = nsetp;
    for (ip = 1; ip <= i__1; ++ip) {
	l = index[ip];
	if (zz[ip] <= 0.) {
	    t = -x[l] / (zz[ip] - x[l]);
	    if (alpha > t) {
		alpha = t;
		jj = ip;
	    }
	}
/* L240: */
    }

/*          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL */
/*          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP. */

    if (alpha == 2.) {
	goto L330;
    }

/*          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO */
/*          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ. */

    i__1 = nsetp;
    for (ip = 1; ip <= i__1; ++ip) {
	l = index[ip];
	x[l] += alpha * (zz[ip] - x[l]);
/* L250: */
    }

/*        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I */
/*        FROM SET P TO SET Z. */

    i__ = index[jj];
L260:
    x[i__] = 0.;

    if (jj != nsetp) {
	++jj;
	i__1 = nsetp;
	for (j = jj; j <= i__1; ++j) {
	    ii = index[j];
	    index[j - 1] = ii;
	    g1_(&a[j - 1 + ii * a_dim1], &a[j + ii * a_dim1], &cc, &ss, &a[j 
		    - 1 + ii * a_dim1]);
	    a[j + ii * a_dim1] = 0.;
	    i__2 = *n;
	    for (l = 1; l <= i__2; ++l) {
		if (l != ii) {

/*                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,
L)) */

		    temp = a[j - 1 + l * a_dim1];
		    a[j - 1 + l * a_dim1] = cc * temp + ss * a[j + l * a_dim1]
			    ;
		    a[j + l * a_dim1] = -ss * temp + cc * a[j + l * a_dim1];
		}
/* L270: */
	    }

/*                 Apply procedure G2 (CC,SS,B(J-1),B(J)) */

	    temp = b[j - 1];
	    b[j - 1] = cc * temp + ss * b[j];
	    b[j] = -ss * temp + cc * b[j];
/* L280: */
	}
    }

    npp1 = nsetp;
    --nsetp;
    --iz1;
    index[iz1] = i__;

/*        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD 
*/
/*        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED. */
/*        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY */
/*        THAT ARE NONPOSITIVE WILL BE SET TO ZERO */
/*        AND MOVED FROM SET P TO SET Z. */

    i__1 = nsetp;
    for (jj = 1; jj <= i__1; ++jj) {
	i__ = index[jj];
	if (x[i__] <= 0.) {
	    goto L260;
	}
/* L300: */
    }

/*         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK. */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L310: */
	zz[i__] = b[i__];
    }
    rtnkey = 2;
    goto L400;
L320:
    goto L210;
/*                      ******  END OF SECONDARY LOOP  ****** */

L330:
    i__1 = nsetp;
    for (ip = 1; ip <= i__1; ++ip) {
	i__ = index[ip];
/* L340: */
	x[i__] = zz[ip];
    }
/*        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING. */
    goto L30;

/*                        ******  END OF MAIN LOOP  ****** */

/*                        COME TO HERE FOR TERMINATION. */
/*                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR. */

L350:
    sm = 0.;
    if (npp1 <= *m) {
	i__1 = *m;
	for (i__ = npp1; i__ <= i__1; ++i__) {
/* L360: */
/* Computing 2nd power */
	    d__1 = b[i__];
	    sm += d__1 * d__1;
	}
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* L380: */
	    w[j] = 0.;
	}
    }
    *rnorm = sqrt(sm);
    return 0;

/*     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE */
/*     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ(). */

L400:
    i__1 = nsetp;
    for (l = 1; l <= i__1; ++l) {
	ip = nsetp + 1 - l;
	if (l != 1) {
	    i__2 = ip;
	    for (ii = 1; ii <= i__2; ++ii) {
		zz[ii] -= a[ii + jj * a_dim1] * zz[ip + 1];
/* L410: */
	    }
	}
	jj = index[ip];
	zz[ip] /= a[ip + jj * a_dim1];
/* L430: */
    }
    switch ((int)rtnkey) {
	case 1:  goto L200;
	case 2:  goto L320;
    }

    /* The next line was added after the f2c translation to keep
       compilers from complaining about a void return from a non-void
       function. */
    return 0;

} /* nnls_ */

/* Subroutine */ int g1_(doublereal *a, doublereal *b, doublereal *cterm, doublereal *sterm, doublereal *sig)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    /* The following line was commented out after the f2c translation */
    /* double sqrt(), d_sign(); */

    /* Local variables */
    /*static*/ doublereal xr, yr;


/*     COMPUTE ORTHOGONAL ROTATION MATRIX.. */

/*  The original version of this code was developed by */
/*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory 
*/
/*  1973 JUN 12, and published in the book */
/*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
/*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */

/*     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2)) */
/*                        (-S,C)         (-S,C)(B)   (   0          ) */
/*     COMPUTE SIG = SQRT(A**2+B**2) */
/*        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT */
/*        SIG MAY BE IN THE SAME LOCATION AS A OR B . */
/*     ------------------------------------------------------------------ 
*/
/*     ------------------------------------------------------------------ 
*/
    if (nnls_abs(*a) > nnls_abs(*b)) {
	xr = *b / *a;
/* Computing 2nd power */
	d__1 = xr;
	yr = sqrt(d__1 * d__1 + 1.);
	d__1 = 1. / yr;
	*cterm = d_sign(&d__1, a);
	*sterm = *cterm * xr;
	*sig = nnls_abs(*a) * yr;
	return 0;
    }
    if (*b != 0.) {
	xr = *a / *b;
/* Computing 2nd power */
	d__1 = xr;
	yr = sqrt(d__1 * d__1 + 1.);
	d__1 = 1. / yr;
	*sterm = d_sign(&d__1, b);
	*cterm = *sterm * xr;
	*sig = nnls_abs(*b) * yr;
	return 0;
    }
    *sig = 0.;
    *cterm = 0.;
    *sterm = 1.;
    return 0;
} /* g1_ */

/*     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV) */

/*  CONSTRUCTION AND/OR APPLICATION OF A SINGLE */
/*  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B */

/*  The original version of this code was developed by */
/*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory */
/*  1973 JUN 12, and published in the book */
/*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
/*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */
/*     ------------------------------------------------------------------ */
/*                     Subroutine Arguments */

/*     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a */
/*            Householder transformation, or Algorithm H2 to apply a */
/*            previously constructed transformation. */
/*     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. */
/*     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO */
/*            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M */
/*            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. */
/*     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot */
/*            vector.  IUE is the storage increment between elements. */
/*            On exit when MODE = 1, U() and UP contain quantities */
/*            defining the vector U of the Householder transformation. */
/*            on entry with MODE = 2, U() and UP should contain */
/*            quantities previously computed with MODE = 1.  These will */
/*            not be modified during the entry with MODE = 2. */
/*     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH */
/*            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE */
/*            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED. */
/*            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS. */
/*     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C(). */
/*     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C(). */
/*     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0 */
/*            NO OPERATIONS WILL BE DONE ON C(). */
/*     ------------------------------------------------------------------ */
/* Subroutine */ int h12_(integer *mode, integer *lpivot, integer *l1, const integer *m, doublereal *u, integer *iue, doublereal *up, doublereal *c__, integer *ice, const integer *icv, integer *ncv)
{
    /* System generated locals */
    integer u_dim1, u_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* The following line was commented out after the f2c translation */
    /* double sqrt(); */

    /* Local variables */
    /*static*/ integer incr;
    /*static*/ doublereal b;
    /*static*/ integer i__, j;
    /*static*/ doublereal clinv;
    /*static*/ integer i2, i3, i4;
    /*static*/ doublereal cl, sm;

/*     ------------------------------------------------------------------ 
*/
/*     double precision U(IUE,M) */
/*     ------------------------------------------------------------------ 
*/
    /* Parameter adjustments */
    u_dim1 = *iue;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    --c__;

    /* Function Body */
    if (0 >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
	return 0;
    }
    cl = (d__1 = u[*lpivot * u_dim1 + 1], nnls_abs(d__1));
    if (*mode == 2) {
	goto L60;
    }
/*                            ****** CONSTRUCT THE TRANSFORMATION. ****** 
*/
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* L10: */
/* Computing MAX */
	d__2 = (d__1 = u[j * u_dim1 + 1], nnls_abs(d__1));
	cl = nnls_max(d__2,cl);
    }
    if (cl <= 0.) {
	goto L130;
    } else {
	goto L20;
    }
L20:
    clinv = 1. / cl;
/* Computing 2nd power */
    d__1 = u[*lpivot * u_dim1 + 1] * clinv;
    sm = d__1 * d__1;
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* L30: */
/* Computing 2nd power */
	d__1 = u[j * u_dim1 + 1] * clinv;
	sm += d__1 * d__1;
    }
    cl *= sqrt(sm);
    if (u[*lpivot * u_dim1 + 1] <= 0.) {
	goto L50;
    } else {
	goto L40;
    }
L40:
    cl = -cl;
L50:
    *up = u[*lpivot * u_dim1 + 1] - cl;
    u[*lpivot * u_dim1 + 1] = cl;
    goto L70;
/*            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ****** 
*/

L60:
    if (cl <= 0.) {
	goto L130;
    } else {
	goto L70;
    }
L70:
    if (*ncv <= 0) {
	return 0;
    }
    b = *up * u[*lpivot * u_dim1 + 1];
/*                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN. 
*/

    if (b >= 0.) {
	goto L130;
    } else {
	goto L80;
    }
L80:
    b = 1. / b;
    i2 = 1 - *icv + *ice * (*lpivot - 1);
    incr = *ice * (*l1 - *lpivot);
    i__1 = *ncv;
    for (j = 1; j <= i__1; ++j) {
	i2 += *icv;
	i3 = i2 + incr;
	i4 = i3;
	sm = c__[i2] * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    sm += c__[i3] * u[i__ * u_dim1 + 1];
/* L90: */
	    i3 += *ice;
	}
	if (sm != 0.) {
	    goto L100;
	} else {
	    goto L120;
	}
L100:
	sm *= b;
	c__[i2] += sm * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    c__[i4] += sm * u[i__ * u_dim1 + 1];
/* L110: */
	    i4 += *ice;
	}
L120:
	;
    }
L130:
    return 0;
} /* h12_ */

doublereal diff_(doublereal* x, doublereal* y)
{
    /* System generated locals */
    doublereal ret_val;


/*  Function used in tests that depend on machine precision. */

/*  The original version of this code was developed by */
/*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory 
*/
/*  1973 JUN 7, and published in the book */
/*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
/*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */

    ret_val = *x - *y;
    return ret_val;
} /* diff_ */


/* The following subroutine was added after the f2c translation */
int nnls_c(double* a, const int* mda, const int* m, const int* n, double* b, 
	 double* x, double* rnorm, double* w, double* zz, int* index, 
	 int* mode)
{
  return (nnls_(a, mda, m, n, b, x, rnorm, w, zz, index, mode));
}
 

/*****************************************************************************/
/* Local function definitions */
int _lss_h12(
  int mode, int lpivot, int l1, int m, double *u, int iue,
  double *up, double *cm, int ice, int icv, int ncv
);
void _lss_g1(
  double a, double b, double *cterm, double *sterm, double *sig
);
/*****************************************************************************/
/*****************************************************************************/
 
int _lss_h12(
  int mode,
  int lpivot,
  int l1,
  int m,
  double *u,
  int u_dim1,
  double *up,
  double *cm,
  int ice,
  int icv,
  int ncv
) {
  /* Check parameters */
  if(mode!=1 && mode!=2) return(1);
  if(m<1 || u==NULL || u_dim1<1 /*|| cm==NULL*/) return(1);
  if(lpivot<0 || lpivot>=l1 || l1>m) return(1);
 
  double cl = fabs(u[lpivot*u_dim1]);
 
  if(mode==2) { /* Apply transformation I+U*(U**T)/B to cm[] */
    if(cl<=0.) return(0);
  } else {   /* Construct the transformation */
  
    /* trying to compensate overflow */
    for(int j=l1; j<m; j++) {  // Computing MAX 
      cl = fmax(fabs(u[j*u_dim1]), cl);
    }
    // zero vector?   
    if(cl<=0.) return(0);
 
    double clinv=1.0/cl;
       
    // cl = sqrt( (u[pivot]*clinv)^2 + sigma(i=l1..m)( (u[i]*clinv)^2 ) )
    double d1=u[lpivot*u_dim1]*clinv; 
    double sm=d1*d1;
    for(int j=l1; j<m; j++) {
      double d2=u[j*u_dim1]*clinv;
      sm+=d2*d2;
    }
    cl*=sqrt(sm);
    if(u[lpivot*u_dim1] > 0.) cl=-cl;
    *up = u[lpivot*u_dim1] - cl; 
    u[lpivot*u_dim1]=cl;
  }
 
  // no vectors where to apply? only change pivot vector! 
  double b=(*up)*u[lpivot*u_dim1];
  
  /* b must be non-positive here; if b>=0., then return */
  if(b>=0.0) return(0); // was if(b==0) before 2013-06-22
  
  // ok, for all vectors we want to apply
  if(cm==NULL) return(2);
  for(int j=0; j<ncv; j++) {
    // take s = c[p,j]*h + sigma(i=l..m){ c[i,j] *v [ i ] }
    double sm = cm[ lpivot*ice + j*icv ] * (*up);
    for(int k=l1; k<m; k++) sm += cm[ k * ice + j*icv ] * u[ k*u_dim1 ]; 
    if(sm!=0.0) {
      sm *= (1.0/b); // was (1/b) before 2013-06-22
      // cm[lpivot, j] = ..
      cm[ lpivot * ice + j*icv] += sm*(*up);
      // for i = l1...m , set c[i,j] = c[i,j] + s*v[i]
      for(int k=l1; k<m; k++) cm[ k*ice + j*icv] += u[k * u_dim1]*sm;
    }
  }
   
  return(0);
} /* _lss_h12 */
/*****************************************************************************/

/*****************************************************************************/
void _lss_g1(double a, double b, double *cterm, double *sterm, double *sig)
{
  double d1, xr, yr;
 
  if(fabs(a)>fabs(b)) {
    xr=b/a; d1=xr; yr=hypot(d1, 1.0); d1=1./yr;
    *cterm=copysign(d1, a);
    *sterm=(*cterm)*xr; *sig=fabs(a)*yr;
  } else if(b!=0.) {
    xr=a/b; d1=xr; yr=hypot(d1, 1.0); d1=1./yr;
    *sterm=copysign(d1, b);
    *cterm=(*sterm)*xr; *sig=fabs(b)*yr;
  } else {
    *sig=0.; *cterm=0.; *sterm=1.;
  }
} /* _lss_g1 */
/*****************************************************************************/
 

//a	On entry, a is M by N matrix A. On exit, a = Q*a, where Q is an m by n orthogonal matrix generated implicitly by this function.
//m	Matrix dimension m
//n	Matrix dimension n
//b	On entry, b is an m-vector B. On exit, b contains Q*B
//x	On exit, x will contain the solution, an n-vector
//rnorm	On exit, rnorm contains the squared Euclidean norm of the residual vector, R^2. If NULL is given, no rnorm is calculated
//wp	An n-array of working space, wp[]. On exit, wp[] will contain the dual solution vector. 
//      wp[i]=0.0 for all i in set p and wp[i]<=0.0 for all i in set z. Can be NULL, which causes this algorithm to allocate memory for it.
//zzp	An m-array of working space, zz[]. Can be NULL, which causes this algorithm to allocate memory for it.
//indexp	An n-array of working space, index[]. Can be NULL, which causes this algorithm to allocate memory for it.
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
) {
  /* Check the parameters and data */
  if(m<=0 || n<=0 || a==NULL || b==NULL || x==NULL) return(2);
  /* Allocate memory for working space, if required */
  int *index=NULL;
  double *w=NULL, *zz=NULL;
  if(wp!=NULL) w=wp; else w=(double*)calloc(n, sizeof(double));
  if(zzp!=NULL) zz=zzp; else zz=(double*)calloc(m, sizeof(double));
  if(indexp!=NULL) index=indexp; else index=(int*)calloc(n, sizeof(int));
  if(w==NULL || zz==NULL || index==NULL) return(2);
 
  /* Initialize the arrays INDEX[] and X[] */
  for(int ni=0; ni<n; ni++) {x[ni]=0.; index[ni]=ni;}
  int iz1=0;
  int iz2=n-1;
  int nsetp=0;
  int npp1=0;
 
  /* Main loop; quit if all coefficients are already in the solution or
     if M cols of A have been triangulated */
  double up=0.0;
  int itmax; if(n<3) itmax=n*3; else itmax=n*n;
  int iter=0; 
  int k, j=0, jj=0;
  while(iz1<=iz2 && nsetp<m) {
    /* Compute components of the dual (negative gradient) vector W[] */
    for(int iz=iz1; iz<=iz2; iz++) {
      int ni=index[iz];
      double sm=0.;
      for(int mi=npp1; mi<m; mi++) sm+=a[ni][mi]*b[mi];
      w[ni]=sm;
    }
 
    double wmax;
    int izmax=0;
    while(1) {
 
      /* Find largest positive W[j] */
      wmax=0.0;
      for(int iz=iz1; iz<=iz2; iz++) {
        int i=index[iz];
        if(w[i]>wmax) {wmax=w[i]; izmax=iz;}
      }
 
      /* Terminate if wmax<=0.; */
      /* it indicates satisfaction of the Kuhn-Tucker conditions */
      if(wmax<=0.0) break;
      j=index[izmax];
 
      /* The sign of W[j] is ok for j to be moved to set P.
         Begin the transformation and check new diagonal element to avoid
         near linear dependence. */
      double asave=a[j][npp1];
      up=0.0;
      _lss_h12(1, npp1, npp1+1, m, &a[j][0], 1, &up, NULL, 1, 1, 0);
      double unorm=0.0;
      if(nsetp!=0) for(int mi=0; mi<nsetp; mi++) unorm+=a[j][mi]*a[j][mi];
      unorm=sqrt(unorm);
      double d=unorm+fabs(a[j][npp1])*0.01;
      if((d-unorm)>0.0) {
        /* Col j is sufficiently independent. Copy B into ZZ, update ZZ
           and solve for ztest ( = proposed new value for X[j] ) */
        for(int mi=0; mi<m; mi++) zz[mi]=b[mi];
        _lss_h12(2, npp1, npp1+1, m, &a[j][0], 1, &up, zz, 1, 1, 1);
        double ztest=zz[npp1]/a[j][npp1];
        /* See if ztest is positive */
        if(ztest>0.) break;
      }
 
      /* Reject j as a candidate to be moved from set Z to set P. Restore
         A[npp1,j], set W[j]=0., and loop back to test dual coefficients again */
      a[j][npp1]=asave; w[j]=0.;
    } /* while(1) */
    if(wmax<=0.0) break;
 
    /* Index j=INDEX[izmax] has been selected to be moved from set Z to set P.
       Update B and indices, apply householder transformations to cols in
       new set Z, zero sub-diagonal elements in col j, set W[j]=0. */
    for(int mi=0; mi<m; mi++) b[mi]=zz[mi];
    index[izmax]=index[iz1]; index[iz1]=j; iz1++; nsetp=npp1+1; npp1++;
    if(iz1<=iz2)
      for(int jz=iz1; jz<=iz2; jz++) {
        jj=index[jz];
        _lss_h12(2, nsetp-1, npp1, m, &a[j][0], 1, &up, &a[jj][0], 1, m, 1);
      }
    if(nsetp!=m) for(int mi=npp1; mi<m; mi++) a[j][mi]=0.;
    w[j]=0.;
    
    /* Solve the triangular system; store the solution temporarily in Z[] */
    for(int mi=0; mi<nsetp; mi++) {
      int ip=nsetp-(mi+1);
      if(mi!=0) for(int ii=0; ii<=ip; ii++) zz[ii]-=a[jj][ii]*zz[ip+1];
      jj=index[ip]; zz[ip]/=a[jj][ip];
    }
 
    /* Secondary loop begins here */
    while(++iter<itmax) {
      /* See if all new constrained coefficients are feasible; if not, compute alpha */
      double alpha=2.0;
      for(int ip=0; ip<nsetp; ip++) {
        int ni=index[ip];
        if(zz[ip]<=0.) {
          double t=-x[ni]/(zz[ip]-x[ni]);
          if(alpha>t) {alpha=t; jj=ip-1;}
        }
      }
 
      /* If all new constrained coefficients are feasible then still alpha==2.
         If so, then exit from the secondary loop to main loop */
      if(alpha==2.0) break;
 
      /* Use alpha (0.<alpha<1.) to interpolate between old X and new ZZ */
      for(int ip=0; ip<nsetp; ip++) {
        int ni=index[ip]; x[ni]+=alpha*(zz[ip]-x[ni]);
      }
 
      /* Modify A and B and the INDEX arrays to move coefficient i from set P to set Z. */
      int pfeas=1;
      k=index[jj+1];
      do {
        x[k]=0.;
        if(jj!=(nsetp-1)) {
          jj++;
          for(int ni=jj+1; ni<nsetp; ni++) {
            int ii=index[ni]; index[ni-1]=ii;
            double ss, cc;
            _lss_g1(a[ii][ni-1], a[ii][ni], &cc, &ss, &a[ii][ni-1]);
            a[ii][ni]=0.0;
            for(int nj=0; nj<n; nj++) if(nj!=ii) {
              /* Apply procedure G2 (CC,SS,A(J-1,L),A(J,L)) */
              double temp=a[nj][ni-1];
              a[nj][ni-1]=cc*temp+ss*a[nj][ni];
              a[nj][ni]=-ss*temp+cc*a[nj][ni];
            }
            /* Apply procedure G2 (CC,SS,B(J-1),B(J)) */
            double temp=b[ni-1]; b[ni-1]=cc*temp+ss*b[ni]; b[ni]=-ss*temp+cc*b[ni];
          }
        }
        npp1=nsetp-1; nsetp--; iz1--; index[iz1]=k;
 
        /* See if the remaining coefficients in set P are feasible; they should be
           because of the way alpha was determined. If any are infeasible 
           it is due to round-off error. Any that are non-positive 
           will be set to zero and moved from set P to set Z. */
        for(jj=0, pfeas=1; jj<nsetp; jj++) {
          k=index[jj]; if(x[k]<=0.) {pfeas=0; break;}
        }
      } while(pfeas==0);
 
      /* Copy B[] into zz[], then solve again and loop back */
      for(int mi=0; mi<m; mi++) zz[mi]=b[mi];
      for(int mi=0; mi<nsetp; mi++) {
        int ip=nsetp-(mi+1);
        if(mi!=0) for(int ii=0; ii<=ip; ii++) zz[ii]-=a[jj][ii]*zz[ip+1];
        jj=index[ip]; zz[ip]/=a[jj][ip];
      }
    } /* end of secondary loop */
 
    if(iter>=itmax) break;
    for(int ip=0; ip<nsetp; ip++) {k=index[ip]; x[k]=zz[ip];}
  } /* end of main loop */
 
  /* Compute the norm of the final residual vector */
  if(rnorm != NULL) {
    double sm=0.0;
    if(npp1<m) for(int mi=npp1; mi<m; mi++) sm+=(b[mi]*b[mi]);
    else for(int ni=0; ni<n; ni++) w[ni]=0.;
    *rnorm=sm; //sqrt(sm);
  }
 
  /* Free working space, if it was allocated here */
  if(wp==NULL) free(w); 
  if(zzp==NULL) free(zz); 
  if(indexp==NULL) free(index);
  if(iter>=itmax) return(1);
  return(0);
} /* nnls */
/*****************************************************************************/
 
/*****************************************************************************/
int nnlsWght(
  int N,
  int M,
  double **A,
  double *b,
  double *weight
) {
  int n, m;
  double *w;
 
  /* Check the arguments */
  if(N<1 || M<1 || A==NULL || b==NULL || weight==NULL) return(1);
 
  /* Allocate memory */
  w=(double*)malloc(M*sizeof(double)); if(w==NULL) return(2);
 
  /* Check that weights are not zero and get the square roots of them to w[] */
  for(m=0; m<M; m++) {
    if(weight[m]<=1.0e-20) w[m]=0.0;
    else w[m]=sqrt(weight[m]);
  }
 
  /* Multiply rows of matrix A and elements of vector b with weights*/
  for(m=0; m<M; m++) {
    for(n=0; n<N; n++) {
      A[n][m]*=w[m];
    }
    b[m]*=w[m];
  }
 
  free(w);
  return(0);
}
/*****************************************************************************/
 
/*****************************************************************************/
int nnlsWghtSquared(
  int N,
  int M,
  double **A,
  double *b,
  double *sweight
) {
  int n, m;
 
  /* Check the arguments */
  if(N<1 || M<1 || A==NULL || b==NULL || sweight==NULL) return(1);
 
  /* Multiply rows of matrix A and elements of vector b with weights*/
  for(m=0; m<M; m++) {
    for(n=0; n<N; n++) {
      A[n][m]*=sweight[m];
    }
    b[m]*=sweight[m];
  }
 
  return(0);
}
/*****************************************************************************/
 
 
 
/*****************************************************************************/
int qrLH(
  const unsigned int m, 
  const unsigned int n,
  double *a,
  double *b,
  double *x,
  double *r2
) {
  /* Check the input */
  if(a==NULL || b==NULL || x==NULL || r2==NULL) return(2);
  if(n<1 || m<n) {*r2=nan(""); return(2);}
 
  /* Initiate output to zeroes, in case of exit because of singularity */
  for(unsigned int ni=0; ni<n; ni++) x[ni]=0.0;
  *r2=0.0;
 
  /* Rotates matrix A into upper triangular form */
  for(unsigned int ni=0; ni<n; ni++) {
    /* Find constants for rotation and diagonal entry */
    double sq=0.0;
    for(unsigned int mi=ni; mi<m; mi++) sq+=a[mi + ni*m]*a[mi + ni*m];
    if(sq==0.0) return(1);
    double qv1=-copysign(sqrt(sq), a[ni + ni*m]);
    double u1=a[ni + ni*m] - qv1;
    a[ni + ni*m]=qv1;
    unsigned int ni1=ni+1;
    /*  Rotate the remaining columns of sub-matrix. */
    for(unsigned int nj=ni1; nj<n; nj++) {
      double dot=u1*a[ni + nj*m];
      for(unsigned int mi=ni1; mi<m; mi++)
        dot+=a[mi + nj*m] * a[mi + ni*m];
      double c=dot/fabs(qv1*u1);
      for(unsigned int mi=ni1; mi<m; mi++)
        a[mi + nj*m]-=c*a[mi + ni*m];
      a[ni + nj*m]-=c*u1;
    }
    /* Rotate vector B */
    double dot=u1*b[ni];
    for(unsigned int mi=ni1; mi<m; mi++)
      dot+=b[mi]*a[mi + ni*m];
    double c=dot/fabs(qv1*u1);
    b[ni]-=c*u1;
    for(unsigned int mi=ni1; mi<m; mi++)
      b[mi]-=c*a[mi + ni*m];
  } // end of rotation loop
 
  /* Solve triangular system by back-substitution. */
  for(unsigned int ni=0; ni<n; ni++) {
    int k=n-ni-1;
    double s=b[k];
    for(unsigned int nj=k+1; nj<n; nj++) 
      s-=a[k + nj*m] * x[nj];
    if(a[k + k*m]==0.0) return(1);
    x[k]=s/a[k + k*m];
  }
 
  /* Calculate the sum of squared residuals. */
  *r2=0.0;
  for(unsigned int mi=n; mi<m; mi++)
    *r2 += b[mi]*b[mi];
 
  return(0);
}
/*****************************************************************************/
 

//key	Enter 0 to solve the problem from scratch, or <>0 to initialize the routine using caller's guess about which parameters are active within their bounds, which are at their lower bounds, and which at their upper bounds, as supplied in the array istate ('warm start'). When key <> 0, the routine initially sets the active components to the averages of their upper and lower bounds.
//m	Number of samples in matrix A and the length of vector b.
//n	Number of parameters in matrix A and the length of vector x. The n must not be > m unless at least n-m variables are non-active (set to their bounds).
//a	Pointer to matrix A; matrix must be given as an n*m array, containing n consecutive m-length vectors. Contents of A are modified in this routine.
//b	Pointer to vector b of length m. Contents of b are modified in this routine.
//bl	Array BL[0..n-1] of lower bounds for parameters.
//bu	Array BU[0..n-1] of upper bounds for parameters.
//x	Pointer to the result vector x of length n.
//w	Pointer to an array of length n, which will be used as working memory, and the minimum 2-norm || a.x-b ||, (R^2), will be written in w[0].
//act	Pointer to an array of length m*(n+2), or m*(m+2) if m<n, which will be used as working memory.
//zz	Pointer to an array of length m, to be used as working memory.
//istate	Pointer to an integer array of length n+1. If parameter key <>0, then the last position istate[n] must contain the total number of components at their bounds (nbound, the ‘bound variables’). The absolute values of the first nbound entries of istate[] are the indices of these ‘bound’ components of x[]. The sign of istate[0.. nbound-1] entries indicates whether x[|istate[i]|] is at its upper (positive) or lower (negative) bound. Entries istate[nbound..n-1] contain the indices of the active components.
//iter	Number of performed iterations is returned in this variable. Maximum nr of iterations is set with this variable, or set it to 0 to use the default maximum, 3*n.
//verbose	Verbose level; if zero, then nothing is printed to stderr or stdout.
/*****************************************************************************/
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
) {
  if(verbose>0) {printf("bvls(%d, %d, %d, ...)\n", key, m, n); fflush(stdout);}
  /* Check the input */
  if(a==NULL || b==NULL || bl==NULL || bu==NULL || x==NULL || w==NULL ||
     act==NULL || zz==NULL || istate==NULL || iter==NULL || n<1 || m<1) 
  {
    if(verbose>0) fprintf(stderr, "Error: invalid input to BVLS.\n");
    return(1);
  }
 
  int maxIter=*iter; if(maxIter<3) maxIter=3*n;
  *iter=0; 
  const double eps=1.0E-13; // stopping rule
 
  /* Step 1. Initialize everything -- active and bound sets, initial values, etc. */
  if(verbose>1) {printf("step 1\n"); fflush(stdout);}
 
  /* Set mm to the smaller of matrix dimensions n and m. */
  int mm; if(m<n) mm=m; else mm=n;
  int aindx=0; // one-based index of X[] that determined the alpha; zero means not set.
  int aindxsign=0; // +1 if zz[aindx]> bu[aindx], -1 if zz[aindx]<bl[aindx]
  /* istateFromStep5 is the one-based index of the parameter that most wants to be active;
     zero value means its is currently not set. */
  int istateFromStep5 = 0;
 
  /* Check the consistency of given bounds bl[] and bu[]. */
  {
    double maxrange=0.0;
    for(int ni=0; ni<n; ni++) {
      double d=bu[ni]-bl[ni];
      if(verbose>3) printf("  bounds[%d]: %g %g\n", 1+ni, bl[ni], bu[ni]);
      if(d<0.0) {
        if(verbose>0) fprintf(stderr, "Error: inconsistent bounds in BVLS.\n"); 
        return(1);
      }
      maxrange=fmax(maxrange, d);
    }
    if(verbose>2) printf("  maxrange := %g\n", maxrange);
    if(maxrange<1.0E-10) {
      if(verbose>0) fprintf(stderr, "Error: no free variables in BVLS.\n");
      return(1);
    }
  }
 
  /* In a fresh initialization (key=0), bind all variables at their lower bounds. 
     If key<>0, use the supplied istate[] array to initialize the variables. */
  int nbound, nact; // number of bound and active parameters, respectively.
  if(key==0) {
    nbound=n;
    /* Write the indices; negative sign indicates lower limit */
    /* Note that these indices must start from 1, because 0 can not have sign */
    for(int ni=0; ni<nbound; ni++) istate[ni]=-(1+ni);
  } else {
    nbound=istate[n];
  }
  nact=n-nbound;
  if(nact>mm) {
    if(verbose>0) fprintf(stderr, "Error: too many active variables in BVLS starting solution.\n");
    return(2);
  }
  for(int ni=0; ni<nbound; ni++) {
    int i=abs(istate[ni])-1;
    if(istate[ni]<0) x[i]=bl[i]; else x[i]=bu[i];
  }
 
  /* In a warm start (key<>0, and nbound<n) initialize the active variables to the mean of
     their bounds. This is needed in case the initial QR results in active variables 
     out-of-bounds and Steps 8-11 get executed the first time through. */
  for(int ni=nbound; ni<n; ni++) {
    int i=abs(istate[ni])-1; // indices of active variables should not have signs, but let's be sure
    x[i]=0.5*(bl[i]+bu[i]);
  }
 
  /* Compute bnorm, the norm of the data vector b, for reference. */
  double bnorm=0.0;
  for(int mi=0; mi<m; mi++) bnorm+=b[mi]*b[mi];
  bnorm=sqrt(bnorm); if(verbose>2) printf("  initial_bnorm := %g\n", bnorm);
 
 
  /*
   *  Main loop
   */
  int skipStep2=0;
  double obj=0.0;
  int iact=0; // The component x[iact] is the one that most wants to become active
 
  for(*iter=1; *iter<=maxIter; ++(*iter)) {
 
    if(verbose>1) {printf("iteration %d\n", *iter); fflush(stdout);}
 
    if(!skipStep2) {
      /* Step 2. */ if(verbose>1) {printf("  step 2\n"); fflush(stdout);}
 
      /*  Initialize the negative gradient vector w(*). */
      for(int ni=0; ni<n; ni++) w[ni]=0.0;
 
      /* Compute the residual vector b-a.x , the negative gradient vector w[*], and 
         the current objective value obj = || a.x - b ||.
         The residual vector is stored in the mm+1'st column of act[*,*]. */
      obj=0.0;
      for(int mi=0; mi<m; mi++) {
        double ri=b[mi];
        for(int ni=0; ni<n; ni++) ri-=a[mi+ni*m]*x[ni];
        obj+=ri*ri;
        for(int ni=0; ni<n; ni++) w[ni]+=a[mi+ni*m]*ri;
        act[mi+mm*m]=ri;
      }
      if(verbose>3) {printf("    obj := %g\n", obj); fflush(stdout);}
 
      /* Converged?  Stop if the misfit << || b ||, or if all components are active 
        (unless this is the first iteration from a 'warm start'). */
      if((sqrt(obj)<=bnorm*eps) || ((*iter)>1 && nbound==0)) {
        if(verbose>1) {printf("bvls converged.\n"); fflush(stdout);}
        istate[n]=nbound;
        w[0]=sqrt(obj);
        return(0);
      }
 
      /* Add the contribution of the active components back into the residual. */
      for(int ni=nbound; ni<n; ni++) {
        int i=abs(istate[ni])-1;
        for(int mi=0; mi<m; mi++) act[mi+mm*m] += a[mi+i*m]*x[i];
      }
      if(verbose>9) {
        printf("Residual vector:\n");
        for(int mi=0; mi<m; mi++) printf("\t%g", act[mi+mm*m]);
        printf("\n");
      }
 
    }
 
 
 
    /* The first iteration in a 'warm start' requires immediate QR in Step 6
       but mostly we want go through Steps 3-5 first . */
    if(key!=0 && (*iter)==1) {
      if(verbose>1) printf("  'warm start' requires immediate QR in Step 6\n");
    } else {
 
      int it; // variable indicating the element in istate that most wants to be active
 
      do {
        /* Steps 3, 4. */ if(verbose>1) {printf("  steps 3 and 4\n"); fflush(stdout);}
        /* Find the bound element that most wants to be active. */
        double worst=0.0;
        it=1;
        for(int ni=0; ni<nbound; ni++) {
          int i=abs(istate[ni])-1;
          double bad; if(istate[ni] < 0) bad=-w[i]; else bad=+w[i];
          if(bad < worst) { it=ni+1; worst=bad; iact=i; }
        }
 
        /* Test whether the Kuhn-Tucker condition is met. */
        if(worst>=0.0) {
          if(verbose>1) {printf("Kuhn-Tucker condition is met.\n"); fflush(stdout);}
          istate[n]=nbound;
          w[0]=sqrt(obj);
          return(0);
        }
 
        /* The component x[iact] is the one that most wants to become active.
           If the last successful change in the active set was to move x[iact] to a bound, 
           don't let x[iact] in now: set the derivative of the misfit with respect to x[iact] 
           to zero and return to the Kuhn-Tucker test. */
        if(iact==(aindx-1)) w[aindx-1]=0.0;
      } while(iact==(aindx-1)); // Step 3 again
 
      /* Step 5. */ if(verbose>1) {printf("  step 5\n"); fflush(stdout);}
 
      /* Undo the effect of the new (potentially) active variable on the residual vector. */
      if(istate[it-1]==0) { // remove this if never happening
        if(verbose>0) fprintf(stderr, "Error: BVLS istate is zero!\n"); 
        return(1);
      }
      {
        double bnd;
        if(istate[it-1]>0) bnd=bu[iact]; else bnd=bl[iact];
        for(int mi=0; mi<m; mi++) act[mi+mm*m]+=bnd*a[mi+iact*m];
      }
 
      /* Set flag istateFromStep5, indicating that Step 6 was entered from Step 5.
         This forms the basis of a test for instability: the gradient calculation shows that x[iact] 
         wants to join the active set; if QR puts x[iact] beyond the bound from which it came, 
         the gradient calculation was in error and the variable should not have been introduced. */
      istateFromStep5=istate[it-1]; 
 
      /* Swap the indices (in istate) of the new active variable and the rightmost bound variable; 
         `unbind' that location by decrementing nbound. */
      istate[it-1]=istate[nbound-1];
      nbound--;  nact++;
      istate[nbound]=1+iact;
      if(mm<nact) {
        if(verbose>0) fprintf(stderr, "Error: too many free variables in BVLS.\n");
        return(2);
      }
 
    } // finalized steps 3-5
 
    do {
 
      skipStep2=0;
 
      /* Step 6. */ if(verbose>1) {printf("  step 6\n"); fflush(stdout);}
 
      /* Load array act with the appropriate columns of A for QR. For added stability, reverse 
         the column ordering so that the most recent addition to the active set is in the last 
         column. Also copy the residual vector from act[., mm] into act[., mm+1]. */
      for(int mi=0; mi<m; mi++) {
        act[mi+(mm+1)*m]=act[mi+mm*m]; // vector b for QR
        for(int ni=nbound; ni<n; ni++) {
          int i=abs(istate[ni])-1;
          act[mi+(nact+nbound-ni-1)*m]=a[mi+i*m];
        }
      }
      if(verbose>9) {
        printf("Matrix A for QR:\n");
        for(int ni=0; ni<nact; ni++) {
          for(int mi=0; mi<m; mi++) printf("\t%g", act[mi+ni*nact]);
          printf("\n");
        }
        printf("Vector B for QR:\n");
        for(int mi=0; mi<m; mi++) printf("\t%g", act[(mm+1)*m + mi]);
        printf("\n");
      }
 
      /* Test for linear dependence in QR, and for an instability that moves the variable 
         just introduced away from the feasible region (rather than into the region or 
         all the way through it). In either case, remove the latest vector introduced from 
         the active set and adjust the residual vector accordingly.
         Set the gradient component (w[iact]) to zero and return to the Kuhn-Tucker test. */
      double r2;
      if(qrLH(m, nact, act, &act[(mm+1)*m], zz, &r2) !=0 || 
         (istateFromStep5>0 && zz[nact-1]>bu[iact]) || 
         (istateFromStep5<0 && zz[nact-1]<bl[iact]) )
      {
        nbound++;
        if(bu[iact]>x[iact]) istate[nbound-1]=-istate[nbound-1];
        nact--;
        for(int mi=0; mi<m; mi++) act[mi+mm*m]-=x[iact]*a[mi+iact*m];
        istateFromStep5 = 0; // not from step 5
        w[iact]=0.0;
        skipStep2=1; // we want to skip Step 2 and go directly to Step 3
        if(verbose>3) {printf("    going from step 6 to step 3\n"); fflush(stdout);}
        break; // go to step 3
      }
 
      /* If Step 6 was entered from Step 5 and we are here, a new variable has been successfully 
         introduced into the active set; the last variable that was fixed at a bound is again 
         permitted to become active. */
      if(istateFromStep5!=0) aindx=0;
      istateFromStep5=0;
 
      /* Step 7. */ if(verbose>1) {printf("  step 7\n"); fflush(stdout);}
      /* Check for strict feasibility of the new QR solution. */
      int qr_solution_feasible=1;
      int indexHolder=0;
      if(verbose>8) printf("    nact=%d  nbound=%d\n", nact, nbound);
      for(int ni=0; ni<nact; ni++) {
        indexHolder=ni; // Loop in step 8 will start from this
        int i=abs(istate[ni+nbound])-1;
        if(verbose>8) {
          printf("      istate[%d]=%d\n", ni+nbound, 1+i);
          printf("      zz[%d]=%g  bl[%d]=%g  bu[%d]=%g\n", nact-ni-1, zz[nact-ni-1], i, bl[i], i, bu[i]);
        }
        if(zz[nact-ni-1]<bl[i] || zz[nact-ni-1]>bu[i]) {
          if(verbose>3) {printf("    new iterate is not feasible\n"); fflush(stdout);}
          qr_solution_feasible=0; break; // go to Step 8
        }
      }
      if(verbose>8) printf("    indexHolder=%d\n", indexHolder);
      if(qr_solution_feasible) {
        if(verbose>3) {printf("    new iterate is feasible\n"); fflush(stdout);}
        for(int ni=0; ni<nact; ni++) {
          int i=abs(istate[ni+nbound])-1;
          x[i]=zz[nact-ni-1];
        }
        /* New iterate is feasible; back to the top. */
        break; // Back to the start of the main loop
      }
 
      { // keep local variables alpha and alf local
        double alpha=2.0;
        /* Steps 8 and 9 */ if(verbose>1) {printf("  steps 8 and 9\n"); fflush(stdout);}
        double alf=alpha;
        for(int ni=indexHolder; ni<nact; ni++) {
          int i=abs(istate[ni+nbound])-1;
          if(zz[nact-ni-1] > bu[i]) alf=(bu[i]-x[i])/(zz[nact-ni-1]-x[i]);
          if(zz[nact-ni-1] < bl[i]) alf=(bl[i]-x[i])/(zz[nact-ni-1]-x[i]);
          if(alf<alpha) {
            alpha=alf;
            aindx=1+i;
            if((zz[nact-ni-1]-bl[i])<0.0) aindxsign=-1; else aindxsign=+1;
          }
        }
        /* Step 10 */ if(verbose>1) {printf("  step 10\n"); fflush(stdout);}
        for(int ni=0; ni<nact; ni++) {
          int i=abs(istate[ni+nbound])-1;
          x[i]+=alpha*(zz[nact-ni-1]-x[i]);
        }
      }
 
      /* Step 11 */ if(verbose>1) {printf("  step 11\n"); fflush(stdout);}
      /* Move the variable that determined alpha to the appropriate bound
         (aindx is its index; sj is + if zz[aindx]> bu[aindx], - if zz[aindx]<bl[aindx] ).
         If any other component of  x  is infeasible at this stage, it must be due to round-off.
         Bind every infeasible component and every component at a bound to the appropriate bound.
         Correct the residual vector for any variables moved to bounds. Since at least one 
         variable is removed from the active set in this step, Loop B
         (Steps 6-11) terminates after at most nact steps. */
      {
        int noldb=nbound;
        for(int ni=0; ni<nact; ni++) {
          int i=abs(istate[ni+noldb])-1;
          if((bu[i]-x[i]<=0.0) || (i==(aindx-1) && aindxsign>0)) {
            /* Move x[i] to its upper bound. */
            x[i]=bu[i];
            istate[ni+noldb]=istate[nbound]; istate[nbound]=+(1+i); nbound++;
            for(int mi=0; mi<m; mi++) act[mi+mm*m]-=bu[i]*a[mi+i*m];
          } else if( ((x[i]-bl[i])<=0.0) || (i==(aindx-1) && aindxsign<0)) {
            /* Move x(j) to its lower bound. */
            x[i]=bl[i];
            istate[ni+noldb]=istate[nbound]; istate[nbound]=-(1+i); nbound++;
            for(int mi=0; mi<m; mi++) act[mi+mm*m]-=bl[i]*a[mi+i*m];
          }
        }
        nact=n-nbound;
      }
      /* If there are still active variables left, repeat the QR; if not, go back to step 6. */
    } while(nact>0);
 
  } // main loop
 
  /* iterMax reached */
  if(verbose>0) fprintf(stderr, "Error: BVLS fails to converge.\n");
  return(-1);
}
/*****************************************************************************/
 
/*****************************************************************************/
int llsqWght(
  int N,
  int M,
  double **A,
  double *a,
  double *b,
  double *weight
) {
  int n, m;
  double *w;
 
  /* Check the arguments */
  if(N<1 || M<1 || (A==NULL && a==NULL) || b==NULL || weight==NULL) return(1);
 
  /* Allocate memory */
  w=(double*)malloc(M*sizeof(double)); if(w==NULL) return(2);
 
  /* Check that weights are not zero and get the square roots of them to w[] */
  for(m=0; m<M; m++) {
    if(weight[m]<=1.0e-20) w[m]=0.0;
    else w[m]=sqrt(weight[m]);
  }
 
  /* Multiply rows of matrix A and elements of vector b with weights*/
  for(m=0; m<M; m++) {
    for(n=0; n<N; n++) {
      if(A!=NULL) A[n][m]*=w[m];
      if(a!=NULL) a[m+n*M]*=w[m];
    }
    b[m]*=w[m];
  }
 
  free(w);
  return(0);
}
/*****************************************************************************/
 
/*****************************************************************************/
int llsqWghtSquared(
  int N,
  int M,
  double **A,
  double *a,
  double *b,
  double *sweight
) {
  int n, m;
 
  /* Check the arguments */
  if(N<1 || M<1 || (A==NULL && a==NULL) || b==NULL || sweight==NULL) return(1);
 
  /* Multiply rows of matrix A and elements of vector b with weights*/
  for(m=0; m<M; m++) {
    for(n=0; n<N; n++) {
      if(A!=NULL) A[n][m]*=sweight[m];
      if(a!=NULL) a[m+n*M]*=sweight[m];
    }
    b[m]*=sweight[m];
  }
 
  return(0);
}
/*****************************************************************************/
 
/*****************************************************************************/




