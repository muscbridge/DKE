/*  -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__8 = 8;
static doublereal c_b41 = 1.;
static integer c__0 = 0;
static doublereal c_b95 = 0.;
static integer c__3 = 3;
static logical c_false = FALSE_;
static integer c_n1 = -1;
static integer c__72 = 72;
static logical c_true = TRUE_;
static integer c__5 = 5;

/* DECK DLSEI */
/* Subroutine */ int dlsei_(doublereal *w, integer *mdw, integer *me, integer 
	*ma, integer *mg, integer *n, doublereal *prgopt, doublereal *x, 
	doublereal *rnorme, doublereal *rnorml, integer *mode, doublereal *ws,
	 integer *ip)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    address a__1[8], a__2[2];
    integer w_dim1, w_offset, i__1, i__2[8], i__3[2], i__4, i__5;
    doublereal d__1, d__2;
    char ch__1[131], ch__2[60], ch__3[61];

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal t;
    static integer n1, n2;
    static doublereal rb, uj, rn, sn, vj, up;
    static integer jp1, np1;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal gam;
    static integer key;
    static doublereal tau;
    static logical cov;
    static integer mep1, lchk, mend;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dlsi_(doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer link, imax, last;
    static doublereal size;
    static integer next, nopt;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static char xern1[8], xern2[8], xern3[8], xern4[8];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer mdeqc;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer nlink;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal enorm, fnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal rnmax, snmax, xnrme;
    extern doublereal d1mach_(integer *);
    static doublereal xnorm;
    static integer mapke1, kranke;
    static doublereal drelpr;
    static integer ntimes;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___7 = { 0, xern2, 0, "(I8)", 8, 1 };
    static icilist io___9 = { 0, xern3, 0, "(I8)", 8, 1 };
    static icilist io___11 = { 0, xern4, 0, "(I8)", 8, 1 };
    static icilist io___13 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___14 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DLSEI */
/* ***PURPOSE  Solve a linearly constrained least squares problem with */
/*            equality and inequality constraints, and optionally compute */
/*            a covariance matrix. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  K1A2A, D9 */
/* ***TYPE      DOUBLE PRECISION (LSEI-S, DLSEI-D) */
/* ***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING, */
/*             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS, */
/*             QUADRATIC PROGRAMMING */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */

/*     This subprogram solves a linearly constrained least squares */
/*     problem with both equality and inequality constraints, and, if the */
/*     user requests, obtains a covariance matrix of the solution */
/*     parameters. */

/*     Suppose there are given matrices E, A and G of respective */
/*     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of */
/*     respective lengths ME, MA and MG.  This subroutine solves the */
/*     linearly constrained least squares problem */

/*                   EX = F, (E ME by N) (equations to be exactly */
/*                                       satisfied) */
/*                   AX = B, (A MA by N) (equations to be */
/*                                       approximately satisfied, */
/*                                       least squares sense) */
/*                   GX .GE. H,(G MG by N) (inequality constraints) */

/*     The inequalities GX .GE. H mean that every component of the */
/*     product GX must be .GE. the corresponding component of H. */

/*     In case the equality constraints cannot be satisfied, a */
/*     generalized inverse solution residual vector length is obtained */
/*     for F-EX.  This is the minimal length possible for F-EX. */

/*     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The */
/*     rank of the matrix E is estimated during the computation.  We call */
/*     this value KRANKE.  It is an output parameter in IP(1) defined */
/*     below.  Using a generalized inverse solution of EX=F, a reduced */
/*     least squares problem with inequality constraints is obtained. */
/*     The tolerances used in these tests for determining the rank */
/*     of E and the rank of the reduced least squares problem are */
/*     given in Sandia Tech. Rept. SAND-78-1290.  They can be */
/*     modified by the user if new values are provided in */
/*     the option list of the array PRGOPT(*). */

/*     The user must dimension all arrays appearing in the call list.. */
/*     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2) */
/*     where K=MAX(MA+MG,N).  This allows for a solution of a range of */
/*     problems in the given working space.  The dimension of WS(*) */
/*     given is a necessary overestimate.  Once a particular problem */
/*     has been run, the output parameter IP(3) gives the actual */
/*     dimension required for that problem. */

/*     The parameters for DLSEI( ) are */

/*     Input.. All TYPE REAL variables are DOUBLE PRECISION */

/*     W(*,*),MDW,   The array W(*,*) is doubly subscripted with */
/*     ME,MA,MG,N    first dimensioning parameter equal to MDW. */
/*                   For this discussion let us call M = ME+MA+MG.  Then */
/*                   MDW must satisfy MDW .GE. M.  The condition */
/*                   MDW .LT. M is an error. */

/*                   The array W(*,*) contains the matrices and vectors */

/*                                  (E  F) */
/*                                  (A  B) */
/*                                  (G  H) */

/*                   in rows and columns 1,...,M and 1,...,N+1 */
/*                   respectively. */

/*                   The integers ME, MA, and MG are the */
/*                   respective matrix row dimensions */
/*                   of E, A and G.  Each matrix has N columns. */

/*     PRGOPT(*)    This real-valued array is the option vector. */
/*                  If the user is satisfied with the nominal */
/*                  subprogram features set */

/*                  PRGOPT(1)=1 (or PRGOPT(1)=1.0) */

/*                  Otherwise PRGOPT(*) is a linked list consisting of */
/*                  groups of data of the following form */

/*                  LINK */
/*                  KEY */
/*                  DATA SET */

/*                  The parameters LINK and KEY are each one word. */
/*                  The DATA SET can be comprised of several words. */
/*                  The number of items depends on the value of KEY. */
/*                  The value of LINK points to the first */
/*                  entry of the next group of data within */
/*                  PRGOPT(*).  The exception is when there are */
/*                  no more options to change.  In that */
/*                  case, LINK=1 and the values KEY and DATA SET */
/*                  are not referenced.  The general layout of */
/*                  PRGOPT(*) is as follows. */

/*               ...PRGOPT(1) = LINK1 (link to first entry of next group) */
/*               .  PRGOPT(2) = KEY1 (key to the option change) */
/*               .  PRGOPT(3) = data value (data value for this change) */
/*               .       . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of */
/*               .                       next group) */
/*               .  PRGOPT(LINK1+1) = KEY2 (key to the option change) */
/*               .  PRGOPT(LINK1+2) = data value */
/*               ...     . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK) = 1 (no more options to change) */

/*                  Values of LINK that are nonpositive are errors. */
/*                  A value of LINK .GT. NLINK=100000 is also an error. */
/*                  This helps prevent using invalid but positive */
/*                  values of LINK that will probably extend */
/*                  beyond the program limits of PRGOPT(*). */
/*                  Unrecognized values of KEY are ignored.  The */
/*                  order of the options is arbitrary and any number */
/*                  of options can be changed with the following */
/*                  restriction.  To prevent cycling in the */
/*                  processing of the option array, a count of the */
/*                  number of options changed is maintained. */
/*                  Whenever this count exceeds NOPT=1000, an error */
/*                  message is printed and the subprogram returns. */

/*                  Options.. */

/*                  KEY=1 */
/*                         Compute in W(*,*) the N by N */
/*                  covariance matrix of the solution variables */
/*                  as an output parameter.  Nominally the */
/*                  covariance matrix will not be computed. */
/*                  (This requires no user input.) */
/*                  The data set for this option is a single value. */
/*                  It must be nonzero when the covariance matrix */
/*                  is desired.  If it is zero, the covariance */
/*                  matrix is not computed.  When the covariance matrix */
/*                  is computed, the first dimensioning parameter */
/*                  of the array W(*,*) must satisfy MDW .GE. MAX(M,N). */

/*                  KEY=10 */
/*                         Suppress scaling of the inverse of the */
/*                  normal matrix by the scale factor RNORM**2/ */
/*                  MAX(1, no. of degrees of freedom).  This option */
/*                  only applies when the option for computing the */
/*                  covariance matrix (KEY=1) is used.  With KEY=1 and */
/*                  KEY=10 used as options the unscaled inverse of the */
/*                  normal matrix is returned in W(*,*). */
/*                  The data set for this option is a single value. */
/*                  When it is nonzero no scaling is done.  When it is */
/*                  zero scaling is done.  The nominal case is to do */
/*                  scaling so if option (KEY=1) is used alone, the */
/*                  matrix will be scaled on output. */

/*                  KEY=2 */
/*                         Scale the nonzero columns of the */
/*                         entire data matrix. */
/*                  (E) */
/*                  (A) */
/*                  (G) */

/*                  to have length one.  The data set for this */
/*                  option is a single value.  It must be */
/*                  nonzero if unit length column scaling */
/*                  is desired. */

/*                  KEY=3 */
/*                         Scale columns of the entire data matrix */
/*                  (E) */
/*                  (A) */
/*                  (G) */

/*                  with a user-provided diagonal matrix. */
/*                  The data set for this option consists */
/*                  of the N diagonal scaling factors, one for */
/*                  each matrix column. */

/*                  KEY=4 */
/*                         Change the rank determination tolerance for */
/*                  the equality constraint equations from */
/*                  the nominal value of SQRT(DRELPR).  This quantity can */
/*                  be no smaller than DRELPR, the arithmetic- */
/*                  storage precision.  The quantity DRELPR is the */
/*                  largest positive number such that T=1.+DRELPR */
/*                  satisfies T .EQ. 1.  The quantity used */
/*                  here is internally restricted to be at */
/*                  least DRELPR.  The data set for this option */
/*                  is the new tolerance. */

/*                  KEY=5 */
/*                         Change the rank determination tolerance for */
/*                  the reduced least squares equations from */
/*                  the nominal value of SQRT(DRELPR).  This quantity can */
/*                  be no smaller than DRELPR, the arithmetic- */
/*                  storage precision.  The quantity used */
/*                  here is internally restricted to be at */
/*                  least DRELPR.  The data set for this option */
/*                  is the new tolerance. */

/*                  For example, suppose we want to change */
/*                  the tolerance for the reduced least squares */
/*                  problem, compute the covariance matrix of */
/*                  the solution parameters, and provide */
/*                  column scaling for the data matrix.  For */
/*                  these options the dimension of PRGOPT(*) */
/*                  must be at least N+9.  The Fortran statements */
/*                  defining these options would be as follows: */

/*                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*)) */
/*                  PRGOPT(2)=1 (covariance matrix key) */
/*                  PRGOPT(3)=1 (covariance matrix wanted) */

/*                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*)) */
/*                  PRGOPT(5)=5 (least squares equas.  tolerance key) */
/*                  PRGOPT(6)=... (new value of the tolerance) */

/*                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*)) */
/*                  PRGOPT(8)=3 (user-provided column scaling key) */

/*                  CALL DCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N */
/*                    scaling factors from the user array D(*) */
/*                    to PRGOPT(9)-PRGOPT(N+8)) */

/*                  PRGOPT(N+9)=1 (no more options to change) */

/*                  The contents of PRGOPT(*) are not modified */
/*                  by the subprogram. */
/*                  The options for WNNLS( ) can also be included */
/*                  in this array.  The values of KEY recognized */
/*                  by WNNLS( ) are 6, 7 and 8.  Their functions */
/*                  are documented in the usage instructions for */
/*                  subroutine WNNLS( ).  Normally these options */
/*                  do not need to be modified when using DLSEI( ). */

/*     IP(1),       The amounts of working storage actually */
/*     IP(2)        allocated for the working arrays WS(*) and */
/*                  IP(*), respectively.  These quantities are */
/*                  compared with the actual amounts of storage */
/*                  needed by DLSEI( ).  Insufficient storage */
/*                  allocated for either WS(*) or IP(*) is an */
/*                  error.  This feature was included in DLSEI( ) */
/*                  because miscalculating the storage formulas */
/*                  for WS(*) and IP(*) might very well lead to */
/*                  subtle and hard-to-find execution errors. */

/*                  The length of WS(*) must be at least */

/*                  LW = 2*(ME+N)+K+(MG+2)*(N+7) */

/*                  where K = max(MA+MG,N) */
/*                  This test will not be made if IP(1).LE.0. */

/*                  The length of IP(*) must be at least */

/*                  LIP = MG+2*N+2 */
/*                  This test will not be made if IP(2).LE.0. */

/*     Output.. All TYPE REAL variables are DOUBLE PRECISION */

/*     X(*),RNORME,  The array X(*) contains the solution parameters */
/*     RNORML        if the integer output flag MODE = 0 or 1. */
/*                   The definition of MODE is given directly below. */
/*                   When MODE = 0 or 1, RNORME and RNORML */
/*                   respectively contain the residual vector */
/*                   Euclidean lengths of F - EX and B - AX.  When */
/*                   MODE=1 the equality constraint equations EX=F */
/*                   are contradictory, so RNORME .NE. 0.  The residual */
/*                   vector F-EX has minimal Euclidean length.  For */
/*                   MODE .GE. 2, none of these parameters is defined. */

/*     MODE          Integer flag that indicates the subprogram */
/*                   status after completion.  If MODE .GE. 2, no */
/*                   solution has been computed. */

/*                   MODE = */

/*                   0  Both equality and inequality constraints */
/*                      are compatible and have been satisfied. */

/*                   1  Equality constraints are contradictory. */
/*                      A generalized inverse solution of EX=F was used */
/*                      to minimize the residual vector length F-EX. */
/*                      In this sense, the solution is still meaningful. */

/*                   2  Inequality constraints are contradictory. */

/*                   3  Both equality and inequality constraints */
/*                      are contradictory. */

/*                   The following interpretation of */
/*                   MODE=1,2 or 3 must be made.  The */
/*                   sets consisting of all solutions */
/*                   of the equality constraints EX=F */
/*                   and all vectors satisfying GX .GE. H */
/*                   have no points in common.  (In */
/*                   particular this does not say that */
/*                   each individual set has no points */
/*                   at all, although this could be the */
/*                   case.) */

/*                   4  Usage error occurred.  The value */
/*                      of MDW is .LT. ME+MA+MG, MDW is */
/*                      .LT. N and a covariance matrix is */
/*                      requested, or the option vector */
/*                      PRGOPT(*) is not properly defined, */
/*                      or the lengths of the working arrays */
/*                      WS(*) and IP(*), when specified in */
/*                      IP(1) and IP(2) respectively, are not */
/*                      long enough. */

/*     W(*,*)        The array W(*,*) contains the N by N symmetric */
/*                   covariance matrix of the solution parameters, */
/*                   provided this was requested on input with */
/*                   the option vector PRGOPT(*) and the output */
/*                   flag is returned with MODE = 0 or 1. */

/*     IP(*)         The integer working array has three entries */
/*                   that provide rank and working array length */
/*                   information after completion. */

/*                      IP(1) = rank of equality constraint */
/*                              matrix.  Define this quantity */
/*                              as KRANKE. */

/*                      IP(2) = rank of reduced least squares */
/*                              problem. */

/*                      IP(3) = the amount of storage in the */
/*                              working array WS(*) that was */
/*                              actually used by the subprogram. */
/*                              The formula given above for the length */
/*                              of WS(*) is a necessary overestimate. */
/*                              If exactly the same problem matrices */
/*                              are used in subsequent executions, */
/*                              the declared dimension of WS(*) can */
/*                              be reduced to this output value. */
/*     User Designated */
/*     Working Arrays.. */

/*     WS(*),IP(*)              These are respectively type real */
/*                              and type integer working arrays. */
/*                              Their required minimal lengths are */
/*                              given above. */

/* ***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for */
/*                 linear least squares problems with equality and */
/*                 nonnegativity constraints, Report SAND77-0552, Sandia */
/*                 Laboratories, June 1978. */
/*               K. H. Haskell and R. J. Hanson, Selected algorithms for */
/*                 the linearly constrained least squares problem - a */
/*                 users guide, Report SAND78-1290, Sandia Laboratories, */
/*                 August 1979. */
/*               K. H. Haskell and R. J. Hanson, An algorithm for */
/*                 linear least squares problems with equality and */
/*                 nonnegativity constraints, Mathematical Programming */
/*                 21 (1981), pp. 98-118. */
/*               R. J. Hanson and K. H. Haskell, Two algorithms for the */
/*                 linearly constrained least squares problem, ACM */
/*                 Transactions on Mathematical Software, September 1982. */
/* ***ROUTINES CALLED  D1MACH, DASUM, DAXPY, DCOPY, DDOT, DH12, DLSI, */
/*                    DNRM2, DSCAL, DSWAP, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and extensively revised (WRB & RWC) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   900604  DP version created from SP version.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DLSEI */



    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --prgopt;
    --x;
    --ws;
    --ip;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DLSEI */

/*     Set the nominal tolerance used in the code for the equality */
/*     constraint equations. */

    if (first) {
	drelpr = d1mach_(&c__4);
    }
    first = FALSE_;
    tau = sqrt(drelpr);

/*     Check that enough storage was allocated in WS(*) and IP(*). */

    *mode = 4;
/* Computing MIN */
    i__1 = min(*n,*me), i__1 = min(i__1,*ma);
    if (min(i__1,*mg) < 0) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&(*me), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___9);
	do_fio(&c__1, (char *)&(*ma), (ftnlen)sizeof(integer));
	e_wsfi();
	s_wsfi(&io___11);
	do_fio(&c__1, (char *)&(*mg), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 78, a__1[0] = "ALL OF THE VARIABLES N, ME, MA, MG MUST BE "
		".GE. 0$$ENTERED ROUTINE WITH$$N  = ";
	i__2[1] = 8, a__1[1] = xern1;
	i__2[2] = 7, a__1[2] = "$$ME = ";
	i__2[3] = 8, a__1[3] = xern2;
	i__2[4] = 7, a__1[4] = "$$MA = ";
	i__2[5] = 8, a__1[5] = xern3;
	i__2[6] = 7, a__1[6] = "$$MG = ";
	i__2[7] = 8, a__1[7] = xern4;
	s_cat(ch__1, a__1, i__2, &c__8, (ftnlen)131);
	xermsg_("SLATEC", "LSEI", ch__1, &c__2, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)131);
	return 0;
    }

    if (ip[1] > 0) {
/* Computing MAX */
	i__1 = *ma + *mg;
	lchk = (*me + *n << 1) + max(i__1,*n) + (*mg + 2) * (*n + 7);
	if (ip[1] < lchk) {
	    s_wsfi(&io___13);
	    do_fio(&c__1, (char *)&lchk, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 52, a__2[0] = "INSUFFICIENT STORAGE ALLOCATED FOR WS(*"
		    "), NEED LW = ";
	    i__3[1] = 8, a__2[1] = xern1;
	    s_cat(ch__2, a__2, i__3, &c__2, (ftnlen)60);
	    xermsg_("SLATEC", "DLSEI", ch__2, &c__2, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)60);
	    return 0;
	}
    }

    if (ip[2] > 0) {
	lchk = *mg + (*n << 1) + 2;
	if (ip[2] < lchk) {
	    s_wsfi(&io___14);
	    do_fio(&c__1, (char *)&lchk, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__3[0] = 53, a__2[0] = "INSUFFICIENT STORAGE ALLOCATED FOR IP(*"
		    "), NEED LIP = ";
	    i__3[1] = 8, a__2[1] = xern1;
	    s_cat(ch__3, a__2, i__3, &c__2, (ftnlen)61);
	    xermsg_("SLATEC", "DLSEI", ch__3, &c__2, &c__1, (ftnlen)6, (
		    ftnlen)5, (ftnlen)61);
	    return 0;
	}
    }

/*     Compute number of possible right multiplying Householder */
/*     transformations. */

    m = *me + *ma + *mg;
    if (*n <= 0 || m <= 0) {
	*mode = 0;
	*rnorme = 0.;
	*rnorml = 0.;
	return 0;
    }

    if (*mdw < m) {
	xermsg_("SLATEC", "DLSEI", "MDW.LT.ME+MA+MG IS AN ERROR", &c__2, &
		c__1, (ftnlen)6, (ftnlen)5, (ftnlen)27);
	return 0;
    }

    np1 = *n + 1;
    kranke = min(*me,*n);
    n1 = (kranke << 1) + 1;
    n2 = n1 + *n;

/*     Set nominal values. */

/*     The nominal column scaling used in the code is */
/*     the identity scaling. */

    dcopy_(n, &c_b41, &c__0, &ws[n1], &c__1);

/*     No covariance matrix is nominally computed. */

    cov = FALSE_;

/*     Process option vector. */
/*     Define bound for number of options to change. */

    nopt = 1000;
    ntimes = 0;

/*     Define bound for positive values of LINK. */

    nlink = 100000;
    last = 1;
    link = (integer) prgopt[1];
    if (link == 0 || link > nlink) {
	xermsg_("SLATEC", "DLSEI", "THE OPTION VECTOR IS UNDEFINED", &c__2, &
		c__1, (ftnlen)6, (ftnlen)5, (ftnlen)30);
	return 0;
    }

L100:
    if (link > 1) {
	++ntimes;
	if (ntimes > nopt) {
	    xermsg_("SLATEC", "DLSEI", "THE LINKS IN THE OPTION VECTOR ARE C"
		    "YCLING.", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)43);
	    return 0;
	}

	key = (integer) prgopt[last + 1];
	if (key == 1) {
	    cov = prgopt[last + 2] != 0.;
	} else if (key == 2 && prgopt[last + 2] != 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		t = dnrm2_(&m, &w[j * w_dim1 + 1], &c__1);
		if (t != 0.) {
		    t = 1. / t;
		}
		ws[j + n1 - 1] = t;
/* L110: */
	    }
	} else if (key == 3) {
	    dcopy_(n, &prgopt[last + 2], &c__1, &ws[n1], &c__1);
	} else if (key == 4) {
/* Computing MAX */
	    d__1 = drelpr, d__2 = prgopt[last + 2];
	    tau = max(d__1,d__2);
	}

	next = (integer) prgopt[link];
	if (next <= 0 || next > nlink) {
	    xermsg_("SLATEC", "DLSEI", "THE OPTION VECTOR IS UNDEFINED", &
		    c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)30);
	    return 0;
	}

	last = link;
	link = next;
	goto L100;
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dscal_(&m, &ws[n1 + j - 1], &w[j * w_dim1 + 1], &c__1);
/* L120: */
    }

    if (cov && *mdw < *n) {
	xermsg_("SLATEC", "DLSEI", "MDW .LT. N WHEN COV MATRIX NEEDED, IS AN"
		" ERROR", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)46);
	return 0;
    }

/*     Problem definition and option vector OK. */

    *mode = 0;

/*     Compute norm of equality constraint matrix and right side. */

    enorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = enorm, d__2 = dasum_(me, &w[j * w_dim1 + 1], &c__1);
	enorm = max(d__1,d__2);
/* L130: */
    }

    fnorm = dasum_(me, &w[np1 * w_dim1 + 1], &c__1);
    snmax = 0.;
    rnmax = 0.;
    i__1 = kranke;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Compute maximum ratio of vector lengths. Partition is at */
/*        column I. */

	i__4 = *me;
	for (k = i__; k <= i__4; ++k) {
	    i__5 = *n - i__ + 1;
	    sn = ddot_(&i__5, &w[k + i__ * w_dim1], mdw, &w[k + i__ * w_dim1],
		     mdw);
	    i__5 = i__ - 1;
	    rn = ddot_(&i__5, &w[k + w_dim1], mdw, &w[k + w_dim1], mdw);
	    if (rn == 0. && sn > snmax) {
		snmax = sn;
		imax = k;
	    } else if (k == i__ || sn * rnmax > rn * snmax) {
		snmax = sn;
		rnmax = rn;
		imax = k;
	    }
/* L140: */
	}

/*        Interchange rows if necessary. */

	if (i__ != imax) {
	    dswap_(&np1, &w[i__ + w_dim1], mdw, &w[imax + w_dim1], mdw);
	}
/* Computing 2nd power */
	d__1 = tau;
	if (snmax > rnmax * (d__1 * d__1)) {

/*        Eliminate elements I+1,...,N in row I. */

	    i__4 = i__ + 1;
	    i__5 = m - i__;
	    dh12_(&c__1, &i__, &i__4, n, &w[i__ + w_dim1], mdw, &ws[i__], &w[
		    i__ + 1 + w_dim1], mdw, &c__1, &i__5);
	} else {
	    kranke = i__ - 1;
	    goto L160;
	}
/* L150: */
    }

/*     Save diagonal terms of lower trapezoidal matrix. */

L160:
    i__1 = *mdw + 1;
    dcopy_(&kranke, &w[w_offset], &i__1, &ws[kranke + 1], &c__1);

/*     Use Householder transformation from left to achieve */
/*     KRANKE by KRANKE upper triangular form. */

    if (kranke < *me) {
	for (k = kranke; k >= 1; --k) {

/*           Apply transformation to matrix cols. 1,...,K-1. */

	    i__1 = kranke + 1;
	    i__4 = k - 1;
	    dh12_(&c__1, &k, &i__1, me, &w[k * w_dim1 + 1], &c__1, &up, &w[
		    w_offset], &c__1, mdw, &i__4);

/*           Apply to rt side vector. */

	    i__1 = kranke + 1;
	    dh12_(&c__2, &k, &i__1, me, &w[k * w_dim1 + 1], &c__1, &up, &w[
		    np1 * w_dim1 + 1], &c__1, &c__1, &c__1);
/* L170: */
	}
    }

/*     Solve for variables 1,...,KRANKE in new coordinates. */

    dcopy_(&kranke, &w[np1 * w_dim1 + 1], &c__1, &x[1], &c__1);
    i__1 = kranke;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__4 = i__ - 1;
	x[i__] = (x[i__] - ddot_(&i__4, &w[i__ + w_dim1], mdw, &x[1], &c__1)) 
		/ w[i__ + i__ * w_dim1];
/* L180: */
    }

/*     Compute residuals for reduced problem. */

    mep1 = *me + 1;
    *rnorml = 0.;
    i__1 = m;
    for (i__ = mep1; i__ <= i__1; ++i__) {
	w[i__ + np1 * w_dim1] -= ddot_(&kranke, &w[i__ + w_dim1], mdw, &x[1], 
		&c__1);
	sn = ddot_(&kranke, &w[i__ + w_dim1], mdw, &w[i__ + w_dim1], mdw);
	i__4 = *n - kranke;
	rn = ddot_(&i__4, &w[i__ + (kranke + 1) * w_dim1], mdw, &w[i__ + (
		kranke + 1) * w_dim1], mdw);
/* Computing 2nd power */
	d__1 = tau;
	if (rn <= sn * (d__1 * d__1) && kranke < *n) {
	    i__4 = *n - kranke;
	    dcopy_(&i__4, &c_b95, &c__0, &w[i__ + (kranke + 1) * w_dim1], mdw)
		    ;
	}
/* L190: */
    }

/*     Compute equality constraint equations residual length. */

    i__1 = *me - kranke;
    *rnorme = dnrm2_(&i__1, &w[kranke + 1 + np1 * w_dim1], &c__1);

/*     Move reduced problem data upward if KRANKE.LT.ME. */

    if (kranke < *me) {
	i__1 = np1;
	for (j = 1; j <= i__1; ++j) {
	    i__4 = m - *me;
	    dcopy_(&i__4, &w[*me + 1 + j * w_dim1], &c__1, &w[kranke + 1 + j *
		     w_dim1], &c__1);
/* L200: */
	}
    }

/*     Compute solution of reduced problem. */

    i__1 = *n - kranke;
    dlsi_(&w[kranke + 1 + (kranke + 1) * w_dim1], mdw, ma, mg, &i__1, &prgopt[
	    1], &x[kranke + 1], rnorml, mode, &ws[n2], &ip[2]);

/*     Test for consistency of equality constraints. */

    if (*me > 0) {
	mdeqc = 0;
	xnrme = dasum_(&kranke, &w[np1 * w_dim1 + 1], &c__1);
	if (*rnorme > tau * (enorm * xnrme + fnorm)) {
	    mdeqc = 1;
	}
	*mode += mdeqc;

/*        Check if solution to equality constraints satisfies inequality */
/*        constraints when there are no degrees of freedom left. */

	if (kranke == *n && *mg > 0) {
	    xnorm = dasum_(n, &x[1], &c__1);
	    mapke1 = *ma + kranke + 1;
	    mend = *ma + kranke + *mg;
	    i__1 = mend;
	    for (i__ = mapke1; i__ <= i__1; ++i__) {
		size = dasum_(n, &w[i__ + w_dim1], mdw) * xnorm + (d__1 = w[
			i__ + np1 * w_dim1], abs(d__1));
		if (w[i__ + np1 * w_dim1] > tau * size) {
		    *mode += 2;
		    goto L290;
		}
/* L210: */
	    }
	}
    }

/*     Replace diagonal terms of lower trapezoidal matrix. */

    if (kranke > 0) {
	i__1 = *mdw + 1;
	dcopy_(&kranke, &ws[kranke + 1], &c__1, &w[w_offset], &i__1);

/*        Reapply transformation to put solution in original coordinates. */

	for (i__ = kranke; i__ >= 1; --i__) {
	    i__1 = i__ + 1;
	    dh12_(&c__2, &i__, &i__1, n, &w[i__ + w_dim1], mdw, &ws[i__], &x[
		    1], &c__1, &c__1, &c__1);
/* L220: */
	}

/*        Compute covariance matrix of equality constrained problem. */

	if (cov) {
/* Computing MIN */
	    i__1 = kranke, i__4 = *n - 1;
	    for (j = min(i__1,i__4); j >= 1; --j) {
		rb = ws[j] * w[j + j * w_dim1];
		if (rb != 0.) {
		    rb = 1. / rb;
		}
		jp1 = j + 1;
		i__1 = *n;
		for (i__ = jp1; i__ <= i__1; ++i__) {
		    i__4 = *n - j;
		    w[i__ + j * w_dim1] = rb * ddot_(&i__4, &w[i__ + jp1 * 
			    w_dim1], mdw, &w[j + jp1 * w_dim1], mdw);
/* L230: */
		}

		i__1 = *n - j;
		gam = rb * .5 * ddot_(&i__1, &w[jp1 + j * w_dim1], &c__1, &w[
			j + jp1 * w_dim1], mdw);
		i__1 = *n - j;
		daxpy_(&i__1, &gam, &w[j + jp1 * w_dim1], mdw, &w[jp1 + j * 
			w_dim1], &c__1);
		i__1 = *n;
		for (i__ = jp1; i__ <= i__1; ++i__) {
		    i__4 = *n;
		    for (k = i__; k <= i__4; ++k) {
			w[i__ + k * w_dim1] = w[i__ + k * w_dim1] + w[j + i__ 
				* w_dim1] * w[k + j * w_dim1] + w[i__ + j * 
				w_dim1] * w[j + k * w_dim1];
			w[k + i__ * w_dim1] = w[i__ + k * w_dim1];
/* L240: */
		    }
/* L250: */
		}
		uj = ws[j];
		vj = gam * uj;
		w[j + j * w_dim1] = uj * vj + uj * vj;
		i__1 = *n;
		for (i__ = jp1; i__ <= i__1; ++i__) {
		    w[j + i__ * w_dim1] = uj * w[i__ + j * w_dim1] + vj * w[j 
			    + i__ * w_dim1];
/* L260: */
		}
		i__1 = *n - j;
		dcopy_(&i__1, &w[j + jp1 * w_dim1], mdw, &w[jp1 + j * w_dim1],
			 &c__1);
/* L270: */
	    }
	}
    }

/*     Apply the scaling to the covariance matrix. */

    if (cov) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dscal_(n, &ws[i__ + n1 - 1], &w[i__ + w_dim1], mdw);
	    dscal_(n, &ws[i__ + n1 - 1], &w[i__ * w_dim1 + 1], &c__1);
/* L280: */
	}
    }

/*     Rescale solution vector. */

L290:
    if (*mode <= 1) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    x[j] *= ws[n1 + j - 1];
/* L300: */
	}
    }

    ip[1] = kranke;
    ip[3] = ip[3] + (kranke << 1) + *n;
    return 0;
} /* dlsei_ */

/* DECK D1MACH */
doublereal d1mach_(integer *i__)
{
    /* System generated locals */
    doublereal ret_val;
    static doublereal equiv_4[6];

    /* Local variables */
#define log10 ((integer *)equiv_4 + 8)
#define dmach (equiv_4)
#define large ((integer *)equiv_4 + 2)
#define small ((integer *)equiv_4)
#define diver ((integer *)equiv_4 + 6)
#define right ((integer *)equiv_4 + 4)
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D1MACH */
/* ***PURPOSE  Return floating point machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   D1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument, and can be referenced as follows: */

/*        D = D1MACH(I) */

/*   where I=1,...,5.  The (output) value of D above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude. */
/*   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   D1MACH( 3) = B**(-T), the smallest relative spacing. */
/*   D1MACH( 4) = B**(1-T), the largest relative spacing. */
/*   D1MACH( 5) = LOG10(B) */

/*   Assume double precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(14) = T, the number of base-B digits. */
/*   I1MACH(15) = EMIN, the smallest exponent E. */
/*   I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890213  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   900911  Added SUN 386i constants.  (WRB) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added CONVEX -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/* ***END PROLOGUE  D1MACH */




/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 / */
/*     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 / */
/*     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 / */
/*     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA SMALL(1) / ZC00800000 / */
/*     DATA SMALL(2) / Z000000000 / */
/*     DATA LARGE(1) / ZDFFFFFFFF / */
/*     DATA LARGE(2) / ZFFFFFFFFF / */
/*     DATA RIGHT(1) / ZCC5800000 / */
/*     DATA RIGHT(2) / Z000000000 / */
/*     DATA DIVER(1) / ZCC6800000 / */
/*     DATA DIVER(2) / Z000000000 / */
/*     DATA LOG10(1) / ZD00E730E7 / */
/*     DATA LOG10(2) / ZC77800DC0 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O0000000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O0007777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O7770000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O7777777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA SMALL(1) / Z"3001800000000000" / */
/*     DATA SMALL(2) / Z"3001000000000000" / */
/*     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" / */
/*     DATA LARGE(2) / Z"4FFE000000000000" / */
/*     DATA RIGHT(1) / Z"3FD2800000000000" / */
/*     DATA RIGHT(2) / Z"3FD2000000000000" / */
/*     DATA DIVER(1) / Z"3FD3800000000000" / */
/*     DATA DIVER(2) / Z"3FD3000000000000" / */
/*     DATA LOG10(1) / Z"3FFF9A209A84FBCF" / */
/*     DATA LOG10(2) / Z"3FFFF7988F8959AC" / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA SMALL(1) / 00564000000000000000B / */
/*     DATA SMALL(2) / 00000000000000000000B / */
/*     DATA LARGE(1) / 37757777777777777777B / */
/*     DATA LARGE(2) / 37157777777777777777B / */
/*     DATA RIGHT(1) / 15624000000000000000B / */
/*     DATA RIGHT(2) / 00000000000000000000B / */
/*     DATA DIVER(1) / 15634000000000000000B / */
/*     DATA DIVER(2) / 00000000000000000000B / */
/*     DATA LOG10(1) / 17164642023241175717B / */
/*     DATA LOG10(2) / 16367571421742254654B / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn OR -pd8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CC0000000000000' / */
/*     DATA DMACH(4) / Z'3CD0000000000000' / */
/*     DATA DMACH(5) / Z'3FF34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F900000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F910000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE CRAY */

/*     DATA SMALL(1) / 201354000000000000000B / */
/*     DATA SMALL(2) / 000000000000000000000B / */
/*     DATA LARGE(1) / 577767777777777777777B / */
/*     DATA LARGE(2) / 000007777777777777774B / */
/*     DATA RIGHT(1) / 376434000000000000000B / */
/*     DATA RIGHT(2) / 000000000000000000000B / */
/*     DATA DIVER(1) / 376444000000000000000B / */
/*     DATA DIVER(2) / 000000000000000000000B / */
/*     DATA LOG10(1) / 377774642023241175717B / */
/*     DATA LOG10(2) / 000007571421742254654B / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */
/*     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD - */
/*     STATIC DMACH(5) */

/*     DATA SMALL /    20K, 3*0 / */
/*     DATA LARGE / 77777K, 3*177777K / */
/*     DATA RIGHT / 31420K, 3*0 / */
/*     DATA DIVER / 32020K, 3*0 / */
/*     DATA LOG10 / 40423K, 42023K, 50237K, 74776K / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA DMACH(1) / '0000000000000010'X / */
/*     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X / */
/*     DATA DMACH(3) / '0000000000003CC0'X / */
/*     DATA DMACH(4) / '0000000000003CD0'X / */
/*     DATA DMACH(5) / '79FF509F44133FF3'X / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FORMAT */

/*     DATA DMACH(1) / '0010000000000000'X / */
/*     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X / */
/*     DATA DMACH(3) / '3CA0000000000000'X / */
/*     DATA DMACH(4) / '3CB0000000000000'X / */
/*     DATA DMACH(5) / '3FD34413509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/ */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/ */
/*     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/ */
/*     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/ */
/*     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/ */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /        128,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /       9344,           0 / */
/*     DATA DIVER(1), DIVER(2) /       9472,           0 / */
/*     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 / */

/*     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /         16,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /      15552,           0 / */
/*     DATA DIVER(1), DIVER(2) /      15568,           0 / */
/*     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 / */

/*     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */
/*     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION) */

/*     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X / */
/*     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X / */
/*     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X / */
/*     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X / */
/*     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA SMALL(1), SMALL(2) / '20000000, '00000201 / */
/*     DATA LARGE(1), LARGE(2) / '37777777, '37777577 / */
/*     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 / */
/*     DATA DIVER(1), DIVER(2) / '20000000, '00000334 / */
/*     DATA LOG10(1), LOG10(2) / '23210115, '10237777 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     THREE WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 / */
/*     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B / */
/*     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B / */
/*     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) /  40000B,       0 / */
/*     DATA SMALL(3), SMALL(4) /       0,       1 / */
/*     DATA LARGE(1), LARGE(2) /  77777B, 177777B / */
/*     DATA LARGE(3), LARGE(4) / 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) /  40000B,       0 / */
/*     DATA RIGHT(3), RIGHT(4) /       0,    225B / */
/*     DATA DIVER(1), DIVER(2) /  40000B,       0 / */
/*     DATA DIVER(3), DIVER(4) /       0,    227B / */
/*     DATA LOG10(1), LOG10(2) /  46420B,  46502B / */
/*     DATA LOG10(3), LOG10(4) /  76747B, 176377B / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B / */
/*     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B / */
/*     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B / */
/*     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B / */
/*     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF / */

/*     MACHINE CONSTANTS FOR THE IBM PC */
/*     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION */
/*     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087. */

/*     DATA SMALL(1) / 2.23D-308  / */
/*     DATA LARGE(1) / 1.79D+308  / */
/*     DATA RIGHT(1) / 1.11D-16   / */
/*     DATA DIVER(1) / 2.22D-16   / */
/*     DATA LOG10(1) / 0.301029995663981195D0 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    8388608,           0 / */
/*     DATA LARGE(1), LARGE(2) / 2147483647,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /  612368384,           0 / */
/*     DATA DIVER(1), DIVER(2) /  620756992,           0 / */
/*     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 / */

/*     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 / */
/*     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 / */
/*     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 / */
/*     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    128,      0 / */
/*     DATA SMALL(3), SMALL(4) /      0,      0 / */
/*     DATA LARGE(1), LARGE(2) /  32767,     -1 / */
/*     DATA LARGE(3), LARGE(4) /     -1,     -1 / */
/*     DATA RIGHT(1), RIGHT(2) /   9344,      0 / */
/*     DATA RIGHT(3), RIGHT(4) /      0,      0 / */
/*     DATA DIVER(1), DIVER(2) /   9472,      0 / */
/*     DATA DIVER(3), DIVER(4) /      0,      0 / */
/*     DATA LOG10(1), LOG10(2) /  16282,   8346 / */
/*     DATA LOG10(3), LOG10(4) / -31493, -12296 / */

/*     DATA SMALL(1), SMALL(2) / O000200, O000000 / */
/*     DATA SMALL(3), SMALL(4) / O000000, O000000 / */
/*     DATA LARGE(1), LARGE(2) / O077777, O177777 / */
/*     DATA LARGE(3), LARGE(4) / O177777, O177777 / */
/*     DATA RIGHT(1), RIGHT(2) / O022200, O000000 / */
/*     DATA RIGHT(3), RIGHT(4) / O000000, O000000 / */
/*     DATA DIVER(1), DIVER(2) / O022400, O000000 / */
/*     DATA DIVER(3), DIVER(4) / O000000, O000000 / */
/*     DATA LOG10(1), LOG10(2) / O037632, O020232 / */
/*     DATA LOG10(3), LOG10(4) / O102373, O147770 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE SUN 386i */

/*     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' / */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' / */
/*     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF' */
/*     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 / */

/* ***FIRST EXECUTABLE STATEMENT  D1MACH */
    if (*i__ < 1 || *i__ > 5) {
	xermsg_("SLATEC", "D1MACH", "I OUT OF BOUNDS", &c__1, &c__2, (ftnlen)
		6, (ftnlen)6, (ftnlen)15);
    }

    ret_val = dmach[*i__ - 1];
    return ret_val;

} /* d1mach_ */

#undef right
#undef diver
#undef small
#undef large
#undef dmach
#undef log10


/* DECK DASUM */
doublereal dasum_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i__, m, ix, mp1;

/* ***BEGIN PROLOGUE  DASUM */
/* ***PURPOSE  Compute the sum of the magnitudes of the elements of a */
/*            vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A3A */
/* ***TYPE      DOUBLE PRECISION (SASUM-S, DASUM-D, SCASUM-C) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */

/*     --Output-- */
/*    DASUM  double precision result (zero if N .LE. 0) */

/*     Returns sum of magnitudes of double precision DX. */
/*     DASUM = sum from 0 to N-1 of ABS(DX(IX+I*INCX)), */
/*     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900821  Modified to correct problem with a negative increment. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DASUM */
/* ***FIRST EXECUTABLE STATEMENT  DASUM */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
	return ret_val;
    }

    if (*incx == 1) {
	goto L20;
    }

/*     Code for increment not equal to 1. */

    ix = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += (d__1 = dx[ix], abs(d__1));
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*     Code for increment equal to 1. */

/*     Clean-up loop so remaining vector length is a multiple of 6. */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += (d__1 = dx[i__], abs(d__1));
/* L30: */
    }
    if (*n < 6) {
	return ret_val;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 6) {
	ret_val = ret_val + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1],
		 abs(d__2)) + (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = dx[
		i__ + 3], abs(d__4)) + (d__5 = dx[i__ + 4], abs(d__5)) + (
		d__6 = dx[i__ + 5], abs(d__6));
/* L50: */
    }
    return ret_val;
} /* dasum_ */

/* DECK DAXPY */
/* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, ix, iy, ns, mp1;

/* ***BEGIN PROLOGUE  DAXPY */
/* ***PURPOSE  Compute a constant times a vector plus a vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A7 */
/* ***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DA  double precision scalar multiplier */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */
/*       DY  double precision vector with N elements */
/*     INCY  storage spacing between elements of DY */

/*     --Output-- */
/*       DY  double precision result (unchanged if N .LE. 0) */

/*     Overwrite double precision DY with double precision DA*DX + DY. */
/*     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) + */
/*       DY(LY+I*INCY), */
/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DAXPY */
/* ***FIRST EXECUTABLE STATEMENT  DAXPY */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0 || *da == 0.) {
	return 0;
    }
    if (*incx == *incy) {
	if ((i__1 = *incx - 1) < 0) {
	    goto L5;
	} else if (i__1 == 0) {
	    goto L20;
	} else {
	    goto L60;
	}
    }

/*     Code for unequal or nonpositive increments. */

L5:
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*     Code for both increments equal to 1. */

/*     Clean-up loop so remaining vector length is a multiple of 4. */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;

/*     Code for equal, positive, non-unit increments. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dy[i__] = *da * dx[i__] + dy[i__];
/* L70: */
    }
    return 0;
} /* daxpy_ */

/* DECK DCOPY */
/* Subroutine */ int dcopy_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, ix, iy, ns, mp1;

/* ***BEGIN PROLOGUE  DCOPY */
/* ***PURPOSE  Copy a vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A5 */
/* ***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I) */
/* ***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */
/*       DY  double precision vector with N elements */
/*     INCY  storage spacing between elements of DY */

/*     --Output-- */
/*       DY  copy of vector DX (unchanged if N .LE. 0) */

/*     Copy double precision DX to double precision DY. */
/*     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY), */
/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DCOPY */
/* ***FIRST EXECUTABLE STATEMENT  DCOPY */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == *incy) {
	if ((i__1 = *incx - 1) < 0) {
	    goto L5;
	} else if (i__1 == 0) {
	    goto L20;
	} else {
	    goto L60;
	}
    }

/*     Code for unequal or nonpositive increments. */

L5:
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*     Code for both increments equal to 1. */

/*     Clean-up loop so remaining vector length is a multiple of 7. */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] = dx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	dy[i__] = dx[i__];
	dy[i__ + 1] = dx[i__ + 1];
	dy[i__ + 2] = dx[i__ + 2];
	dy[i__ + 3] = dx[i__ + 3];
	dy[i__ + 4] = dx[i__ + 4];
	dy[i__ + 5] = dx[i__ + 5];
	dy[i__ + 6] = dx[i__ + 6];
/* L50: */
    }
    return 0;

/*     Code for equal, positive, non-unit increments. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dy[i__] = dx[i__];
/* L70: */
    }
    return 0;
} /* dcopy_ */

/* DECK DDOT */
doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, ns, mp1;

/* ***BEGIN PROLOGUE  DDOT */
/* ***PURPOSE  Compute the inner product of two vectors. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A4 */
/* ***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C) */
/* ***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */
/*       DY  double precision vector with N elements */
/*     INCY  storage spacing between elements of DY */

/*     --Output-- */
/*     DDOT  double precision dot product (zero if N .LE. 0) */

/*     Returns the dot product of double precision DX and DY. */
/*     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY), */
/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DDOT */
/* ***FIRST EXECUTABLE STATEMENT  DDOT */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == *incy) {
	if ((i__1 = *incx - 1) < 0) {
	    goto L5;
	} else if (i__1 == 0) {
	    goto L20;
	} else {
	    goto L60;
	}
    }

/*     Code for unequal or nonpositive increments. */

L5:
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return ret_val;

/*     Code for both increments equal to 1. */

/*     Clean-up loop so remaining vector length is a multiple of 5. */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	return ret_val;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	ret_val = ret_val + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + 
		dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[
		i__ + 4] * dy[i__ + 4];
/* L50: */
    }
    return ret_val;

/*     Code for equal, positive, non-unit increments. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	ret_val += dx[i__] * dy[i__];
/* L70: */
    }
    return ret_val;
} /* ddot_ */

/* DECK DH12 */
/* Subroutine */ int dh12_(integer *mode, integer *lpivot, integer *l1, 
	integer *m, doublereal *u, integer *iue, doublereal *up, doublereal *
	c__, integer *ice, integer *icv, integer *ncv)
{
    /* System generated locals */
    integer u_dim1, u_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal b;
    static integer i__, j, i2, i3, i4;
    static doublereal cl, sm;
    static integer kl1, kl2, l1m1;
    static doublereal one;
    static integer klp;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer incr;
    static doublereal ul1m1;
    static integer mml1p2;
    static doublereal clinv;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DH12 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DHFTI, DLSEI and DWNNLS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (H12-S, DH12-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*      *** DOUBLE PRECISION VERSION OF H12 ****** */

/*     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12 */
/*     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974 */

/*     Construction and/or application of a single */
/*     Householder transformation..     Q = I + U*(U**T)/B */

/*     MODE    = 1 or 2   to select algorithm  H1  or  H2 . */
/*     LPIVOT is the index of the pivot element. */
/*     L1,M   If L1 .LE. M   the transformation will be constructed to */
/*            zero elements indexed from L1 through M.   If L1 GT. M */
/*            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION. */
/*     U(),IUE,UP    On entry to H1 U() contains the pivot vector. */
/*                   IUE is the storage increment between elements. */
/*                                       On exit from H1 U() and UP */
/*                   contain quantities defining the vector U of the */
/*                   Householder transformation.   On entry to H2 U() */
/*                   and UP should contain quantities previously computed */
/*                   by H1.  These will not be modified by H2. */
/*     C()    On entry to H1 or H2 C() contains a matrix which will be */
/*            regarded as a set of vectors to which the Householder */
/*            transformation is to be applied.  On exit C() contains the */
/*            set of transformed vectors. */
/*     ICE    Storage increment between elements of vectors in C(). */
/*     ICV    Storage increment between vectors in C(). */
/*     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0 */
/*            no operations will be done on C(). */

/* ***SEE ALSO  DHFTI, DLSEI, DWNNLS */
/* ***ROUTINES CALLED  DAXPY, DDOT, DSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900911  Added DDOT to DOUBLE PRECISION statement.  (WRB) */
/* ***END PROLOGUE  DH12 */
/*     BEGIN BLOCK PERMITTING ...EXITS TO 140 */
/* ***FIRST EXECUTABLE STATEMENT  DH12 */
    /* Parameter adjustments */
    u_dim1 = *iue;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --c__;

    /* Function Body */
    one = 1.;

/*     ...EXIT */
    if (0 >= *lpivot || *lpivot >= *l1 || *l1 > *m) {
	goto L140;
    }
    cl = (d__1 = u[*lpivot * u_dim1 + 1], abs(d__1));
    if (*mode == 2) {
	goto L40;
    }
/*           ****** CONSTRUCT THE TRANSFORMATION. ****** */
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = (d__1 = u[j * u_dim1 + 1], abs(d__1));
	cl = max(d__2,cl);
/* L10: */
    }
    if (cl > 0.) {
	goto L20;
    }
/*     .........EXIT */
    goto L140;
L20:
    clinv = one / cl;
/* Computing 2nd power */
    d__1 = u[*lpivot * u_dim1 + 1] * clinv;
    sm = d__1 * d__1;
    i__1 = *m;
    for (j = *l1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = u[j * u_dim1 + 1] * clinv;
	sm += d__1 * d__1;
/* L30: */
    }
    cl *= sqrt(sm);
    if (u[*lpivot * u_dim1 + 1] > 0.) {
	cl = -cl;
    }
    *up = u[*lpivot * u_dim1 + 1] - cl;
    u[*lpivot * u_dim1 + 1] = cl;
    goto L50;
L40:
/*        ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ****** */

    if (cl > 0.) {
	goto L50;
    }
/*     ......EXIT */
    goto L140;
L50:
/*     ...EXIT */
    if (*ncv <= 0) {
	goto L140;
    }
    b = *up * u[*lpivot * u_dim1 + 1];
/*        B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN. */

    if (b < 0.) {
	goto L60;
    }
/*     ......EXIT */
    goto L140;
L60:
    b = one / b;
    mml1p2 = *m - *l1 + 2;
    if (mml1p2 <= 20) {
	goto L80;
    }
    l1m1 = *l1 - 1;
    kl1 = (l1m1 - 1) * *ice + 1;
    kl2 = kl1;
    klp = (*lpivot - 1) * *ice + 1;
    ul1m1 = u[l1m1 * u_dim1 + 1];
    u[l1m1 * u_dim1 + 1] = *up;
    if (*lpivot != l1m1) {
	dswap_(ncv, &c__[kl1], icv, &c__[klp], icv);
    }
    i__1 = *ncv;
    for (j = 1; j <= i__1; ++j) {
	sm = ddot_(&mml1p2, &u[l1m1 * u_dim1 + 1], iue, &c__[kl1], ice);
	sm *= b;
	daxpy_(&mml1p2, &sm, &u[l1m1 * u_dim1 + 1], iue, &c__[kl1], ice);
	kl1 += *icv;
/* L70: */
    }
    u[l1m1 * u_dim1 + 1] = ul1m1;
/*     ......EXIT */
    if (*lpivot == l1m1) {
	goto L140;
    }
    kl1 = kl2;
    dswap_(ncv, &c__[kl1], icv, &c__[klp], icv);
    goto L130;
L80:
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
	    i3 += *ice;
/* L90: */
	}
	if (sm == 0.) {
	    goto L110;
	}
	sm *= b;
	c__[i2] += sm * *up;
	i__2 = *m;
	for (i__ = *l1; i__ <= i__2; ++i__) {
	    c__[i4] += sm * u[i__ * u_dim1 + 1];
	    i4 += *ice;
/* L100: */
	}
L110:
/* L120: */
	;
    }
L130:
L140:
    return 0;
} /* dh12_ */

/* DECK DHFTI */
/* Subroutine */ int dhfti_(doublereal *a, integer *mda, integer *m, integer *
	n, doublereal *b, integer *mdb, integer *nb, doublereal *tau, integer 
	*krank, doublereal *rnorm, doublereal *h__, doublereal *g, integer *
	ip)
{
    /* Initialized data */

    static doublereal releps = 0.;

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, jb, ii, jj;
    static doublereal sm;
    static integer ip1, kp1;
    static doublereal sm1;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal tmp, hmax;
    static integer lmax, nerr, iopt, ldiag;
    static doublereal dzero;
    extern doublereal d1mach_(integer *);
    static doublereal szero, factor;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DHFTI */
/* ***PURPOSE  Solve a least squares problem for banded matrices using */
/*            sequential accumulation of rows of the data matrix. */
/*            Exactly one right-hand side vector is permitted. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D9 */
/* ***TYPE      DOUBLE PRECISION (HFTI-S, DHFTI-D) */
/* ***KEYWORDS  CURVE FITTING, LEAST SQUARES */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N) */

/*     This subroutine solves a linear least squares problem or a set of */
/*     linear least squares problems having the same matrix but different */
/*     right-side vectors.  The problem data consists of an M by N matrix */
/*     A, an M by NB matrix B, and an absolute tolerance parameter TAU */
/*     whose usage is described below.  The NB column vectors of B */
/*     represent right-side vectors for NB distinct linear least squares */
/*     problems. */

/*     This set of problems can also be written as the matrix least */
/*     squares problem */

/*                       AX = B, */

/*     where X is the N by NB solution matrix. */

/*     Note that if B is the M by M identity matrix, then X will be the */
/*     pseudo-inverse of A. */

/*     This subroutine first transforms the augmented matrix (A B) to a */
/*     matrix (R C) using premultiplying Householder transformations with */
/*     column interchanges.  All subdiagonal elements in the matrix R are */
/*     zero and its diagonal elements satisfy */

/*                       ABS(R(I,I)).GE.ABS(R(I+1,I+1)), */

/*                       I = 1,...,L-1, where */

/*                       L = MIN(M,N). */

/*     The subroutine will compute an integer, KRANK, equal to the number */
/*     of diagonal terms of R that exceed TAU in magnitude. Then a */
/*     solution of minimum Euclidean length is computed using the first */
/*     KRANK rows of (R C). */

/*     To be specific we suggest that the user consider an easily */
/*     computable matrix norm, such as, the maximum of all column sums of */
/*     magnitudes. */

/*     Now if the relative uncertainty of B is EPS, (norm of uncertainty/ */
/*     norm of B), it is suggested that TAU be set approximately equal to */
/*     EPS*(norm of A). */

/*     The user must dimension all arrays appearing in the call list.. */
/*     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This */
/*     permits the solution of a range of problems in the same array */
/*     space. */

/*     The entire set of parameters for DHFTI are */

/*     INPUT.. All TYPE REAL variables are DOUBLE PRECISION */

/*     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N */
/*                       matrix A of the least squares problem AX = B. */
/*                       The first dimensioning parameter of the array */
/*                       A(*,*) is MDA, which must satisfy MDA.GE.M */
/*                       Either M.GE.N or M.LT.N is permitted.  There */
/*                       is no restriction on the rank of A.  The */
/*                       condition MDA.LT.M is considered an error. */

/*     B(*),MDB,NB       If NB = 0 the subroutine will perform the */
/*                       orthogonal decomposition but will make no */
/*                       references to the array B(*).  If NB.GT.0 */
/*                       the array B(*) must initially contain the M by */
/*                       NB matrix B of the least squares problem AX = */
/*                       B.  If NB.GE.2 the array B(*) must be doubly */
/*                       subscripted with first dimensioning parameter */
/*                       MDB.GE.MAX(M,N).  If NB = 1 the array B(*) may */
/*                       be either doubly or singly subscripted.  In */
/*                       the latter case the value of MDB is arbitrary */
/*                       but it should be set to some valid integer */
/*                       value such as MDB = M. */

/*                       The condition of NB.GT.1.AND.MDB.LT. MAX(M,N) */
/*                       is considered an error. */

/*     TAU               Absolute tolerance parameter provided by user */
/*                       for pseudorank determination. */

/*     H(*),G(*),IP(*)   Arrays of working space used by DHFTI. */

/*     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION */

/*     A(*,*)            The contents of the array A(*,*) will be */
/*                       modified by the subroutine. These contents */
/*                       are not generally required by the user. */

/*     B(*)              On return the array B(*) will contain the N by */
/*                       NB solution matrix X. */

/*     KRANK             Set by the subroutine to indicate the */
/*                       pseudorank of A. */

/*     RNORM(*)          On return, RNORM(J) will contain the Euclidean */
/*                       norm of the residual vector for the problem */
/*                       defined by the J-th column vector of the array */
/*                       B(*,*) for J = 1,...,NB. */

/*     H(*),G(*)         On return these arrays respectively contain */
/*                       elements of the pre- and post-multiplying */
/*                       Householder transformations used to compute */
/*                       the minimum Euclidean length solution. */

/*     IP(*)             Array in which the subroutine records indices */
/*                       describing the permutation of column vectors. */
/*                       The contents of arrays H(*),G(*) and IP(*) */
/*                       are not generally required by the user. */

/* ***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares */
/*                 Problems, Prentice-Hall, Inc., 1974, Chapter 14. */
/* ***ROUTINES CALLED  D1MACH, DH12, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   901005  Replace usage of DDIFF with usage of D1MACH.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DHFTI */
    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *mdb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --rnorm;
    --h__;
    --g;
    --ip;

    /* Function Body */
/*     BEGIN BLOCK PERMITTING ...EXITS TO 360 */
/* ***FIRST EXECUTABLE STATEMENT  DHFTI */
    if (releps == 0.) {
	releps = d1mach_(&c__4);
    }
    szero = 0.;
    dzero = 0.;
    factor = .001;

    k = 0;
    ldiag = min(*m,*n);
    if (ldiag <= 0) {
	goto L350;
    }
/*           BEGIN BLOCK PERMITTING ...EXITS TO 130 */
/*              BEGIN BLOCK PERMITTING ...EXITS TO 120 */
    if (*mda >= *m) {
	goto L10;
    }
    nerr = 1;
    iopt = 2;
    xermsg_("SLATEC", "DHFTI", "MDA.LT.M, PROBABLE ERROR.", &nerr, &iopt, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
/*     ...............EXIT */
    goto L360;
L10:

    if (*nb <= 1 || max(*m,*n) <= *mdb) {
	goto L20;
    }
    nerr = 2;
    iopt = 2;
    xermsg_("SLATEC", "DHFTI", "MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERROR.",
	     &nerr, &iopt, (ftnlen)6, (ftnlen)5, (ftnlen)44);
/*     ...............EXIT */
    goto L360;
L20:

    i__1 = ldiag;
    for (j = 1; j <= i__1; ++j) {
/*                    BEGIN BLOCK PERMITTING ...EXITS TO 70 */
	if (j == 1) {
	    goto L40;
	}

/*                           UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX */
/*                          .. */
	lmax = j;
	i__2 = *n;
	for (l = j; l <= i__2; ++l) {
/* Computing 2nd power */
	    d__1 = a[j - 1 + l * a_dim1];
	    h__[l] -= d__1 * d__1;
	    if (h__[l] > h__[lmax]) {
		lmax = l;
	    }
/* L30: */
	}
/*                    ......EXIT */
	if (factor * h__[lmax] > hmax * releps) {
	    goto L70;
	}
L40:

/*                        COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX */
/*                       .. */
	lmax = j;
	i__2 = *n;
	for (l = j; l <= i__2; ++l) {
	    h__[l] = 0.;
	    i__3 = *m;
	    for (i__ = j; i__ <= i__3; ++i__) {
/* Computing 2nd power */
		d__1 = a[i__ + l * a_dim1];
		h__[l] += d__1 * d__1;
/* L50: */
	    }
	    if (h__[l] > h__[lmax]) {
		lmax = l;
	    }
/* L60: */
	}
	hmax = h__[lmax];
L70:
/*                    .. */
/*                     LMAX HAS BEEN DETERMINED */

/*                     DO COLUMN INTERCHANGES IF NEEDED. */
/*                    .. */
	ip[j] = lmax;
	if (ip[j] == j) {
	    goto L90;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tmp = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = a[i__ + lmax * a_dim1];
	    a[i__ + lmax * a_dim1] = tmp;
/* L80: */
	}
	h__[lmax] = h__[j];
L90:

/*                     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A */
/*                     AND B. */
/*                    .. */
	i__2 = j + 1;
	i__3 = *n - j;
	dh12_(&c__1, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &a[(j 
		+ 1) * a_dim1 + 1], &c__1, mda, &i__3);
	i__2 = j + 1;
	dh12_(&c__2, &j, &i__2, m, &a[j * a_dim1 + 1], &c__1, &h__[j], &b[
		b_offset], &c__1, mdb, nb);
/* L100: */
    }

/*                  DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, */
/*                  TAU. */
/*                 .. */
    i__1 = ldiag;
    for (j = 1; j <= i__1; ++j) {
/*              ......EXIT */
	if ((d__1 = a[j + j * a_dim1], abs(d__1)) <= *tau) {
	    goto L120;
	}
/* L110: */
    }
    k = ldiag;
/*           ......EXIT */
    goto L130;
L120:
    k = j - 1;
L130:
    kp1 = k + 1;

/*           COMPUTE THE NORMS OF THE RESIDUAL VECTORS. */

    if (*nb < 1) {
	goto L170;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
	tmp = szero;
	if (*m < kp1) {
	    goto L150;
	}
	i__2 = *m;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    d__1 = b[i__ + jb * b_dim1];
	    tmp += d__1 * d__1;
/* L140: */
	}
L150:
	rnorm[jb] = sqrt(tmp);
/* L160: */
    }
L170:
/*           SPECIAL FOR PSEUDORANK = 0 */
    if (k > 0) {
	goto L210;
    }
    if (*nb < 1) {
	goto L200;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    b[i__ + jb * b_dim1] = szero;
/* L180: */
	}
/* L190: */
    }
L200:
    goto L340;
L210:

/*               IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER */
/*               DECOMPOSITION OF FIRST K ROWS. */
/*              .. */
    if (k == *n) {
	goto L230;
    }
    i__1 = k;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = kp1 - ii;
	i__2 = i__ - 1;
	dh12_(&c__1, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &a[
		a_offset], mda, &c__1, &i__2);
/* L220: */
    }
L230:


    if (*nb < 1) {
	goto L330;
    }
    i__1 = *nb;
    for (jb = 1; jb <= i__1; ++jb) {

/*                  SOLVE THE K BY K TRIANGULAR SYSTEM. */
/*                 .. */
	i__2 = k;
	for (l = 1; l <= i__2; ++l) {
	    sm = dzero;
	    i__ = kp1 - l;
	    ip1 = i__ + 1;
	    if (k < ip1) {
		goto L250;
	    }
	    i__3 = k;
	    for (j = ip1; j <= i__3; ++j) {
		sm += a[i__ + j * a_dim1] * b[j + jb * b_dim1];
/* L240: */
	    }
L250:
	    sm1 = sm;
	    b[i__ + jb * b_dim1] = (b[i__ + jb * b_dim1] - sm1) / a[i__ + i__ 
		    * a_dim1];
/* L260: */
	}

/*                  COMPLETE COMPUTATION OF SOLUTION VECTOR. */
/*                 .. */
	if (k == *n) {
	    goto L290;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    b[j + jb * b_dim1] = szero;
/* L270: */
	}
	i__2 = k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dh12_(&c__2, &i__, &kp1, n, &a[i__ + a_dim1], mda, &g[i__], &b[jb 
		    * b_dim1 + 1], &c__1, mdb, &c__1);
/* L280: */
	}
L290:

/*                   RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE */
/*                   COLUMN INTERCHANGES. */
/*                 .. */
	i__2 = ldiag;
	for (jj = 1; jj <= i__2; ++jj) {
	    j = ldiag + 1 - jj;
	    if (ip[j] == j) {
		goto L300;
	    }
	    l = ip[j];
	    tmp = b[l + jb * b_dim1];
	    b[l + jb * b_dim1] = b[j + jb * b_dim1];
	    b[j + jb * b_dim1] = tmp;
L300:
/* L310: */
	    ;
	}
/* L320: */
    }
L330:
L340:
L350:
/*        .. */
/*         THE SOLUTION VECTORS, X, ARE NOW */
/*         IN THE FIRST  N  ROWS OF THE ARRAY B(,). */

    *krank = k;
L360:
    return 0;
} /* dhfti_ */

/* DECK DLPDP */
/* Subroutine */ int dlpdp_(doublereal *a, integer *mda, integer *m, integer *
	n1, integer *n2, doublereal *prgopt, doublereal *x, doublereal *wnorm,
	 integer *mode, doublereal *ws, integer *is)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal fac = .1;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l, n;
    static doublereal sc;
    static integer iw, ix, np1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer modew;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal rnorm, ynorm;
    extern /* Subroutine */ int dwnnls_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *);

/* ***BEGIN PROLOGUE  DLPDP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DLSEI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LPDP-S, DLPDP-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*  **** Double Precision version of LPDP **** */
/*     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1), */
/*     where N=N1+N2.  This is a slight overestimate for WS(*). */

/*     Determine an N1-vector W, and */
/*               an N2-vector Z */
/*     which minimizes the Euclidean length of W */
/*     subject to G*W+H*Z .GE. Y. */
/*     This is the least projected distance problem, LPDP. */
/*     The matrices G and H are of respective */
/*     dimensions M by N1 and M by N2. */

/*     Called by subprogram DLSI( ). */

/*     The matrix */
/*                (G H Y) */

/*     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*). */

/*     The solution (W) is returned in X(*). */
/*                  (Z) */

/*     The value of MODE indicates the status of */
/*     the computation after returning to the user. */

/*          MODE=1  The solution was successfully obtained. */

/*          MODE=2  The inequalities are inconsistent. */

/* ***SEE ALSO  DLSEI */
/* ***ROUTINES CALLED  DCOPY, DDOT, DNRM2, DSCAL, DWNNLS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  DLPDP */

    /* Parameter adjustments */
    a_dim1 = *mda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --prgopt;
    --x;
    --ws;
    --is;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DLPDP */
    n = *n1 + *n2;
    *mode = 1;
    if (*m > 0) {
	goto L20;
    }
    if (n <= 0) {
	goto L10;
    }
    x[1] = zero;
    dcopy_(&n, &x[1], &c__0, &x[1], &c__1);
L10:
    *wnorm = zero;
    goto L200;
L20:
/*        BEGIN BLOCK PERMITTING ...EXITS TO 190 */
    np1 = n + 1;

/*           SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sc = dnrm2_(&n, &a[i__ + a_dim1], mda);
	if (sc == zero) {
	    goto L30;
	}
	sc = one / sc;
	dscal_(&np1, &sc, &a[i__ + a_dim1], mda);
L30:
/* L40: */
	;
    }

/*           SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO). */
    ynorm = dnrm2_(m, &a[np1 * a_dim1 + 1], &c__1);
    if (ynorm == zero) {
	goto L50;
    }
    sc = one / ynorm;
    dscal_(m, &sc, &a[np1 * a_dim1 + 1], &c__1);
L50:

/*           SCALE COLS OF MATRIX H. */
    j = *n1 + 1;
L60:
    if (j > n) {
	goto L70;
    }
    sc = dnrm2_(m, &a[j * a_dim1 + 1], &c__1);
    if (sc != zero) {
	sc = one / sc;
    }
    dscal_(m, &sc, &a[j * a_dim1 + 1], &c__1);
    x[j] = sc;
    ++j;
    goto L60;
L70:
    if (*n1 <= 0) {
	goto L130;
    }

/*              COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*). */
    iw = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*                 MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY. */
	dcopy_(n2, &a[i__ + (*n1 + 1) * a_dim1], mda, &ws[iw + 1], &c__1);
	iw += *n2;

/*                 MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY. */
	dcopy_(n1, &a[i__ + a_dim1], mda, &ws[iw + 1], &c__1);
	iw += *n1;

/*                 MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY. */
	ws[iw + 1] = a[i__ + np1 * a_dim1];
	++iw;
/* L80: */
    }
    ws[iw + 1] = zero;
    dcopy_(&n, &ws[iw + 1], &c__0, &ws[iw + 1], &c__1);
    iw += n;
    ws[iw + 1] = one;
    ++iw;

/*              SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE */
/*              MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR */
/*              F = TRANSPOSE OF (0,...,0,1). */
    ix = iw + 1;
    iw += *m;

/*              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF */
/*              DWNNLS( ). */
    is[1] = 0;
    is[2] = 0;
    i__1 = np1 - *n2;
    dwnnls_(&ws[1], &np1, n2, &i__1, m, &c__0, &prgopt[1], &ws[ix], &rnorm, &
	    modew, &is[1], &ws[iw + 1]);

/*              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W. */
    sc = one - ddot_(m, &a[np1 * a_dim1 + 1], &c__1, &ws[ix], &c__1);
    if (one + fac * abs(sc) == one || rnorm <= zero) {
	goto L110;
    }
    sc = one / sc;
    i__1 = *n1;
    for (j = 1; j <= i__1; ++j) {
	x[j] = sc * ddot_(m, &a[j * a_dim1 + 1], &c__1, &ws[ix], &c__1);
/* L90: */
    }

/*                 COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS */
/*                 VECTOR. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ + np1 * a_dim1] -= ddot_(n1, &a[i__ + a_dim1], mda, &x[1], &
		c__1);
/* L100: */
    }
    goto L120;
L110:
    *mode = 2;
/*        .........EXIT */
    goto L190;
L120:
L130:
    if (*n2 <= 0) {
	goto L180;
    }

/*              COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*). */
    iw = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dcopy_(n2, &a[i__ + (*n1 + 1) * a_dim1], mda, &ws[iw + 1], &c__1);
	iw += *n2;
	ws[iw + 1] = a[i__ + np1 * a_dim1];
	++iw;
/* L140: */
    }
    ws[iw + 1] = zero;
    dcopy_(n2, &ws[iw + 1], &c__0, &ws[iw + 1], &c__1);
    iw += *n2;
    ws[iw + 1] = one;
    ++iw;
    ix = iw + 1;
    iw += *m;

/*              SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE */
/*              OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE */
/*              OF (0,...,0,1)). */

/*              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF */
/*              DWNNLS( ). */
    is[1] = 0;
    is[2] = 0;
    i__1 = *n2 + 1;
    i__2 = *n2 + 1;
    dwnnls_(&ws[1], &i__1, &c__0, &i__2, m, &c__0, &prgopt[1], &ws[ix], &
	    rnorm, &modew, &is[1], &ws[iw + 1]);

/*              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z. */
    sc = one - ddot_(m, &a[np1 * a_dim1 + 1], &c__1, &ws[ix], &c__1);
    if (one + fac * abs(sc) == one || rnorm <= zero) {
	goto L160;
    }
    sc = one / sc;
    i__1 = *n2;
    for (j = 1; j <= i__1; ++j) {
	l = *n1 + j;
	x[l] = sc * ddot_(m, &a[l * a_dim1 + 1], &c__1, &ws[ix], &c__1) * x[l]
		;
/* L150: */
    }
    goto L170;
L160:
    *mode = 2;
/*        .........EXIT */
    goto L190;
L170:
L180:

/*           ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION. */
    dscal_(&n, &ynorm, &x[1], &c__1);
    *wnorm = dnrm2_(n1, &x[1], &c__1);
L190:
L200:
    return 0;
} /* dlpdp_ */

/* DECK DLSI */
/* Subroutine */ int dlsi_(doublereal *w, integer *mdw, integer *ma, integer *
	mg, integer *n, doublereal *prgopt, doublereal *x, doublereal *rnorm, 
	integer *mode, doublereal *ws, integer *ip)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, m, n1, n2, n3;
    static doublereal rb;
    static integer np1;
    static doublereal fac;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal gam;
    static integer key;
    static doublereal tau;
    static logical cov;
    static doublereal tol;
    static integer map1, krm1, krp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer link, last, next;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dhfti_(doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *), dlpdp_(
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *);
    static integer krank;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static doublereal anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    extern doublereal d1mach_(integer *);
    static doublereal xnorm;
    static integer minman, mdlpdp;
    static doublereal drelpr;
    static logical sclcov;

/* ***BEGIN PROLOGUE  DLSI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DLSEI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LSI-S, DLSI-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to DLSEI.  The documentation for */
/*     DLSEI has complete usage instructions. */

/*     Solve.. */
/*              AX = B,  A  MA by N  (least squares equations) */
/*     subject to.. */

/*              GX.GE.H, G  MG by N  (inequality constraints) */

/*     Input.. */

/*      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1. */
/*                       (G H) */

/*     MDW,MA,MG,N */
/*              contain (resp) var. dimension of W(*,*), */
/*              and matrix dimensions. */

/*     PRGOPT(*), */
/*              Program option vector. */

/*     OUTPUT.. */

/*      X(*),RNORM */

/*              Solution vector(unless MODE=2), length of AX-B. */

/*      MODE */
/*              =0   Inequality constraints are compatible. */
/*              =2   Inequality constraints contradictory. */

/*      WS(*), */
/*              Working storage of dimension K+N+(MG+2)*(N+7), */
/*              where K=MAX(MA+MG,N). */
/*      IP(MG+2*N+1) */
/*              Integer working storage */

/* ***ROUTINES CALLED  D1MACH, DASUM, DAXPY, DCOPY, DDOT, DH12, DHFTI, */
/*                    DLPDP, DSCAL, DSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and extensively revised (WRB & RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900604  DP version created from SP version.  (RWC) */
/*   920422  Changed CALL to DHFTI to include variable MA.  (WRB) */
/* ***END PROLOGUE  DLSI */



    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --prgopt;
    --x;
    --ws;
    --ip;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  DLSI */

/*     Set the nominal tolerance used in the code. */

    if (first) {
	drelpr = d1mach_(&c__4);
    }
    first = FALSE_;
    tol = sqrt(drelpr);

    *mode = 0;
    *rnorm = 0.;
    m = *ma + *mg;
    np1 = *n + 1;
    krank = 0;
    if (*n <= 0 || m <= 0) {
	goto L370;
    }

/*     To process option vector. */

    cov = FALSE_;
    sclcov = TRUE_;
    last = 1;
    link = (integer) prgopt[1];

L100:
    if (link > 1) {
	key = (integer) prgopt[last + 1];
	if (key == 1) {
	    cov = prgopt[last + 2] != 0.;
	}
	if (key == 10) {
	    sclcov = prgopt[last + 2] == 0.;
	}
	if (key == 5) {
/* Computing MAX */
	    d__1 = drelpr, d__2 = prgopt[last + 2];
	    tol = max(d__1,d__2);
	}
	next = (integer) prgopt[link];
	last = link;
	link = next;
	goto L100;
    }

/*     Compute matrix norm of least squares equations. */

    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = dasum_(ma, &w[j * w_dim1 + 1], &c__1);
	anorm = max(d__1,d__2);
/* L110: */
    }

/*     Set tolerance for DHFTI( ) rank test. */

    tau = tol * anorm;

/*     Compute Householder orthogonal decomposition of matrix. */

    dcopy_(n, &c_b95, &c__0, &ws[1], &c__1);
    dcopy_(ma, &w[np1 * w_dim1 + 1], &c__1, &ws[1], &c__1);
    k = max(m,*n);
    minman = min(*ma,*n);
    n1 = k + 1;
    n2 = n1 + *n;
    dhfti_(&w[w_offset], mdw, ma, n, &ws[1], ma, &c__1, &tau, &krank, rnorm, &
	    ws[n2], &ws[n1], &ip[1]);
    fac = 1.;
    gam = (doublereal) (*ma - krank);
    if (krank < *ma && sclcov) {
/* Computing 2nd power */
	d__1 = *rnorm;
	fac = d__1 * d__1 / gam;
    }

/*     Reduce to DLPDP and solve. */

    map1 = *ma + 1;

/*     Compute inequality rt-hand side for DLPDP. */

    if (*ma < m) {
	if (minman > 0) {
	    i__1 = m;
	    for (i__ = map1; i__ <= i__1; ++i__) {
		w[i__ + np1 * w_dim1] -= ddot_(n, &w[i__ + w_dim1], mdw, &ws[
			1], &c__1);
/* L120: */
	    }

/*           Apply permutations to col. of inequality constraint matrix. */

	    i__1 = minman;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dswap_(mg, &w[map1 + i__ * w_dim1], &c__1, &w[map1 + ip[i__] *
			 w_dim1], &c__1);
/* L130: */
	    }

/*           Apply Householder transformations to constraint matrix. */

	    if (krank > 0 && krank < *n) {
		for (i__ = krank; i__ >= 1; --i__) {
		    i__1 = krank + 1;
		    dh12_(&c__2, &i__, &i__1, n, &w[i__ + w_dim1], mdw, &ws[
			    n1 + i__ - 1], &w[map1 + w_dim1], mdw, &c__1, mg);
/* L140: */
		}
	    }

/*           Compute permuted inequality constraint matrix times r-inv. */

	    i__1 = m;
	    for (i__ = map1; i__ <= i__1; ++i__) {
		i__2 = krank;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = j - 1;
		    w[i__ + j * w_dim1] = (w[i__ + j * w_dim1] - ddot_(&i__3, 
			    &w[j * w_dim1 + 1], &c__1, &w[i__ + w_dim1], mdw))
			     / w[j + j * w_dim1];
/* L150: */
		}
/* L160: */
	    }
	}

/*        Solve the reduced problem with DLPDP algorithm, */
/*        the least projected distance problem. */

	i__1 = *n - krank;
	dlpdp_(&w[map1 + w_dim1], mdw, mg, &krank, &i__1, &prgopt[1], &x[1], &
		xnorm, &mdlpdp, &ws[n2], &ip[*n + 1]);

/*        Compute solution in original coordinates. */

	if (mdlpdp == 1) {
	    for (i__ = krank; i__ >= 1; --i__) {
		i__1 = krank - i__;
		x[i__] = (x[i__] - ddot_(&i__1, &w[i__ + (i__ + 1) * w_dim1], 
			mdw, &x[i__ + 1], &c__1)) / w[i__ + i__ * w_dim1];
/* L170: */
	    }

/*           Apply Householder transformation to solution vector. */

	    if (krank < *n) {
		i__1 = krank;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = krank + 1;
		    dh12_(&c__2, &i__, &i__2, n, &w[i__ + w_dim1], mdw, &ws[
			    n1 + i__ - 1], &x[1], &c__1, &c__1, &c__1);
/* L180: */
		}
	    }

/*           Repermute variables to their input order. */

	    if (minman > 0) {
		for (i__ = minman; i__ >= 1; --i__) {
		    dswap_(&c__1, &x[i__], &c__1, &x[ip[i__]], &c__1);
/* L190: */
		}

/*              Variables are now in original coordinates. */
/*              Add solution of unconstrained problem. */

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    x[i__] += ws[i__];
/* L200: */
		}

/*              Compute the residual vector norm. */

/* Computing 2nd power */
		d__1 = *rnorm;
/* Computing 2nd power */
		d__2 = xnorm;
		*rnorm = sqrt(d__1 * d__1 + d__2 * d__2);
	    }
	} else {
	    *mode = 2;
	}
    } else {
	dcopy_(n, &ws[1], &c__1, &x[1], &c__1);
    }

/*     Compute covariance matrix based on the orthogonal decomposition */
/*     from DHFTI( ). */

    if (! cov || krank <= 0) {
	goto L370;
    }
    krm1 = krank - 1;
    krp1 = krank + 1;

/*     Copy diagonal terms to working array. */

    i__1 = *mdw + 1;
    dcopy_(&krank, &w[w_offset], &i__1, &ws[n2], &c__1);

/*     Reciprocate diagonal terms. */

    i__1 = krank;
    for (j = 1; j <= i__1; ++j) {
	w[j + j * w_dim1] = 1. / w[j + j * w_dim1];
/* L210: */
    }

/*     Invert the upper triangular QR factor on itself. */

    if (krank > 1) {
	i__1 = krm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = krank;
	    for (j = i__ + 1; j <= i__2; ++j) {
		i__3 = j - i__;
		w[i__ + j * w_dim1] = -ddot_(&i__3, &w[i__ + i__ * w_dim1], 
			mdw, &w[i__ + j * w_dim1], &c__1) * w[j + j * w_dim1];
/* L220: */
	    }
/* L230: */
	}
    }

/*     Compute the inverted factor times its transpose. */

    i__1 = krank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = krank;
	for (j = i__; j <= i__2; ++j) {
	    i__3 = krank + 1 - j;
	    w[i__ + j * w_dim1] = ddot_(&i__3, &w[i__ + j * w_dim1], mdw, &w[
		    j + j * w_dim1], mdw);
/* L240: */
	}
/* L250: */
    }

/*     Zero out lower trapezoidal part. */
/*     Copy upper triangular to lower triangular part. */

    if (krank < *n) {
	i__1 = krank;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(&j, &w[j * w_dim1 + 1], &c__1, &w[j + w_dim1], mdw);
/* L260: */
	}

	i__1 = *n;
	for (i__ = krp1; i__ <= i__1; ++i__) {
	    dcopy_(&i__, &c_b95, &c__0, &w[i__ + w_dim1], mdw);
/* L270: */
	}

/*        Apply right side transformations to lower triangle. */

	n3 = n2 + krp1;
	i__1 = krank;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    l = n1 + i__;
	    k = n2 + i__;
	    rb = ws[l - 1] * ws[k - 1];

/*           If RB.GE.0.D0, transformation can be regarded as zero. */

	    if (rb < 0.) {
		rb = 1. / rb;

/*              Store unscaled rank one Householder update in work array. */

		dcopy_(n, &c_b95, &c__0, &ws[n3], &c__1);
		l = n1 + i__;
		k = n3 + i__;
		ws[k - 1] = ws[l - 1];

		i__2 = *n;
		for (j = krp1; j <= i__2; ++j) {
		    ws[n3 + j - 1] = w[i__ + j * w_dim1];
/* L280: */
		}

		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = j - i__;
		    i__4 = *n - j + 1;
		    ws[j] = rb * (ddot_(&i__3, &w[j + i__ * w_dim1], mdw, &ws[
			    n3 + i__ - 1], &c__1) + ddot_(&i__4, &w[j + j * 
			    w_dim1], &c__1, &ws[n3 + j - 1], &c__1));
/* L290: */
		}

		l = n3 + i__;
		i__2 = *n - i__ + 1;
		gam = rb * .5 * ddot_(&i__2, &ws[l - 1], &c__1, &ws[i__], &
			c__1);
		i__2 = *n - i__ + 1;
		daxpy_(&i__2, &gam, &ws[l - 1], &c__1, &ws[i__], &c__1);
		i__2 = *n;
		for (j = i__; j <= i__2; ++j) {
		    i__3 = i__ - 1;
		    for (l = 1; l <= i__3; ++l) {
			w[j + l * w_dim1] += ws[n3 + j - 1] * ws[l];
/* L300: */
		    }

		    i__3 = j;
		    for (l = i__; l <= i__3; ++l) {
			w[j + l * w_dim1] = w[j + l * w_dim1] + ws[j] * ws[n3 
				+ l - 1] + ws[l] * ws[n3 + j - 1];
/* L310: */
		    }
/* L320: */
		}
	    }
/* L330: */
	}

/*        Copy lower triangle to upper triangle to symmetrize the */
/*        covariance matrix. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dcopy_(&i__, &w[i__ + w_dim1], mdw, &w[i__ * w_dim1 + 1], &c__1);
/* L340: */
	}
    }

/*     Repermute rows and columns. */

    for (i__ = minman; i__ >= 1; --i__) {
	k = ip[i__];
	if (i__ != k) {
	    dswap_(&c__1, &w[i__ + i__ * w_dim1], &c__1, &w[k + k * w_dim1], &
		    c__1);
	    i__1 = i__ - 1;
	    dswap_(&i__1, &w[i__ * w_dim1 + 1], &c__1, &w[k * w_dim1 + 1], &
		    c__1);
	    i__1 = k - i__ - 1;
	    dswap_(&i__1, &w[i__ + (i__ + 1) * w_dim1], mdw, &w[i__ + 1 + k * 
		    w_dim1], &c__1);
	    i__1 = *n - k;
	    dswap_(&i__1, &w[i__ + (k + 1) * w_dim1], mdw, &w[k + (k + 1) * 
		    w_dim1], mdw);
	}
/* L350: */
    }

/*     Put in normalized residual sum of squares scale factor */
/*     and symmetrize the resulting covariance matrix. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dscal_(&j, &fac, &w[j * w_dim1 + 1], &c__1);
	dcopy_(&j, &w[j * w_dim1 + 1], &c__1, &w[j + w_dim1], mdw);
/* L360: */
    }

L370:
    ip[1] = krank;
    ip[2] = *n + max(m,*n) + (*mg + 2) * (*n + 7);
    return 0;
} /* dlsi_ */

/* DECK DNRM2 */
doublereal dnrm2_(integer *n, doublereal *dx, integer *incx)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal cutlo = 8.232e-11;
    static doublereal cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nn;
    static doublereal sum, xmax;
    static integer next;
    static doublereal hitest;

    /* Assigned format variables */
    static char *next_fmt;

/* ***BEGIN PROLOGUE  DNRM2 */
/* ***PURPOSE  Compute the Euclidean length (L2 norm) of a vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A3B */
/* ***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C) */
/* ***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2, */
/*             LINEAR ALGEBRA, UNITARY, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */

/*     --Output-- */
/*    DNRM2  double precision result (zero if N .LE. 0) */

/*     Euclidean norm of the N-vector stored in DX with storage */
/*     increment INCX. */
/*     If N .LE. 0, return with result = 0. */
/*     If N .GE. 1, then INCX must be .GE. 1 */

/*     Four phase method using two built-in constants that are */
/*     hopefully applicable to all machines. */
/*         CUTLO = maximum of  SQRT(U/EPS)  over all known machines. */
/*         CUTHI = minimum of  SQRT(V)      over all known machines. */
/*     where */
/*         EPS = smallest no. such that EPS + 1. .GT. 1. */
/*         U   = smallest positive no.   (underflow limit) */
/*         V   = largest  no.            (overflow  limit) */

/*     Brief outline of algorithm. */

/*     Phase 1 scans zero components. */
/*     move to phase 2 when a component is nonzero and .LE. CUTLO */
/*     move to phase 3 when a component is .GT. CUTLO */
/*     move to phase 4 when a component is .GE. CUTHI/M */
/*     where M = N for X() real and M = 2*N for complex. */

/*     Values for CUTLO and CUTHI. */
/*     From the environmental parameters listed in the IMSL converter */
/*     document the limiting values are as follows: */
/*     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are */
/*                   Univac and DEC at 2**(-103) */
/*                   Thus CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC. */
/*                   Thus CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC. */
/*                   Thus CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/ */
/*     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/ */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DNRM2 */
    /* Parameter adjustments */
    --dx;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  DNRM2 */
    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;

/*                                                 BEGIN MAIN LOOP */

    i__ = 1;
L20:
    switch (next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        PHASE 1.  SUM IS ZERO */

L50:
    if (dx[i__] == zero) {
	goto L200;
    }
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L85;
    }

/*                                PREPARE FOR PHASE 2. */

    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                PREPARE FOR PHASE 4. */

L100:
    i__ = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / dx[i__] / dx[i__];
L105:
    xmax = (d__1 = dx[i__], abs(d__1));
    goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

L70:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L75;
    }

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. */

L110:
    if ((d__1 = dx[i__], abs(d__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i__];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i__], abs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i__] / xmax;
    sum += d__1 * d__1;
    goto L200;

/*                  PREPARE FOR PHASE 3. */

L75:
    sum = sum * xmax * xmax;

/*     FOR REAL OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

L85:
    hitest = cuthi / *n;

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	if ((d__1 = dx[j], abs(d__1)) >= hitest) {
	    goto L100;
	}
/* L95: */
/* Computing 2nd power */
	d__1 = dx[j];
	sum += d__1 * d__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i__ += *incx;
    if (i__ <= nn) {
	goto L20;
    }

/*              END OF MAIN LOOP. */

/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* dnrm2_ */

/* DECK DROTM */
/* Subroutine */ int drotm_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *dparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal w, z__;
    static integer kx, ky;
    static doublereal dh11, dh12, dh22, dh21, dflag;
    static integer nsteps;

/* ***BEGIN PROLOGUE  DROTM */
/* ***PURPOSE  Apply a modified Givens transformation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A8 */
/* ***TYPE      DOUBLE PRECISION (SROTM-S, DROTM-D) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */
/*       DY  double precision vector with N elements */
/*     INCY  storage spacing between elements of DY */
/*   DPARAM  5-element D.P. vector.  DPARAM(1) is DFLAG described below. */
/*           Locations 2-5 of SPARAM contain elements of the */
/*           transformation matrix H described below. */

/*     --Output-- */
/*       DX  rotated vector (unchanged if N .LE. 0) */
/*       DY  rotated vector (unchanged if N .LE. 0) */

/*     Apply the modified Givens transformation, H, to the 2 by N matrix */
/*     (DX**T) */
/*     (DY**T) , where **T indicates transpose.  The elements of DX are */
/*     in DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else */
/*     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY. */

/*     With DPARAM(1)=DFLAG, H has one of the following forms: */

/*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */

/*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */

/*     See DROTMG for a description of data storage in DPARAM. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DROTM */
    /* Parameter adjustments */
    --dparam;
    --dy;
    --dx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DROTM */
    dflag = dparam[1];
    if (*n <= 0 || dflag + two == zero) {
	goto L140;
    }
    if (! (*incx == *incy && *incx > 0)) {
	goto L70;
    }

    nsteps = *n * *incx;
    if (dflag < 0.) {
	goto L50;
    } else if (dflag == 0) {
	goto L10;
    } else {
	goto L30;
    }
L10:
    dh12 = dparam[4];
    dh21 = dparam[3];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = dx[i__];
	z__ = dy[i__];
	dx[i__] = w + z__ * dh12;
	dy[i__] = w * dh21 + z__;
/* L20: */
    }
    goto L140;
L30:
    dh11 = dparam[2];
    dh22 = dparam[5];
    i__2 = nsteps;
    i__1 = *incx;
    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
	w = dx[i__];
	z__ = dy[i__];
	dx[i__] = w * dh11 + z__;
	dy[i__] = -w + dh22 * z__;
/* L40: */
    }
    goto L140;
L50:
    dh11 = dparam[2];
    dh12 = dparam[4];
    dh21 = dparam[3];
    dh22 = dparam[5];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = dx[i__];
	z__ = dy[i__];
	dx[i__] = w * dh11 + z__ * dh12;
	dy[i__] = w * dh21 + z__ * dh22;
/* L60: */
    }
    goto L140;
L70:
    kx = 1;
    ky = 1;
    if (*incx < 0) {
	kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0) {
	ky = (1 - *n) * *incy + 1;
    }

    if (dflag < 0.) {
	goto L120;
    } else if (dflag == 0) {
	goto L80;
    } else {
	goto L100;
    }
L80:
    dh12 = dparam[4];
    dh21 = dparam[3];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = dx[kx];
	z__ = dy[ky];
	dx[kx] = w + z__ * dh12;
	dy[ky] = w * dh21 + z__;
	kx += *incx;
	ky += *incy;
/* L90: */
    }
    goto L140;
L100:
    dh11 = dparam[2];
    dh22 = dparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = dx[kx];
	z__ = dy[ky];
	dx[kx] = w * dh11 + z__;
	dy[ky] = -w + dh22 * z__;
	kx += *incx;
	ky += *incy;
/* L110: */
    }
    goto L140;
L120:
    dh11 = dparam[2];
    dh12 = dparam[4];
    dh21 = dparam[3];
    dh22 = dparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = dx[kx];
	z__ = dy[ky];
	dx[kx] = w * dh11 + z__ * dh12;
	dy[ky] = w * dh21 + z__ * dh22;
	kx += *incx;
	ky += *incy;
/* L130: */
    }
L140:
    return 0;
} /* drotm_ */

/* DECK DROTMG */
/* Subroutine */ int drotmg_(doublereal *dd1, doublereal *dd2, doublereal *
	dx1, doublereal *dy1, doublereal *dparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal gam = 4096.;
    static doublereal gamsq = 16777216.;
    static doublereal rgamsq = 5.9604645e-8;

    /* Format strings */
    static char fmt_120[] = "";
    static char fmt_150[] = "";
    static char fmt_180[] = "";
    static char fmt_210[] = "";

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal du, dp1, dp2, dq1, dq2, dh11, dh12, dh21, dh22;
    static integer igo;
    static doublereal dflag, dtemp;

    /* Assigned format variables */
    static char *igo_fmt;

/* ***BEGIN PROLOGUE  DROTMG */
/* ***PURPOSE  Construct a modified Givens transformation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1B10 */
/* ***TYPE      DOUBLE PRECISION (SROTMG-S, DROTMG-D) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*      DD1  double precision scalar */
/*      DD2  double precision scalar */
/*      DX1  double precision scalar */
/*      DX2  double precision scalar */
/*   DPARAM  D.P. 5-vector. DPARAM(1)=DFLAG defined below. */
/*           Locations 2-5 contain the rotation matrix. */

/*     --Output-- */
/*      DD1  changed to represent the effect of the transformation */
/*      DD2  changed to represent the effect of the transformation */
/*      DX1  changed to represent the effect of the transformation */
/*      DX2  unchanged */

/*     Construct the modified Givens transformation matrix H which zeros */
/*     the second component of the 2-vector  (SQRT(DD1)*DX1,SQRT(DD2)* */
/*     DY2)**T. */
/*     With DPARAM(1)=DFLAG, H has one of the following forms: */

/*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */

/*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */

/*     Locations 2-5 of DPARAM contain DH11, DH21, DH12, and DH22, */
/*     respectively.  (Values of 1.D0, -1.D0, or 0.D0 implied by the */
/*     value of DPARAM(1) are not stored in DPARAM.) */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920316  Prologue corrected.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DROTMG */
    /* Parameter adjustments */
    --dparam;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DROTMG */
    if (! (*dd1 < zero)) {
	goto L10;
    }
/*       GO ZERO-H-D-AND-DX1.. */
    goto L60;
L10:
/*     CASE-DD1-NONNEGATIVE */
    dp2 = *dd2 * *dy1;
    if (! (dp2 == zero)) {
	goto L20;
    }
    dflag = -two;
    goto L260;
/*     REGULAR-CASE.. */
L20:
    dp1 = *dd1 * *dx1;
    dq2 = dp2 * *dy1;
    dq1 = dp1 * *dx1;

    if (! (abs(dq1) > abs(dq2))) {
	goto L40;
    }
    dh21 = -(*dy1) / *dx1;
    dh12 = dp2 / dp1;

    du = one - dh12 * dh21;

    if (! (du <= zero)) {
	goto L30;
    }
/*         GO ZERO-H-D-AND-DX1.. */
    goto L60;
L30:
    dflag = zero;
    *dd1 /= du;
    *dd2 /= du;
    *dx1 *= du;
/*         GO SCALE-CHECK.. */
    goto L100;
L40:
    if (! (dq2 < zero)) {
	goto L50;
    }
/*         GO ZERO-H-D-AND-DX1.. */
    goto L60;
L50:
    dflag = one;
    dh11 = dp1 / dp2;
    dh22 = *dx1 / *dy1;
    du = one + dh11 * dh22;
    dtemp = *dd2 / du;
    *dd2 = *dd1 / du;
    *dd1 = dtemp;
    *dx1 = *dy1 * du;
/*         GO SCALE-CHECK */
    goto L100;
/*     PROCEDURE..ZERO-H-D-AND-DX1.. */
L60:
    dflag = -one;
    dh11 = zero;
    dh12 = zero;
    dh21 = zero;
    dh22 = zero;

    *dd1 = zero;
    *dd2 = zero;
    *dx1 = zero;
/*         RETURN.. */
    goto L220;
/*     PROCEDURE..FIX-H.. */
L70:
    if (! (dflag >= zero)) {
	goto L90;
    }

    if (! (dflag == zero)) {
	goto L80;
    }
    dh11 = one;
    dh22 = one;
    dflag = -one;
    goto L90;
L80:
    dh21 = -one;
    dh12 = one;
    dflag = -one;
L90:
    switch (igo) {
	case 0: goto L120;
	case 1: goto L150;
	case 2: goto L180;
	case 3: goto L210;
    }
/*     PROCEDURE..SCALE-CHECK */
L100:
L110:
    if (! (*dd1 <= rgamsq)) {
	goto L130;
    }
    if (*dd1 == zero) {
	goto L160;
    }
    igo = 0;
    igo_fmt = fmt_120;
/*              FIX-H.. */
    goto L70;
L120:
/* Computing 2nd power */
    d__1 = gam;
    *dd1 *= d__1 * d__1;
    *dx1 /= gam;
    dh11 /= gam;
    dh12 /= gam;
    goto L110;
L130:
L140:
    if (! (*dd1 >= gamsq)) {
	goto L160;
    }
    igo = 1;
    igo_fmt = fmt_150;
/*              FIX-H.. */
    goto L70;
L150:
/* Computing 2nd power */
    d__1 = gam;
    *dd1 /= d__1 * d__1;
    *dx1 *= gam;
    dh11 *= gam;
    dh12 *= gam;
    goto L140;
L160:
L170:
    if (! (abs(*dd2) <= rgamsq)) {
	goto L190;
    }
    if (*dd2 == zero) {
	goto L220;
    }
    igo = 2;
    igo_fmt = fmt_180;
/*              FIX-H.. */
    goto L70;
L180:
/* Computing 2nd power */
    d__1 = gam;
    *dd2 *= d__1 * d__1;
    dh21 /= gam;
    dh22 /= gam;
    goto L170;
L190:
L200:
    if (! (abs(*dd2) >= gamsq)) {
	goto L220;
    }
    igo = 3;
    igo_fmt = fmt_210;
/*              FIX-H.. */
    goto L70;
L210:
/* Computing 2nd power */
    d__1 = gam;
    *dd2 /= d__1 * d__1;
    dh21 *= gam;
    dh22 *= gam;
    goto L200;
L220:
    if (dflag < 0.) {
	goto L250;
    } else if (dflag == 0) {
	goto L230;
    } else {
	goto L240;
    }
L230:
    dparam[3] = dh21;
    dparam[4] = dh12;
    goto L260;
L240:
    dparam[2] = dh11;
    dparam[5] = dh22;
    goto L260;
L250:
    dparam[2] = dh11;
    dparam[3] = dh21;
    dparam[4] = dh12;
    dparam[5] = dh22;
L260:
    dparam[1] = dflag;
    return 0;
} /* drotmg_ */

/* DECK DSCAL */
/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, mp1;

/* ***BEGIN PROLOGUE  DSCAL */
/* ***PURPOSE  Multiply a vector by a constant. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A6 */
/* ***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DA  double precision scale factor */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */

/*     --Output-- */
/*       DX  double precision result (unchanged if N.LE.0) */

/*     Replace double precision DX by double precision DA*DX. */
/*     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX), */
/*     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900821  Modified to correct problem with a negative increment. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DSCAL */
/* ***FIRST EXECUTABLE STATEMENT  DSCAL */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*     Code for increment not equal to 1. */

    ix = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dx[ix] = *da * dx[ix];
	ix += *incx;
/* L10: */
    }
    return 0;

/*     Code for increment equal to 1. */

/*     Clean-up loop so remaining vector length is a multiple of 5. */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dx[i__] = *da * dx[i__];
	dx[i__ + 1] = *da * dx[i__ + 1];
	dx[i__ + 2] = *da * dx[i__ + 2];
	dx[i__ + 3] = *da * dx[i__ + 3];
	dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* DECK DSWAP */
/* Subroutine */ int dswap_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, ix, iy, ns, mp1;
    static doublereal dtemp1, dtemp2, dtemp3;

/* ***BEGIN PROLOGUE  DSWAP */
/* ***PURPOSE  Interchange two vectors. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A5 */
/* ***TYPE      DOUBLE PRECISION (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I) */
/* ***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */
/*       DY  double precision vector with N elements */
/*     INCY  storage spacing between elements of DY */

/*     --Output-- */
/*       DX  input vector DY (unchanged if N .LE. 0) */
/*       DY  input vector DX (unchanged if N .LE. 0) */

/*     Interchange double precision DX and double precision DY. */
/*     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY), */
/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DSWAP */
/* ***FIRST EXECUTABLE STATEMENT  DSWAP */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == *incy) {
	if ((i__1 = *incx - 1) < 0) {
	    goto L5;
	} else if (i__1 == 0) {
	    goto L20;
	} else {
	    goto L60;
	}
    }

/*     Code for unequal or nonpositive increments. */

L5:
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp1 = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp1;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*     Code for both increments equal to 1. */

/*     Clean-up loop so remaining vector length is a multiple of 3. */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp1 = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp1;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
	dtemp1 = dx[i__];
	dtemp2 = dx[i__ + 1];
	dtemp3 = dx[i__ + 2];
	dx[i__] = dy[i__];
	dx[i__ + 1] = dy[i__ + 1];
	dx[i__ + 2] = dy[i__ + 2];
	dy[i__] = dtemp1;
	dy[i__ + 1] = dtemp2;
	dy[i__ + 2] = dtemp3;
/* L50: */
    }
    return 0;

/*     Code for equal, positive, non-unit increments. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dtemp1 = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp1;
/* L70: */
    }
    return 0;
} /* dswap_ */

/* DECK DWNLIT */
/* Subroutine */ int dwnlit_(doublereal *w, integer *mdw, integer *m, integer 
	*n, integer *l, integer *ipivot, integer *itype, doublereal *h__, 
	doublereal *scale, doublereal *rnorm, integer *idope, doublereal *
	dope, logical *done)
{
    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer i1, j1, l1, lb, me, jj, jp, ir;
    static doublereal rn;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static doublereal tau;
    static integer niv;
    static doublereal hbar;
    static integer lend, mend;
    static doublereal amax;
    static integer imax;
    static doublereal alsq;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical indep;
    static integer krank;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), drotm_(integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer nsoln;
    extern /* Subroutine */ int dwnlt1_(integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern logical dwnlt2_(integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int dwnlt3_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *);
    static logical recalc;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal factor, eanorm, sparam[5];
    extern /* Subroutine */ int drotmg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

/* ***BEGIN PROLOGUE  DWNLIT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DWNNLS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLIT-S, DWNLIT-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to DWNNLS( ). */
/*     The documentation for DWNNLS( ) has complete usage instructions. */

/*     Note  The M by (N+1) matrix W( , ) contains the rt. hand side */
/*           B as the (N+1)st col. */

/*     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with */
/*     col interchanges. */

/* ***SEE ALSO  DWNNLS */
/* ***ROUTINES CALLED  DCOPY, DH12, DROTM, DROTMG, DSCAL, DSWAP, DWNLT1, */
/*                    DWNLT2, DWNLT3, IDAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and revised.  (WRB & RWC) */
/*   890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900604  DP version created from SP version. .  (RWC) */
/* ***END PROLOGUE  DWNLIT */



/* ***FIRST EXECUTABLE STATEMENT  DWNLIT */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --ipivot;
    --itype;
    --h__;
    --scale;
    --idope;
    --dope;

    /* Function Body */
    me = idope[1];
    nsoln = idope[2];
    l1 = idope[3];

    alsq = dope[1];
    eanorm = dope[2];
    tau = dope[3];

/* Computing MIN */
    i__1 = *m - 1;
    lb = min(i__1,*l);
    recalc = TRUE_;
    *rnorm = 0.;
    krank = 0;

/*     We set FACTOR=1.0 so that the heavy weight ALAMDA will be */
/*     included in the test for column independence. */

    factor = 1.;
    lend = *l;
    i__1 = lb;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Set IR to point to the I-th row. */

	ir = i__;
	mend = *m;
	dwnlt1_(&i__, &lend, m, &ir, mdw, &recalc, &imax, &hbar, &h__[1], &
		scale[1], &w[w_offset]);

/*        Update column SS and find pivot column. */

	dwnlt3_(&i__, &imax, m, mdw, &ipivot[1], &h__[1], &w[w_offset]);

/*        Perform column interchange. */
/*        Test independence of incoming column. */

L130:
	if (dwnlt2_(&me, &mend, &ir, &factor, &tau, &scale[1], &w[i__ * 
		w_dim1 + 1])) {

/*           Eliminate I-th column below diagonal using modified Givens */
/*           transformations applied to (A B). */

/*           When operating near the ME line, use the largest element */
/*           above it as the pivot. */

	    i__2 = i__ + 1;
	    for (j = *m; j >= i__2; --j) {
		jp = j - 1;
		if (j == me + 1) {
		    imax = me;
/* Computing 2nd power */
		    d__1 = w[me + i__ * w_dim1];
		    amax = scale[me] * (d__1 * d__1);
		    i__3 = i__;
		    for (jp = j - 1; jp >= i__3; --jp) {
/* Computing 2nd power */
			d__1 = w[jp + i__ * w_dim1];
			t = scale[jp] * (d__1 * d__1);
			if (t > amax) {
			    imax = jp;
			    amax = t;
			}
/* L150: */
		    }
		    jp = imax;
		}

		if (w[j + i__ * w_dim1] != 0.) {
		    drotmg_(&scale[jp], &scale[j], &w[jp + i__ * w_dim1], &w[
			    j + i__ * w_dim1], sparam);
		    w[j + i__ * w_dim1] = 0.;
		    i__3 = *n + 1 - i__;
		    drotm_(&i__3, &w[jp + (i__ + 1) * w_dim1], mdw, &w[j + (
			    i__ + 1) * w_dim1], mdw, sparam);
		}
/* L160: */
	    }
	} else if (lend > i__) {

/*           Column I is dependent.  Swap with column LEND. */
/*           Perform column interchange, */
/*           and find column in remaining set with largest SS. */

	    dwnlt3_(&i__, &lend, m, mdw, &ipivot[1], &h__[1], &w[w_offset]);
	    --lend;
	    i__2 = lend - i__ + 1;
	    imax = idamax_(&i__2, &h__[i__], &c__1) + i__ - 1;
	    hbar = h__[imax];
	    goto L130;
	} else {
	    krank = i__ - 1;
	    goto L190;
	}
/* L180: */
    }
    krank = l1;

L190:
    if (krank < me) {
	factor = alsq;
	i__1 = me;
	for (i__ = krank + 1; i__ <= i__1; ++i__) {
	    dcopy_(l, &c_b95, &c__0, &w[i__ + w_dim1], mdw);
/* L200: */
	}

/*        Determine the rank of the remaining equality constraint */
/*        equations by eliminating within the block of constrained */
/*        variables.  Remove any redundant constraints. */

	recalc = TRUE_;
/* Computing MIN */
	i__1 = *l + me - krank;
	lb = min(i__1,*n);
	i__1 = lb;
	for (i__ = *l + 1; i__ <= i__1; ++i__) {
	    ir = krank + i__ - *l;
	    lend = *n;
	    mend = me;
	    dwnlt1_(&i__, &lend, &me, &ir, mdw, &recalc, &imax, &hbar, &h__[1]
		    , &scale[1], &w[w_offset]);

/*           Update col ss and find pivot col */

	    dwnlt3_(&i__, &imax, m, mdw, &ipivot[1], &h__[1], &w[w_offset]);

/*           Perform column interchange */
/*           Eliminate elements in the I-th col. */

	    i__2 = ir + 1;
	    for (j = me; j >= i__2; --j) {
		if (w[j + i__ * w_dim1] != 0.) {
		    drotmg_(&scale[j - 1], &scale[j], &w[j - 1 + i__ * w_dim1]
			    , &w[j + i__ * w_dim1], sparam);
		    w[j + i__ * w_dim1] = 0.;
		    i__3 = *n + 1 - i__;
		    drotm_(&i__3, &w[j - 1 + (i__ + 1) * w_dim1], mdw, &w[j + 
			    (i__ + 1) * w_dim1], mdw, sparam);
		}
/* L240: */
	    }

/*           I=column being eliminated. */
/*           Test independence of incoming column. */
/*           Remove any redundant or dependent equality constraints. */

	    if (! dwnlt2_(&me, &mend, &ir, &factor, &tau, &scale[1], &w[i__ * 
		    w_dim1 + 1])) {
		jj = ir;
		i__2 = me;
		for (ir = jj; ir <= i__2; ++ir) {
		    dcopy_(n, &c_b95, &c__0, &w[ir + w_dim1], mdw);
		    *rnorm += scale[ir] * w[ir + (*n + 1) * w_dim1] / alsq * 
			    w[ir + (*n + 1) * w_dim1];
		    w[ir + (*n + 1) * w_dim1] = 0.;
		    scale[ir] = 1.;

/*                 Reclassify the zeroed row as a least squares equation. */

		    itype[ir] = 1;
/* L260: */
		}

/*              Reduce ME to reflect any discovered dependent equality */
/*              constraints. */

		me = jj - 1;
		goto L280;
	    }
/* L270: */
	}
    }

/*     Try to determine the variables KRANK+1 through L1 from the */
/*     least squares equations.  Continue the triangularization with */
/*     pivot element W(ME+1,I). */

L280:
    if (krank < l1) {
	recalc = TRUE_;

/*        Set FACTOR=ALSQ to remove effect of heavy weight from */
/*        test for column independence. */

	factor = alsq;
	i__1 = l1;
	for (i__ = krank + 1; i__ <= i__1; ++i__) {

/*           Set IR to point to the ME+1-st row. */

	    ir = me + 1;
	    lend = *l;
	    mend = *m;
	    dwnlt1_(&i__, l, m, &ir, mdw, &recalc, &imax, &hbar, &h__[1], &
		    scale[1], &w[w_offset]);

/*           Update column SS and find pivot column. */

	    dwnlt3_(&i__, &imax, m, mdw, &ipivot[1], &h__[1], &w[w_offset]);

/*           Perform column interchange. */
/*           Eliminate I-th column below the IR-th element. */

	    i__2 = ir + 1;
	    for (j = *m; j >= i__2; --j) {
		if (w[j + i__ * w_dim1] != 0.) {
		    drotmg_(&scale[j - 1], &scale[j], &w[j - 1 + i__ * w_dim1]
			    , &w[j + i__ * w_dim1], sparam);
		    w[j + i__ * w_dim1] = 0.;
		    i__3 = *n + 1 - i__;
		    drotm_(&i__3, &w[j - 1 + (i__ + 1) * w_dim1], mdw, &w[j + 
			    (i__ + 1) * w_dim1], mdw, sparam);
		}
/* L320: */
	    }

/*           Test if new pivot element is near zero. */
/*           If so, the column is dependent. */
/*           Then check row norm test to be classified as independent. */

/* Computing 2nd power */
	    d__1 = w[ir + i__ * w_dim1];
	    t = scale[ir] * (d__1 * d__1);
/* Computing 2nd power */
	    d__1 = tau * eanorm;
	    indep = t > d__1 * d__1;
	    if (indep) {
		rn = 0.;
		i__2 = *m;
		for (i1 = ir; i1 <= i__2; ++i1) {
		    i__3 = *n;
		    for (j1 = i__ + 1; j1 <= i__3; ++j1) {
/* Computing MAX */
/* Computing 2nd power */
			d__3 = w[i1 + j1 * w_dim1];
			d__1 = rn, d__2 = scale[i1] * (d__3 * d__3);
			rn = max(d__1,d__2);
/* L330: */
		    }
/* L340: */
		}
/* Computing 2nd power */
		d__1 = tau;
		indep = t > rn * (d__1 * d__1);
	    }

/*           If independent, swap the IR-th and KRANK+1-th rows to */
/*           maintain the triangular form.  Update the rank indicator */
/*           KRANK and the equality constraint pointer ME. */

	    if (! indep) {
		goto L360;
	    }
	    i__2 = *n + 1;
	    dswap_(&i__2, &w[krank + 1 + w_dim1], mdw, &w[ir + w_dim1], mdw);
	    dswap_(&c__1, &scale[krank + 1], &c__1, &scale[ir], &c__1);

/*           Reclassify the least square equation as an equality */
/*           constraint and rescale it. */

	    itype[ir] = 0;
	    t = sqrt(scale[krank + 1]);
	    i__2 = *n + 1;
	    dscal_(&i__2, &t, &w[krank + 1 + w_dim1], mdw);
	    scale[krank + 1] = alsq;
	    ++me;
	    ++krank;
/* L350: */
	}
    }

/*     If pseudorank is less than L, apply Householder transformation. */
/*     from right. */

L360:
    if (krank < *l) {
	for (j = krank; j >= 1; --j) {
	    i__1 = krank + 1;
	    i__2 = j - 1;
	    dh12_(&c__1, &j, &i__1, l, &w[j + w_dim1], mdw, &h__[j], &w[
		    w_offset], mdw, &c__1, &i__2);
/* L370: */
	}
    }

    niv = krank + nsoln - *l;
    if (*l == *n) {
	*done = TRUE_;
    }

/*     End of initial triangularization. */

    idope[1] = me;
    idope[2] = krank;
    idope[3] = niv;
    return 0;
} /* dwnlit_ */

/* DECK DWNLSM */
/* Subroutine */ int dwnlsm_(doublereal *w, integer *mdw, integer *mme, 
	integer *ma, integer *n, integer *l, doublereal *prgopt, doublereal *
	x, doublereal *rnorm, integer *mode, integer *ipivot, integer *itype, 
	doublereal *wd, doublereal *h__, doublereal *scale, doublereal *z__, 
	doublereal *temp, doublereal *d__)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, m;
    static doublereal t;
    static integer l1;
    static doublereal z2;
    static integer me, jp;
    static doublereal sm, zz, fac;
    extern /* Subroutine */ int dh12_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *);
    static integer key;
    static doublereal tau;
    static integer niv;
    static logical pos, done;
    static doublereal amax, dope[3];
    static integer jcon, link, imax;
    static doublereal alsq;
    static integer iter, last, isol;
    static doublereal wmax;
    static integer next, nopt;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer idope[3], krank;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static integer nlink;
    static doublereal bnorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer itemp, itmax;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), drotm_(integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *);
    static integer iwmax, nsoln;
    extern doublereal d1mach_(integer *);
    static doublereal alamda;
    static logical feasbl;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal eanorm, sparam[5];
    static logical hitcon;
    static doublereal drelpr;
    extern /* Subroutine */ int drotmg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer ntimes;
    extern /* Subroutine */ int dwnlit_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, logical *), 
	    xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    static doublereal blowup;

/* ***BEGIN PROLOGUE  DWNLSM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DWNNLS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLSM-S, DWNLSM-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     This is a companion subprogram to DWNNLS. */
/*     The documentation for DWNNLS has complete usage instructions. */

/*     In addition to the parameters discussed in the prologue to */
/*     subroutine DWNNLS, the following work arrays are used in */
/*     subroutine DWNLSM  (they are passed through the calling */
/*     sequence from DWNNLS for purposes of variable dimensioning). */
/*     Their contents will in general be of no interest to the user. */

/*     Variables of type REAL are DOUBLE PRECISION. */

/*         IPIVOT(*) */
/*            An array of length N.  Upon completion it contains the */
/*         pivoting information for the cols of W(*,*). */

/*         ITYPE(*) */
/*            An array of length M which is used to keep track */
/*         of the classification of the equations.  ITYPE(I)=0 */
/*         denotes equation I as an equality constraint. */
/*         ITYPE(I)=1 denotes equation I as a least squares */
/*         equation. */

/*         WD(*) */
/*            An array of length N.  Upon completion it contains the */
/*         dual solution vector. */

/*         H(*) */
/*            An array of length N.  Upon completion it contains the */
/*         pivot scalars of the Householder transformations performed */
/*         in the case KRANK.LT.L. */

/*         SCALE(*) */
/*            An array of length M which is used by the subroutine */
/*         to store the diagonal matrix of weights. */
/*         These are used to apply the modified Givens */
/*         transformations. */

/*         Z(*),TEMP(*) */
/*            Working arrays of length N. */

/*         D(*) */
/*            An array of length N that contains the */
/*         column scaling for the matrix (E). */
/*                                       (A) */

/* ***SEE ALSO  DWNNLS */
/* ***ROUTINES CALLED  D1MACH, DASUM, DAXPY, DCOPY, DH12, DNRM2, DROTM, */
/*                    DROTMG, DSCAL, DSWAP, DWNLIT, IDAMAX, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and revised.  (WRB & RWC) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900328  Added TYPE section.  (WRB) */
/*   900510  Fixed an error message.  (RWC) */
/*   900604  DP version created from SP version.  (RWC) */
/*   900911  Restriction on value of ALAMDA included.  (WRB) */
/* ***END PROLOGUE  DWNLSM */



    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --prgopt;
    --x;
    --ipivot;
    --itype;
    --wd;
    --h__;
    --scale;
    --z__;
    --temp;
    --d__;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DWNLSM */

/*     Initialize variables. */
/*     DRELPR is the precision for the particular machine */
/*     being used.  This logic avoids resetting it every entry. */

    if (first) {
	drelpr = d1mach_(&c__4);
    }
    first = FALSE_;

/*     Set the nominal tolerance used in the code. */

    tau = sqrt(drelpr);

    m = *ma + *mme;
    me = *mme;
    *mode = 2;

/*     To process option vector */

    fac = 1e-4;

/*     Set the nominal blow up factor used in the code. */

    blowup = tau;

/*     The nominal column scaling used in the code is */
/*     the identity scaling. */

    dcopy_(n, &c_b41, &c__0, &d__[1], &c__1);

/*     Define bound for number of options to change. */

    nopt = 1000;

/*     Define bound for positive value of LINK. */

    nlink = 100000;
    ntimes = 0;
    last = 1;
    link = (integer) prgopt[1];
    if (link <= 0 || link > nlink) {
	xermsg_("SLATEC", "DWNLSM", "IN DWNNLS, THE OPTION VECTOR IS UNDEFIN"
		"ED", &c__3, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)41);
	return 0;
    }

L100:
    if (link > 1) {
	++ntimes;
	if (ntimes > nopt) {
	    xermsg_("SLATEC", "DWNLSM", "IN DWNNLS, THE LINKS IN THE OPTION "
		    "VECTOR ARE CYCLING.", &c__3, &c__1, (ftnlen)6, (ftnlen)6, 
		    (ftnlen)54);
	    return 0;
	}

	key = (integer) prgopt[last + 1];
	if (key == 6 && prgopt[last + 2] != 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		t = dnrm2_(&m, &w[j * w_dim1 + 1], &c__1);
		if (t != 0.) {
		    t = 1. / t;
		}
		d__[j] = t;
/* L110: */
	    }
	}

	if (key == 7) {
	    dcopy_(n, &prgopt[last + 2], &c__1, &d__[1], &c__1);
	}
	if (key == 8) {
/* Computing MAX */
	    d__1 = drelpr, d__2 = prgopt[last + 2];
	    tau = max(d__1,d__2);
	}
	if (key == 9) {
/* Computing MAX */
	    d__1 = drelpr, d__2 = prgopt[last + 2];
	    blowup = max(d__1,d__2);
	}

	next = (integer) prgopt[link];
	if (next <= 0 || next > nlink) {
	    xermsg_("SLATEC", "DWNLSM", "IN DWNNLS, THE OPTION VECTOR IS UND"
		    "EFINED", &c__3, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)41);
	    return 0;
	}

	last = link;
	link = next;
	goto L100;
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dscal_(&m, &d__[j], &w[j * w_dim1 + 1], &c__1);
/* L120: */
    }

/*     Process option vector */

    done = FALSE_;
    iter = 0;
    itmax = (*n - *l) * 3;
    *mode = 0;
    nsoln = *l;
    l1 = min(m,*l);

/*     Compute scale factor to apply to equality constraint equations. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wd[j] = dasum_(&m, &w[j * w_dim1 + 1], &c__1);
/* L130: */
    }

    imax = idamax_(n, &wd[1], &c__1);
    eanorm = wd[imax];
    bnorm = dasum_(&m, &w[(*n + 1) * w_dim1 + 1], &c__1);
    alamda = eanorm / (drelpr * fac);

/*     On machines, such as the VAXes using D floating, with a very */
/*     limited exponent range for double precision values, the previously */
/*     computed value of ALAMDA may cause an overflow condition. */
/*     Therefore, this code further limits the value of ALAMDA. */

/* Computing MIN */
    d__1 = alamda, d__2 = sqrt(d1mach_(&c__2));
    alamda = min(d__1,d__2);

/*     Define scaling diagonal matrix for modified Givens usage and */
/*     classify equation types. */

/* Computing 2nd power */
    d__1 = alamda;
    alsq = d__1 * d__1;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        When equation I is heavily weighted ITYPE(I)=0, */
/*        else ITYPE(I)=1. */

	if (i__ <= me) {
	    t = alsq;
	    itemp = 0;
	} else {
	    t = 1.;
	    itemp = 1;
	}
	scale[i__] = t;
	itype[i__] = itemp;
/* L140: */
    }

/*     Set the solution vector X(*) to zero and the column interchange */
/*     matrix to the identity. */

    dcopy_(n, &c_b95, &c__0, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipivot[i__] = i__;
/* L150: */
    }

/*     Perform initial triangularization in the submatrix */
/*     corresponding to the unconstrained variables. */
/*     Set first L components of dual vector to zero because */
/*     these correspond to the unconstrained variables. */

    dcopy_(l, &c_b95, &c__0, &wd[1], &c__1);

/*     The arrays IDOPE(*) and DOPE(*) are used to pass */
/*     information to DWNLIT().  This was done to avoid */
/*     a long calling sequence or the use of COMMON. */

    idope[0] = me;
    idope[1] = nsoln;
    idope[2] = l1;

    dope[0] = alsq;
    dope[1] = eanorm;
    dope[2] = tau;
    dwnlit_(&w[w_offset], mdw, &m, n, l, &ipivot[1], &itype[1], &h__[1], &
	    scale[1], rnorm, idope, dope, &done);
    me = idope[0];
    krank = idope[1];
    niv = idope[2];

/*     Perform WNNLS algorithm using the following steps. */

/*     Until(DONE) */
/*        compute search direction and feasible point */
/*        when (HITCON) add constraints */
/*        else perform multiplier test and drop a constraint */
/*        fin */
/*     Compute-Final-Solution */

/*     To compute search direction and feasible point, */
/*     solve the triangular system of currently non-active */
/*     variables and store the solution in Z(*). */

/*     To solve system */
/*     Copy right hand side into TEMP vector to use overwriting method. */

L160:
    if (done) {
	goto L330;
    }
    isol = *l + 1;
    if (nsoln >= isol) {
	dcopy_(&niv, &w[(*n + 1) * w_dim1 + 1], &c__1, &temp[1], &c__1);
	i__1 = isol;
	for (j = nsoln; j >= i__1; --j) {
	    if (j > krank) {
		i__ = niv - nsoln + j;
	    } else {
		i__ = j;
	    }

	    if (j > krank && j <= *l) {
		z__[j] = 0.;
	    } else {
		z__[j] = temp[i__] / w[i__ + j * w_dim1];
		i__2 = i__ - 1;
		d__1 = -z__[j];
		daxpy_(&i__2, &d__1, &w[j * w_dim1 + 1], &c__1, &temp[1], &
			c__1);
	    }
/* L170: */
	}
    }

/*     Increment iteration counter and check against maximum number */
/*     of iterations. */

    ++iter;
    if (iter > itmax) {
	*mode = 1;
	done = TRUE_;
    }

/*     Check to see if any constraints have become active. */
/*     If so, calculate an interpolation factor so that all */
/*     active constraints are removed from the basis. */

    alpha = 2.;
    hitcon = FALSE_;
    i__1 = nsoln;
    for (j = *l + 1; j <= i__1; ++j) {
	zz = z__[j];
	if (zz <= 0.) {
	    t = x[j] / (x[j] - zz);
	    if (t < alpha) {
		alpha = t;
		jcon = j;
	    }
	    hitcon = TRUE_;
	}
/* L180: */
    }

/*     Compute search direction and feasible point */

    if (hitcon) {

/*        To add constraints, use computed ALPHA to interpolate between */
/*        last feasible solution X(*) and current unconstrained (and */
/*        infeasible) solution Z(*). */

	i__1 = nsoln;
	for (j = *l + 1; j <= i__1; ++j) {
	    x[j] += alpha * (z__[j] - x[j]);
/* L190: */
	}
	feasbl = FALSE_;

/*        Remove column JCON and shift columns JCON+1 through N to the */
/*        left.  Swap column JCON into the N th position.  This achieves */
/*        upper Hessenberg form for the nonactive constraints and */
/*        leaves an upper Hessenberg matrix to retriangularize. */

L200:
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    t = w[i__ + jcon * w_dim1];
	    i__2 = *n - jcon;
	    dcopy_(&i__2, &w[i__ + (jcon + 1) * w_dim1], mdw, &w[i__ + jcon * 
		    w_dim1], mdw);
	    w[i__ + *n * w_dim1] = t;
/* L210: */
	}

/*        Update permuted index vector to reflect this shift and swap. */

	itemp = ipivot[jcon];
	i__1 = *n - 1;
	for (i__ = jcon; i__ <= i__1; ++i__) {
	    ipivot[i__] = ipivot[i__ + 1];
/* L220: */
	}
	ipivot[*n] = itemp;

/*        Similarly permute X(*) vector. */

	i__1 = *n - jcon;
	dcopy_(&i__1, &x[jcon + 1], &c__1, &x[jcon], &c__1);
	x[*n] = 0.;
	--nsoln;
	--niv;

/*        Retriangularize upper Hessenberg matrix after adding */
/*        constraints. */

	i__ = krank + jcon - *l;
	i__1 = nsoln;
	for (j = jcon; j <= i__1; ++j) {
	    if (itype[i__] == 0 && itype[i__ + 1] == 0) {

/*              Zero IP1 to I in column J */

		if (w[i__ + 1 + j * w_dim1] != 0.) {
		    drotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * w_dim1]
			    , &w[i__ + 1 + j * w_dim1], sparam);
		    w[i__ + 1 + j * w_dim1] = 0.;
		    i__2 = *n + 1 - j;
		    drotm_(&i__2, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ + 1 
			    + (j + 1) * w_dim1], mdw, sparam);
		}
	    } else if (itype[i__] == 1 && itype[i__ + 1] == 1) {

/*              Zero IP1 to I in column J */

		if (w[i__ + 1 + j * w_dim1] != 0.) {
		    drotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * w_dim1]
			    , &w[i__ + 1 + j * w_dim1], sparam);
		    w[i__ + 1 + j * w_dim1] = 0.;
		    i__2 = *n + 1 - j;
		    drotm_(&i__2, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ + 1 
			    + (j + 1) * w_dim1], mdw, sparam);
		}
	    } else if (itype[i__] == 1 && itype[i__ + 1] == 0) {
		i__2 = *n + 1;
		dswap_(&i__2, &w[i__ + w_dim1], mdw, &w[i__ + 1 + w_dim1], 
			mdw);
		dswap_(&c__1, &scale[i__], &c__1, &scale[i__ + 1], &c__1);
		itemp = itype[i__ + 1];
		itype[i__ + 1] = itype[i__];
		itype[i__] = itemp;

/*              Swapped row was formerly a pivot element, so it will */
/*              be large enough to perform elimination. */
/*              Zero IP1 to I in column J. */

		if (w[i__ + 1 + j * w_dim1] != 0.) {
		    drotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * w_dim1]
			    , &w[i__ + 1 + j * w_dim1], sparam);
		    w[i__ + 1 + j * w_dim1] = 0.;
		    i__2 = *n + 1 - j;
		    drotm_(&i__2, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ + 1 
			    + (j + 1) * w_dim1], mdw, sparam);
		}
	    } else if (itype[i__] == 0 && itype[i__ + 1] == 1) {
/* Computing 2nd power */
		d__1 = w[i__ + j * w_dim1];
/* Computing 2nd power */
		d__2 = tau * eanorm;
		if (scale[i__] * (d__1 * d__1) / alsq > d__2 * d__2) {

/*                 Zero IP1 to I in column J */

		    if (w[i__ + 1 + j * w_dim1] != 0.) {
			drotmg_(&scale[i__], &scale[i__ + 1], &w[i__ + j * 
				w_dim1], &w[i__ + 1 + j * w_dim1], sparam);
			w[i__ + 1 + j * w_dim1] = 0.;
			i__2 = *n + 1 - j;
			drotm_(&i__2, &w[i__ + (j + 1) * w_dim1], mdw, &w[i__ 
				+ 1 + (j + 1) * w_dim1], mdw, sparam);
		    }
		} else {
		    i__2 = *n + 1;
		    dswap_(&i__2, &w[i__ + w_dim1], mdw, &w[i__ + 1 + w_dim1],
			     mdw);
		    dswap_(&c__1, &scale[i__], &c__1, &scale[i__ + 1], &c__1);
		    itemp = itype[i__ + 1];
		    itype[i__ + 1] = itype[i__];
		    itype[i__] = itemp;
		    w[i__ + 1 + j * w_dim1] = 0.;
		}
	    }
	    ++i__;
/* L230: */
	}

/*        See if the remaining coefficients in the solution set are */
/*        feasible.  They should be because of the way ALPHA was */
/*        determined.  If any are infeasible, it is due to roundoff */
/*        error.  Any that are non-positive will be set to zero and */
/*        removed from the solution set. */

	i__1 = nsoln;
	for (jcon = *l + 1; jcon <= i__1; ++jcon) {
	    if (x[jcon] <= 0.) {
		goto L250;
	    }
/* L240: */
	}
	feasbl = TRUE_;
L250:
	if (! feasbl) {
	    goto L200;
	}
    } else {

/*        To perform multiplier test and drop a constraint. */

	dcopy_(&nsoln, &z__[1], &c__1, &x[1], &c__1);
	if (nsoln < *n) {
	    i__1 = *n - nsoln;
	    dcopy_(&i__1, &c_b95, &c__0, &x[nsoln + 1], &c__1);
	}

/*        Reclassify least squares equations as equalities as necessary. */

	i__ = niv + 1;
L260:
	if (i__ <= me) {
	    if (itype[i__] == 0) {
		++i__;
	    } else {
		i__1 = *n + 1;
		dswap_(&i__1, &w[i__ + w_dim1], mdw, &w[me + w_dim1], mdw);
		dswap_(&c__1, &scale[i__], &c__1, &scale[me], &c__1);
		itemp = itype[i__];
		itype[i__] = itype[me];
		itype[me] = itemp;
		--me;
	    }
	    goto L260;
	}

/*        Form inner product vector WD(*) of dual coefficients. */

	i__1 = *n;
	for (j = nsoln + 1; j <= i__1; ++j) {
	    sm = 0.;
	    i__2 = m;
	    for (i__ = nsoln + 1; i__ <= i__2; ++i__) {
		sm += scale[i__] * w[i__ + j * w_dim1] * w[i__ + (*n + 1) * 
			w_dim1];
/* L270: */
	    }
	    wd[j] = sm;
/* L280: */
	}

/*        Find J such that WD(J)=WMAX is maximum.  This determines */
/*        that the incoming column J will reduce the residual vector */
/*        and be positive. */

L290:
	wmax = 0.;
	iwmax = nsoln + 1;
	i__1 = *n;
	for (j = nsoln + 1; j <= i__1; ++j) {
	    if (wd[j] > wmax) {
		wmax = wd[j];
		iwmax = j;
	    }
/* L300: */
	}
	if (wmax <= 0.) {
	    goto L330;
	}

/*        Set dual coefficients to zero for incoming column. */

	wd[iwmax] = 0.;

/*        WMAX .GT. 0.D0, so okay to move column IWMAX to solution set. */
/*        Perform transformation to retriangularize, and test for near */
/*        linear dependence. */

/*        Swap column IWMAX into NSOLN-th position to maintain upper */
/*        Hessenberg form of adjacent columns, and add new column to */
/*        triangular decomposition. */

	++nsoln;
	++niv;
	if (nsoln != iwmax) {
	    dswap_(&m, &w[nsoln * w_dim1 + 1], &c__1, &w[iwmax * w_dim1 + 1], 
		    &c__1);
	    wd[iwmax] = wd[nsoln];
	    wd[nsoln] = 0.;
	    itemp = ipivot[nsoln];
	    ipivot[nsoln] = ipivot[iwmax];
	    ipivot[iwmax] = itemp;
	}

/*        Reduce column NSOLN so that the matrix of nonactive constraints */
/*        variables is triangular. */

	i__1 = niv + 1;
	for (j = m; j >= i__1; --j) {
	    jp = j - 1;

/*           When operating near the ME line, test to see if the pivot */
/*           element is near zero.  If so, use the largest element above */
/*           it as the pivot.  This is to maintain the sharp interface */
/*           between weighted and non-weighted rows in all cases. */

	    if (j == me + 1) {
		imax = me;
/* Computing 2nd power */
		d__1 = w[me + nsoln * w_dim1];
		amax = scale[me] * (d__1 * d__1);
		i__2 = niv;
		for (jp = j - 1; jp >= i__2; --jp) {
/* Computing 2nd power */
		    d__1 = w[jp + nsoln * w_dim1];
		    t = scale[jp] * (d__1 * d__1);
		    if (t > amax) {
			imax = jp;
			amax = t;
		    }
/* L310: */
		}
		jp = imax;
	    }

	    if (w[j + nsoln * w_dim1] != 0.) {
		drotmg_(&scale[jp], &scale[j], &w[jp + nsoln * w_dim1], &w[j 
			+ nsoln * w_dim1], sparam);
		w[j + nsoln * w_dim1] = 0.;
		i__2 = *n + 1 - nsoln;
		drotm_(&i__2, &w[jp + (nsoln + 1) * w_dim1], mdw, &w[j + (
			nsoln + 1) * w_dim1], mdw, sparam);
	    }
/* L320: */
	}

/*        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if */
/*        this is nonpositive or too large.  If this was true or if the */
/*        pivot term was zero, reject the column as dependent. */

	if (w[niv + nsoln * w_dim1] != 0.) {
	    isol = niv;
	    z2 = w[isol + (*n + 1) * w_dim1] / w[isol + nsoln * w_dim1];
	    z__[nsoln] = z2;
	    pos = z2 > 0.;
	    if (z2 * eanorm >= bnorm && pos) {
		pos = ! (blowup * z2 * eanorm >= bnorm);
	    }

/*           Try to add row ME+1 as an additional equality constraint. */
/*           Check size of proposed new solution component. */
/*           Reject it if it is too large. */

	} else if (niv <= me && w[me + 1 + nsoln * w_dim1] != 0.) {
	    isol = me + 1;
	    if (pos) {

/*              Swap rows ME+1 and NIV, and scale factors for these rows. */

		i__1 = *n + 1;
		dswap_(&i__1, &w[me + 1 + w_dim1], mdw, &w[niv + w_dim1], mdw)
			;
		dswap_(&c__1, &scale[me + 1], &c__1, &scale[niv], &c__1);
		itemp = itype[me + 1];
		itype[me + 1] = itype[niv];
		itype[niv] = itemp;
		++me;
	    }
	} else {
	    pos = FALSE_;
	}

	if (! pos) {
	    --nsoln;
	    --niv;
	}
	if (! (pos || done)) {
	    goto L290;
	}
    }
    goto L160;

/*     Else perform multiplier test and drop a constraint.  To compute */
/*     final solution.  Solve system, store results in X(*). */

/*     Copy right hand side into TEMP vector to use overwriting method. */

L330:
    isol = 1;
    if (nsoln >= isol) {
	dcopy_(&niv, &w[(*n + 1) * w_dim1 + 1], &c__1, &temp[1], &c__1);
	i__1 = isol;
	for (j = nsoln; j >= i__1; --j) {
	    if (j > krank) {
		i__ = niv - nsoln + j;
	    } else {
		i__ = j;
	    }

	    if (j > krank && j <= *l) {
		z__[j] = 0.;
	    } else {
		z__[j] = temp[i__] / w[i__ + j * w_dim1];
		i__2 = i__ - 1;
		d__1 = -z__[j];
		daxpy_(&i__2, &d__1, &w[j * w_dim1 + 1], &c__1, &temp[1], &
			c__1);
	    }
/* L340: */
	}
    }

/*     Solve system. */

    dcopy_(&nsoln, &z__[1], &c__1, &x[1], &c__1);

/*     Apply Householder transformations to X(*) if KRANK.LT.L */

    if (krank < *l) {
	i__1 = krank;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = krank + 1;
	    dh12_(&c__2, &i__, &i__2, l, &w[i__ + w_dim1], mdw, &h__[i__], &x[
		    1], &c__1, &c__1, &c__1);
/* L350: */
	}
    }

/*     Fill in trailing zeroes for constrained variables not in solution. */

    if (nsoln < *n) {
	i__1 = *n - nsoln;
	dcopy_(&i__1, &c_b95, &c__0, &x[nsoln + 1], &c__1);
    }

/*     Permute solution vector to natural order. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__;
L360:
	if (ipivot[j] == i__) {
	    goto L370;
	}
	++j;
	goto L360;

L370:
	ipivot[j] = ipivot[i__];
	ipivot[i__] = j;
	dswap_(&c__1, &x[j], &c__1, &x[i__], &c__1);
/* L380: */
    }

/*     Rescale the solution using the column scaling. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] *= d__[j];
/* L390: */
    }

    i__1 = m;
    for (i__ = nsoln + 1; i__ <= i__1; ++i__) {
	t = w[i__ + (*n + 1) * w_dim1];
	if (i__ <= me) {
	    t /= alamda;
	}
	t = scale[i__] * t * t;
	*rnorm += t;
/* L400: */
    }

    *rnorm = sqrt(*rnorm);
    return 0;
} /* dwnlsm_ */

/* DECK DWNLT1 */
/* Subroutine */ int dwnlt1_(integer *i__, integer *lend, integer *mend, 
	integer *ir, integer *mdw, logical *recalc, integer *imax, doublereal 
	*hbar, doublereal *h__, doublereal *scale, doublereal *w)
{
    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer j, k;
    extern integer idamax_(integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DWNLT1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to WNLIT */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLT1-S, DWNLT1-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     To update the column Sum Of Squares and find the pivot column. */
/*     The column Sum of Squares Vector will be updated at each step. */
/*     When numerically necessary, these values will be recomputed. */

/* ***SEE ALSO  DWNLIT */
/* ***ROUTINES CALLED  IDAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890620  Code extracted from WNLIT and made a subroutine.  (RWC)) */
/*   900604  DP version created from SP version.  (RWC) */
/* ***END PROLOGUE  DWNLT1 */



/* ***FIRST EXECUTABLE STATEMENT  DWNLT1 */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --h__;
    --scale;

    /* Function Body */
    if (*ir != 1 && ! (*recalc)) {

/*        Update column SS=sum of squares. */

	i__1 = *lend;
	for (j = *i__; j <= i__1; ++j) {
/* Computing 2nd power */
	    d__1 = w[*ir - 1 + j * w_dim1];
	    h__[j] -= scale[*ir - 1] * (d__1 * d__1);
/* L10: */
	}

/*        Test for numerical accuracy. */

	i__1 = *lend - *i__ + 1;
	*imax = idamax_(&i__1, &h__[*i__], &c__1) + *i__ - 1;
	*recalc = *hbar + h__[*imax] * .001f == *hbar;
    }

/*     If required, recalculate column SS, using rows IR through MEND. */

    if (*recalc) {
	i__1 = *lend;
	for (j = *i__; j <= i__1; ++j) {
	    h__[j] = 0.;
	    i__2 = *mend;
	    for (k = *ir; k <= i__2; ++k) {
/* Computing 2nd power */
		d__1 = w[k + j * w_dim1];
		h__[j] += scale[k] * (d__1 * d__1);
/* L20: */
	    }
/* L30: */
	}

/*        Find column with largest SS. */

	i__1 = *lend - *i__ + 1;
	*imax = idamax_(&i__1, &h__[*i__], &c__1) + *i__ - 1;
	*hbar = h__[*imax];
    }
    return 0;
} /* dwnlt1_ */

/* DECK DWNLT2 */
logical dwnlt2_(integer *me, integer *mend, integer *ir, doublereal *factor, 
	doublereal *tau, doublereal *scale, doublereal *wic)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    logical ret_val;

    /* Local variables */
    static integer j;
    static doublereal t, rn, sn;

/* ***BEGIN PROLOGUE  DWNLT2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to WNLIT */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLT2-S, DWNLT2-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     To test independence of incoming column. */

/*     Test the column IC to determine if it is linearly independent */
/*     of the columns already in the basis.  In the initial tri. step, */
/*     we usually want the heavy weight ALAMDA to be included in the */
/*     test for independence.  In this case, the value of FACTOR will */
/*     have been set to 1.E0 before this procedure is invoked. */
/*     In the potentially rank deficient problem, the value of FACTOR */
/*     will have been set to ALSQ=ALAMDA**2 to remove the effect of the */
/*     heavy weight from the test for independence. */

/*     Write new column as partitioned vector */
/*           (A1)  number of components in solution so far = NIV */
/*           (A2)  M-NIV components */
/*     And compute  SN = inverse weighted length of A1 */
/*                  RN = inverse weighted length of A2 */
/*     Call the column independent when RN .GT. TAU*SN */

/* ***SEE ALSO  DWNLIT */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890620  Code extracted from WNLIT and made a subroutine.  (RWC)) */
/*   900604  DP version created from SP version.  (RWC) */
/* ***END PROLOGUE  DWNLT2 */


/* ***FIRST EXECUTABLE STATEMENT  DWNLT2 */
    /* Parameter adjustments */
    --wic;
    --scale;

    /* Function Body */
    sn = 0.f;
    rn = 0.f;
    i__1 = *mend;
    for (j = 1; j <= i__1; ++j) {
	t = scale[j];
	if (j <= *me) {
	    t /= *factor;
	}
/* Computing 2nd power */
	d__1 = wic[j];
	t *= d__1 * d__1;

	if (j < *ir) {
	    sn += t;
	} else {
	    rn += t;
	}
/* L10: */
    }
/* Computing 2nd power */
    d__1 = *tau;
    ret_val = rn > sn * (d__1 * d__1);
    return ret_val;
} /* dwnlt2_ */

/* DECK DWNLT3 */
/* Subroutine */ int dwnlt3_(integer *i__, integer *imax, integer *m, integer 
	*mdw, integer *ipivot, doublereal *h__, doublereal *w)
{
    /* System generated locals */
    integer w_dim1, w_offset;

    /* Local variables */
    static doublereal t;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer itemp;

/* ***BEGIN PROLOGUE  DWNLT3 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to WNLIT */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     Perform column interchange. */
/*     Exchange elements of permuted index vector and perform column */
/*     interchanges. */

/* ***SEE ALSO  DWNLIT */
/* ***ROUTINES CALLED  DSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890620  Code extracted from WNLIT and made a subroutine.  (RWC)) */
/*   900604  DP version created from SP version.  (RWC) */
/* ***END PROLOGUE  DWNLT3 */



/* ***FIRST EXECUTABLE STATEMENT  DWNLT3 */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --ipivot;
    --h__;

    /* Function Body */
    if (*imax != *i__) {
	itemp = ipivot[*i__];
	ipivot[*i__] = ipivot[*imax];
	ipivot[*imax] = itemp;

	dswap_(m, &w[*imax * w_dim1 + 1], &c__1, &w[*i__ * w_dim1 + 1], &c__1)
		;

	t = h__[*imax];
	h__[*imax] = h__[*i__];
	h__[*i__] = t;
    }
    return 0;
} /* dwnlt3_ */

/* DECK DWNNLS */
/* Subroutine */ int dwnnls_(doublereal *w, integer *mdw, integer *me, 
	integer *ma, integer *n, integer *l, doublereal *prgopt, doublereal *
	x, doublereal *rnorm, integer *mode, integer *iwork, doublereal *work)
{
    /* System generated locals */
    address a__1[2];
    integer w_dim1, w_offset, i__1[2];
    char ch__1[62], ch__2[64];

    /* Builtin functions */
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer l1, l2, l3, l4, l5, lw, liw;
    static char xern1[8];
    extern /* Subroutine */ int dwnlsm_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xermsg_(char *, char *, char *, integer *, integer 
	    *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___300 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___302 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DWNNLS */
/* ***PURPOSE  Solve a linearly constrained least squares problem with */
/*            equality constraints and nonnegativity constraints on */
/*            selected variables. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  K1A2A */
/* ***TYPE      DOUBLE PRECISION (WNNLS-S, DWNNLS-D) */
/* ***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING, */
/*             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS, */
/*             NONNEGATIVITY CONSTRAINTS, QUADRATIC PROGRAMMING */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */

/*     This subprogram solves a linearly constrained least squares */
/*     problem.  Suppose there are given matrices E and A of */
/*     respective dimensions ME by N and MA by N, and vectors F */
/*     and B of respective lengths ME and MA.  This subroutine */
/*     solves the problem */

/*               EX = F, (equations to be exactly satisfied) */

/*               AX = B, (equations to be approximately satisfied, */
/*                        in the least squares sense) */

/*               subject to components L+1,...,N nonnegative */

/*     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted. */

/*     The problem is reposed as problem DWNNLS */

/*               (WT*E)X = (WT*F) */
/*               (   A)    (   B), (least squares) */
/*               subject to components L+1,...,N nonnegative. */

/*     The subprogram chooses the heavy weight (or penalty parameter) WT. */

/*     The parameters for DWNNLS are */

/*     INPUT.. All TYPE REAL variables are DOUBLE PRECISION */

/*     W(*,*),MDW,  The array W(*,*) is double subscripted with first */
/*     ME,MA,N,L    dimensioning parameter equal to MDW.  For this */
/*                  discussion let us call M = ME + MA.  Then MDW */
/*                  must satisfy MDW.GE.M.  The condition MDW.LT.M */
/*                  is an error. */

/*                  The array W(*,*) contains the matrices and vectors */

/*                       (E  F) */
/*                       (A  B) */

/*                  in rows and columns 1,...,M and 1,...,N+1 */
/*                  respectively.  Columns 1,...,L correspond to */
/*                  unconstrained variables X(1),...,X(L).  The */
/*                  remaining variables are constrained to be */
/*                  nonnegative. The condition L.LT.0 or L.GT.N is */
/*                  an error. */

/*     PRGOPT(*)    This double precision array is the option vector. */
/*                  If the user is satisfied with the nominal */
/*                  subprogram features set */

/*                  PRGOPT(1)=1 (or PRGOPT(1)=1.0) */

/*                  Otherwise PRGOPT(*) is a linked list consisting of */
/*                  groups of data of the following form */

/*                  LINK */
/*                  KEY */
/*                  DATA SET */

/*                  The parameters LINK and KEY are each one word. */
/*                  The DATA SET can be comprised of several words. */
/*                  The number of items depends on the value of KEY. */
/*                  The value of LINK points to the first */
/*                  entry of the next group of data within */
/*                  PRGOPT(*).  The exception is when there are */
/*                  no more options to change.  In that */
/*                  case LINK=1 and the values KEY and DATA SET */
/*                  are not referenced. The general layout of */
/*                  PRGOPT(*) is as follows. */

/*               ...PRGOPT(1)=LINK1 (link to first entry of next group) */
/*               .  PRGOPT(2)=KEY1 (key to the option change) */
/*               .  PRGOPT(3)=DATA VALUE (data value for this change) */
/*               .       . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK1)=LINK2 (link to the first entry of */
/*               .                       next group) */
/*               .  PRGOPT(LINK1+1)=KEY2 (key to the option change) */
/*               .  PRGOPT(LINK1+2)=DATA VALUE */
/*               ...     . */
/*               .       . */
/*               .       . */
/*               ...PRGOPT(LINK)=1 (no more options to change) */

/*                  Values of LINK that are nonpositive are errors. */
/*                  A value of LINK.GT.NLINK=100000 is also an error. */
/*                  This helps prevent using invalid but positive */
/*                  values of LINK that will probably extend */
/*                  beyond the program limits of PRGOPT(*). */
/*                  Unrecognized values of KEY are ignored.  The */
/*                  order of the options is arbitrary and any number */
/*                  of options can be changed with the following */
/*                  restriction.  To prevent cycling in the */
/*                  processing of the option array a count of the */
/*                  number of options changed is maintained. */
/*                  Whenever this count exceeds NOPT=1000 an error */
/*                  message is printed and the subprogram returns. */

/*                  OPTIONS.. */

/*                  KEY=6 */
/*                         Scale the nonzero columns of the */
/*                  entire data matrix */
/*                  (E) */
/*                  (A) */
/*                  to have length one. The DATA SET for */
/*                  this option is a single value.  It must */
/*                  be nonzero if unit length column scaling is */
/*                  desired. */

/*                  KEY=7 */
/*                         Scale columns of the entire data matrix */
/*                  (E) */
/*                  (A) */
/*                  with a user-provided diagonal matrix. */
/*                  The DATA SET for this option consists */
/*                  of the N diagonal scaling factors, one for */
/*                  each matrix column. */

/*                  KEY=8 */
/*                         Change the rank determination tolerance from */
/*                  the nominal value of SQRT(SRELPR).  This quantity */
/*                  can be no smaller than SRELPR, The arithmetic- */
/*                  storage precision.  The quantity used */
/*                  here is internally restricted to be at */
/*                  least SRELPR.  The DATA SET for this option */
/*                  is the new tolerance. */

/*                  KEY=9 */
/*                         Change the blow-up parameter from the */
/*                  nominal value of SQRT(SRELPR).  The reciprocal of */
/*                  this parameter is used in rejecting solution */
/*                  components as too large when a variable is */
/*                  first brought into the active set.  Too large */
/*                  means that the proposed component times the */
/*                  reciprocal of the parameter is not less than */
/*                  the ratio of the norms of the right-side */
/*                  vector and the data matrix. */
/*                  This parameter can be no smaller than SRELPR, */
/*                  the arithmetic-storage precision. */

/*                  For example, suppose we want to provide */
/*                  a diagonal matrix to scale the problem */
/*                  matrix and change the tolerance used for */
/*                  determining linear dependence of dropped col */
/*                  vectors.  For these options the dimensions of */
/*                  PRGOPT(*) must be at least N+6.  The FORTRAN */
/*                  statements defining these options would */
/*                  be as follows. */

/*                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*)) */
/*                  PRGOPT(2)=7 (user-provided scaling key) */

/*                  CALL DCOPY(N,D,1,PRGOPT(3),1) (copy the N */
/*                  scaling factors from a user array called D(*) */
/*                  into PRGOPT(3)-PRGOPT(N+2)) */

/*                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*)) */
/*                  PRGOPT(N+4)=8 (linear dependence tolerance key) */
/*                  PRGOPT(N+5)=... (new value of the tolerance) */

/*                  PRGOPT(N+6)=1 (no more options to change) */


/*     IWORK(1),    The amounts of working storage actually allocated */
/*     IWORK(2)     for the working arrays WORK(*) and IWORK(*), */
/*                  respectively.  These quantities are compared with */
/*                  the actual amounts of storage needed for DWNNLS( ). */
/*                  Insufficient storage allocated for either WORK(*) */
/*                  or IWORK(*) is considered an error.  This feature */
/*                  was included in DWNNLS( ) because miscalculating */
/*                  the storage formulas for WORK(*) and IWORK(*) */
/*                  might very well lead to subtle and hard-to-find */
/*                  execution errors. */

/*                  The length of WORK(*) must be at least */

/*                  LW = ME+MA+5*N */
/*                  This test will not be made if IWORK(1).LE.0. */

/*                  The length of IWORK(*) must be at least */

/*                  LIW = ME+MA+N */
/*                  This test will not be made if IWORK(2).LE.0. */

/*     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION */

/*     X(*)         An array dimensioned at least N, which will */
/*                  contain the N components of the solution vector */
/*                  on output. */

/*     RNORM        The residual norm of the solution.  The value of */
/*                  RNORM contains the residual vector length of the */
/*                  equality constraints and least squares equations. */

/*     MODE         The value of MODE indicates the success or failure */
/*                  of the subprogram. */

/*                  MODE = 0  Subprogram completed successfully. */

/*                       = 1  Max. number of iterations (equal to */
/*                            3*(N-L)) exceeded. Nearly all problems */
/*                            should complete in fewer than this */
/*                            number of iterations. An approximate */
/*                            solution and its corresponding residual */
/*                            vector length are in X(*) and RNORM. */

/*                       = 2  Usage error occurred.  The offending */
/*                            condition is noted with the error */
/*                            processing subprogram, XERMSG( ). */

/*     User-designated */
/*     Working arrays.. */

/*     WORK(*)      A double precision working array of length at least */
/*                  M + 5*N. */

/*     IWORK(*)     An integer-valued working array of length at least */
/*                  M+N. */

/* ***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for */
/*                 linear least squares problems with equality and */
/*                 nonnegativity constraints, Report SAND77-0552, Sandia */
/*                 Laboratories, June 1978. */
/*               K. H. Haskell and R. J. Hanson, Selected algorithms for */
/*                 the linearly constrained least squares problem - a */
/*                 users guide, Report SAND78-1290, Sandia Laboratories, */
/*                 August 1979. */
/*               K. H. Haskell and R. J. Hanson, An algorithm for */
/*                 linear least squares problems with equality and */
/*                 nonnegativity constraints, Mathematical Programming */
/*                 21 (1981), pp. 98-118. */
/*               R. J. Hanson and K. H. Haskell, Two algorithms for the */
/*                 linearly constrained least squares problem, ACM */
/*                 Transactions on Mathematical Software, September 1982. */
/*               C. L. Lawson and R. J. Hanson, Solving Least Squares */
/*                 Problems, Prentice-Hall, Inc., 1974. */
/* ***ROUTINES CALLED  DWNLSM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890618  Completely restructured and revised.  (WRB & RWC) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Convert XERRWV calls to XERMSG calls, change Prologue */
/*           comments to agree with WNNLS.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DWNNLS */
/* ***FIRST EXECUTABLE STATEMENT  DWNNLS */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --prgopt;
    --x;
    --iwork;
    --work;

    /* Function Body */
    *mode = 0;
    if (*ma + *me <= 0 || *n <= 0) {
	return 0;
    }

    if (iwork[1] > 0) {
	lw = *me + *ma + *n * 5;
	if (iwork[1] < lw) {
	    s_wsfi(&io___300);
	    do_fio(&c__1, (char *)&lw, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 54, a__1[0] = "INSUFFICIENT STORAGE ALLOCATED FOR WORK"
		    "(*), NEED LW = ";
	    i__1[1] = 8, a__1[1] = xern1;
	    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)62);
	    xermsg_("SLATEC", "DWNNLS", ch__1, &c__2, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)62);
	    *mode = 2;
	    return 0;
	}
    }

    if (iwork[2] > 0) {
	liw = *me + *ma + *n;
	if (iwork[2] < liw) {
	    s_wsfi(&io___302);
	    do_fio(&c__1, (char *)&liw, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 56, a__1[0] = "INSUFFICIENT STORAGE ALLOCATED FOR IWOR"
		    "K(*), NEED LIW = ";
	    i__1[1] = 8, a__1[1] = xern1;
	    s_cat(ch__2, a__1, i__1, &c__2, (ftnlen)64);
	    xermsg_("SLATEC", "DWNNLS", ch__2, &c__2, &c__1, (ftnlen)6, (
		    ftnlen)6, (ftnlen)64);
	    *mode = 2;
	    return 0;
	}
    }

    if (*mdw < *me + *ma) {
	xermsg_("SLATEC", "DWNNLS", "THE VALUE MDW.LT.ME+MA IS AN ERROR", &
		c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)34);
	*mode = 2;
	return 0;
    }

    if (*l < 0 || *l > *n) {
	xermsg_("SLATEC", "DWNNLS", "L.GE.0 .AND. L.LE.N IS REQUIRED", &c__2, 
		&c__1, (ftnlen)6, (ftnlen)6, (ftnlen)31);
	*mode = 2;
	return 0;
    }

/*     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS */
/*     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS */
/*     REQUIRED BY THE MAIN SUBROUTINE DWNLSM( ). */

    l1 = *n + 1;
    l2 = l1 + *n;
    l3 = l2 + *me + *ma;
    l4 = l3 + *n;
    l5 = l4 + *n;

    dwnlsm_(&w[w_offset], mdw, me, ma, n, l, &prgopt[1], &x[1], rnorm, mode, &
	    iwork[1], &iwork[l1], &work[1], &work[l1], &work[l2], &work[l3], &
	    work[l4], &work[l5]);
    return 0;
} /* dwnnls_ */

/* DECK FDUMP */
/* Subroutine */ int fdump_(void)
{
/* ***BEGIN PROLOGUE  FDUMP */
/* ***PURPOSE  Symbolic dump (should be locally written). */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3 */
/* ***TYPE      ALL (FDUMP-A) */
/* ***KEYWORDS  ERROR, XERMSG */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*        ***Note*** Machine Dependent Routine */
/*        FDUMP is intended to be replaced by a locally written */
/*        version which produces a symbolic dump.  Failing this, */
/*        it should be replaced by a version which prints the */
/*        subprogram nesting list.  Note that this dump must be */
/*        printed on each of up to five files, as indicated by the */
/*        XGETUA routine.  See XSETUA and XGETUA for details. */

/*     Written by Ron Jones, with SLATEC Common Math Library Subcommittee */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  FDUMP */
/* ***FIRST EXECUTABLE STATEMENT  FDUMP */
    return 0;
} /* fdump_ */

/* DECK I1MACH */
integer i1mach_(integer *i__)
{
    /* Format strings */
    static char fmt_9000[] = "(\0021ERROR    1 IN I1MACH - I OUT OF BOUND"
	    "S\002)";

    /* System generated locals */
    integer ret_val;
    static integer equiv_0[16];

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
#define imach (equiv_0)
#define output (equiv_0 + 3)

    /* Fortran I/O blocks */
    static cilist io___310 = { 0, 0, 0, fmt_9000, 0 };


/* ***BEGIN PROLOGUE  I1MACH */
/* ***PURPOSE  Return integer machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      INTEGER (I1MACH-I) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   I1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument and can be referenced as follows: */

/*        K = I1MACH(I) */

/*   where I=1,...,16.  The (output) value of K above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   I/O unit numbers: */
/*     I1MACH( 1) = the standard input unit. */
/*     I1MACH( 2) = the standard output unit. */
/*     I1MACH( 3) = the standard punch unit. */
/*     I1MACH( 4) = the standard error message unit. */

/*   Words: */
/*     I1MACH( 5) = the number of bits per integer storage unit. */
/*     I1MACH( 6) = the number of characters per integer storage unit. */

/*   Integers: */
/*     assume integers are represented in the S-digit, base-A form */

/*                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) */

/*                where 0 .LE. X(I) .LT. A for I=0,...,S-1. */
/*     I1MACH( 7) = A, the base. */
/*     I1MACH( 8) = S, the number of base-A digits. */
/*     I1MACH( 9) = A**S - 1, the largest magnitude. */

/*   Floating-Point Numbers: */
/*     Assume floating-point numbers are represented in the T-digit, */
/*     base-B form */
/*                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*                where 0 .LE. X(I) .LT. B for I=1,...,T, */
/*                0 .LT. X(1), and EMIN .LE. E .LE. EMAX. */
/*     I1MACH(10) = B, the base. */

/*   Single-Precision: */
/*     I1MACH(11) = T, the number of base-B digits. */
/*     I1MACH(12) = EMIN, the smallest exponent E. */
/*     I1MACH(13) = EMAX, the largest exponent E. */

/*   Double-Precision: */
/*     I1MACH(14) = T, the number of base-B digits. */
/*     I1MACH(15) = EMIN, the smallest exponent E. */
/*     I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for */
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   891012  Added VAX G-floating constants.  (WRB) */
/*   891012  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16. */
/*           (RWC) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added Convex -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/*   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler */
/*           options.  (DWL, RWC and WRB). */
/* ***END PROLOGUE  I1MACH */


/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        129 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1025 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA IMACH( 1) /          7 / */
/*     DATA IMACH( 2) /          2 / */
/*     DATA IMACH( 3) /          2 / */
/*     DATA IMACH( 4) /          2 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         33 / */
/*     DATA IMACH( 9) / Z1FFFFFFFF / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -256 / */
/*     DATA IMACH(13) /        255 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /       -256 / */
/*     DATA IMACH(16) /        255 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /        -50 / */
/*     DATA IMACH(16) /         76 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /     -32754 / */
/*     DATA IMACH(16) /      32780 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -4095 / */
/*     DATA IMACH(13) /       4094 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -4095 / */
/*     DATA IMACH(16) /       4094 / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /    6LOUTPUT/ */
/*     DATA IMACH( 5) /         60 / */
/*     DATA IMACH( 6) /         10 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         48 / */
/*     DATA IMACH( 9) / 00007777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /       -929 / */
/*     DATA IMACH(13) /       1070 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /       -929 / */
/*     DATA IMACH(16) /       1069 / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / Z'7FFFFFFF' / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1023 / */
/*     DATA IMACH(13) /       1023 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16383 / */
/*     DATA IMACH(16) /      16383 / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -pd8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1023 / */
/*     DATA IMACH(13) /       1023 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 46 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         46 / */
/*     DATA IMACH( 9) / 1777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 64 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 777777777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */

/*     DATA IMACH( 1) /         11 / */
/*     DATA IMACH( 2) /         12 / */
/*     DATA IMACH( 3) /          8 / */
/*     DATA IMACH( 4) /         10 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         24 / */
/*     DATA IMACH( 6) /          3 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         23 / */
/*     DATA IMACH( 9) /    8388607 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         38 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /         43 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         63 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         39 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         55 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          7 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1015 / */
/*     DATA IMACH(16) /       1017 / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) /  Z7FFFFFFF / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE IBM PC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         54 / */
/*     DATA IMACH(15) /       -101 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         62 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1021 / */
/*     DATA IMACH(13) /       1024 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16381 / */
/*     DATA IMACH(16) /      16384 / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          1 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /      -1024 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR */

/*     DATA IMACH( 1) /          1 / */
/*     DATA IMACH( 2) /          1 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/* ***FIRST EXECUTABLE STATEMENT  I1MACH */
    if (*i__ < 1 || *i__ > 16) {
	goto L10;
    }

    ret_val = imach[*i__ - 1];
    return ret_val;

L10:
    io___310.ciunit = *output;
    s_wsfe(&io___310);
    e_wsfe();

/*     CALL FDUMP */

    s_stop("", (ftnlen)0);
    return ret_val;
} /* i1mach_ */

#undef output
#undef imach


/* DECK IDAMAX */
integer idamax_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, ix;
    static doublereal dmax__, xmag;

/* ***BEGIN PROLOGUE  IDAMAX */
/* ***PURPOSE  Find the smallest index of that component of a vector */
/*            having the maximum magnitude. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A2 */
/* ***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */

/*     --Output-- */
/*   IDAMAX  smallest index (zero if N .LE. 0) */

/*     Find smallest index of maximum magnitude of double precision DX. */
/*     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)), */
/*     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900821  Modified to correct problem with a negative increment. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  IDAMAX */
/* ***FIRST EXECUTABLE STATEMENT  IDAMAX */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }

    if (*incx == 1) {
	goto L20;
    }

/*     Code for increments not equal to 1. */

    ix = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    dmax__ = (d__1 = dx[ix], abs(d__1));
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	xmag = (d__1 = dx[ix], abs(d__1));
	if (xmag > dmax__) {
	    ret_val = i__;
	    dmax__ = xmag;
	}
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*     Code for increments equal to 1. */

L20:
    dmax__ = abs(dx[1]);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	xmag = (d__1 = dx[i__], abs(d__1));
	if (xmag > dmax__) {
	    ret_val = i__;
	    dmax__ = xmag;
	}
/* L30: */
    }
    return ret_val;
} /* idamax_ */

/* DECK J4SAVE */
integer j4save_(integer *iwhich, integer *ivalue, logical *iset)
{
    /* Initialized data */

    static integer iparam[9] = { 0,2,0,10,1,0,0,0,0 };

    /* System generated locals */
    integer ret_val;

/* ***BEGIN PROLOGUE  J4SAVE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Save or recall global variables needed by error */
/*            handling routines. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***TYPE      INTEGER (J4SAVE-I) */
/* ***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        J4SAVE saves and recalls several global variables needed */
/*        by the library error handling routines. */

/*     Description of Parameters */
/*      --Input-- */
/*        IWHICH - Index of item desired. */
/*                = 1 Refers to current error number. */
/*                = 2 Refers to current error control flag. */
/*                = 3 Refers to current unit number to which error */
/*                    messages are to be sent.  (0 means use standard.) */
/*                = 4 Refers to the maximum number of times any */
/*                     message is to be printed (as set by XERMAX). */
/*                = 5 Refers to the total number of units to which */
/*                     each error message is to be written. */
/*                = 6 Refers to the 2nd unit for error messages */
/*                = 7 Refers to the 3rd unit for error messages */
/*                = 8 Refers to the 4th unit for error messages */
/*                = 9 Refers to the 5th unit for error messages */
/*        IVALUE - The value to be set for the IWHICH-th parameter, */
/*                 if ISET is .TRUE. . */
/*        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE */
/*                 given the value, IVALUE.  If ISET=.FALSE., the */
/*                 IWHICH-th parameter will be unchanged, and IVALUE */
/*                 is a dummy parameter. */
/*      --Output-- */
/*        The (old) value of the IWHICH-th parameter will be returned */
/*        in the function value, J4SAVE. */

/* ***SEE ALSO  XERMSG */
/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900205  Minor modifications to prologue.  (WRB) */
/*   900402  Added TYPE section.  (WRB) */
/*   910411  Added KEYWORDS section.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  J4SAVE */
/* ***FIRST EXECUTABLE STATEMENT  J4SAVE */
    ret_val = iparam[(0 + (0 + (*iwhich - 1 << 2))) / 4];
    if (*iset) {
	iparam[*iwhich - 1] = *ivalue;
    }
    return ret_val;
} /* j4save_ */

/* DECK XERCNT */
/* Subroutine */ int xercnt_(char *librar, char *subrou, char *messg, integer 
	*nerr, integer *level, integer *kontrl, ftnlen librar_len, ftnlen 
	subrou_len, ftnlen messg_len)
{
/* ***BEGIN PROLOGUE  XERCNT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Allow user control over handling of errors. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERCNT-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        Allows user control over handling of individual errors. */
/*        Just after each message is recorded, but before it is */
/*        processed any further (i.e., before it is printed or */
/*        a decision to abort is made), a call is made to XERCNT. */
/*        If the user has provided his own version of XERCNT, he */
/*        can then override the value of KONTROL used in processing */
/*        this message by redefining its value. */
/*        KONTRL may be set to any value from -2 to 2. */
/*        The meanings for KONTRL are the same as in XSETF, except */
/*        that the value of KONTRL changes only for this message. */
/*        If KONTRL is set to a value outside the range from -2 to 2, */
/*        it will be moved back into that range. */

/*     Description of Parameters */

/*      --Input-- */
/*        LIBRAR - the library that the routine is in. */
/*        SUBROU - the subroutine that XERMSG is being called from */
/*        MESSG  - the first 20 characters of the error message. */
/*        NERR   - same as in the call to XERMSG. */
/*        LEVEL  - same as in the call to XERMSG. */
/*        KONTRL - the current value of the control flag as set */
/*                 by a call to XSETF. */

/*      --Output-- */
/*        KONTRL - the new value of KONTRL.  If KONTRL is not */
/*                 defined, it will remain at its original value. */
/*                 This changed value of control affects only */
/*                 the current occurrence of the current message. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900206  Routine changed from user-callable to subsidiary.  (WRB) */
/*   900510  Changed calling sequence to include LIBRARY and SUBROUTINE */
/*           names, changed routine name from XERCTL to XERCNT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERCNT */
/* ***FIRST EXECUTABLE STATEMENT  XERCNT */
    return 0;
} /* xercnt_ */

/* DECK XERHLT */
/* Subroutine */ int xerhlt_(char *messg, ftnlen messg_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

/* ***BEGIN PROLOGUE  XERHLT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Abort program execution and print error message. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERHLT-A) */
/* ***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        ***Note*** machine dependent routine */
/*        XERHLT aborts the execution of the program. */
/*        The error message causing the abort is given in the calling */
/*        sequence, in case one needs it for printing on a dayfile, */
/*        for example. */

/*     Description of Parameters */
/*        MESSG is as in XERMSG. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900206  Routine changed from user-callable to subsidiary.  (WRB) */
/*   900510  Changed calling sequence to delete length of character */
/*           and changed routine name from XERABT to XERHLT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERHLT */
/* ***FIRST EXECUTABLE STATEMENT  XERHLT */
    s_stop("", (ftnlen)0);
    return 0;
} /* xerhlt_ */

/* DECK XERMSG */
/* Subroutine */ int xermsg_(char *librar, char *subrou, char *messg, integer 
	*nerr, integer *level, ftnlen librar_len, ftnlen subrou_len, ftnlen 
	messg_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    char ch__1[87];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_len(char *, ftnlen), s_wsfi(icilist *), do_fio(integer *, char *
	    , ftnlen), e_wsfi(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, lerr;
    static char temp[72];
    extern /* Subroutine */ int fdump_(void);
    static char xlibr[8];
    static integer ltemp, kount;
    static char xsubr[8];
    extern integer j4save_(integer *, integer *, logical *);
    static integer llevel, maxmes;
    static char lfirst[20];
    extern /* Subroutine */ int xercnt_(char *, char *, char *, integer *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer lkntrl, kdummy;
    extern /* Subroutine */ int xerhlt_(char *, ftnlen);
    static integer mkntrl;
    extern /* Subroutine */ int xersve_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen), xerprn_(
	    char *, integer *, char *, integer *, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___329 = { 0, temp, 0, "('ERROR NUMBER = ', I8)", 72, 1 }
	    ;


/* ***BEGIN PROLOGUE  XERMSG */
/* ***PURPOSE  Process error messages for SLATEC and other libraries. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERMSG-A) */
/* ***KEYWORDS  ERROR MESSAGE, XERROR */
/* ***AUTHOR  Fong, Kirby, (NMFECC at LLNL) */
/* ***DESCRIPTION */

/*   XERMSG processes a diagnostic message in a manner determined by the */
/*   value of LEVEL and the current value of the library error control */
/*   flag, KONTRL.  See subroutine XSETF for details. */

/*    LIBRAR   A character constant (or character variable) with the name */
/*             of the library.  This will be 'SLATEC' for the SLATEC */
/*             Common Math Library.  The error handling package is */
/*             general enough to be used by many libraries */
/*             simultaneously, so it is desirable for the routine that */
/*             detects and reports an error to identify the library name */
/*             as well as the routine name. */

/*    SUBROU   A character constant (or character variable) with the name */
/*             of the routine that detected the error.  Usually it is the */
/*             name of the routine that is calling XERMSG.  There are */
/*             some instances where a user callable library routine calls */
/*             lower level subsidiary routines where the error is */
/*             detected.  In such cases it may be more informative to */
/*             supply the name of the routine the user called rather than */
/*             the name of the subsidiary routine that detected the */
/*             error. */

/*    MESSG    A character constant (or character variable) with the text */
/*             of the error or warning message.  In the example below, */
/*             the message is a character constant that contains a */
/*             generic message. */

/*                   CALL XERMSG ('SLATEC', 'MMPY', */
/*                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION', */
/*                  *3, 1) */

/*             It is possible (and is sometimes desirable) to generate a */
/*             specific message--e.g., one that contains actual numeric */
/*             values.  Specific numeric values can be converted into */
/*             character strings using formatted WRITE statements into */
/*             character variables.  This is called standard Fortran */
/*             internal file I/O and is exemplified in the first three */
/*             lines of the following example.  You can also catenate */
/*             substrings of characters to construct the error message. */
/*             Here is an example showing the use of both writing to */
/*             an internal file and catenating character strings. */

/*                   CHARACTER*5 CHARN, CHARL */
/*                   WRITE (CHARN,10) N */
/*                   WRITE (CHARL,10) LDA */
/*                10 FORMAT(I5) */
/*                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN// */
/*                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'// */
/*                  *   CHARL, 3, 1) */

/*             There are two subtleties worth mentioning.  One is that */
/*             the // for character catenation is used to construct the */
/*             error message so that no single character constant is */
/*             continued to the next line.  This avoids confusion as to */
/*             whether there are trailing blanks at the end of the line. */
/*             The second is that by catenating the parts of the message */
/*             as an actual argument rather than encoding the entire */
/*             message into one large character variable, we avoid */
/*             having to know how long the message will be in order to */
/*             declare an adequate length for that large character */
/*             variable.  XERMSG calls XERPRN to print the message using */
/*             multiple lines if necessary.  If the message is very long, */
/*             XERPRN will break it into pieces of 72 characters (as */
/*             requested by XERMSG) for printing on multiple lines. */
/*             Also, XERMSG asks XERPRN to prefix each line with ' *  ' */
/*             so that the total line length could be 76 characters. */
/*             Note also that XERPRN scans the error message backwards */
/*             to ignore trailing blanks.  Another feature is that */
/*             the substring '$$' is treated as a new line sentinel */
/*             by XERPRN.  If you want to construct a multiline */
/*             message without having to count out multiples of 72 */
/*             characters, just use '$$' as a separator.  '$$' */
/*             obviously must occur within 72 characters of the */
/*             start of each line to have its intended effect since */
/*             XERPRN is asked to wrap around at 72 characters in */
/*             addition to looking for '$$'. */

/*    NERR     An integer value that is chosen by the library routine's */
/*             author.  It must be in the range -99 to 999 (three */
/*             printable digits).  Each distinct error should have its */
/*             own error number.  These error numbers should be described */
/*             in the machine readable documentation for the routine. */
/*             The error numbers need be unique only within each routine, */
/*             so it is reasonable for each routine to start enumerating */
/*             errors from 1 and proceeding to the next integer. */

/*    LEVEL    An integer value in the range 0 to 2 that indicates the */
/*             level (severity) of the error.  Their meanings are */

/*            -1  A warning message.  This is used if it is not clear */
/*                that there really is an error, but the user's attention */
/*                may be needed.  An attempt is made to only print this */
/*                message once. */

/*             0  A warning message.  This is used if it is not clear */
/*                that there really is an error, but the user's attention */
/*                may be needed. */

/*             1  A recoverable error.  This is used even if the error is */
/*                so serious that the routine cannot return any useful */
/*                answer.  If the user has told the error package to */
/*                return after recoverable errors, then XERMSG will */
/*                return to the Library routine which can then return to */
/*                the user's routine.  The user may also permit the error */
/*                package to terminate the program upon encountering a */
/*                recoverable error. */

/*             2  A fatal error.  XERMSG will not return to its caller */
/*                after it receives a fatal error.  This level should */
/*                hardly ever be used; it is much better to allow the */
/*                user a chance to recover.  An example of one of the few */
/*                cases in which it is permissible to declare a level 2 */
/*                error is a reverse communication Library routine that */
/*                is likely to be called repeatedly until it integrates */
/*                across some interval.  If there is a serious error in */
/*                the input such that another step cannot be taken and */
/*                the Library routine is called again without the input */
/*                error having been corrected by the caller, the Library */
/*                routine will probably be called forever with improper */
/*                input.  In this case, it is reasonable to declare the */
/*                error to be fatal. */

/*    Each of the arguments to XERMSG is input; none will be modified by */
/*    XERMSG.  A routine may make multiple calls to XERMSG with warning */
/*    level messages; however, after a call to XERMSG with a recoverable */
/*    error, the routine should return to the user.  Do not try to call */
/*    XERMSG with a second recoverable error after the first recoverable */
/*    error because the error package saves the error number.  The user */
/*    can retrieve this error number by calling another entry point in */
/*    the error handling package and then clear the error number when */
/*    recovering from the error.  Calling XERMSG in succession causes the */
/*    old error number to be overwritten by the latest error number. */
/*    This is considered harmless for error numbers associated with */
/*    warning messages but must not be done for error numbers of serious */
/*    errors.  After a call to XERMSG with a recoverable error, the user */
/*    must be given a chance to call NUMXER or XERCLR to retrieve or */
/*    clear the error number. */
/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   880101  DATE WRITTEN */
/*   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988. */
/*           THERE ARE TWO BASIC CHANGES. */
/*           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO */
/*               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES */
/*               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS */
/*               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE */
/*               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER */
/*               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY */
/*               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE */
/*               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76. */
/*           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE */
/*               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE */
/*               OF LOWER CASE. */
/*   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30. */
/*           THE PRINCIPAL CHANGES ARE */
/*           1.  CLARIFY COMMENTS IN THE PROLOGUES */
/*           2.  RENAME XRPRNT TO XERPRN */
/*           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES */
/*               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE / */
/*               CHARACTER FOR NEW RECORDS. */
/*   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO */
/*           CLEAN UP THE CODING. */
/*   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN */
/*           PREFIX. */
/*   891013  REVISED TO CORRECT COMMENTS. */
/*   891214  Prologue converted to Version 4.0 format.  (WRB) */
/*   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but */
/*           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added */
/*           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and */
/*           XERCTL to XERCNT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERMSG */
/* ***FIRST EXECUTABLE STATEMENT  XERMSG */
    lkntrl = j4save_(&c__2, &c__0, &c_false);
    maxmes = j4save_(&c__4, &c__0, &c_false);

/*       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL. */
/*       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE */
/*          SHOULD BE PRINTED. */

/*       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN */
/*          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE, */
/*          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2. */

    if (*nerr < -9999999 || *nerr > 99999999 || *nerr == 0 || *level < -1 || *
	    level > 2) {
	xerprn_(" ***", &c_n1, "FATAL ERROR IN...$$ XERMSG -- INVALID ERROR "
		"NUMBER OR LEVEL$$ JOB ABORT DUE TO FATAL ERROR.", &c__72, (
		ftnlen)4, (ftnlen)91);
	xersve_(" ", " ", " ", &c__0, &c__0, &c__0, &kdummy, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
	xerhlt_(" ***XERMSG -- INVALID INPUT", (ftnlen)27);
	return 0;
    }

/*       RECORD THE MESSAGE. */

    i__ = j4save_(&c__1, nerr, &c_true);
    xersve_(librar, subrou, messg, &c__1, nerr, level, &kount, librar_len, 
	    subrou_len, messg_len);

/*       HANDLE PRINT-ONCE WARNING MESSAGES. */

    if (*level == -1 && kount > 1) {
	return 0;
    }

/*       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG. */

    s_copy(xlibr, librar, (ftnlen)8, librar_len);
    s_copy(xsubr, subrou, (ftnlen)8, subrou_len);
    s_copy(lfirst, messg, (ftnlen)20, messg_len);
    lerr = *nerr;
    llevel = *level;
    xercnt_(xlibr, xsubr, lfirst, &lerr, &llevel, &lkntrl, (ftnlen)8, (ftnlen)
	    8, (ftnlen)20);

/* Computing MAX */
    i__1 = -2, i__2 = min(2,lkntrl);
    lkntrl = max(i__1,i__2);
    mkntrl = abs(lkntrl);

/*       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS */
/*       ZERO AND THE ERROR IS NOT FATAL. */

    if (*level < 2 && lkntrl == 0) {
	goto L30;
    }
    if (*level == 0 && kount > maxmes) {
	goto L30;
    }
    if (*level == 1 && kount > maxmes && mkntrl == 1) {
	goto L30;
    }
    if (*level == 2 && kount > max(1,maxmes)) {
	goto L30;
    }

/*       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A */
/*       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS) */
/*       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG */
/*       IS NOT ZERO. */

    if (lkntrl != 0) {
	s_copy(temp, "MESSAGE FROM ROUTINE ", (ftnlen)21, (ftnlen)21);
/* Computing MIN */
	i__1 = i_len(subrou, subrou_len);
	i__ = min(i__1,16);
	s_copy(temp + 21, subrou, i__, i__);
	i__1 = i__ + 21;
	s_copy(temp + i__1, " IN LIBRARY ", i__ + 33 - i__1, (ftnlen)12);
	ltemp = i__ + 33;
/* Computing MIN */
	i__1 = i_len(librar, librar_len);
	i__ = min(i__1,16);
	i__1 = ltemp;
	s_copy(temp + i__1, librar, ltemp + i__ - i__1, i__);
	i__1 = ltemp + i__;
	s_copy(temp + i__1, ".", ltemp + i__ + 1 - i__1, (ftnlen)1);
	ltemp = ltemp + i__ + 1;
	xerprn_(" ***", &c_n1, temp, &c__72, (ftnlen)4, ltemp);
    }

/*       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE */
/*       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE */
/*       FROM EACH OF THE FOLLOWING THREE OPTIONS. */
/*       1.  LEVEL OF THE MESSAGE */
/*              'INFORMATIVE MESSAGE' */
/*              'POTENTIALLY RECOVERABLE ERROR' */
/*              'FATAL ERROR' */
/*       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE */
/*              'PROG CONTINUES' */
/*              'PROG ABORTED' */
/*       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK */
/*           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS */
/*           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.) */
/*              'TRACEBACK REQUESTED' */
/*              'TRACEBACK NOT REQUESTED' */
/*       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT */
/*       EXCEED 74 CHARACTERS. */
/*       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED. */

    if (lkntrl > 0) {

/*       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL. */

	if (*level <= 0) {
	    s_copy(temp, "INFORMATIVE MESSAGE,", (ftnlen)20, (ftnlen)20);
	    ltemp = 20;
	} else if (*level == 1) {
	    s_copy(temp, "POTENTIALLY RECOVERABLE ERROR,", (ftnlen)30, (
		    ftnlen)30);
	    ltemp = 30;
	} else {
	    s_copy(temp, "FATAL ERROR,", (ftnlen)12, (ftnlen)12);
	    ltemp = 12;
	}

/*       THEN WHETHER THE PROGRAM WILL CONTINUE. */

	if (mkntrl == 2 && *level >= 1 || mkntrl == 1 && *level == 2) {
	    i__1 = ltemp;
	    s_copy(temp + i__1, " PROG ABORTED,", ltemp + 14 - i__1, (ftnlen)
		    14);
	    ltemp += 14;
	} else {
	    i__1 = ltemp;
	    s_copy(temp + i__1, " PROG CONTINUES,", ltemp + 16 - i__1, (
		    ftnlen)16);
	    ltemp += 16;
	}

/*       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK. */

	if (lkntrl > 0) {
	    i__1 = ltemp;
	    s_copy(temp + i__1, " TRACEBACK REQUESTED", ltemp + 20 - i__1, (
		    ftnlen)20);
	    ltemp += 20;
	} else {
	    i__1 = ltemp;
	    s_copy(temp + i__1, " TRACEBACK NOT REQUESTED", ltemp + 24 - i__1,
		     (ftnlen)24);
	    ltemp += 24;
	}
	xerprn_(" ***", &c_n1, temp, &c__72, (ftnlen)4, ltemp);
    }

/*       NOW SEND OUT THE MESSAGE. */

    xerprn_(" *  ", &c_n1, messg, &c__72, (ftnlen)4, messg_len);

/*       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A */
/*          TRACEBACK. */

    if (lkntrl > 0) {
	s_wsfi(&io___329);
	do_fio(&c__1, (char *)&(*nerr), (ftnlen)sizeof(integer));
	e_wsfi();
	for (i__ = 16; i__ <= 22; ++i__) {
	    if (*(unsigned char *)&temp[i__ - 1] != ' ') {
		goto L20;
	    }
/* L10: */
	}

L20:
/* Writing concatenation */
	i__3[0] = 15, a__1[0] = temp;
	i__3[1] = 23 - (i__ - 1), a__1[1] = temp + (i__ - 1);
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)87);
	xerprn_(" *  ", &c_n1, ch__1, &c__72, (ftnlen)4, 23 - (i__ - 1) + 15);
	fdump_();
    }

/*       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE. */

    if (lkntrl != 0) {
	xerprn_(" *  ", &c_n1, " ", &c__72, (ftnlen)4, (ftnlen)1);
	xerprn_(" ***", &c_n1, "END OF MESSAGE", &c__72, (ftnlen)4, (ftnlen)
		14);
	xerprn_("    ", &c__0, " ", &c__72, (ftnlen)4, (ftnlen)1);
    }

/*       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE */
/*       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN. */

L30:
    if (*level <= 0 || *level == 1 && mkntrl <= 1) {
	return 0;
    }

/*       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A */
/*       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR */
/*       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT. */

    if (lkntrl > 0 && kount < max(1,maxmes)) {
	if (*level == 1) {
	    xerprn_(" ***", &c_n1, "JOB ABORT DUE TO UNRECOVERED ERROR.", &
		    c__72, (ftnlen)4, (ftnlen)35);
	} else {
	    xerprn_(" ***", &c_n1, "JOB ABORT DUE TO FATAL ERROR.", &c__72, (
		    ftnlen)4, (ftnlen)29);
	}
	xersve_(" ", " ", " ", &c_n1, &c__0, &c__0, &kdummy, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
	xerhlt_(" ", (ftnlen)1);
    } else {
	xerhlt_(messg, messg_len);
    }
    return 0;
} /* xermsg_ */

/* DECK XERPRN */
/* Subroutine */ int xerprn_(char *prefix, integer *npref, char *messg, 
	integer *nwrap, ftnlen prefix_len, ftnlen messg_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     i_indx(char *, char *, ftnlen, ftnlen), s_cmp(char *, char *, 
	    ftnlen, ftnlen);

    /* Local variables */
    static integer i__, n, iu[5];
    static char cbuff[148];
    static integer lpref, nextc, lwrap, nunit;
    extern integer i1mach_(integer *);
    static integer lpiece, idelta, lenmsg;
    extern /* Subroutine */ int xgetua_(integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___338 = { 0, 0, 0, "(A)", 0 };
    static cilist io___342 = { 0, 0, 0, "(A)", 0 };


/* ***BEGIN PROLOGUE  XERPRN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Print error messages processed by XERMSG. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERPRN-A) */
/* ***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR */
/* ***AUTHOR  Fong, Kirby, (NMFECC at LLNL) */
/* ***DESCRIPTION */

/* This routine sends one or more lines to each of the (up to five) */
/* logical units to which error messages are to be sent.  This routine */
/* is called several times by XERMSG, sometimes with a single line to */
/* print and sometimes with a (potentially very long) message that may */
/* wrap around into multiple lines. */

/* PREFIX  Input argument of type CHARACTER.  This argument contains */
/*         characters to be put at the beginning of each line before */
/*         the body of the message.  No more than 16 characters of */
/*         PREFIX will be used. */

/* NPREF   Input argument of type INTEGER.  This argument is the number */
/*         of characters to use from PREFIX.  If it is negative, the */
/*         intrinsic function LEN is used to determine its length.  If */
/*         it is zero, PREFIX is not used.  If it exceeds 16 or if */
/*         LEN(PREFIX) exceeds 16, only the first 16 characters will be */
/*         used.  If NPREF is positive and the length of PREFIX is less */
/*         than NPREF, a copy of PREFIX extended with blanks to length */
/*         NPREF will be used. */

/* MESSG   Input argument of type CHARACTER.  This is the text of a */
/*         message to be printed.  If it is a long message, it will be */
/*         broken into pieces for printing on multiple lines.  Each line */
/*         will start with the appropriate prefix and be followed by a */
/*         piece of the message.  NWRAP is the number of characters per */
/*         piece; that is, after each NWRAP characters, we break and */
/*         start a new line.  In addition the characters '$$' embedded */
/*         in MESSG are a sentinel for a new line.  The counting of */
/*         characters up to NWRAP starts over for each new line.  The */
/*         value of NWRAP typically used by XERMSG is 72 since many */
/*         older error messages in the SLATEC Library are laid out to */
/*         rely on wrap-around every 72 characters. */

/* NWRAP   Input argument of type INTEGER.  This gives the maximum size */
/*         piece into which to break MESSG for printing on multiple */
/*         lines.  An embedded '$$' ends a line, and the count restarts */
/*         at the following character.  If a line break does not occur */
/*         on a blank (it would split a word) that word is moved to the */
/*         next line.  Values of NWRAP less than 16 will be treated as */
/*         16.  Values of NWRAP greater than 132 will be treated as 132. */
/*         The actual line length will be NPREF + NWRAP after NPREF has */
/*         been adjusted to fall between 0 and 16 and NWRAP has been */
/*         adjusted to fall between 16 and 132. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  I1MACH, XGETUA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   880621  DATE WRITTEN */
/*   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF */
/*           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK */
/*           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE */
/*           SLASH CHARACTER IN FORMAT STATEMENTS. */
/*   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO */
/*           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK */
/*           LINES TO BE PRINTED. */
/*   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF */
/*           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH. */
/*   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH. */
/*   891214  Prologue converted to Version 4.0 format.  (WRB) */
/*   900510  Added code to break messages between words.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERPRN */
/* ***FIRST EXECUTABLE STATEMENT  XERPRN */
    xgetua_(iu, &nunit);

/*       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD */
/*       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD */
/*       ERROR MESSAGE UNIT. */

    n = i1mach_(&c__4);
    i__1 = nunit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iu[i__ - 1] == 0) {
	    iu[i__ - 1] = n;
	}
/* L10: */
    }

/*       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE */
/*       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING */
/*       THE REST OF THIS ROUTINE. */

    if (*npref < 0) {
	lpref = i_len(prefix, prefix_len);
    } else {
	lpref = *npref;
    }
    lpref = min(16,lpref);
    if (lpref != 0) {
	s_copy(cbuff, prefix, lpref, prefix_len);
    }

/*       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE */
/*       TIME FROM MESSG TO PRINT ON ONE LINE. */

/* Computing MAX */
    i__1 = 16, i__2 = min(132,*nwrap);
    lwrap = max(i__1,i__2);

/*       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS. */

    lenmsg = i_len(messg, messg_len);
    n = lenmsg;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&messg[lenmsg - 1] != ' ') {
	    goto L30;
	}
	--lenmsg;
/* L20: */
    }
L30:

/*       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE. */

    if (lenmsg == 0) {
	i__1 = lpref;
	s_copy(cbuff + i__1, " ", lpref + 1 - i__1, (ftnlen)1);
	i__1 = nunit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___338.ciunit = iu[i__ - 1];
	    s_wsfe(&io___338);
	    do_fio(&c__1, cbuff, lpref + 1);
	    e_wsfe();
/* L40: */
	}
	return 0;
    }

/*       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING */
/*       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL. */
/*       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT. */
/*       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED. */

/*       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE */
/*       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE */
/*       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH */
/*       OF THE SECOND ARGUMENT. */

/*       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE */
/*       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER */
/*       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT */
/*       POSITION NEXTC. */

/*       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE */
/*                       REMAINDER OF THE CHARACTER STRING.  LPIECE */
/*                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC, */
/*                       WHICHEVER IS LESS. */

/*       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC: */
/*                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE */
/*                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY */
/*                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION */
/*                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF */
/*                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE */
/*                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC */
/*                       SHOULD BE INCREMENTED BY 2. */

/*       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP. */

/*       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1 */
/*                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS */
/*                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ. */
/*                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY */
/*                       AT THE END OF A LINE. */

    nextc = 1;
L50:
    lpiece = i_indx(messg + (nextc - 1), "$$", lenmsg - (nextc - 1), (ftnlen)
	    2);
    if (lpiece == 0) {

/*       THERE WAS NO NEW LINE SENTINEL FOUND. */

	idelta = 0;
/* Computing MIN */
	i__1 = lwrap, i__2 = lenmsg + 1 - nextc;
	lpiece = min(i__1,i__2);
	if (lpiece < lenmsg + 1 - nextc) {
	    for (i__ = lpiece + 1; i__ >= 2; --i__) {
		i__1 = nextc + i__ - 2;
		if (s_cmp(messg + i__1, " ", nextc + i__ - 1 - i__1, (ftnlen)
			1) == 0) {
		    lpiece = i__ - 1;
		    idelta = 1;
		    goto L54;
		}
/* L52: */
	    }
	}
L54:
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
	nextc = nextc + lpiece + idelta;
    } else if (lpiece == 1) {

/*       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1). */
/*       DON'T PRINT A BLANK LINE. */

	nextc += 2;
	goto L50;
    } else if (lpiece > lwrap + 1) {

/*       LPIECE SHOULD BE SET DOWN TO LWRAP. */

	idelta = 0;
	lpiece = lwrap;
	for (i__ = lpiece + 1; i__ >= 2; --i__) {
	    i__1 = nextc + i__ - 2;
	    if (s_cmp(messg + i__1, " ", nextc + i__ - 1 - i__1, (ftnlen)1) ==
		     0) {
		lpiece = i__ - 1;
		idelta = 1;
		goto L58;
	    }
/* L56: */
	}
L58:
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
	nextc = nextc + lpiece + idelta;
    } else {

/*       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1. */
/*       WE SHOULD DECREMENT LPIECE BY ONE. */

	--lpiece;
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
	nextc = nextc + lpiece + 2;
    }

/*       PRINT */

    i__1 = nunit;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___342.ciunit = iu[i__ - 1];
	s_wsfe(&io___342);
	do_fio(&c__1, cbuff, lpref + lpiece);
	e_wsfe();
/* L60: */
    }

    if (nextc <= lenmsg) {
	goto L50;
    }
    return 0;
} /* xerprn_ */

/* DECK XERSVE */
/* Subroutine */ int xersve_(char *librar, char *subrou, char *messg, integer 
	*kflag, integer *nerr, integer *level, integer *icount, ftnlen 
	librar_len, ftnlen subrou_len, ftnlen messg_len)
{
    /* Initialized data */

    static integer kountx = 0;
    static integer nmsg = 0;

    /* Format strings */
    static char fmt_9000[] = "(\0020          ERROR MESSAGE SUMMARY\002/\002"
	    " LIBRARY    SUBROUTINE MESSAGE START             NERR\002,\002  "
	    "   LEVEL     COUNT\002)";
    static char fmt_9010[] = "(1x,a,3x,a,3x,a,3i10)";
    static char fmt_9020[] = "(\0020OTHER ERRORS NOT INDIVIDUALLY TABULATED "
	    "= \002,i10)";
    static char fmt_9030[] = "(1x)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static char lib[8], mes[20], sub[8];
    static integer lun[5], iunit, kunit, nunit, kount[10];
    extern integer i1mach_(integer *);
    static char libtab[8*10], mestab[20*10];
    static integer nertab[10], levtab[10];
    static char subtab[8*10];
    extern /* Subroutine */ int xgetua_(integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___349 = { 0, 0, 0, fmt_9000, 0 };
    static cilist io___351 = { 0, 0, 0, fmt_9010, 0 };
    static cilist io___358 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io___359 = { 0, 0, 0, fmt_9030, 0 };


/* ***BEGIN PROLOGUE  XERSVE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Record that an error has occurred. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3 */
/* ***TYPE      ALL (XERSVE-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/* *Usage: */

/*        INTEGER  KFLAG, NERR, LEVEL, ICOUNT */
/*        CHARACTER * (len) LIBRAR, SUBROU, MESSG */

/*        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT) */

/* *Arguments: */

/*        LIBRAR :IN    is the library that the message is from. */
/*        SUBROU :IN    is the subroutine that the message is from. */
/*        MESSG  :IN    is the message to be saved. */
/*        KFLAG  :IN    indicates the action to be performed. */
/*                      when KFLAG > 0, the message in MESSG is saved. */
/*                      when KFLAG=0 the tables will be dumped and */
/*                      cleared. */
/*                      when KFLAG < 0, the tables will be dumped and */
/*                      not cleared. */
/*        NERR   :IN    is the error number. */
/*        LEVEL  :IN    is the error severity. */
/*        ICOUNT :OUT   the number of times this message has been seen, */
/*                      or zero if the table has overflowed and does not */
/*                      contain this message specifically.  When KFLAG=0, */
/*                      ICOUNT will not be altered. */

/* *Description: */

/*   Record that this error occurred and possibly dump and clear the */
/*   tables. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  I1MACH, XGETUA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800319  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900413  Routine modified to remove reference to KFLAG.  (WRB) */
/*   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling */
/*           sequence, use IF-THEN-ELSE, make number of saved entries */
/*           easily changeable, changed routine name from XERSAV to */
/*           XERSVE.  (RWC) */
/*   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERSVE */
/* ***FIRST EXECUTABLE STATEMENT  XERSVE */

    if (*kflag <= 0) {

/*        Dump the table. */

	if (nmsg == 0) {
	    return 0;
	}

/*        Print to each unit. */

	xgetua_(lun, &nunit);
	i__1 = nunit;
	for (kunit = 1; kunit <= i__1; ++kunit) {
	    iunit = lun[kunit - 1];
	    if (iunit == 0) {
		iunit = i1mach_(&c__4);
	    }

/*           Print the table header. */

	    io___349.ciunit = iunit;
	    s_wsfe(&io___349);
	    e_wsfe();

/*           Print body of table. */

	    i__2 = nmsg;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		io___351.ciunit = iunit;
		s_wsfe(&io___351);
		do_fio(&c__1, libtab + (i__ - 1 << 3), (ftnlen)8);
		do_fio(&c__1, subtab + (i__ - 1 << 3), (ftnlen)8);
		do_fio(&c__1, mestab + (i__ - 1) * 20, (ftnlen)20);
		do_fio(&c__1, (char *)&nertab[i__ - 1], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&levtab[i__ - 1], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&kount[i__ - 1], (ftnlen)sizeof(integer)
			);
		e_wsfe();
/* L10: */
	    }

/*           Print number of other errors. */

	    if (kountx != 0) {
		io___358.ciunit = iunit;
		s_wsfe(&io___358);
		do_fio(&c__1, (char *)&kountx, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    io___359.ciunit = iunit;
	    s_wsfe(&io___359);
	    e_wsfe();
/* L20: */
	}

/*        Clear the error tables. */

	if (*kflag == 0) {
	    nmsg = 0;
	    kountx = 0;
	}
    } else {

/*        PROCESS A MESSAGE... */
/*        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG, */
/*        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL. */

	s_copy(lib, librar, (ftnlen)8, librar_len);
	s_copy(sub, subrou, (ftnlen)8, subrou_len);
	s_copy(mes, messg, (ftnlen)20, messg_len);
	i__1 = nmsg;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (s_cmp(lib, libtab + (i__ - 1 << 3), (ftnlen)8, (ftnlen)8) == 
		    0 && s_cmp(sub, subtab + (i__ - 1 << 3), (ftnlen)8, (
		    ftnlen)8) == 0 && s_cmp(mes, mestab + (i__ - 1) * 20, (
		    ftnlen)20, (ftnlen)20) == 0 && *nerr == nertab[i__ - 1] &&
		     *level == levtab[i__ - 1]) {
		++kount[i__ - 1];
		*icount = kount[i__ - 1];
		return 0;
	    }
/* L30: */
	}

	if (nmsg < 10) {

/*           Empty slot found for new message. */

	    ++nmsg;
	    s_copy(libtab + (i__ - 1 << 3), lib, (ftnlen)8, (ftnlen)8);
	    s_copy(subtab + (i__ - 1 << 3), sub, (ftnlen)8, (ftnlen)8);
	    s_copy(mestab + (i__ - 1) * 20, mes, (ftnlen)20, (ftnlen)20);
	    nertab[i__ - 1] = *nerr;
	    levtab[i__ - 1] = *level;
	    kount[i__ - 1] = 1;
	    *icount = 1;
	} else {

/*           Table is full. */

	    ++kountx;
	    *icount = 0;
	}
    }
    return 0;

/*     Formats. */

} /* xersve_ */

/* DECK XGETUA */
/* Subroutine */ int xgetua_(integer *iunita, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, index;
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XGETUA */
/* ***PURPOSE  Return unit number(s) to which error messages are being */
/*            sent. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XGETUA-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XGETUA may be called to determine the unit number or numbers */
/*        to which error messages are being sent. */
/*        These unit numbers may have been set by a call to XSETUN, */
/*        or a call to XSETUA, or may be a default value. */

/*     Description of Parameters */
/*      --Output-- */
/*        IUNIT - an array of one to five unit numbers, depending */
/*                on the value of N.  A value of zero refers to the */
/*                default unit, as defined by the I1MACH machine */
/*                constant routine.  Only IUNIT(1),...,IUNIT(N) are */
/*                defined by XGETUA.  The values of IUNIT(N+1),..., */
/*                IUNIT(5) are not defined (for N .LT. 5) or altered */
/*                in any way by XGETUA. */
/*        N     - the number of units to which copies of the */
/*                error messages are being sent.  N will be in the */
/*                range from 1 to 5. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XGETUA */
/* ***FIRST EXECUTABLE STATEMENT  XGETUA */
    /* Parameter adjustments */
    --iunita;

    /* Function Body */
    *n = j4save_(&c__5, &c__0, &c_false);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	index = i__ + 4;
	if (i__ == 1) {
	    index = 3;
	}
	iunita[i__] = j4save_(&index, &c__0, &c_false);
/* L30: */
    }
    return 0;
} /* xgetua_ */

