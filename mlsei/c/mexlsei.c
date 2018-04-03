#include "mex.h"

/*
 * mexlsei.c - gateway to dlsei 
 *
 * gateway to dlsei fortan program converted to c source using f2c
 *
 * This is a MEX-file for MATLAB.
 *
 */

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *W, *PRGOPT, *WS, *IP;
  long MDW, ME, MA, MG, N, MODE;
  double  *X, *RNORME, *RNORML;
  mwSize matm, matn;

/*
  printf("long int size: %zu", sizeof(long));
*/
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs != 9) 
    mexErrMsgTxt("Nine inputs required.");
  if(nlhs != 4) 
    mexErrMsgTxt("Four outputs required.");
  
  /* check the input arguments */
  if( mxIsComplex(prhs[1]) || mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) {
    mexErrMsgTxt("Input mdw must be a scalar.");
  }
  
  if( mxIsComplex(prhs[2]) || mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
    mexErrMsgTxt("Input me must be a scalar.");
  }
  
  if( mxIsComplex(prhs[3]) || mxGetN(prhs[3])*mxGetM(prhs[3])!=1 ) {
    mexErrMsgTxt("Input ma must be a scalar.");
  }
  
  if( mxIsComplex(prhs[4]) || mxGetN(prhs[4])*mxGetM(prhs[4])!=1 ) {
    mexErrMsgTxt("Input mg must be a scalar.");
  }
  
  if( mxIsComplex(prhs[5]) || mxGetN(prhs[5])*mxGetM(prhs[5])!=1 ) {
    mexErrMsgTxt("Input n must be a scalar.");
  }

  
  /*  create a pointer to the input vector W */
  W = mxGetPr(prhs[0]);
  
  /*  get the scalar input MDW */
  MDW = (long) mxGetScalar(prhs[1]);
  
    /*  get the scalar input ME */
  ME = (long) mxGetScalar(prhs[2]);
  
    /*  get the scalar input MA */
  MA = (long) mxGetScalar(prhs[3]);
  
    /*  get the scalar input MG */
  MG = (long) mxGetScalar(prhs[4]);
  
    /*  get the scalar input N */
  N = (long) mxGetScalar(prhs[5]);

    /*  get the vector input PRGOPT */
  PRGOPT = mxGetPr(prhs[6]);
  
  MODE = 0;

    /*  get the vector input WS */  
  WS = mxGetPr(prhs[7]);
  
    /*  get the vector input IP */
  IP = mxGetPr(prhs[8]);

  plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[1] = mxCreateDoubleScalar(0.0);
  plhs[2] = mxCreateDoubleScalar(0.0);

  X = mxGetPr(plhs[0]);
  RNORME = mxGetPr(plhs[1]);
  RNORML = mxGetPr(plhs[2]);

  /*  call the C subroutine */
  dlsei_(W, &MDW, &ME, &MA, &MG, &N, PRGOPT, X, RNORME, RNORML, &MODE, WS, IP);
  
  plhs[3] = mxCreateDoubleScalar(MODE);

}
