#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif

#include "mex.h"
#include <cmath>
#include <cstring>
/* mxGetPr() used to get DOUBLE arrays... mxGetData() used for ints shorts*/

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  double *Gt, *Ht, *Gi, *Hi;
  int bins;
  double *t_diff, *i_diff;
  
  int mrows_Gt,ncols_Gt, mrows_Ht,ncols_Ht, \
    mrows_Gi,ncols_Gi, mrows_Hi,ncols_Hi;
  int returnsizedims[] = {1,1};

  /*  Check for proper number of arguments. */
  if (nrhs != 5) 
    mexErrMsgIdAndTxt("distcalcmh:NotEnoughInputs","distcalcmh - Five inputs required.");
  if (nlhs != 2) 
    mexErrMsgIdAndTxt("distcalcmh:NotEnoughOutputs","distcalcmh - Two outputs required.");
  
  /* Create a pointer to the input matrix data. */
  Gt = mxGetPr(prhs[0]);
  Ht = mxGetPr(prhs[1]);
  Gi = mxGetPr(prhs[2]);
  Hi = mxGetPr(prhs[3]);
  /* Get the scalar input */
  bins = mxGetScalar(prhs[4]);

  /* Get the dimensions of the matrix input. */
  mrows_Gt = mxGetM(prhs[0]); // number of pixels in 1
  ncols_Gt = mxGetN(prhs[0]); // number of features
  mrows_Ht = mxGetM(prhs[1]); // number of pixels in 2
  ncols_Ht = mxGetN(prhs[1]);
  mrows_Gi = mxGetM(prhs[2]); // number of pixels in 1
  ncols_Gi = mxGetN(prhs[2]); // number of channels
  mrows_Hi = mxGetM(prhs[3]); // number of pixels in 2
  ncols_Hi = mxGetN(prhs[3]);
  
  plhs[0] = mxCreateNumericArray(2, returnsizedims , mxDOUBLE_CLASS,mxREAL);
  plhs[1] = mxCreateNumericArray(2, returnsizedims , mxDOUBLE_CLASS,mxREAL);
  t_diff = mxGetPr(plhs[0]);
  i_diff = mxGetPr(plhs[1]);
  
  *t_diff = 0.0; *i_diff = 0.0;
  
  // create histograms
  double *h1 = new double [ncols_Gt * bins];
  memset( (void *)h1, 0, ncols_Gt * bins * sizeof(double) );
  double *h2 = new double [ncols_Gt * bins];
  memset( (void *)h2, 0, ncols_Gt * bins * sizeof(double) );
  double *c_h1 = new double [ncols_Gt * bins];
  memset( (void *)c_h1, 0, ncols_Gt * bins * sizeof(double) );
  double *c_h2 = new double [ncols_Gt * bins];
  memset( (void *)c_h2, 0, ncols_Gt * bins * sizeof(double) );

  double *c_h_diff = new double [ncols_Gt * bins];
  memset( (void *)c_h_diff, 0, ncols_Gt * bins * sizeof(double) );
  double *d_vec = new double [ncols_Gt];

  int k,n,b;
  
  // create histo of region 1 wavelet features
  for(k=0; k < ncols_Gt; k++ ) 
  {
      for(int n=0; n < mrows_Gt; n++ )
      {
          //histObj1(r1(n,k),k) = histObj1(r1(n,k),k) +1;
          h1[ (k*bins) + (int)( Gt[ (k*mrows_Gt) + n ] ) - 1 ] += 1.0;
      }
  }
  // create histo of region 2 wavelet features
  for(k=0; k < ncols_Ht; k++ ) 
  {
      for(n=0; n < mrows_Ht; n++ ) 
      {
          h2[ (k*bins) + (int)( Ht[ (k*mrows_Ht) + n ] ) - 1 ] += 1.0;
      }
  }  
  // normalise
  for(k=0; k < ncols_Gt; k++ ) 
  {
      for(b=0; b < bins; b++ )
      {
          h1[ (k*bins) + b ] /= mrows_Gt;
          h2[ (k*bins) + b ] /= mrows_Ht;
      }
  }
  
  // cumulative sum histos (sum down rows) and do EMD
  double m=0.0;
  for(k=0; k < ncols_Gt; k++ ) 
  {
      // init first value
      c_h1[(k*bins)] = h1[(k*bins)];
      c_h2[(k*bins)] = h2[(k*bins)];
      c_h_diff[(k*bins)] = abs(h1[(k*bins)] - h2[(k*bins)]);

      d_vec[k] = 0.0;
      
      for(b=1; b < bins; b++ )
      {
          //histObj1 = cumsum(histObj1,1);
          //histObj2 = cumsum(histObj2,1);
          c_h1[(k*bins) + b] = c_h1[(k*bins) + (b-1)] + h1[(k*bins) + b];
          c_h2[(k*bins) + b] = c_h2[(k*bins) + (b-1)] + h2[(k*bins) + b];
          
          //histObj1 = abs(histObj1-histObj2);
          c_h_diff[(k*bins) + b] = abs(c_h1[(k*bins) + b] - c_h2[(k*bins) + b]);
          
          //dvec = sum(histObj1,1)./(binsperdim-1);
          d_vec[k] += c_h_diff[(k*bins) + b];
      }
      d_vec[k] /= (bins - 1);
      
      // max(dvec)
      if( d_vec[k] > m )
          m = d_vec[k];
  }
  m *= 2.0;
  *t_diff = (m < 1.0) ? (m) : (1.0);
  
  
  // delete
  delete h1;
  delete h2;  
  delete c_h1;
  delete c_h2;  
  delete c_h_diff;
  delete d_vec;
  
  // ************************
  // Repeat for intensity data
  // create histograms
  h1 = new double [ncols_Gi * bins];
  memset( (void *)h1, 0, ncols_Gi * bins * sizeof(double) );
  h2 = new double [ncols_Gi * bins];
  memset( (void *)h2, 0, ncols_Gi * bins * sizeof(double) );
  c_h1 = new double [ncols_Gi * bins];
  memset( (void *)c_h1, 0, ncols_Gi * bins * sizeof(double) );
  c_h2 = new double [ncols_Gi * bins];
  memset( (void *)c_h2, 0, ncols_Gi * bins * sizeof(double) );

  c_h_diff = new double [ncols_Gi * bins];
  memset( (void *)c_h_diff, 0, ncols_Gi * bins * sizeof(double) );
  d_vec = new double [ncols_Gi];

  
  // create histo of region 1 wavelet features
  for(k=0; k < ncols_Gi; k++ ) 
  {
      for(n=0; n < mrows_Gi; n++ )
      {
          //histObj1(r1(n,k),k) = histObj1(r1(n,k),k) +1;
          h1[ (k*bins) + (int)( Gi[ (k*mrows_Gi) + n ] ) - 1 ] += 1.0;
      }
  }
  // create histo of region 2 wavelet features
  for(k=0; k < ncols_Hi; k++ ) 
  {
      for(n=0; n < mrows_Hi; n++ ) 
      {
          h2[ (k*bins) + (int)( Hi[ (k*mrows_Hi) + n ] ) - 1 ] += 1.0;
      }
  }  
  // normalise
  for(k=0; k < ncols_Gi; k++ ) 
  {
      for(b=0; b < bins; b++ )
      {
          h1[ (k*bins) + b ] /= mrows_Gi;
          h2[ (k*bins) + b ] /= mrows_Hi;
      }
  }
  
  // cumulative sum histos (sum down rows) and do EMD
  m=0.0;
  for(k=0; k < ncols_Gi; k++ ) 
  {
      // init first value
      c_h1[(k*bins)] = h1[(k*bins)];
      c_h2[(k*bins)] = h2[(k*bins)];
      c_h_diff[(k*bins)] = abs(h1[(k*bins)] - h2[(k*bins)]);

      d_vec[k] = 0.0;
      
      for(b=1; b < bins; b++ )
      {
          //histObj1 = cumsum(histObj1,1);
          //histObj2 = cumsum(histObj2,1);
          c_h1[(k*bins) + b] = c_h1[(k*bins) + (b-1)] + h1[(k*bins) + b];
          c_h2[(k*bins) + b] = c_h2[(k*bins) + (b-1)] + h2[(k*bins) + b];
          
          //histObj1 = abs(histObj1-histObj2);
          c_h_diff[(k*bins) + b] = abs(c_h1[(k*bins) + b] - c_h2[(k*bins) + b]);
          
          //dvec = sum(histObj1,1)./(binsperdim-1);
          d_vec[k] += c_h_diff[(k*bins) + b];
      }
      d_vec[k] /= (bins - 1);
      
      // max(dvec)
      if( d_vec[k] > m )
          m = d_vec[k];
  }
  m *= 2.0;
  *i_diff = (m < 1.0) ? (m) : (1.0);
  
  
  // delete
  delete h1;
  delete h2;  
  delete c_h1;
  delete c_h2;  
  delete c_h_diff;
  delete d_vec;
  

}  
