
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <complex>
#include <unistd.h>
#include "fftw3.h"


#include <stdio.h>
#include <string>
#include <iostream>
#include <sys/time.h>
#include <iomanip> 

#include <string.h>
#define WISDOMNAME_SINGLE "myfftwWisdom"
#define WISDOMNAME_DOUBLE "myfftwWisdom_double"

// #define _MYFFTW_VERBOSE

#define REAL double


inline double
gettime(void)
{
  static struct timeval timeS;
  gettimeofday( &timeS, NULL);
  return (timeS.tv_sec + timeS.tv_usec/1000000.0);
}




inline void report_timing( std::string message)
{
  static double startTime = -1.0;
  static double lastTime = -1.0;
  if( startTime < 0.0) {
    startTime = gettime();
    lastTime = startTime;
  }
  double currentTime = gettime();

  std::cerr.setf(std::ios::fixed, std::ios::floatfield);


  std::cerr << std::setw(7)
          << std::setprecision(3) << currentTime - startTime << "s ("
          << std::setw(6) << std::setprecision(4)
          << currentTime - lastTime << "s)  "
          << message << std::endl;
  lastTime = currentTime;
}



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 6 ) {
	printf("\nUsage:\n");
    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}


    ///////////////////////////////////////////// fetching data 

    int pcnt = 0;  
    const mxArray *Measurement;
    Measurement = prhs[pcnt++];       
    std::complex<double> *meas = ( std::complex<double> *) mxGetData(Measurement);
    

    
    const mxArray *Sens;
    Sens = prhs[pcnt++];       
    const int *dims_sens = mxGetDimensions(Sens);
    const int numdim_sens = mxGetNumberOfDimensions(Sens);
    std::complex<double> *sens = ( std::complex<double> *) mxGetData(Sens);
    
    int numsens;
    if (numdim_sens == 4)
        numsens = 1;
    else
        numsens = dims_sens[4];
        
    const int numdim =4;
    const int dims_sz[] = {2, dims_sens[1], dims_sens[2], dims_sens[3] };
    int w = (int)dims_sz[1];
    int h = (int)dims_sz[2];
    int d = (int)dims_sz[3];
    int totsz = w*h*d;
    
 
    const mxArray *Ipk_index;
    Ipk_index = prhs[pcnt++];       
    const int *dims_ipk = mxGetDimensions(Ipk_index);
    double *ipk_index = (double*) mxGetData(Ipk_index);

    const mxArray *Ipk_we;
    Ipk_we = prhs[pcnt++];       
    std::complex<double> *ipk_we = (std::complex<double>*) mxGetData(Ipk_we);
  
    int numP = dims_ipk[0];
    int numK = dims_ipk[1];
    
    int the_index[numP*numK];
    for(int i = 0; i < numP*numK; i++)
        the_index[i] = (int)(ipk_index[i]-1);
    

  
    const mxArray *Dims_pad;
    Dims_pad = prhs[pcnt++];       
    double *dims_pad_d = (double*) mxGetData(Dims_pad);
    int w_pad = (int)dims_pad_d[0];
    int h_pad = (int)dims_pad_d[1];
    int d_pad = (int)dims_pad_d[2];
    int totsz_pad  = w_pad*h_pad*d_pad;
    const int dims_pad[] = {2, w_pad, h_pad, d_pad};
    
 
    
    const mxArray *Params;
    Params = prhs[pcnt++];       
    double *params = (double*) mxGetData(Params);
    int numit = (int) params[0];
    double lambda = params[1];
    double tol = params[2];
    int verbose = (int) params[3];
 
    
    
    
    ///////////////////////////////////////////// init FFTW

    
    int dims_sz_[] = {2, w_pad, h_pad, d_pad};
 
    int rank = 3;

    int fftdim[numdim];
    for (int i = 0; i < numdim;i++)
	fftdim[i] = 0;
    for (int i = 0; i < rank; i++)
	fftdim[ i+1] = 1;

    fftwf_iodim *dims = new fftwf_iodim[rank];
    fftwf_iodim *howmany_dims = new fftwf_iodim[numdim-rank];



    int howmany_rank = 0;
    rank = 0;
    int stride = 1;
    dims_sz_[0] /= 2;
    for (int i = 0; i < numdim; i++)
    {
	if (fftdim[i] == 1)
	{
		dims[rank].n = dims_sz_[i];
		dims[rank].is = stride;
		dims[rank].os = dims[rank].is;
		rank++;
	}
	else
	{
		howmany_dims[howmany_rank].n = dims_sz_[i];
		howmany_dims[howmany_rank].is = stride;
		howmany_dims[howmany_rank].os = howmany_dims[howmany_rank].is;
		howmany_rank++;
	}
	stride *= dims_sz_[i];
	
    }



    unsigned flags =FFTW_MEASURE;



	fftw_forget_wisdom();
	mxArray *wisdom = mexGetVariable("global",WISDOMNAME_DOUBLE);
	if (wisdom != NULL)
	{
		int buflen = mxGetN(wisdom)*sizeof(mxChar)+1;
		char fftwwisdom[buflen];
		mxGetString(wisdom,fftwwisdom,buflen);
		int ret = fftw_import_wisdom_from_string(fftwwisdom);
	}

    
    /////////////////////////////////////// MALLOCs
    
    
    
	plhs[0]             =  mxCreateNumericArray(numdim,dims_sz,mxGetClassID(Sens),mxREAL);
    const mxArray *Tmp1 =  mxCreateNumericArray(numdim,dims_pad,mxGetClassID(Sens),mxREAL);
    const mxArray *Tmp2 =  mxCreateNumericArray(numdim,dims_pad,mxGetClassID(Sens),mxREAL);
    const mxArray *Tmp3 =  mxCreateNumericArray(numdim,dims_sz,mxGetClassID(Sens),mxREAL);
    const mxArray *_D =  mxCreateNumericArray(numdim,dims_sz,mxGetClassID(Sens),mxREAL);
    const mxArray *_Z =  mxCreateNumericArray(numdim,dims_sz,mxGetClassID(Sens),mxREAL);
    const mxArray *_R =  mxCreateNumericArray(numdim,dims_sz,mxGetClassID(Sens),mxREAL);
   
    std::complex<double> *tmp1 = (std::complex<double>*) mxGetData(Tmp1);	
    std::complex<double> *tmp2 = (std::complex<double>*) mxGetData(Tmp2);	
    std::complex<double> *tmp3 = (std::complex<double>*) mxGetData(Tmp3);	
    std::complex<double> *_d = (std::complex<double>*) mxGetData(_D);	
    std::complex<double> *_z = (std::complex<double>*) mxGetData(plhs[0]);	
    std::complex<double> *_r = (std::complex<double>*) mxGetData(_R);	
    
    const fftw_plan plan_fw = fftw_plan_guru_dft(rank,(const fftwf_iodim*) dims, howmany_rank, (const fftwf_iodim*)howmany_dims,(fftw_complex*) tmp1,(fftw_complex*)tmp2, -1, flags| FFTW_PRESERVE_INPUT);
    const fftw_plan plan_bw = fftw_plan_guru_dft(rank,(const fftwf_iodim*) dims, howmany_rank, (const fftwf_iodim*)howmany_dims,(fftw_complex*) tmp1,(fftw_complex*)tmp2, 1, flags| FFTW_PRESERVE_INPUT);

    
    
    
    /////////////////////////////////////////////////////// init CG
    
    
    
     double normrr = 0;
     double dAAd = 0;
     double alpha = 0;
     double normrr2 = 0;
     double beta = 0;
    
                   
    // backproject measurement -- x=A'b
    for (int i = 0; i < numsens; i++)
    {
        memset(tmp1,0, sizeof(std::complex<double>)*totsz_pad); 
       
        for (int k = 0; k < numK; k++)
        {
            std::complex<double> val = meas[i*numK+k];
            for (int p = 0; p < numP; p++)
            {
                 int idx = the_index[numP*k + p];
                 tmp1[idx] += val*std::conj(ipk_we[numP*k + p]) ;

            }
        }
        fftw_execute_dft(plan_bw,(fftw_complex*)tmp1,(fftw_complex*)tmp2);

        int k = 0;
        for (int z = 0; z < d; z++)
            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x ++)
                {
                     _r[k] += tmp2[z*w_pad*h_pad+y*w_pad+x]*std::conj(sens[k+totsz*i]);         
                    k++;
                }

     }


    
        
     for (int k = 0 ; k < totsz;k++)
          _d[k] = _r[k];

     
     normrr = 0;
     for (int k = 0 ; k < totsz;k++)
          normrr += real(_r[k])*real(_r[k]) + imag(_r[k])*imag(_r[k]);            
  
     double normrr0 = normrr;
  //   mexPrintf("first residual: %f\n",normrr);
    
    
    ////////////////////////////////////////////////////////////// start CG
     
    for (int it = 0; it < numit; it++)
    {
         
 //   Qd = (A'*(A*d)) + l2*d;
        
        for (int k = 0 ; k < totsz;k++)
                tmp3[k] = lambda*_d[k];

        for (int i = 0; i < numsens; i++)
        {

            memset(tmp1,0, sizeof(std::complex<double>)*totsz_pad);
            int k = 0;
            for (int z = 0; z < d; z++)
                for (int y = 0; y < h; y++)
                    for (int x = 0; x < w; x ++)
                    {
                        tmp1[z*w_pad*h_pad+y*w_pad+x] = _d[k]*sens[k+totsz*i];
                        k++;
                    }
            
            fftw_execute_dft(plan_fw,(fftw_complex*)tmp1,(fftw_complex*)tmp2);
            
            memset(tmp1,0, sizeof(std::complex<double>)*totsz_pad);
            for (int k = 0; k < numK; k++)
            {
                std::complex<double> val = 0;
                for (int p = 0; p < numP; p++)
                {
                    int idx = the_index[numP*k + p];
                    val += tmp2[idx]*ipk_we[numP*k + p];
                }
                for (int p = 0; p < numP; p++)
                {
                     int idx = the_index[numP*k + p];
                     tmp1[idx] += val*std::conj(ipk_we[numP*k + p]) ;

                }
            }
            fftw_execute_dft(plan_bw,(fftw_complex*)tmp1,(fftw_complex*)tmp2);
            
            k = 0;
            for (int z = 0; z < d; z++)
                for (int y = 0; y < h; y++)
                    for (int x = 0; x < w; x ++)
                    {
                        tmp3[k] += tmp2[z*w_pad*h_pad+y*w_pad+x]*std::conj(sens[k+totsz*i]);        
                        k++;
                    }
            
        }
    
   //    dAAd = d(:)'*Qd(:);
  
        dAAd = 0;
        for (int k = 0 ; k < totsz;k++)
            dAAd += real(_d[k])*real(tmp3[k]) + imag(_d[k])*imag(tmp3[k]);
            
//        alpha = normrr/dAAd;
        alpha = normrr/dAAd;

//    z = z + alpha*d;
//    r = r - alpha*Qd;
        
        for (int k = 0 ; k < totsz;k++)
            _z[k] += alpha*_d[k];
        for (int k = 0 ; k < totsz;k++)
            _r[k] -= alpha*tmp3[k];                
   
//    normrr2 = r(:)'*r(:);
        normrr2 = 0;
        for (int k = 0 ; k < totsz;k++)
              normrr2 += real(_r[k])*real(_r[k]) + imag(_r[k])*imag(_r[k]);            
    	
        beta = normrr2/normrr;
        normrr = normrr2;
   
        for (int k = 0 ; k < totsz;k++)    
          	_d[k] = _r[k] + beta*_d[k];
 
        if (sqrt(normrr/normrr0) < tol)
            break;       
    
        if (verbose ==  1)
            mexPrintf("tol : %f\n",sqrt(normrr/normrr0));
        
        mexPrintf("."); mexEvalString("drawnow");
    
    }
    
    
    
    
    fftw_destroy_plan(plan_fw);
    fftw_destroy_plan(plan_bw);
    
    
    
    mxArray *newwisdom = mxCreateString(fftw_export_wisdom_to_string());
    mexPutVariable("global",WISDOMNAME_DOUBLE,newwisdom);


}
