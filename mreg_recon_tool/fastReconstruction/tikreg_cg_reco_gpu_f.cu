
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <complex>
#include <unistd.h>
#include "fftw3.h"



#include "cufft.h"
#include "cuda_runtime.h"
#include <cuda.h> 
#include <cublas.h>


#include <stdio.h>
#include <string>
#include <iostream>
#include <sys/time.h>


#include <string.h>
#include <sys/time.h>

#define GET_TIME(now) { \
   struct timeval t; \
   gettimeofday(&t, NULL); \
   now = t.tv_sec + t.tv_usec/1000000.0; \
}


#define MAX_BLOCK_SZ 512


#include "tikreg_CG_kernels.cu"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 9 ) {
	printf("\nUsage:\n");
    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}

    //////////////////////////////////// fetching data from MATLAB

     int pcnt = 0;  
    const mxArray *Measurement;
    Measurement = prhs[pcnt++];       
    std::complex<float> *meas = ( std::complex<float> *) mxGetData(Measurement);
    

    
    const mxArray *Sens;
    Sens = prhs[pcnt++];       
    const int *dims_sens = mxGetDimensions(Sens);
    const int numdim_sens = mxGetNumberOfDimensions(Sens);
    std::complex<float> *sens = ( std::complex<float> *) mxGetData(Sens);
    
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
    float *ipk_index = (float*) mxGetData(Ipk_index);

    const mxArray *Ipk_we;
    Ipk_we = prhs[pcnt++];       
    std::complex<float> *ipk_we = (std::complex<float>*) mxGetData(Ipk_we);
  
    int numP = dims_ipk[0];
    int numK = dims_ipk[1];
    
    int the_index[numP*numK];
    for(int i = 0; i < numP*numK; i++)
        the_index[i] = (int)(ipk_index[i]-1);
    
    

  
    const mxArray *Dims_pad;
    Dims_pad = prhs[pcnt++];       
    float *dims_pad_d = (float*) mxGetData(Dims_pad);
    int w_pad = (int)dims_pad_d[0];
    int h_pad = (int)dims_pad_d[1];
    int d_pad = (int)dims_pad_d[2];
    int totsz_pad  = w_pad*h_pad*d_pad;

    
    const mxArray *BPidx;
    BPidx = prhs[pcnt++];       
    int numVox= mxGetM(BPidx);
    int * bpidx = (int*) mxGetData(BPidx);
    
    const mxArray *BPmidx;
    BPmidx = prhs[pcnt++];       
   
    const mxArray *BPweight;
    BPweight = prhs[pcnt++];       
   
    
    const mxArray *Params;
    Params = prhs[pcnt++];       
    float *params = (float*) mxGetData(Params);
    int numit = (int) params[0];
    float lambda = params[1];
    int device_num = (int) params[2];
    float tol = params[3];
    int VERBOSE = (int) params[4];
    
    if (VERBOSE == 1)  
        mexPrintf("gpuDevice: %i  lambda^2: %f\n",device_num,lambda);

   /**************** Init Cuda *****************/
    
    cudaError_t rv; 
    CUdevice dev; 
    
    if (cuCtxGetDevice(&dev) == CUDA_SUCCESS)
    {
    //   CUcontext  pctx ;
    //   cuCtxPopCurrent(&pctx);	      
    }   
    
    mexPrintf("dev:%i\n",dev);
    
//     
//    cudaSetDevice(device_num); 
// 
//    rv = cudaSetDeviceFlags(cudaDeviceMapHost);
//    if (rv != cudaSuccess )
//    {
//       mexPrintf("Call to cudaSetDeviceFlags failed\n");
//       return;
//    } 
//     
//     /******** Allocate mapped tmps for dot product calc **********/
//    
//    int dot_threads = 128;
//    int dot_blocks = w*h*d/dot_threads;
//    
//    float *dot_z_h,*dot_z_d;
//    rv = cudaHostAlloc(&dot_z_h, dot_blocks*sizeof(float), cudaHostAllocMapped);
// 
//    if (rv != cudaSuccess) {
//       mexPrintf("Call to cudaHostAlloc failed: %i\n",rv);
//       return;
//    } 
//    cudaHostGetDevicePointer(&dot_z_d, dot_z_h, 0);
//    cudaMemset(dot_z_d,0,dot_blocks*sizeof(float));
//    
   
    /////////////////////////////////////// MALLOCs
    
    double start,finish;
     
    GET_TIME(start);
    
    cufftComplex *tmp1,*tmp2,*tmp3, *_r, *_d, *_z , *_meas,*_sens, *_ipk_we, *tmpsens;
    int *_the_index;
    cufftHandle            plan;
    
	plhs[0]             =  mxCreateNumericArray(numdim,dims_sz,mxGetClassID(Sens),mxREAL);
     
    std::complex<float> *res = (std::complex<float> *) mxGetData(plhs[0]);
   
    cudaMalloc( (void **) &tmp1,sizeof(cufftComplex)*totsz_pad);
    cudaMalloc( (void **) &tmp2,sizeof(cufftComplex)*totsz_pad);
    cudaMalloc( (void **) &tmp3,sizeof(cufftComplex)*totsz);
    cudaMalloc( (void **) &_r,sizeof(cufftComplex)*totsz);
    cudaMalloc( (void **) &_d,sizeof(cufftComplex)*totsz);
    cudaMalloc( (void **) &_z,sizeof(cufftComplex)*totsz);
    cudaMalloc( (void **) &_meas,sizeof(cufftComplex)*numsens*numK);

    cudaMalloc( (void **) &_sens,sizeof(cufftComplex)*numsens*totsz);
    cudaMalloc( (void **) &_ipk_we,sizeof(cufftComplex)*numP*numK);
    cudaMalloc( (void **) &_the_index,sizeof(int)*numP*numK);
    cudaMalloc( (void **) &tmpsens,sizeof(cufftComplex)*numK);

    cudaThreadSynchronize();
   
    cudaMemset( tmp1,0,sizeof(cufftComplex)*totsz_pad);
    cudaMemset( tmp2,0,sizeof(cufftComplex)*totsz_pad);
    cudaMemset(  tmp3,0,sizeof(cufftComplex)*totsz);
    cudaMemset(  _r,0,sizeof(cufftComplex)*totsz);
    cudaMemset(  _d,0,sizeof(cufftComplex)*totsz);
    cudaMemset(  _z,0,sizeof(cufftComplex)*totsz);
 
     cudaThreadSynchronize();
 
  
     /************** copy data on device **********************/

     cudaMemcpy( _meas, meas, sizeof(cufftComplex)*numsens*numK, cudaMemcpyHostToDevice);
     cudaMemcpy( _sens, sens, sizeof(cufftComplex)*numsens*totsz, cudaMemcpyHostToDevice);
     cudaMemcpy( _ipk_we, ipk_we, sizeof(cufftComplex)*numP*numK, cudaMemcpyHostToDevice);
     cudaMemcpy( _the_index, the_index, sizeof(int)*numP*numK, cudaMemcpyHostToDevice);
   
     cudaMemcpy( ipk_we, _ipk_we, sizeof(cufftComplex)*numP*numK, cudaMemcpyDeviceToHost);
     cudaMemcpy( the_index, _the_index, sizeof(int)*numP*numK, cudaMemcpyDeviceToHost);
     
 
     cudaThreadSynchronize();
    
    if (VERBOSE == 1) 
        mexPrintf("numP: %i  numK: %i whd %i %i %i pad %i %i %i numsens: %i\n",numP,numK,w,h,d,w_pad,h_pad,d_pad,numsens);
            
      
    /************** copy bpidx on device **********************/
    int *_bpmidx;
    cufftComplex *_bpweight;
    int *bpsize = (int*) malloc(sizeof(int)*numVox);
    int *bponset  = (int*) malloc(sizeof(int)*(numVox+1));
    int *_bpsize, *_bponset, *_bpidx;
    bponset[0] = 0;
    for (int j = 0; j < numVox;j++)
    {
        mxArray *Midx = mxGetCell(BPmidx,j);
        bpsize[j] = mxGetM(Midx);
        bponset[j+1] = bponset[j] + bpsize[j];
    }
    
    int *tmp_bpmidx;
    cufftComplex *tmp_bpweight;
    tmp_bpmidx = (int*) malloc(sizeof(int)*bponset[numVox]);
    tmp_bpweight = (cufftComplex*) malloc(sizeof(cufftComplex)*bponset[numVox]);
    if (tmp_bpmidx == 0)
    {
        mexPrintf("out of mem (host)\n");
        return;
    }
    if (tmp_bpweight == 0)
    {
        mexPrintf("out of mem (host)\n");
        return;
    }
    
    for (int j = 0; j < numVox;j++)
    {
        mxArray *Midx = mxGetCell(BPmidx,j);
        mxArray *Weight = mxGetCell(BPweight,j);
        int *midx = (int*)  mxGetData(Midx);
        cufftComplex *bpwei = (cufftComplex*) mxGetData(Weight);
        memcpy(tmp_bpmidx + bponset[j] , midx, sizeof(int)* bpsize[j]);
        memcpy(tmp_bpweight + bponset[j] , bpwei, sizeof(cufftComplex)* bpsize[j]);    
    }
    
    cudaMalloc( (void **) &_bpmidx,sizeof(int)* bponset[numVox]);
    cudaMalloc( (void **) &_bpweight,sizeof(cufftComplex)* bponset[numVox]);
      
    cudaMemcpy(_bpmidx,tmp_bpmidx,sizeof(int)*bponset[numVox], cudaMemcpyHostToDevice);
    cudaMemcpy(_bpweight,tmp_bpweight,sizeof(cufftComplex)*bponset[numVox], cudaMemcpyHostToDevice);
 
    free(tmp_bpmidx);
    free(tmp_bpweight);

    cudaMalloc( (void **) &_bpsize,sizeof(int)* numVox);   
    cudaMalloc( (void **) &_bpidx,sizeof(int)* numVox);
    cudaMalloc( (void **) &_bponset,sizeof(int)* numVox+1);    
    cudaMemcpy(_bpsize,bpsize,sizeof(int)* numVox, cudaMemcpyHostToDevice);
    cudaMemcpy(_bpidx,bpidx,sizeof(int)* numVox, cudaMemcpyHostToDevice);
    cudaMemcpy(_bponset,bponset,sizeof(int)* numVox+1, cudaMemcpyHostToDevice);


            
    GET_TIME(finish);

    
    if (VERBOSE == 1) {
        mexPrintf("num active Vox: %i\n",numVox);    
        mexPrintf("alloc/copy time: %f\n",finish-start);
    }
        
    cufftPlan3d(&plan, d_pad, h_pad, w_pad, CUFFT_C2C) ;
        
  
     
    // thread managements 
    int vx_block = 128;
    dim3 dimBlock_vx(vx_block,1);
    dim3 dimGrid_vx (numVox/vx_block + 1,1);
 
    dim3 dimBlock_dw(d,1);
    dim3 dimGrid_dw (w,h);

    dim3 dimBlock_sq(d,1);
    dim3 dimGrid_sq (w*h,1);
  
    // for sensing 
    int sens_block = 256;
    dim3 dimBlock_se(sens_block,1);
    dim3 dimGrid_se (numK/sens_block + 1,1);

 
    
     double AA_time = 0;
     double cg_time = 0;
     
      int err;
    
   
     float normrr = 0;
     float dAAd = 0;
     float alpha = 0;
     float normrr2 = 0;
     float beta = 0;
      
    /////////////////////////////////////////////////////// init CG
    

     // we need this because first fft fails
     int _res = cufftExecC2C(plan, tmp1, tmp2, CUFFT_FORWARD);
    if (VERBOSE == 1)
      mexPrintf("first fft call ret: %i\n", _res);

     
    cudaMemset(_r,0, sizeof(cufftComplex)*totsz);
    cudaMemset(tmp3,0, sizeof(cufftComplex)*totsz);
    cudaMemset(_z,0, sizeof(cufftComplex)*totsz);   
    cudaMemset( tmp2,0,sizeof(cufftComplex)*totsz_pad);
   
                   
    // backproject measurement -- x=A'b
    for (int i = 0; i < numsens; i++)
    { 
        cudaMemset(tmp1,0, sizeof(cufftComplex)*totsz_pad);
        backprojVX<<<dimGrid_vx,dimBlock_vx>>>(_bpidx,_bponset,_bpweight,_bpmidx,_bpsize,_meas + i*numK, tmp1,numVox);
        if (err=cufftExecC2C(plan, tmp1, tmp2, CUFFT_INVERSE) != CUFFT_SUCCESS)
        {
            mexPrintf("cufft has failed with err %i \n",err);
            return;
        }    
        downwind<<<dimGrid_dw,dimBlock_dw>>>(_r,tmp2, _sens + i*totsz, w, h, d, w_pad, h_pad, d_pad);
     }
  
     cudaMemcpy( res, _r, sizeof(cufftComplex)*totsz,cudaMemcpyDeviceToHost);    
     cudaMemcpy( _d, _r, sizeof(cufftComplex)*totsz,cudaMemcpyDeviceToDevice);
     
    
     //normrr = Dot_wrapper(_r, _r, dot_z_d, dot_z_h, totsz, dot_blocks, dot_threads) ;
     
     normrr = cublasCdotc(totsz,(cuComplex*)_r,1,(cuComplex*)_r,1).x;           
     
     float normrr0 = normrr;
     if (VERBOSE == 1)
        mexPrintf("first residual: %f\n",normrr0);
   
    
     
    ////////////////////////////////////////////////////////////// start CG
    GET_TIME(start);
   
    for (int it = 0; it < numit; it++)
    {
     
        GET_TIME(start);
        scmult<<<dimGrid_sq,dimBlock_sq>>>(tmp3,_d,lambda,totsz);
        GET_TIME(finish); cg_time += finish-start;
        
        for (int i = 0; i < numsens; i++)
        {
            
            cudaMemset(tmp1,0, sizeof(cufftComplex)*totsz_pad);           
            upwind<<<dimGrid_dw,dimBlock_dw>>>(tmp1,_d, _sens + i*totsz, w, h, d, w_pad, h_pad, d_pad);            
            if (err=cufftExecC2C(plan, tmp1, tmp2, CUFFT_FORWARD) != CUFFT_SUCCESS)
            {
                mexPrintf("cufft has failed with err %i \n",err);
                return;
            }
            cudaMemset(tmpsens,0, sizeof(cufftComplex)*numK);
            
            
            dosens<<<dimGrid_se,dimBlock_se>>>(tmpsens,tmp2,_ipk_we,_the_index,numP,numK);
            cudaMemset(tmp1,0, sizeof(cufftComplex)*totsz_pad);            
            
            
            backprojVX<<<dimGrid_vx,dimBlock_vx>>>(_bpidx,_bponset,_bpweight,_bpmidx,_bpsize,tmpsens, tmp1,numVox);
              
                        
            if (err=cufftExecC2C(plan, tmp1, tmp2, CUFFT_INVERSE) != CUFFT_SUCCESS)
            {
                mexPrintf("cufft has failed with err %i \n",err);
                return;
            }
            downwind<<<dimGrid_dw,dimBlock_dw>>>(tmp3,tmp2, _sens + i*totsz, w, h, d, w_pad, h_pad, d_pad);                
         
        }
        
        cudaThreadSynchronize();
       
        GET_TIME(finish); AA_time += finish-start;
               
        GET_TIME(start);    
        //dAAd = Dot_wrapper(_d, tmp3, dot_z_d, dot_z_h, totsz, dot_blocks, dot_threads) ;        
        dAAd = cublasCdotc(totsz,(cuComplex*)_d,1,(cuComplex*)tmp3,1).x;           
     
        alpha = normrr/dAAd;        
        scpm<<<dimGrid_sq,dimBlock_sq>>>(_z,_d,alpha,totsz);
        scpm<<<dimGrid_sq,dimBlock_sq>>>(_r,tmp3,-alpha,totsz);                       	
        //normrr2 = Dot_wrapper(_r, _r, dot_z_d, dot_z_h, totsz, dot_blocks, dot_threads) ;   
        normrr2 = cublasCdotc(totsz,(cuComplex*)_r,1,(cuComplex*)_r,1).x;     
        beta = normrr2/normrr;
        normrr = normrr2;        
        scmultplus<<<dimGrid_sq,dimBlock_sq>>>(_d,_r,beta,totsz);     
        
        cudaThreadSynchronize();
        
        GET_TIME(finish);  cg_time += finish-start;
       
        if (sqrt(normrr/normrr0) < tol)
            break;       
        
            
        if (VERBOSE == 1)  
            mexPrintf("tol: %f\n",sqrt(normrr/normrr0));

        mexPrintf("."); mexEvalString("drawnow");
    
    }
         
   
    if (VERBOSE == 1)
    {
        mexPrintf("\n");        
        mexPrintf(" AA time: %f \n",AA_time);
        mexPrintf(" cg  time: %f \n",cg_time);
    }
        
    cudaMemcpy( res, _z, sizeof(cufftComplex)*totsz,cudaMemcpyDeviceToHost);
   

//    cudaFreeHost(dot_z_h);     
    cudaFree(tmp1);
    cudaFree(tmp2);
    cudaFree(tmp3);
    cudaFree(_r); 
    cudaFree(_d);
    cudaFree(_z);
    cudaFree(_meas);
    cudaFree(_sens);
    cudaFree(_ipk_we);
    cudaFree(_the_index);
    cudaFree(tmpsens);
    cudaFree(_bpmidx);
    cudaFree(_bpweight);
    cudaFree(_bpsize);
    cudaFree(_bpidx);
    cudaFree(_bponset);    
    
    cufftDestroy(plan);
    free(bpsize);
    free(bponset);
 

     CUcontext  pctx ;
     cuCtxPopCurrent(&pctx);	
    
}













