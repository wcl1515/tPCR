

static __global__ void  backproj(cufftDoubleComplex *val,
                                    cufftDoubleComplex *tmp1,
                                     cufftDoubleComplex *_ipk_we,
                                     int *_the_index)
{
      int p   = threadIdx.x;
      int idx = _the_index[p];
      tmp1[idx].x +=   val[0].x * _ipk_we[p].x +  val[0].y*_ipk_we[p].y;
      tmp1[idx].y += - val[0].x*_ipk_we[p].y +  val[0].y*_ipk_we[p].x ;  
}



static __global__ void  backprojWS(cufftDoubleComplex *val,
                                    cufftDoubleComplex *tmp1,
                                     cufftDoubleComplex *_ipk_we,
                                     int *_the_index, int numP, int *ws_indices, int ws_size)
{
      int p = threadIdx.x;
      int k = blockIdx.x; 
      if (k < ws_size)
      {
          int j = ws_indices[k];     
          int q = p + numP*j;

          int idx = _the_index[q];
          tmp1[idx].x +=   val[j].x * _ipk_we[q].x +  val[j].y*_ipk_we[q].y;
          tmp1[idx].y += - val[j].x *_ipk_we[q].y  +  val[j].y*_ipk_we[q].x;  
      }
}
     

static __global__ void  backprojVX(int *vxIdx,
                                    int *onset,
                                    cufftDoubleComplex *we, 
                                    int *id,
                                    int *sz, 
                                    cufftDoubleComplex *val, cufftDoubleComplex *tmp1, int numVox)
                                  
{
      int t = blockIdx.x*blockDim.x + threadIdx.x;
      if (t < numVox)
      {
          int idx = vxIdx[t];
          int ons = onset[t];
          int size = sz[t];
          
          for (int k = 0; k < size; k++)
          {          
              int j = id[ons+k];
              tmp1[idx].x +=   val[j].x * we[ons+k].x +  val[j].y*we[ons+k].y;
              tmp1[idx].y +=   - val[j].x * we[ons+k].y  +  val[j].y*we[ons+k].x;            
          }
      }
}




static __global__ void  dosens(cufftDoubleComplex *val,
                                    cufftDoubleComplex *tmp2,
                                     cufftDoubleComplex *_ipk_we,
                                     int *_the_index,int numP, int numK)
{     
      int k = blockDim.x * blockIdx.x + threadIdx.x;
      if (k < numK)
      {
          val[k].x = 0;  val[k].y = 0;          
          for (int p = 0; p < numP; p++)
          { 
              int idx = _the_index[numP*k + p];          
              val[k].x += tmp2[idx].x*_ipk_we[numP*k + p].x - tmp2[idx].y*_ipk_we[numP*k + p].y;
              val[k].y += tmp2[idx].x*_ipk_we[numP*k + p].y + tmp2[idx].y*_ipk_we[numP*k + p].x;
          }
      }
}




static __global__ void  downwind(cufftDoubleComplex *_r,
                                    cufftDoubleComplex *tmp2,
                                     cufftDoubleComplex *_sens, int w, int h, int d, int w_pad, int h_pad, int d_pad)
{
    int z = threadIdx.x;
    int y = blockIdx.y;
    int x = blockIdx.x;
    int idx_pad = z*w_pad*h_pad+y*w_pad+x;
    int idx = z*w*h + y*w + x;
    
    _r[idx].x +=  tmp2[idx_pad].x*_sens[idx].x + tmp2[idx_pad].y*_sens[idx].y;
    _r[idx].y +=  - tmp2[idx_pad].x*_sens[idx].y + tmp2[idx_pad].y*_sens[idx].x;

}

static __global__ void  upwind(cufftDoubleComplex *_r,
                                    cufftDoubleComplex *tmp2,
                                     cufftDoubleComplex *_sens, int w, int h, int d, int w_pad, int h_pad, int d_pad)
{
    int z = threadIdx.x;
    int y = blockIdx.y;
    int x = blockIdx.x;
    int idx_pad = z*w_pad*h_pad+y*w_pad+x;
    int idx = z*w*h + y*w + x;
    
    _r[idx_pad].x +=  tmp2[idx].x*_sens[idx].x - tmp2[idx].y*_sens[idx].y;
    _r[idx_pad].y +=  tmp2[idx].x*_sens[idx].y + tmp2[idx].y*_sens[idx].x;

}



static __global__ void  scmult(cufftDoubleComplex *_a,cufftDoubleComplex *_b, double alpha, int n)
{
    int t = blockDim.x * blockIdx.x + threadIdx.x;
    if (t < n)
    {        
         _a[t].x = _b[t].x * alpha;
         _a[t].y = _b[t].y * alpha;
    }
}



static __global__ void  scmultplus(cufftDoubleComplex *_a,cufftDoubleComplex *_b, double alpha, int n)
{
    int t = blockDim.x * blockIdx.x + threadIdx.x;
    if (t < n)
    {        
         _a[t].x = _a[t].x * alpha + _b[t].x;
         _a[t].y = _a[t].y * alpha + _b[t].y;
    }
}



static __global__ void  scpm(cufftDoubleComplex *_a,cufftDoubleComplex *_b, double alpha, int n)
{
    int t = blockDim.x * blockIdx.x + threadIdx.x;
    if (t < n)
    {        
         _a[t].x += _b[t].x * alpha;
         _a[t].y += _b[t].y * alpha;
    }    
}

__global__ void Dev_dot(cufftDoubleComplex x[], cufftDoubleComplex y[], double z[], int n) {
   /* Use tmp to store products of vector components in each block */
   /* Can't use variable dimension here                            */
   __shared__ double tmp[MAX_BLOCK_SZ];
   int t = blockDim.x * blockIdx.x + threadIdx.x;
   int loc_t = threadIdx.x;


   if (t < n) 
   {
       tmp[loc_t] = x[t].x*y[t].x + x[t].y*y[t].y;
   }
   __syncthreads();
   
   /* This uses a tree structure to do the addtions */
   for (int stride = blockDim.x/2; stride >  0; stride /= 2) {
      if (loc_t < stride)
      {
         tmp[loc_t] += tmp[loc_t + stride];
      }
      __syncthreads();
   }

   /* Store the result from this cache block in z[blockIdx.x] */
   if (threadIdx.x == 0) {
      z[blockIdx.x] = tmp[0];
   }
}  /* Dev_dot */    

 
double Dot_wrapper(cufftDoubleComplex x_d[], cufftDoubleComplex y_d[], double z_d[], double z_h[],
      int n, int blocks, int threads) { 
   int i;
   double dot = 0;

   /* Invoke kernel */
   Dev_dot<<<blocks, threads>>>(x_d, y_d, z_d, n);
   cudaThreadSynchronize();

   /* Note that we don't need to copy z_d back to host */
   for (i = 0; i < blocks; i++)
   {
      dot += z_h[i];
   }
   return dot;
} 


