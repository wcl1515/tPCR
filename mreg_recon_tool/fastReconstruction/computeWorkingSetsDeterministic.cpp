
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <complex>
#include <unistd.h>
#include "fftw3.h"

#include <vector>

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


using namespace std;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 2 ) {
	printf("\nUsage:\n");
    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}

    ///////////////////////////////////////////// fetching data 

    int pcnt = 0;  
    const mxArray *Traj;
    Traj = prhs[pcnt++];       
    double *traj = ( double *) mxGetData(Traj);

    const mxArray *Params;
    Params = prhs[pcnt++];       
    double *params = ( double *) mxGetData(Params);
    
    
    int MAX_D1 = (int) params[0];
    int MAX_D2 = (int) params[1];
    int MAX_D3 = (int) params[2];
    int N1 = (int) params[3];
    int N2 = (int) params[4];
    int N3 = (int) params[5];
    
 
    int numK = mxGetN(Traj);
    
    
    int MAX_WS_SIZE = 256;
    
    int tr[numK*3];   
    int idx[numK];
    for (int i = 0; i < numK*3 ; i++)
        tr[i] = int(traj[i]);
     for (int i = 0; i < numK ; i++)
        idx[i] = i;
    
    int start = 0;
    int end = 1;
    int candidate = 1;
    int k,tmp;
    vector<int> start_p;
    vector<int> end_p;    
    while(true)
    {
        candidate++;
        if (candidate >= numK || end-start > MAX_WS_SIZE) // no candidates anymore
        {
            // save WS 
       
            start_p.push_back(start);
            end_p.push_back(end);
            start = end;
            end = start+1;            
            candidate = end;
        }
        
        if (end >= numK)
        {
            start_p.push_back(start);
            end_p.push_back(end);
            break;
        }
            
        
        for (k = start; k < end; k++)
        {
            if (abs(tr[candidate*3] - tr[k*3])  < MAX_D1)
                break;
            if (abs(tr[candidate*3+1] - tr[k*3+1])  < MAX_D2)
                break;
            if (abs(tr[candidate*3+2] - tr[k*3+2])  < MAX_D3)
                break;                    

            if (N1-abs(tr[candidate*3] - tr[k*3])  < MAX_D1)
                break;
            if (N2-abs(tr[candidate*3+1] - tr[k*3+1])  < MAX_D2)
                break;
            if (N3-abs(tr[candidate*3+2] - tr[k*3+2])  < MAX_D3)
                break;                    
        }
        if (k == end) // no collision
        {
            for (int d = 0; d < 3; d++)
            {
                tmp = tr[candidate*3+d]; tr[candidate*3+d] = tr[end*3+d]; tr[end*3+d] = tmp;
            }
            tmp = idx[candidate]; idx[candidate] = idx[end]; idx[end] = tmp;
            end++;
            candidate = end;            
        }
        
                
    }
    
    int dims = start_p.size();
    plhs[0] =  mxCreateCellArray(1, &dims);
    for (int i = 0; i < start_p.size(); i++)
    {
        int nd = end_p[i]-start_p[i];
        mxArray *mxt = mxCreateNumericArray(1,&nd,mxINT32_CLASS,mxREAL);
        int *t = (int*) mxGetData(mxt);
        for (int j = 0; j < end_p[i]-start_p[i];j++)
        {
            t[j] = idx[start_p[i] + j];
        }
        
        mxSetCell(plhs[0],i,mxt);
    }


    
    
}    





