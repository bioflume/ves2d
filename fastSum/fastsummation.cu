/* CUDA MATLAB FAST SUMMATION FUNCTION
The input signature for this program is PointX, PointY, fx, fy.
The output signature for this program is PotentialX, PotentialY
Matlab command line:
[PotentialX,PotentialY]=fastsummation(PointX,PointY,fx,fy)

Sample run in Matlab:
f=0:1:63
[gx,gy]=fastsummation(f,f,f,f)

The dimensions of PointX, PointY, fx, fy should be the multiples of 64
PointX, PointY,fx,fy is 1D array in which number of columns equal to number of targets, (the number of rows is 1)



*/
#include <stdlib.h>

#include "mex.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "cufft.h"
#include "driver_types.h"
#include </opt/cuda_sdk/common/inc/cutil.h>

#define BLOCK_SIZE 64
//#define CPU
#define GPU

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }
}


//Device potential calculation function called Potential
__global__ void Potential(int numberOfPoints,float *pointX, float *pointY, float *densityX, float *densityY, float *potentialX, float *potentialY){

          __shared__ float ptX[BLOCK_SIZE];
          __shared__ float ptY[BLOCK_SIZE];
          __shared__ float fx[BLOCK_SIZE];
          __shared__ float fy[BLOCK_SIZE];



   int i;


   float fastSummation[2],rx,ry,dist;
   //Block index

   //Thread index
   int tx =threadIdx.x;

   int tid=blockIdx.x*blockDim.x + tx;

   float mypositionX=pointX[tid];
   float mypositionY=pointY[tid];

   fastSummation[0]=0;
   fastSummation[1]=0;

   for (int a=0; a<numberOfPoints;a+=BLOCK_SIZE){

      ptX[tx]=pointX[a+tx];
      ptY[tx]=pointY[a+tx];
      fx[tx]=densityX[a+tx];
      fy[tx]=densityY[a+tx];
      __syncthreads();

      for (i=0;i<BLOCK_SIZE;i++){

        if (tid!=i+a){
         rx=ptX[i]-mypositionX;
         ry=ptY[i]-mypositionY;
//         dist = sqrtf(rx*rx+ry*ry);
//         fastSummation[0]+=-logf(dist)*fx[i]+(rx*fx[i]+ry*fy[i])*rx/dist;
//         fastSummation[1]+=-logf(dist)*fy[i]+(rx*fx[i]+ry*fy[i])*ry/dist;

         dist = rx*rx+ry*ry;
         fastSummation[0]+=-logf(dist)*fx[i]/2+(rx*fx[i]+ry*fy[i])*rx/dist;
         fastSummation[1]+=-logf(dist)*fy[i]/2+(rx*fx[i]+ry*fy[i])*ry/dist;

                  }
                }
      __syncthreads();
        }
             potentialX[tid]= fastSummation[0];
             potentialY[tid]= fastSummation[1];

}

int devcheck(int gpudevice)
{
    int device_count=0;
    int device;  // used with cudaGetDevice() to verify cudaSetDevice()

    // get the number of non-emulation devices detected
    cudaGetDeviceCount( &device_count);
    if (gpudevice >= device_count)
    {
        printf("gpudevice = %d , valid devices = [ ", gpudevice);
   for (int i=0; i<device_count; i++)
      printf("%d ", i);
   printf("] ... exiting \n");
        exit(1);
    }
    cudaError_t cudareturn;
    cudaDeviceProp deviceProp;

    // cudaGetDeviceProperties() is also demonstrated in the deviceQuery/ 
    // of the sdk projects directory
    cudaGetDeviceProperties(&deviceProp, gpudevice);
    printf("[deviceProp.major.deviceProp.minor] = [%d.%d]\n",
     deviceProp.major, deviceProp.minor);

    if (deviceProp.major > 999)
    {
        printf("warning, CUDA Device Emulation (CPU) detected, exiting\n");
                exit(1);
    }

    // choose a cuda device for kernel execution
    cudareturn=cudaSetDevice(gpudevice);
    if (cudareturn == cudaErrorInvalidDevice)
    {
        perror("cudaSetDevice returned cudaErrorInvalidDevice");
    }
    else
    {
        // double check that device was properly selected
        cudaGetDevice(&device);
        printf("cudaGetDevice()=%d\n",device);
    }
    return(0);
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){


//Checking input
   if (nrhs!=4) {
      printf("Number of inputs must be 4. The signature must be PointX, PointY, fx, fy \n");
      return;
   }

   if ((mxGetN(prhs[0])!=mxGetN(prhs[1]))||
      (mxGetN(prhs[0])!=mxGetN(prhs[2])) ||
      (mxGetN(prhs[0])!=mxGetN(prhs[3]))){
      printf("Dimensions of all inputs must be the same\n");
      exit(-1);
    }


//Starting ...
   int i,j;

  int numberOfPoints=mxGetN(prhs[0]); //=(int) *mxGetPr(prhs[0]);


   double *pointX,*pointY,*fx,*fy,*potentialX,*potentialY;
   float *pointXf,*pointYf,*fxf,*fyf,*potentialXf,*potentialYf;


   float* cudapointX,*cudapointY,*cudafx,*cudafy,*cudaPotentialX,*cudaPotentialY;
   int size= numberOfPoints*sizeof(float);

   mxClassID category;

   
   cudaMalloc((void **) &cudapointX,size);
   cudaMalloc((void **) &cudapointY,size);
   cudaMalloc((void **) &cudafx,size);
   cudaMalloc((void **) &cudafy,size);
   cudaMalloc((void **) &cudaPotentialX,size);
   cudaMalloc((void **) &cudaPotentialY,size);
 
   plhs[0]=mxCreateDoubleMatrix(1,numberOfPoints,mxREAL);
   plhs[1]=mxCreateDoubleMatrix(1,numberOfPoints,mxREAL);

   potentialXf= (float *) mxMalloc(sizeof(float)*numberOfPoints);
   potentialYf= (float *) mxMalloc(sizeof(float)*numberOfPoints);


   pointX = mxGetPr(prhs[0]);
   pointY = mxGetPr(prhs[1]);
   fx = mxGetPr(prhs[2]);
   fy= mxGetPr(prhs[3]);

#ifdef CPU
   float rx,ry,dist,fastSummation[2];

 for (i=0;i<numberOfPoints;++i){


      fastSummation[0]=0;
      fastSummation[1]=0;
            for (j=0;j<numberOfPoints;j++){

      if (i!=j){
        rx=pointX[j]-pointX[i];
        ry=pointY[j]-pointY[i];
//        dist = sqrtf(rx*rx+ry*ry);

//        fastSummation[0]+=-logf(dist)*fx[j]+(rx*fx[j]+ry*fy[j])*rx/dist;
//        fastSummation[1]+=-logf(dist)*fy[j]+(rx*fx[j]+ry*fy[j])*ry/dist;


        dist = rx*rx+ry*ry;
        fastSummation[0]+=-logf(dist)*fx[i]/2+(rx*fx[i]+ry*fy[i])*rx/dist;
        fastSummation[1]+=-logf(dist)*fy[i]/2+(rx*fx[i]+ry*fy[i])*ry/dist;

      }

        }

        potentialXf[i]= fastSummation[0];
        potentialYf[i]= fastSummation[1];

      }

   
#endif


#ifdef GPU
   category = mxGetClassID(prhs[0]);


 if( category == mxSINGLE_CLASS){
      cudaMemcpy(cudapointX,pointX,size,cudaMemcpyHostToDevice);
      cudaMemcpy(cudapointY,pointY,size,cudaMemcpyHostToDevice);
      cudaMemcpy(cudafx,fx,size,cudaMemcpyHostToDevice);
      cudaMemcpy(cudafy,fy,size,cudaMemcpyHostToDevice);
   }
    if( category == mxDOUBLE_CLASS)
    {
      // The input array is in double precision, it needs to be converted t
//floats before being sent to the card 
      pointXf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
      pointYf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
      fxf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
      fyf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
      for (i=0;i<numberOfPoints;i++){
         pointXf[i]=pointX[i];
         pointYf[i]=pointY[i];
         fxf[i]=fx[i];
         fyf[i]=fy[i];
      }

      cudaMemcpy(cudapointX,pointXf,size,cudaMemcpyHostToDevice);
      cudaMemcpy(cudapointY,pointYf,size,cudaMemcpyHostToDevice);
      cudaMemcpy(cudafx,fxf,size,cudaMemcpyHostToDevice);
      cudaMemcpy(cudafy,fyf,size,cudaMemcpyHostToDevice);
    }

   dim3 dimBlock(BLOCK_SIZE,1);
   dim3 dimGrid(numberOfPoints/BLOCK_SIZE,1);

  
   Potential<<<dimGrid,dimBlock>>>(numberOfPoints,cudapointX,cudapointY,cudafx,cudafy,cudaPotentialX,cudaPotentialY);

   checkCUDAError("kernel execution");

   potentialXf= (float *) mxMalloc(sizeof(float)*numberOfPoints);
   potentialYf= (float *) mxMalloc(sizeof(float)*numberOfPoints);


//Read potential from the CUDA device  
   cudaMemcpy(potentialXf,cudaPotentialX,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(potentialYf,cudaPotentialY,size,cudaMemcpyDeviceToHost);

#endif

  potentialX=mxGetPr(plhs[0]);
   potentialY=mxGetPr(plhs[1]);


/* Convert from single to double before returning */
   for (j = 0; j < numberOfPoints; j++)
   {
        potentialX[j] = (double) potentialXf[j];
        potentialY[j] = (double) potentialYf[j];

   }


//Free device memory
   cudaFree(cudapointX);
   cudaFree(cudapointY);
   cudaFree(cudafx);
   cudaFree(cudafy);
   cudaFree(cudaPotentialX);
   cudaFree(cudaPotentialY);

}
