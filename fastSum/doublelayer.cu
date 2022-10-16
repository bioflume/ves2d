/* CUDA MATLAB FAST SUMMATION FUNCTION
   The input signature for this program is PointX, PointY, fx, fy.
   The output signature for this program is PotentialX, PotentialY
   Matlab command line:
   [PotentialX,PotentialY]=doublelayer(xs,ys,fx,fy,nx,ny,xt,yt)

   Sample run in Matlab:
   f=0:1:127
   g=0:1:64
   [gx,gy]=doublelayer(f,f,f,f,f,f,g,g)

   The dimensions of all of the input, xs,ys, fx, fy,nx, ny, xt, yt should be the multiples of 64
   All inputs are 1D arrays in which number of columns equal to number of targets, (the number of rows is 1)



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


__global__ void Potential(int numberOfPoints,float *pointX, float *pointY, float *densityX, float *densityY,float *cudanx,float *cudany,float *targetX, float *targetY, float *potentialX, float *potentialY){

    __shared__ float ptX[BLOCK_SIZE];
    __shared__ float ptY[BLOCK_SIZE];
    __shared__ float fx[BLOCK_SIZE];
    __shared__ float fy[BLOCK_SIZE];
    __shared__ float nx[BLOCK_SIZE];
    __shared__ float ny[BLOCK_SIZE];
         


    int i;


    float fastSummation[2],rx,ry,dist,c;
    //Block index

    //Thread index
    int tx =threadIdx.x;

    int tid=blockIdx.x*blockDim.x + tx;

    float mypositionX=targetX[tid];
    float mypositionY=targetY[tid];

    fastSummation[0]=0;
    fastSummation[1]=0;

    for (int a=0; a<numberOfPoints;a+=BLOCK_SIZE){

        ptX[tx]=pointX[a+tx];
        ptY[tx]=pointY[a+tx];
        fx[tx]=densityX[a+tx];
        fy[tx]=densityY[a+tx];
        nx[tx]=cudanx[a+tx];
        ny[tx]=cudany[a+tx];
        __syncthreads();

        for (i=0;i<BLOCK_SIZE;i++){

            rx=mypositionX-ptX[i];
            ry=mypositionY-ptY[i];

            dist = rx*rx+ry*ry;
            c =(rx*fx[i]+ry*fy[i])*(rx*nx[i]+ry*ny[i] )/dist/dist;
            fastSummation[0]+=c*rx;
            fastSummation[1]+=c*ry;

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
    if (nrhs!=8) {
        printf("Number of inputs must be 8. The signature must be PointX, PointY, fx, fy \n");
        return;
    }

    if ((mxGetN(prhs[0])!=mxGetN(prhs[1]))||
        (mxGetN(prhs[0])!=mxGetN(prhs[2])) ||
        (mxGetN(prhs[0])!=mxGetN(prhs[3])) ||
        (mxGetN(prhs[0])!=mxGetN(prhs[4]))||
        (mxGetN(prhs[0])!=mxGetN(prhs[5])) ){
        printf("The signature should be (xs,ys,fx,fy,nx,ny,xt,yt)");
        printf("Dimensions of source inputs xs,ys,fx,fy,nx,ny must be the same\n");
        exit(-1);
    }

    if ((mxGetN(prhs[6])!=mxGetN(prhs[7]))) {
        printf("Dimensions of target inputs xt,yt must be the same \n");
        return;
    }



    //Starting ...
    int i,j;

    int numberOfPoints=mxGetN(prhs[0]); //=(int) *mxGetPr(prhs[0]);
    int numberOfTargetPoints=mxGetN(prhs[6]);

    double *pointX,*pointY,*fx,*fy,*nx,*ny,*xt,*yt,*potentialX,*potentialY;
    float *pointXf,*pointYf,*fxf,*fyf,*nxf,*nyf,*xtf,*ytf,*potentialXf,*potentialYf;


    float* cudapointX,*cudapointY,*cudafx,*cudafy,*cudanx,*cudany,*cudaxt,*cudayt,*cudaPotentialX,*cudaPotentialY;
    int size= numberOfPoints*sizeof(float);
    int sizetarget = numberOfTargetPoints*sizeof(float);
    mxClassID category;

    //  printf("number of target is %d\n",numberOfTargetPoints); 
    cudaMalloc((void **) &cudapointX,size);
    cudaMalloc((void **) &cudapointY,size);
    cudaMalloc((void **) &cudafx,size);
    cudaMalloc((void **) &cudafy,size);
    cudaMalloc((void **) &cudanx,size);
    cudaMalloc((void **) &cudany,size);

    cudaMalloc((void **) &cudaxt,sizetarget);
    cudaMalloc((void **) &cudayt,sizetarget);


    cudaMalloc((void **) &cudaPotentialX,sizetarget);
    cudaMalloc((void **) &cudaPotentialY,sizetarget);

 
    plhs[0]=mxCreateDoubleMatrix(1,numberOfTargetPoints,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(1,numberOfTargetPoints,mxREAL);

    potentialXf= (float *) mxMalloc(sizeof(float)*numberOfTargetPoints);
    potentialYf= (float *) mxMalloc(sizeof(float)*numberOfTargetPoints);


    pointX = mxGetPr(prhs[0]);
    pointY = mxGetPr(prhs[1]);
    fx = mxGetPr(prhs[2]);
    fy= mxGetPr(prhs[3]);
    nx = mxGetPr(prhs[4]);
    ny= mxGetPr(prhs[5]);

    xt = mxGetPr(prhs[6]);
    yt = mxGetPr(prhs[7]);


#ifdef CPU
    float rx,ry,dist,fastSummation[2],c;

    for (i=0;i<numberOfTargetPoints;++i){


        fastSummation[0]=0;
        fastSummation[1]=0;
        for (j=0;j<numberOfPoints;j++){

            rx=-pointX[j]+xt[i];
            ry=-pointY[j]+yt[i];
            dist = rx*rx+ry*ry;
            c =(rx*fx[j]+ry*fy[j])*(rx*nx[j]+ry*ny[j] )/dist/dist;
            fastSummation[0]+=c*rx;
            fastSummation[1]+=c*ry;

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
        cudaMemcpy(cudanx,nx,size,cudaMemcpyHostToDevice);
        cudaMemcpy(cudany,ny,size,cudaMemcpyHostToDevice);
        cudaMemcpy(cudaxt,xt,sizetarget,cudaMemcpyHostToDevice);
        cudaMemcpy(cudayt,yt,sizetarget,cudaMemcpyHostToDevice);

    }
    if( category == mxDOUBLE_CLASS)
    {
        // The input array is in double precision, it needs to be converted t
        //floats before being sent to the card 
        pointXf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
        pointYf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
        fxf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
        fyf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
        nxf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
        nyf=(float *) mxMalloc(sizeof(float)*numberOfPoints);
        xtf=(float *) mxMalloc(sizeof(float)*numberOfTargetPoints);
        ytf=(float *) mxMalloc(sizeof(float)*numberOfTargetPoints);


        for (i=0;i<numberOfPoints;i++){
            pointXf[i]=pointX[i];
            pointYf[i]=pointY[i];
            fxf[i]=fx[i];
            fyf[i]=fy[i];
            nxf[i]=nx[i];
            nyf[i]=ny[i];

        }

        for (i=0;i<numberOfTargetPoints;i++){

            xtf[i]=xt[i];
            ytf[i]=yt[i];

        }
        cudaMemcpy(cudapointX,pointXf,size,cudaMemcpyHostToDevice);
        cudaMemcpy(cudapointY,pointYf,size,cudaMemcpyHostToDevice);
        cudaMemcpy(cudafx,fxf,size,cudaMemcpyHostToDevice);
        cudaMemcpy(cudafy,fyf,size,cudaMemcpyHostToDevice);
        cudaMemcpy(cudanx,nxf,size,cudaMemcpyHostToDevice);
        cudaMemcpy(cudany,nyf,size,cudaMemcpyHostToDevice);
        cudaMemcpy(cudaxt,xtf,sizetarget,cudaMemcpyHostToDevice);
        cudaMemcpy(cudayt,ytf,sizetarget,cudaMemcpyHostToDevice);

    }

    dim3 dimBlock(BLOCK_SIZE,1);
    dim3 dimGrid(numberOfTargetPoints/BLOCK_SIZE,1);

    Potential<<<dimGrid,dimBlock>>>(numberOfPoints,cudapointX,cudapointY,cudafx,cudafy,cudanx,cudany,cudaxt,cudayt,cudaPotentialX,cudaPotentialY);

    checkCUDAError("kernel execution");

    potentialXf= (float *) mxMalloc(sizeof(float)*numberOfTargetPoints);
    potentialYf= (float *) mxMalloc(sizeof(float)*numberOfTargetPoints);


    cudaMemcpy(potentialXf,cudaPotentialX,sizetarget,cudaMemcpyDeviceToHost);
    cudaMemcpy(potentialYf,cudaPotentialY,sizetarget,cudaMemcpyDeviceToHost);

#endif

    potentialX=mxGetPr(plhs[0]);
    potentialY=mxGetPr(plhs[1]);


    /* Convert from single to double before returning */
    for (j = 0; j < numberOfTargetPoints; j++)
    {
        potentialX[j] = (double) potentialXf[j];
        potentialY[j] = (double) potentialYf[j];

    }


    //Free device memory
    cudaFree(cudapointX);
    cudaFree(cudapointY);
    cudaFree(cudafx);
    cudaFree(cudafy);
    cudaFree(cudanx);
    cudaFree(cudany);
    cudaFree(cudaxt);
    cudaFree(cudayt);

    cudaFree(cudaPotentialX);
    cudaFree(cudaPotentialY);



}
