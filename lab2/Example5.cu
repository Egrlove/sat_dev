
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BLOCK_SIZE 16

__global__ void kernel_shared(float * a, float * b,	int n, float * c)
{
	int bx = blockIdx.x, by = blockIdx.y;
	int tx = threadIdx.x, ty = threadIdx.y;
	int aBegin = n * BLOCK_SIZE * by;
	int aEnd = aBegin + n - 1;
	int bBegin = BLOCK_SIZE * bx;
	int aStep = BLOCK_SIZE, bStep = BLOCK_SIZE * n;
	float sum = 0.0f;
	__shared__ float as[BLOCK_SIZE][BLOCK_SIZE + 1];
	__shared__ float bs[BLOCK_SIZE][BLOCK_SIZE + 1];

	for (int ia = aBegin, ib = bBegin; ia <= aEnd; ia += aStep, ib += bStep)
	{
		as[tx][ty] = a[ia + n * ty + tx];
		bs[tx][ty] = b[ib + n * ty + tx];

		__syncthreads(); 	// Synchronize to make sure the matrices are loaded 
		for (int k = 0; k < BLOCK_SIZE; k++) sum += as[k][ty] * bs[tx][k];

		__syncthreads(); 	// Synchronize to make sure submatrices not needed
	}
	c[n * BLOCK_SIZE * by + BLOCK_SIZE * bx + n * ty + tx] = sum;
}

__global__ void kernel_shared_1(float * a, float * b, int n, float * c)
{
	int bx = blockIdx.x, by = blockIdx.y;
	int tx = threadIdx.x, ty = threadIdx.y;
	int aBegin = n * BLOCK_SIZE * by;
	int aEnd = aBegin + n - 1;
	int bBegin = BLOCK_SIZE * bx;
	int aStep = BLOCK_SIZE, bStep = BLOCK_SIZE * n;
	float sum = 0.0f;
	for (int ia = aBegin, ib = bBegin; ia <= aEnd; ia += aStep, ib += bStep)
	{
		__shared__ float as[BLOCK_SIZE][BLOCK_SIZE];
		__shared__ float bs[BLOCK_SIZE][BLOCK_SIZE];
		as[ty][tx] = a[ia + n * ty + tx];
		bs[ty][tx] = b[ib + n * ty + tx];

		__syncthreads(); 	// Synchronize to make sure the matrices are loaded 
		for (int k = 0; k < BLOCK_SIZE; k++) sum += as[ty][k] * bs[k][tx];

		__syncthreads(); 	// Synchronize to make sure submatrices not needed
	}
	c[n * BLOCK_SIZE * by + BLOCK_SIZE * bx + n * ty + tx] = sum;
}

__global__ void kernel_global(float * a, float * b, int n, float * c)
{
	int   bx = blockIdx.x;
	int   by = blockIdx.y;
	int   tx = threadIdx.x;
	int   ty = threadIdx.y;
	float sum = 0.0f;
	int   ia = n * BLOCK_SIZE * by + n * ty;
	int   ib = BLOCK_SIZE * bx + tx;

	int   ic = n * BLOCK_SIZE * by + BLOCK_SIZE * bx;

	for (int k = 0; k < n; k++) sum += a[ia + k] * b[ib + k*n];

	c[ic + n * ty + tx] = sum;
}

int main()
{
	int N = 1024;
	int m, n, k;

	float CPUstart, CPUstop;
	float timerValueGPU, timerValueCPU;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	int numBytes = N*N*sizeof(float);
	float *devA, *devB, *devC, *a, *b, *c, *cc, *bT, *aT;

	a =  (float*)malloc(numBytes);
	b =  (float*)malloc(numBytes);
	bT = (float*)malloc(numBytes);
	aT = (float*)malloc(numBytes);
	c =  (float*)malloc(numBytes);
	cc = (float*)malloc(numBytes);

	for (n = 0; n<N; n++)
	{
		for (m = 0; m<N; m++)
		{
			a[m + n*N] = 2.0f*m + n;
			b[m + n*N] = m - n;
			aT[m + n*N] = m + n*2.0f;
			bT[m + n*N] = n - m;
		}
	}

	cudaMalloc((void**)&devA, numBytes);	// allocate DRAM
	cudaMalloc((void**)&devB, numBytes); // allocate DRAM
	cudaMalloc((void**)&devC, numBytes); // allocate DRAM

	dim3 threads(BLOCK_SIZE, BLOCK_SIZE);
	dim3 blocks(N / threads.x, N / threads.y);
	
	// DEVICE ------------------------------------------------------
	
	cudaEventRecord(start, 0);
	cudaMemcpy(devA, a, numBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(devB, b, numBytes, cudaMemcpyHostToDevice);

	kernel_shared <<<blocks, threads >>> (devA, devB, N, devC);
	//kernel_shared_1 << <blocks, threads >> > (devA, devB, N, devC);
	//kernel_global <<< blocks, threads >>> ( devA, devB, N, devC ); 
	cudaMemcpy(c, devC, numBytes, cudaMemcpyDeviceToHost);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&timerValueGPU, start, stop);
	printf("\n GPU calculation time %f msec\n", timerValueGPU);	
	//---------------------------------------------------------------

	// HOST ---------------------------------------------------------
	CPUstart = clock();
	for (n = 0; n<N; n++)
	{
		for (m = 0; m<N; m++)
		{
			cc[m + n*N] = 0.f;
			for (k = 0; k<N; k++) cc[m + n*N] += a[k + n*N] * bT[k + m*N]; // T
			//  for(k=0;k<N;k++) cc[m+n*N]+=a[k+n*N]*b[m+k*N]; // 
		}
	}
	CPUstop = clock();
	timerValueCPU = 1000.*(CPUstop - CPUstart) / CLOCKS_PER_SEC;
	printf("CPU time : %.3f ms\n", timerValueCPU);

	printf("Rate : %.3f \n", timerValueCPU / timerValueGPU);	
	//---------------------------------------------------------------

	cudaFree(devA);
	cudaFree(devB);
	cudaFree(devC);
	free(a);
	free(b);
	free(bT);
	free(aT);
	free(c);
	free(cc);
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	return 0;
}