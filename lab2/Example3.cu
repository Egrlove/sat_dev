
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (1024)

__global__ void mult(float *A, float *B, float *C) {
	unsigned int idx_X = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int idx_Y = threadIdx.y + blockIdx.y * blockDim.y;
	float sum = 0.;

	if ((idx_X < N) && (idx_Y < N)) {
		for (int i = 0; i < N; i++) {
			sum += A[idx_X*N + i] * B[idx_Y + i*N];
		}
		C[idx_X*N + idx_Y] = sum;
	}
}

int main(void) {
	cudaEvent_t GPUstart, GPUstop;
	float GPUtime = 0.0f;

	float *hostA, *hostB;
	float *hostC;

	float *devA, *devB;
	float *devC;

	size_t mem_size = N*N*sizeof(float);

	hostA = (float *)malloc(mem_size);
	hostB = (float *)malloc(mem_size);
	hostC = (float *)malloc(mem_size);

	cudaMalloc((void**)&devA, mem_size);
	cudaMalloc((void**)&devB, mem_size);
	cudaMalloc((void**)&devC, mem_size);

	for (int i = 0; i < N*N; i++) {
		hostA[i] = sqrtf(i);
		hostB[i] = sinf(i);
		hostC[i] = 0.;
	}

	int N_Threads = 8;
	int N_Blocks = 0;

	if (((N) % N_Threads) == 0) {
		N_Blocks = ((N) / N_Threads);
	}
	else {
		N_Blocks = ((N) / N_Threads) + 1;
	}
	dim3 Threads(N_Threads,N_Threads);
	dim3 Blocks(N_Blocks, N_Blocks);

	cudaEventCreate(&GPUstart);
	cudaEventCreate(&GPUstop);

	cudaEventRecord(GPUstart, 0);

	cudaMemcpy(devA, hostA, mem_size, cudaMemcpyHostToDevice);
	cudaMemcpy(devB, hostB, mem_size, cudaMemcpyHostToDevice);
	cudaMemset(devC, 0, mem_size);

	mult <<< Blocks, Threads >>> (devA, devB, devC);

	cudaMemcpy(hostC, devC, mem_size, cudaMemcpyDeviceToHost);

	cudaEventRecord(GPUstop, 0);
	cudaEventSynchronize(GPUstop);

	cudaEventElapsedTime(&GPUtime, GPUstart, GPUstop);
	printf("GPU time : %.3f ms\n", GPUtime);

	cudaFree(devA);
	cudaFree(devB);
	cudaFree(devC);
	free(hostA);
	free(hostB);
	free(hostC);

	return 0;
}