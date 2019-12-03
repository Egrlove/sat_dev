
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N_STREAMS 4

static void HandleError(cudaError_t err,
	const char *file,
	int line)
{
	if (err != cudaSuccess) {
		printf("%s in %s at line %d\n", cudaGetErrorString(err),
			file, line);
		exit(EXIT_FAILURE);
	}
}
#define HANDLE_ERROR( error ) (HandleError( error, __FILE__, __LINE__ ))

__global__ void addKernel(const float *a, const float *b, float *c, const int size) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < size) {
		c[i] = a[i] + b[i];
	}
}

__global__ void multKernel(const float *a, const float *b, float *c, const int size) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < size) {
		c[i] = a[i] * b[i];
	}
}


__global__ void func1Kernel(const float *a, const float *b, float *c, const int size) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i < size) {
		for (int j = 0; j < 100; j++) {
			c[i] += cosf(sqrtf(a[i]) * tanhf(b[i])) * sqrtf(j);
		}
	}
}

void initialization(float *hostA, float *hostB, const int size) {
	for (int i = 0; i < size; i++) {
		hostA[i] = sqrtf(i);
		hostB[i] = 2.*i;
	}
}

void workFunction() {
	float *hostA, *hostB, *hostC;
	float *devA, *devB, *devC;
	int arraySize = ( 512 * 5000  ) / N_STREAMS;

	cudaEvent_t GPUstart, GPUstop;
	float CPUstart, CPUstop;

	float GPUtime = 0.0f;
	float CPUtime = 0.0f;

	int N_threads = 512;
	int N_blocks;

	cudaStream_t stream[N_STREAMS];

	for (int i = 0; i < N_STREAMS; i++) {
		cudaStreamCreate(&stream[i]);
	}

	size_t mem_size = sizeof(float)*arraySize;

	HANDLE_ERROR(cudaMallocHost((void**)&hostA, mem_size*N_STREAMS));
	HANDLE_ERROR(cudaMallocHost((void**)&hostB, mem_size*N_STREAMS));
	HANDLE_ERROR(cudaMallocHost((void**)&hostC, mem_size*N_STREAMS));
	
	HANDLE_ERROR(cudaMalloc((void**)&devA, mem_size*N_STREAMS));
	HANDLE_ERROR(cudaMalloc((void**)&devB, mem_size*N_STREAMS));
	HANDLE_ERROR(cudaMalloc((void**)&devC, mem_size*N_STREAMS));
	

	initialization(hostA, hostB, arraySize);

	if ((arraySize % N_threads) == 0) {
		N_blocks = (arraySize / N_threads);
	}
	else {
		N_blocks = (arraySize / N_threads) + 1;
	}

	dim3 Threads(N_threads);
	dim3 Blocks(N_blocks);

	cudaMemset(devC, 0, mem_size);

	for (int i = 0; i < N_STREAMS; ++i) {
		HANDLE_ERROR(cudaMemcpyAsync(devA + i*arraySize / N_STREAMS, hostA + i*arraySize / N_STREAMS, mem_size,
			cudaMemcpyHostToDevice, stream[i]) );
		HANDLE_ERROR(cudaMemcpyAsync(devB + i*arraySize / N_STREAMS, hostB + i*arraySize / N_STREAMS, mem_size,
			cudaMemcpyHostToDevice, stream[i]) );
	}

	cudaEventCreate(&GPUstart);
	cudaEventCreate(&GPUstop);

	cudaEventRecord(GPUstart, 0);
	
	for (int i = 0; i < N_STREAMS; i++) {
		func1Kernel << < Blocks, Threads, 0, stream[i] >> > (devA + i*arraySize / N_STREAMS,
															 devB + i*arraySize / N_STREAMS,
															 devC + i*arraySize / N_STREAMS,
															 arraySize);
	}

	HANDLE_ERROR(cudaGetLastError());

	cudaEventRecord(GPUstop, 0);
	cudaEventSynchronize(GPUstop);

	cudaEventElapsedTime(&GPUtime, GPUstart, GPUstop);
	printf("GPU time : %.3f ms\n", GPUtime);

	for (int i = 0; i < N_STREAMS; i++) {
		cudaMemcpyAsync(hostC + i*mem_size, devC + i*mem_size, mem_size,
			cudaMemcpyDeviceToHost, stream[i]);
	}

	CPUstart = clock();

	for (int i = 0; i < arraySize; i++) {
		//hostC[i] = hostA[i] + hostB[i];

		for (int j = 0; j < 100; j++) {
			hostC[i] += cosf(sqrtf(hostA[i]) * tanf(hostB[i])) * sqrtf(j);
		}
	}

	CPUstop = clock();
	CPUtime = 1000.*(CPUstop - CPUstart) / CLOCKS_PER_SEC;
	printf("CPU time : %.3f ms\n", CPUtime);

	printf("Rate : %.3f \n", CPUtime / GPUtime);
	
	HANDLE_ERROR(cudaFreeHost(hostA));
	HANDLE_ERROR(cudaFreeHost(hostB));
	HANDLE_ERROR(cudaFreeHost(hostC));

	HANDLE_ERROR(cudaFree(devA));
	HANDLE_ERROR(cudaFree(devB));
	HANDLE_ERROR(cudaFree(devC));

	HANDLE_ERROR(cudaEventDestroy(GPUstart));
	HANDLE_ERROR(cudaEventDestroy(GPUstop));
}

int main() {
	workFunction();

	return 0;
}