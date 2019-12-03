
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <stdio.h>

__global__ void Kernel(float *X, float *Y, float *Z) {
	unsigned int idx_X = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int idx_Y = threadIdx.y + blockIdx.y * blockDim.y;

	Z[idx_X*dim + idx_Y] = 20. + X[idx_X]*X[idx_X] +  Y[idx_Y]*Y[idx_Y] - 10.*(cosf(2.*3.14*X[idx_X]) + cosf(2.*3.14*X[idx_X]));
}

void initialization(const float leftB, const float rightB, float *X, const unsigned int dim) {
	try {
		float step = (rightB - leftB)/(float)dim;
		if (step <= 0) throw "error";

		X[0] = leftB;
		for( unsigned int i=1; i<dim ; i++) {
			X[i] = X[i-1] + step;
		}
	}
	catch (...) {
		fprintf(stderr, "step failed!");
        exit(1);
	}
}

int main()
{
	int dim = 2048;
	size_t mem_size = sizeof(float)*dim;
	cudaError cudaStatus;

	float *hostX, *hostY, *hostZ;
	float *devX, *devY, *devZ;

	float rightB, leftB;
	leftB = -5;
	rightB = 5;

	hostX = (float*)malloc(mem_size);
	hostY = (float*)malloc(mem_size);
	hostZ = (float*)malloc(mem_size*mem_size);

	initialization(leftB, rightB, hostX, dim);

	//memcpy(hostY,hostX,mem_size);

	cudaMalloc((void**)&devX, mem_size);
	cudaMalloc((void**)&devY, mem_size);
	cudaMalloc((void**)&devZ, mem_size*mem_size);

	cudaMemcpy(devX, hostX, mem_size, cudaMemcpyHostToDevice);
	cudaMemcpy(devY, devX, mem_size, cudaMemcpyDeviceToDevice);

	dim3 N_Grid  (dim/32,dim/32,1);
	dim3 N_Block (48,48,1);
	
	Kernel <<< N_Grid, N_Block >>> (devX,devY,devZ);
	cudaStatus = cudaGetLastError();

	if(cudaStatus != cudaSuccess) {
		printf("Last error: %s\n", cudaGetErrorString(cudaStatus));
		return 0;
	}

	cudaMemcpy(hostZ, devZ, mem_size*mem_size, cudaMemcpyDeviceToHost);

	for(unsigned int i=0; i<dim*dim; i++) {
		std::cout << "i: " << hostZ[i] << std::endl;
	}

	cudaFree(devX);
    cudaFree(devY);
    cudaFree(devZ);
	free(hostX);
	free(hostY);
	free(hostZ);

    return 0;
}