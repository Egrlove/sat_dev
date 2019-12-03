#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void initialization(float *hostA, float *hostB, const int size) {
	for (int i = 0; i < size; i++) {
		hostA[i] = sqrtf(i);
		hostB[i] = 2.*i;
	}
}

void workFunction() {
	float *hostA, *hostB, *hostC;
	float *devA, *devB, *devC;
	int arraySize = 512 * 50000;

	float CPUstart, CPUstop;

	float CPUtime = 0.0f;

	size_t mem_size = sizeof(float)* arraySize;

	hostA = (float*)malloc(mem_size);
	hostB = (float*)malloc(mem_size);
	hostC = (float*)malloc(mem_size);

	initialization(hostA, hostB, arraySize);
	
	CPUstart = clock();

#pragma omp parallel num_threads(8) 
	{
#pragma omp for
		for (int i = 0; i < arraySize; i++) {
			//hostC[i] = hostA[i] + hostB[i];
			
			for (int j = 0; j < 100; j++) {
				hostC[i] += cosf(sqrtf(hostA[i]) * tanhf(hostB[i])) * sqrtf(j);
			}
			
		}
	}

	CPUstop = clock();
	CPUtime = 1000.*(CPUstop - CPUstart) / CLOCKS_PER_SEC;
	printf("CPU time : %.3f ms\n", CPUtime);

	free(hostA);
	free(hostB);
	free(hostC);
}

int main() {
	workFunction();

	return 0;
}
