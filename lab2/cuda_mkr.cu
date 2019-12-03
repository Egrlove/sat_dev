#define LENX 40
#define LENY 25
#define LENE 25


// #define LENX 8
// #define LENY 5
// #define LENE 5
#define T0 0
#define T1 200
#define T2 50

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"


static void HandleError(cudaError_t err,const char *file,int line)
{
  if (err != cudaSuccess)
  {
    printf("%s in %s at line %d\n", cudaGetErrorString(err),
           file, line);
    exit(EXIT_FAILURE);
  }
}

#define HANDLE_ERROR(error) (HandleError(error, __FILE__, __LINE__))

using namespace std;

void mtrx_print();
void sltn_print();

long int Nm, Nb, Nth, modTime;
long int cntB = 0;
long double *mA, *vX, *vB;
long double L, H, deltaX, deltaT, startT, Tle, Tre;
long double a = 1.0;
int leftEdge, rightEdge;
long int nodePerL, NoX, NoY, edge;

clock_t time_start() { return clock(); }
clock_t time_stop() { return clock(); }

/* Gauss */

void frw_one_th() {
  long int i, j, k;
  long double dgE;

  for (k = 0; k < Nm; k++) {
    dgE = mA[Nm * k + k];
    for (j = k; j < Nm; j++)
      mA[Nm * k + j] /= dgE;
    vB[k] /= dgE;
    // cout << k << endl;
    for (i = k + 1; i < Nm; i++) {
      dgE = mA[Nm * i + k];
      for (j = k; j < Nm; j++)
        mA[Nm * i + j] -= mA[Nm * k + j] * dgE;
      vB[i] -= vB[k] * dgE;
      // cout << vB[i] << endl;
    }
    // cout << endl;
  }
}


__global__ void cuda_gauss_step(long double *dev_mA, long double *dev_vB, int Nm, int k)
{
  double dgE;
  int  j;
  int cuda_i = threadIdx.x + blockIdx.x * blockDim.x;
  if (cuda_i < Nm && cuda_i > k )
  {
    dgE =dev_mA[Nm * (cuda_i) + k] / dev_mA[Nm * (k) + k];
    for (j = k; j < Nm; j++)
      dev_mA[Nm * (cuda_i) + j] -= dev_mA[Nm * k + j] * dgE;
    dev_vB[cuda_i] -= dev_vB[k] * dgE;
  }

}

void cuda_frw_one_th()
{
  long int k;
  // long double dgE;

  long double *dev_mA;
  long double *dev_vB;

  int N_threads = 680;
  int N_blocks;
  if ((Nm % N_threads) == 0)
  {
    N_blocks = (Nm / N_threads);
  }
  else
  {
    N_blocks = (Nm / N_threads) + 1;
  }
  dim3 Threads(N_threads);
  dim3 Blocks(N_blocks);

  HANDLE_ERROR(cudaMalloc((void **)&dev_mA, sizeof(long double) * Nm * Nm));
  HANDLE_ERROR(cudaMalloc((void **)&dev_vB, sizeof(long double) * Nm));
  
  cudaMemcpy(dev_mA, mA, sizeof(long double) * Nm * Nm, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_vB, vB, sizeof(long double) * Nm, cudaMemcpyHostToDevice);
  for (k = 0; k < Nm; k++)
  {
    // dgE = mA[Nm * k + k];
    // for (j = k; j < Nm; j++)
    //   mA[Nm * k + j] /= dgE;
    // vB[k] /= dgE;
    cuda_gauss_step << < Blocks, Threads >> > (dev_mA, dev_vB, Nm, k);
    
  }
  cudaMemcpy(mA, dev_mA, sizeof(long double) * Nm* Nm, cudaMemcpyDeviceToHost);
  cudaMemcpy(vB, dev_vB, sizeof(long double) * Nm, cudaMemcpyDeviceToHost);
	HANDLE_ERROR(cudaFree(dev_mA));
	HANDLE_ERROR(cudaFree(dev_vB));
  
}



void bck_one_th()
{
  long int i, j;

  // vX[Nm - 1] = vB[Nm - 1] / mA[Nm * (Nm - 1) + Nm - 1];
  // for (i = Nm - 2; i >= 0; i--) {
  //   vX[i] = vB[i];
  //   for (j = i + 1; j < Nm; j++)
  //     vX[i] -= mA[Nm * i + j] * vX[j];
  // }

  for(i = Nm - 1; i >= 0; i--){
      long double coeff = mA[i * (Nm) + i];
      vX[i] = vB[i] / coeff;

      for(j = i - 1; j >= 0; j-- ){
          vB[j] -= vB[i] * mA[Nm * j + i] / coeff;
          mA[Nm * j + i] = 0;
      }
  }
}

/* generating */

void all_gen() {
  long int i;
  mA = (long double *)calloc(Nm * Nm, sizeof(long double));
  vX = (long double *)calloc(Nm, sizeof(long double));
  vB = (long double *)calloc(Nm, sizeof(long double));

  for (i = 0; i < Nm; i++) {
    vX[i] = startT;
  }
}

void regen_m_x() {
  long i, j, corr = 0;
  int flag, flagl = 1, crr = 1;

  for (i = 0; i < Nm * Nm; i++)
    mA[i] = 0.0;
  for (i = 0; i < Nm; i++)
    vB[i] = 0.0;
  //
  edge = LENE * nodePerL;
  flag = NoX - edge;
  for (i = 0; i < Nm; i++) {
    vB[i] = vX[i] / deltaT;
    if (i == flagl - 1)
      vB[i] += Tle * a / deltaX * deltaX; // Левая граница
    if (i == flag - 1) {
      vB[i] += Tre * a / deltaX * deltaX; // Правая граница
      flagl = flag + 1;
      flag += NoX - (edge--) + crr;
    }
    if (edge <= 0) {
      edge = 0;
      crr = 0;
    }
  }
  int cnt = 0;
  edge = LENE * nodePerL;
  for (j = 0; j < NoY; j++) {
    for (i = 0; i < NoX - edge; i++) {

      if ((j >= NoY / 4 && j < NoY * 3 / 4 - 1) && (i == NoX / 8)) {
        mA[Nm * (corr + i) + corr + i - 1] = 1;
        mA[Nm * (corr + i) + corr + i] = -(1 + deltaX);
        vB[cnt] = 0.0;
      } else if ((j >= NoY / 4 && j < NoY * 3 / 4 - 1) &&
                 (i == NoX - LENE * nodePerL)) {
        mA[Nm * (corr + i) + corr + i + 1] = 1;
        mA[Nm * (corr + i) + corr + i] = -(1 + deltaX);
        vB[cnt] = 0.0;
      } else {
        mA[Nm * (corr + i) + corr + i] =
            1.0 / deltaT + 2.0 * a / deltaX * deltaX;
        if ((NoX + i) % (NoX) != 0)
          mA[Nm * (corr + i) + corr + i - 1] = -a / deltaX * deltaX;
        if ((NoX - edge + i + 1) % (NoX - edge) != 0)
          mA[Nm * (corr + i) + corr + i + 1] = -a / deltaX * deltaX;
      }

      cnt++;
    }

    corr += NoX - edge;
    if (--edge < 0)
      edge = 0;
  }
}

void regen_m_y() {
  long i, j, corr = 0;
  int cnt = 0, localCnt = 0, cntR = 0;

  for (i = 0; i < Nm * Nm; i++)
    mA[i] = 0.0;
  for (i = 0; i < Nm; i++)
    vB[i] = 0.0;
  for (i = 0; i < Nm; i++) {
    vB[i] = vX[i] / deltaT;
  }
  edge = LENE * nodePerL;
  int oldEdge = edge;
  for (j = 0; j < NoY; j++) {
    for (i = 0; i < NoX - edge; i++) {
      if (((j > 0) && (cnt != NoX - oldEdge + localCnt)) || cntR >= NoX) {
        mA[Nm * (corr + i) + corr + i - (NoX - oldEdge)] = -a / deltaX * deltaX;

        ;
      } else {
        cntR++;
        vB[cnt] += Tre * a / deltaX * deltaX;
      }
      // mA[Nm * (corr + i) + corr + i + (NoX - edge)] = -a / deltaX * deltaX;

      if ((i > NoX / 8 && i < NoX - LENE * nodePerL) &&
          (j >= NoY / 4 && j < NoY * 3 / 4 - 1)) {
        vB[cnt] = 0;
      }
      mA[Nm * (corr + i) + corr + i] = 1.0 / deltaT + 2.0 * a / deltaX * deltaX;

      if (j < NoY - 1) {
        if ((j == NoY / 4 - 1) &&
            (i > NoX / 8 && (i < NoX - LENE * nodePerL))) {
          mA[Nm * (corr + i) + corr + i] = -(1 + deltaX);
          mA[Nm * (corr + i) + corr + i + NoX - edge] = 0;
          mA[Nm * (corr + i) + corr + i - NoX + oldEdge] = 1;
          vB[cnt] = 0;
        } else if ((j == NoY * 3 / 4) &&
                   (i > NoX / 8 && (i < NoX - LENE * nodePerL))) {
          mA[Nm * (corr + i) + corr + i] = 1;
          mA[Nm * (corr + i) + corr + i + NoX - edge] = 0;
          mA[Nm * (corr + i) + corr + i - NoX + oldEdge] = -(1 + deltaX);
          vB[cnt] = 0;
        }
      }

      if (j == NoY - 1) // Теплоизолированность внизу
      {
        vB[cnt] = 0;
        mA[Nm * (corr + i) + corr + NoX - edge + i - (NoX - edge)] = 1;
        if (i == NoX - 1 - edge) {
          mA[Nm * (corr + i) + corr + i - 1 - (NoX - oldEdge)] = -1;
        }
      }

      cnt++;
    }
    corr += NoX - edge;
    localCnt += NoX - edge;
    if (--edge < 0)
      edge = 0;
    if (j > 0)
      oldEdge = edge + 1;
  }
}

/* print */
void mtrx_print() {
  long int i, j;

  for (i = 0; i < Nm; i++) {
    for (j = 0; j < Nm; j++) {
      printf("%.3f ", (double)mA[Nm * i + j]);
      // if (j % NoX == 0 && j != 0)
      //   printf("\n");
    }

    printf("= %.1f", (double)vB[i]);
    printf("\n\n");
  }
}

void sltn_print() {
  long int i;

  for (i = 0; i < Nm; i++) {
    printf("%.3f\t", (double)vX[i]);
  }
  printf("\n");
}

void sltn_pl_print() {
  long int i, j, cnt = 0;

  edge = LENE * nodePerL;
  for (j = 0; j < NoY; j++) {
    for (i = 0; i < NoX - edge; i++) {
      printf("%.1f ", (double)vX[cnt++]);
    }
    printf("\n");
    edge--;
    if (edge < 0)
      edge = 0;
  }
}


/* main */
int main(int argc, char *argv[]) {
  long int i, j, k, ret;
  double cpu_time = 0;
  double gpu_time = 0;

  char string[80];

  ofstream myfile;
  myfile.open("result.txt");

  nodePerL =1;
  modTime = 10;

  Nth = 1; // number of blocks

  L = 1 * LENX;
  H = 1 * LENY;
  Nm = LENY * LENX * nodePerL * nodePerL / 1;
  NoX = LENX * nodePerL;
  NoY = LENY * nodePerL;
  edge = LENE * nodePerL;

  for (i = 1; i <= LENE * nodePerL; i++) {
    Nm -= i;
  }
  // deltaX = H * L / Nm;
  deltaX = 1;

  printf("Step %lf \n", (double)deltaX);

  deltaT = 1;
  startT = T0;
  leftEdge = 1;
  Tle = T1;
  rightEdge = 1;
  Tre = T2;


  printf("Nm: %ld\n", Nm);
  all_gen();

  for (i = 0; i < modTime / deltaT; i++) // 10*
  {
    regen_m_x();
    clock_t start_time = time_start();
    cuda_frw_one_th();
    bck_one_th();
    gpu_time += time_stop() - start_time;

    regen_m_y();

    start_time = time_start();
    cuda_frw_one_th();
    bck_one_th();
    gpu_time += time_stop() - start_time;

    edge = LENE * nodePerL;
    ret = 0;
    for (j = 0; j < NoX - edge + 1; j++) {
      string[0] = '\0';
      sprintf(string, "%.3f ", (double)Tre); //Верхняя строка
      // fwrite(string, 1, strlen(string), fds);
      myfile << string ;
    }
    for (j = NoX - edge; j < NoX + 1; j++) {
      // fwrite("-nan ", 1, 5, fds);
      myfile << "-nan ";
    }
    // fwrite("\n", 1, 1, fds);
    myfile << "\n";

    for (j = 0; j < NoY; j++) {
      string[0] = '\0';
      sprintf(string, "%.3f ", (double)Tle); // Левая граница
      // fwrite(string, 1, strlen(string), fds);
      myfile << string;

      for (k = 0; k < NoX - edge; k++) {
        string[0] = '\0';
        if ((k > NoX / 8 && k < NoX - LENE * nodePerL) &&
            (j >= NoY / 4 && j < NoY * 3 / 4 - 1)) {
          // fwrite("-nan ", 1, 5, fds);
          myfile << "-nan ";

          ret++;
        } else {
          sprintf(string, "%.3f ", (double)vX[ret++]);
          // fwrite(string, 1, strlen(string), fds);
          myfile << string;
        }
      }

      string[0] = '\0';
      sprintf(string, "%.3f ", (double)Tre);
      // fwrite(string, 1, strlen(string), fds);
      myfile << string;

      for (k = NoX - edge; k < NoX; k++) {
        // fwrite("-nan ", 1, 5, fds);
        myfile << "-nan ";
      }
      // fwrite("\n", 1, 1, fds);
      myfile << "\n";

      edge--;
      if (edge < 0)
        edge = 0;
    }
    // fwrite("\n\n", 1, 2, fds);
    myfile << "\n\n";
  }
  printf("GPU solution time: %.3f\n", gpu_time / 1000.0);
  for (i = 0; i < modTime / deltaT; i++) // 10*
  {
    regen_m_x();

    clock_t start_time = time_start();
    frw_one_th();
    bck_one_th();
    cpu_time += time_stop() - start_time;

    regen_m_y();

    start_time = time_start();
    frw_one_th();
    bck_one_th();
    cpu_time += time_stop() - start_time;

    // time += time_stop();

    edge = LENE * nodePerL;
    ret = 0;
    for (j = 0; j < NoX - edge + 1; j++) {
      string[0] = '\0';
      sprintf(string, "%.3f ", (double)Tre); //Верхняя строка
      // fwrite(string, 1, strlen(string), fds);
      myfile << string ;
    }
    for (j = NoX - edge; j < NoX + 1; j++) {
      // fwrite("-nan ", 1, 5, fds);
      myfile << "-nan ";
    }
    // fwrite("\n", 1, 1, fds);
    myfile << "\n";

    for (j = 0; j < NoY; j++) {
      string[0] = '\0';
      sprintf(string, "%.3f ", (double)Tle); // Левая граница
      // fwrite(string, 1, strlen(string), fds);
      myfile << string;

      for (k = 0; k < NoX - edge; k++) {
        string[0] = '\0';
        if ((k > NoX / 8 && k < NoX - LENE * nodePerL) &&
            (j >= NoY / 4 && j < NoY * 3 / 4 - 1)) {
          // fwrite("-nan ", 1, 5, fds);
          myfile << "-nan ";

          ret++;
        } else {
          sprintf(string, "%.3f ", (double)vX[ret++]);
          // fwrite(string, 1, strlen(string), fds);
          myfile << string;
        }
      }

      string[0] = '\0';
      sprintf(string, "%.3f ", (double)Tre);
      // fwrite(string, 1, strlen(string), fds);
      myfile << string;

      for (k = NoX - edge; k < NoX; k++) {
        // fwrite("-nan ", 1, 5, fds);
        myfile << "-nan ";
      }
      // fwrite("\n", 1, 1, fds);
      myfile << "\n";

      edge--;
      if (edge < 0)
        edge = 0;
    }
    // fwrite("\n\n", 1, 2, fds);
    myfile << "\n\n";
  }
  printf("CPU solution time: %.3f\n", cpu_time / 1000.0);

  printf("TIME RATE: %.3f\n", cpu_time / gpu_time);

  myfile.close();
  free(mA);
  free(vX);
  free(vB);
  return (0);
}
