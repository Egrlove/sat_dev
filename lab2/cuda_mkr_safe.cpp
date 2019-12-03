#define LENX 40
#define LENY 25
#define LENE 25
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
// #include <unistd.h>

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
// struct timeval tv1, tv2, dtv;
// struct timezone tz;
int ij;

// void time_start() { gettimeofday(&tv1, &tz); }
clock_t time_start() { return clock(); }
clock_t time_stop() { return clock(); }

// long time_stop() {
//   gettimeofday(&tv2, &tz);
//   dtv.tv_sec = tv2.tv_sec - tv1.tv_sec;
//   dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
//   if (dtv.tv_usec < 0) {
//     dtv.tv_sec--;
//     dtv.tv_usec += 1000000;
//   }
//   return dtv.tv_sec * 1000 + dtv.tv_usec / 1000;
// }

/* Gauss */

void frw_one_th() {
  long int i, j, k;
  long double dgE;

  for (k = 0; k < Nm; k++) {
    dgE = mA[Nm * k + k];
    for (j = k; j < Nm; j++)
      mA[Nm * k + j] /= dgE;
    vB[k] /= dgE;

    for (i = k + 1; i < Nm; i++) {
      dgE = mA[Nm * i + k];
      for (j = k; j < Nm; j++)
        mA[Nm * i + j] -= mA[Nm * k + j] * dgE;
      vB[i] -= vB[k] * dgE;
    }
  }
}

void bck_one_th() {
  long int i, j;

  vX[Nm - 1] = vB[Nm - 1];

  for (i = Nm - 2; i >= 0; i--) {
    vX[i] = vB[i];
    for (j = i + 1; j < Nm; j++)
      vX[i] -= mA[Nm * i + j] * vX[j];
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
  long int i, j, k, ret, time = 0;
  char string[80];

  ofstream myfile;
  myfile.open("result.txt");

  nodePerL = 1;
  modTime = 3;

  Nth = 1; // number of blocks

  L = 1 * LENX;
  H = 1 * LENY;
  Nm = LENY * LENX * nodePerL * nodePerL / 1;
  NoX = LENX * nodePerL;
  NoY = LENY * nodePerL;
  edge = LENE * nodePerL;
  int setY = edge + 2;

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

    clock_t start_time = time_start();
    regen_m_x();
    frw_one_th();
    bck_one_th();
    time += time_stop() - start_time;

    time_start();
    regen_m_y();
    frw_one_th();
    bck_one_th();
    time += time_stop() - start_time;

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

  // sltn_pl_print();

  printf("Solution time: %.3f\n", time / 1000.0);
  myfile.close();
  free(mA);
  free(vX);
  free(vB);
  return (0);
}
