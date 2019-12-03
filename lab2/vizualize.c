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
// #include <unistd.h>


void visualize() {
  int i;
  int num_time = 3;
  FILE *gnu = popen("gnuplot -persist", "w");
  // fprintf(gnu, "set view map\n");
  // fprintf(gnu, "set palette rgbformulae 22,13,-31\n");

  fprintf(gnu, "set pm3d map\n");
  // fprintf(gnu, "set palette rgbformulae 21,22,23\n");
  fprintf(gnu, "set palette rgbformulae 22,13,-31\n");

  fprintf(gnu, "set size 0.8,1.0\n");

  // fprintf(gnu, "set yrange [17:0]\n", Yrange);
  fprintf(gnu, "set yrange [%d:0]\n", LENY);

  fprintf(gnu, "set xrange [0:%d]\n", LENX - 1);

  for (i = 0; i < num_time; i++) {
    fprintf(gnu, "splot \"result.txt\" index %d matrix \n", i);
    fprintf(gnu, "pause 0.2\n");
  }

  fflush(gnu);
}

int main(){
    visualize();
}