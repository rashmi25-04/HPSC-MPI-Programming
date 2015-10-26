#include <stdio.h>
#include <sys/time.h>
#define SIZE 2

int main(int argc, char* argv[]) 
{
  int c, d, k, sum = 0;
  int first[SIZE][SIZE], second[SIZE][SIZE], multiply[SIZE][SIZE];
 int i,j;
 struct timeval  tv1, tv2;
gettimeofday(&tv1, NULL);
  
      static int n = 0;
        for (i = 0; i < SIZE; i++){
            for (j = 0; j < SIZE; j++){
                first[i][j] = n++;
                second[i][j] = n++;
            }
        }
 
    for (c = 0; c < SIZE; c++) {
      for (d = 0; d < SIZE; d++) {
        for (k = 0; k < SIZE; k++) {
          sum = sum + first[c][k]*second[k][d];
        }
 
        multiply[c][d] = sum;
        sum = 0;
      }
    }
 
    printf("Product of entered matrices:-\n");
 
    for (c = 0; c < SIZE; c++) {
      for (d = 0; d < SIZE; d++)
        printf("%d\t", multiply[c][d]);
 
      printf("\n");
    }
  gettimeofday(&tv2, NULL);

printf ("Total time = %f seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
 
  return 0;
}