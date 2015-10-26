#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
double f(double x);

  
int main(int argc, char** argv) {
  int flag,k=1;                         
  double  sum;
  double f(double x)     ;  
  double a = 0.00;
  double b = 2.00;
  int n;
  float h;
  n = atoi(argv[1]) ; 
  h = (b-a)/n;
  sum=f(a); 
  printf("f(local_a)%f\n",sum)   ;                  
  flag=4;
    while (k <= n-1)               
    {
                        
      sum = sum + flag*f(a+k*h); 
      flag=6-flag;    
      k++;                      
    }                   
    

  sum = ( sum + f(b) )*h/3 ;    
  printf("final value of integral %f\n",sum) ;                        
  return sum;

 }

double f(double x)
    {
     double function;
     function = x*x;
     return(function);
    }

