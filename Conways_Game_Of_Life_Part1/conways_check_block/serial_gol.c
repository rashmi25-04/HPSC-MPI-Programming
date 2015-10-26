/***********************

Conway Game of Life

serial version

************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define DEAD 0
#define ALIVE 1
#define TWO 2
#define THREE 3

#define NSTEPS 2000   /* number of time steps */

int main(int argc, char *argv[]) {

  int i, j, n, im, ip, jm, jp, ni, nj, neighbors, isum;
  int **old, **new;  
  float x;
  char header[10];
    int width, height, depth;

  
  /**************************** Reading from pgm file *************************************/
    char *filename="life.pgm";
    FILE *fp = fopen( filename, "r" );
    if( !fp )
    {
    printf( "Error: The file '%s' could not be opened.\n", filename );
    return 0;
    }
    
    int rv = fscanf( fp, "%6s\n%i %i\n%i\n", header, &width, &height, &depth );
    if( rv != 4 )
    {
      printf( "Error: The file '%s' did not have a valid PGM header\n", filename );
       return 0;
    }
    printf( "%s: %s %i %i %i\n", filename, header, width, height, depth );
     // Make sure the header is valid
    //if( strcmp( header, "P1") )
    if( strcmp( header, "P5") )
    {
      printf( "Error: PGM file is not a valid P5 pixmap.\n" );
      return 0;
    }
    if( depth != 255)
    {
      printf( "Error: PGM file has depth=%i, require depth=255 \n", depth );
      return 0;
    }
    printf("width = %d height = %d depth = %d\n",width,height,depth);
    int b;
    ni = height + 2;  /* add 2 for left and right ghost cells */
      nj = width + 2;
      old = malloc(ni*sizeof(int*));
      new = malloc(ni*sizeof(int*));

      for(i=0; i<ni; i++){
        old[i] = malloc(nj*sizeof(int));
        new[i] = malloc(nj*sizeof(int));
      }
    for( int y=1; y<=height; y++ )
    {
    for( int x=1; x<=width; x++ )
    {
      // Read the next character

      b = fgetc( fp );
      //printf("initial b  = %d\t", b);
      if( b == EOF )
      {
        printf( "Error: Encountered EOF at [%i,%i]\n", y,x );
        return 0;
      }
      
      // From the PGM, black cells (b=0) are bugs, all other 
      // cells are background 
      /*if( b==48 || b==49 )
      {
        b=1;
        //printf("[%d][%d]b = %d\n",y,x,b );
        old[y][x] = b;
      }
      else
      {
        b=0;
        //printf("[%d][%d]b = %d\n",y,x,b );
        old[y][x] = b;
    
      }*/
      if( b==0 )
      {
        b=1;
        old[y][x] = b;
      }
      else
      {
        b=0;
        old[y][x] = b;
    
      }
    }
  }
  
  printf("width = %d height = %d depth = %d\n",width,height,depth);
  for(int p=1; p<=height; p++){
    for(int q=1;q<=width;q++){
      //printf("old[%d][%d] %d\n",p,q,old[p][q] );
    }
  }
  

  

  /*****************************************************************/
  /*  time steps */
  for(n=0; n<NSTEPS; n++){

    
    /* Make extreme left and right column ghost cells DEAD */

    for(i=0; i<=height+1; i++){
      old[i][0] = DEAD;
      old[i][width+1] = DEAD;
    }

    /* Make top and bottom cells ghost cells DEAD  */

    for(j=0; j<=width+1; j++){
      old[0][j] = DEAD;
      old[height+1][j] = DEAD;
    }
  for(i=0; i<=height+1; i++){
        for(j=0; j<=width+1; j++){
      //printf("old_old[%d][%d] = %d\t",i,j,old[i][j] );
          //printf("%d %d %d\t",i,j,old[i][j] );
        }
      }
    for(i=1; i<=width; i++){
      for(j=1; j<=height; j++){
  

    neighbors = old[i-1][j-1] + old[i-1][j] + old[i-1][j+1] + old[i][j-1] + old[i][j+1] + old[i+1][j-1] + old[i+1][j] + old[i+1][j+1];
    
  
  
    //printf("neighbors for old[%d][%d] == %d\n", i,j,neighbors);
    if (old[i][j] == ALIVE && (neighbors < TWO || neighbors > THREE))
        new[i][j] = DEAD;
    else if(old[i][j] == ALIVE && (neighbors == TWO || neighbors == THREE))
        new[i][j] = ALIVE;
    else if (old[i][j] == DEAD && neighbors == THREE)
        new[i][j] = old[i][j]+1;
    else
        new[i][j] = old[i][j];

  
    }
  }

    /* copy new state into old state */

    for(i=1; i<=height; i++){
      for(j=1; j<=width; j++){
  old[i][j] = new[i][j];
      }
    }
    for(i=1; i<=height; i++){
      for(j=1; j<=width; j++){
  //printf("new_old[%d][%d] = %d\n",i,j,old[i][j] );
      }
    }
  
}

  
  isum = 0;
  for(i=1; i<=height; i++){
    for(j=1; j<=width; j++){
      if(new[i][j] == 1)
        isum = isum + 1;
    }
  }
  printf("\nNumber of live cells = %d\n", isum);

  return 0;
}