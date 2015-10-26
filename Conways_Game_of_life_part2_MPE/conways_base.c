// Conway's Game of Life
// Main Executable Program
//
// CSCI 4576/5576 High Performance Scientific Computing
// Michael Oberg, modified from code supplied by Dr. Matthew Woitaszek

// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

// Include global variables. Only this file needs the #define
#define __MAIN 
#include "globals.h"
#undef __MAIN

// User includes
#include "pprintf.h"
#include "pgm.h"
#define DEAD 0
#define ALIVE 1
#define TWO 2
#define THREE 3

#define NSTEPS 1000 /* number of time steps */

int main(int argc, char* argv[]) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Get the communicator and process information
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  // Print rank and hostname
  MPI_Get_processor_name(my_name, &my_name_len);
  printf("Rank %i is running on %s\n", rank, my_name );

  // Initialize the pretty printer
  init_pprintf( rank );
  pp_set_banner( "main" );
  if( rank==0 )
  pprintf( "Welcome to Conway's Game of Life!\n" );

  
  // Determine the partitioning
  nrows = np;
  ncols = 1;

  
  int steps;
  
  MPI_Status status;
  if( np != nrows * ncols )
  {
    if( rank==0 )
      pprintf("Error: %ix%i partitioning requires %i np (%i provided)\n", 
        nrows, ncols, nrows * ncols, np );
    MPI_Finalize();
    return 1;
  }
  my_col = 0;
  my_row = rank;

  if(!readpgm("life.pgm"))
  {
    if( rank==0 )
    pprintf( "An error occured while reading the pgm file\n" );
    MPI_Finalize();
    return 1;
  }
  int j=0;
    for( int y=1; y<local_height+1; y++ )
    {
      for( int x=1; x<local_width+1; x++ )
      {
        if( field_a[ y * field_width + x ] )
        {
          j++;
        }
      }
    }
   // pprintf( "%i local buggies\n", i );
    int total_initial;
    MPI_Allreduce( &j, &total_initial, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    if( rank==0 )
      pprintf( "%i total initial buggies\n", total_initial );
  // we have field_a and field_b with ghost nodes
  /* set ghost nodes to be DEAD */
      //printf("Rank: %d : : field_width = %d field_height = %d local_width = %d local_height = %d\n",rank,field_width,field_height,local_width,local_height );
  
     //printf("field_a[%d] = %d\n",7,field_a[12] );
  for(steps = 0; steps<NSTEPS; steps++){
    for( int x=0; x<field_width; x++ )
      {
        field_a[x] = DEAD;
        field_a[(field_height-1) * field_width + x] = DEAD;
      }
   

      
      for(int y=0; y<field_height; y++)
      {
        field_a[y * field_width] = DEAD;
      }
      for(int y=0; y<field_height; y++)
      {
        field_a[y * field_width + (field_width - 1) ] = DEAD;
      }

       for( int x=0; x<field_width * field_height; x++ ){
        //printf("Rank: %d field_a[%d] = %d\n",rank,x,field_a[x] );
      }
     //while(1){}

      
        if(rank < np-1){
          //printf("inside mpi send\n");
          MPI_Send(&field_a[(field_height-2)*field_width],field_width, MPI_INT, rank+1, 0 ,MPI_COMM_WORLD);
        }
        if(rank > 0){
          MPI_Send(&field_a[1 * field_width],field_width, MPI_INT, rank-1, 1 ,MPI_COMM_WORLD);
        }
        if(rank >0){
          //dummy_field_a = (int *)malloc( field_width * sizeof(int));
          MPI_Recv(&field_a[0], field_width, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
          //int count = 0;
        
        }
        
        
        if(rank < np-1){
          //dummy_field_a1 = (int *)malloc( field_width * sizeof(int));
          MPI_Recv(&field_a[((field_height-1) * field_width)], field_width, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
            
        }
        //while(1){}
        int x,y;
        for(y=0; y<field_height; y++ )
        {
        for(x=0; x<field_width; x++ )
          {
            
          //printf("Rank:%d field_a[%d]=%d\t", rank, (y*field_width + x),field_a[y*field_width + x]);
          }
          //printf("\n");
        }
        //while(1){}

    for( int y=1; y<field_height-1; y++ )
    {
      for( int x=1; x<field_width-1; x++ )
      {
          //int neighbors = 0;
          int ll = y * field_width + x;
          //neighbors[ll] = 0;
          int lx_N = (x);
          int ly_N = (y-1);
          int ll_N = (ly_N * field_width + lx_N );
          //printf("Rank : %d ll_N = %d field_a[ll_N] = %d\n",rank,ll_N,field_a[ll_N] );
          
          int lx_S = (x);
          int ly_S = (y+1);
          int ll_S = (ly_S * field_width + lx_S );
         // printf("Rank : %d ll_S = %d field_a[ll_S] = %d\n",rank,ll_S,field_a[ll_S] );
          

          int lx_E = (x+1);
          int ly_E = (y);
          int ll_E = (ly_E * field_width + lx_E );
          //printf("Rank : %d ll_E = %d field_a[ll_E] = %d\n",rank,ll_E ,field_a[ll_E]);
          
          
          int lx_W = (x-1);
          int ly_W = (y);
          int ll_W = (ly_W * field_width + lx_W );
          //printf("Rank : %d ll_W = %d field_a[ll_W] = %d\n",rank,ll_W ,field_a[ll_W]);
          
          int lx_NW = (x-1);
          int ly_NW = (y-1);
          int ll_NW = (ly_NW * field_width + lx_NW );
          //printf("Rank : %d ll_NW = %d field_a[ll_NW] = %d\n",rank,ll_NW,field_a[ll_NW] );
          
          int lx_NE = (x+1);
          int ly_NE = (y-1);
          int ll_NE = (ly_NE * field_width + lx_NE );
          //printf("Rank : %d ll_NE = %d field_a[ll_NE] = %d\n",rank,ll_NE ,field_a[ll_NE]);
          

          int lx_SW = (x-1);
          int ly_SW = (y+1);
          int ll_SW = (ly_SW * field_width + lx_SW );
          //printf("Rank : %d ll_SW = %d field_a[ll_SW] = %d\n",rank,ll_SW,field_a[ll_SW] );
          //neighbors[ll] = neighbors[ll] + field_a[ll_SW];

          
          int lx_SE = (x+1);
          int ly_SE = (y+1);
          int ll_SE = (ly_SE * field_width + lx_SE );
          //printf("Rank : %d ll_SE = %d field_a[ll_SE] = %d\n",rank,ll_SE,field_a[ll_SW] );
          
          
          int neighbors = field_a[ll_N] + field_a[ll_S] + field_a[ll_E] + field_a[ll_W] + field_a[ll_NW] + field_a[ll_NE] + field_a[ll_SW] + field_a[ll_SE] ;

          //printf("Rank : %d neighbors = %d of ll %d\n",rank, neighbors ,ll);

          
          if (field_a[ll] == ALIVE && (neighbors < TWO || neighbors > THREE)){
            field_b[ll]= DEAD;
            //printf("rank %d field_b[%d] = %d\n",rank,ll,field_b[ll] );
          }
          else if(field_a[ll] == ALIVE && (neighbors == TWO || neighbors == THREE)){
            field_b[ll] = ALIVE;
            //printf("rank %d field_b[%d] = %d\n",rank,ll,field_b[ll] );
          }
          else if (field_a[ll] == DEAD && neighbors == THREE){
            field_b[ll] = field_a[ll]+1;
            //printf("rank %d field_b[%d] = %d\n",rank,ll,field_b[ll] );
          }
          else{
            field_b[ll] = field_a[ll];
            //printf("rank %d field_b[%d] = %d\n",rank,ll,field_b[ll] );
          }
         
          
      }
    }
   
    
    for(int p=1; p<field_height-1; p++){
      for(int q=1; q<field_width-1; q++){
        int ll = p*field_width + q;
        field_a[ll] = field_b[ll];
      }
    }
    

    // Count the life forms. Note that we count from [1,1] - [height+1,width+1];
    // we need to ignore the ghost row!
    int i=0;
    for( int y=1; y<local_height+1; y++ )
    {
      for( int x=1; x<local_width+1; x++ )
      {
        if( field_a[ y * field_width + x ] )
        {
          i++;
        }
      }
    }
    int *temp_field_a = (int *)malloc( local_width * local_height * sizeof(int));
    int temp = 0;
    for(int j=1; j<local_height+1; j++){
      for(int k=1;k<local_width+1; k++){
        temp_field_a[temp] = field_a[j * field_width + k];
        //printf("Rank : %d temp_field_a[%d] = %d\n",rank,temp,temp_field_a[temp] );
        temp ++;
        
      }
    }

   // pprintf( "%i local buggies\n", i );
    int total;
    MPI_Allreduce( &i, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    if( rank==0 )
      pprintf( "%i total buggies\n", total );

    
    int height_main = local_height * np;
    int width_main = local_width * 1;
    int *anim = (int *)malloc( width_main * height_main * sizeof(int));
    //printf("height_main = %d and width_main = %d\n",height_main,width_main );

    
    /**************** writing into .pgm file and create animation *******************/
    MPI_Gather(temp_field_a,local_height*local_width,MPI_INT,anim,local_height*local_width,MPI_INT,0,MPI_COMM_WORLD); //Gather local matrices
    /********************************************************************************/
    if(rank == 0){
      for(int k=0;k< width_main*height_main; k++){
        //printf("anim[%d] = %d\n",k,anim[k] );
      }
    }
    /********** write into .pgm file *********************/
    if(rank ==0){
    char filename[1000];
    FILE *pgmfile;
    sprintf(filename, "test%d.pbm", steps);
    pgmfile = fopen(filename, "wb");
    if (pgmfile == NULL) {
        perror("cannot open file to write");
        exit(EXIT_FAILURE);
    }
 
    fprintf(pgmfile, "P1 ");
    fprintf(pgmfile, "%d %d ", width_main, height_main);
    fprintf(pgmfile, "%d ", 1);
    for(int y=0;y<height_main;y++)
          {
            for(int x=0;x<width_main;x++)
            {
              fputc(anim[y*width_main +x]+48, pgmfile );   
                          
            }
            
          }
        fclose(pgmfile);
    /*****************************************************/
      }
  }

    
    
  // Free the fields
  if( field_a != NULL ) free( field_a );
  if( field_b != NULL ) free( field_b );

  // Finalize MPI and terminate
  if( rank==0 )
    pprintf( "Terminating normally\n" );
  MPI_Finalize();
  return 0;
} 

