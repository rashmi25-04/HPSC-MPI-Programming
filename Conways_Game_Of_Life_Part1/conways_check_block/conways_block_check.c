// Conway's Game of Life
// Main Executable Program
//
// CSCI 4576/5576 High Performance Scientific Computing
// Michael Oberg, modified from code supplied by Dr. Matthew Woitaszek

// Rashmi Oak
/* program accepts 4 command line arguments
* 1 -> -block or -check for block and checkerboard implementation respectively
* 2 -> Number of iterations
* 3 -> .pgm input file
* 4 -> If and when count total bugs 0 if you dont want bugs to be counted and <Number> if you want bugs to be counted
* Note: for this implementation 2nd Argument has to be greater than 4th argument
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
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


int main(int argc, char* argv[]) {
  // Initialize MPI
  double start,end,elapsed_time,local_elapsed;
  char *pgm;
  int count;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  int total_steps;
  total_steps = atoi(argv[2]);
  if(rank == 0)
  printf("Number of steps = %d\n",total_steps );
  //take innput .pgm file as command line srgument
  pgm = argv[3];
  
  //if and when to count
  if(atoi(argv[4]) == 0){
    if(rank ==0)
    printf("Counter is turned off : Dont count the bugs\n");
  }
  else{
    count = atoi(argv[4]);
    if(count > total_steps){
        if(rank == 0)
        printf("Count should be less than STEPS, Enter Again\n");
        MPI_Finalize();
        return 1;
      }
    if(rank ==0){

    printf("Count No of bugs for %d\n", count);
    }
  }

  if(strcmp(argv[1],"-block") == 0){
  // Get the communicator and process information
  
  // Print rank and hostname
  MPI_Get_processor_name(my_name, &my_name_len);
  //printf("Rank %i is running on %s\n", rank, my_name );
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
  if(!readpgm(pgm))
  {
    if( rank==0 )
    pprintf( "An error occured while reading the pgm file\n" );
    MPI_Finalize();
    return 1;
  }
  int j=0;
    for( int y=1; y<local_height+1; y++ ){
      for( int x=1; x<local_width+1; x++ ){
        if( field_a[ y * field_width + x ] ){
          j++;
        }
      }
    }
    int total_initial;
    MPI_Allreduce( &j, &total_initial, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    if( rank==0 )
      printf( "%d Block Implementation : total initial buggies\n", total_initial );
    int flag;
  for(steps = 0; steps<total_steps; steps++){
    flag = flag +1;
    for( int x=0; x<field_width; x++ )
      {
        field_a[x] = DEAD;
        field_a[(field_height-1) * field_width + x] = DEAD;
      }
      for(int y=0; y<field_height; y++){
        field_a[y * field_width] = DEAD;
      }
      for(int y=0; y<field_height; y++){
        field_a[y * field_width + (field_width - 1) ] = DEAD;
      }
      for( int x=0; x<field_width * field_height; x++ ){
        //printf("Rank: %d field_a[%d] = %d\n",rank,x,field_a[x] );
      }
      if(rank < np-1){
        MPI_Send(&field_a[(field_height-2)*field_width],field_width, MPI_INT, rank+1, 0 ,MPI_COMM_WORLD);
      }
      if(rank > 0){
        MPI_Send(&field_a[1 * field_width],field_width, MPI_INT, rank-1, 1 ,MPI_COMM_WORLD);
      }
      if(rank >0){
        MPI_Recv(&field_a[0], field_width, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
      }
        
        if(rank < np-1){
          MPI_Recv(&field_a[((field_height-1) * field_width)], field_width, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
        }
        int x,y;
        for(y=0; y<field_height; y++ ){
          for(x=0; x<field_width; x++ ){
            //printf("Rank:%d field_a[%d]=%d\t", rank, (y*field_width + x),field_a[y*field_width + x]);
          }
            //printf("\n");
        }
        
    for( int y=1; y<field_height-1; y++ ){
      for( int x=1; x<field_width-1; x++ ){
          
          int ll = y * field_width + x;
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
   // pprintf( "%i local buggies\n", i );
    int total;
    MPI_Allreduce( &i, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    if( rank==0 && flag == count)
      printf( "%d Block implementation : total buggies and count = %d\n", total,flag );
    
    }
       
  // Free the fields
  if( field_a != NULL ) free( field_a );
  if( field_b != NULL ) free( field_b );
  // Finalize MPI and terminate
  if( rank==0 )
    pprintf( "Block Implementation Terminating normally\n" );
  }


  else{
  
  MPI_Get_processor_name(my_name, &my_name_len);
  //printf("Rank %i is running on %s\n", rank, my_name );
  // Initialize the pretty printer
  init_pprintf( rank );
  pp_set_banner( "main" );
  if( rank==0 )
  pprintf( "Welcome to Conway's Game of Life!\n" );
  // Determine the partitioning
  nrows = sqrt(np);
  
  ncols = sqrt(np);
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
  /********************** set up grid **************************************/
  int old_rank;
  int dimensions[2];
  int wrap_around[2];
  int coordinates[2];
  MPI_Comm comm;

  /* Set up Global Grid Information */
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

  /* We assume p is a perfect square */
  int q = (int) sqrt((double)np);
  dimensions[0] = dimensions[1] = q;

  /* We want a circular shift in second dimension. */
  /* Don't care about first                        */
  wrap_around[0] = wrap_around[1] = 1;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Cart_coords(comm,rank, 2,coordinates);
  my_row = coordinates[0];
  my_col = coordinates[1];
  //printf("Rank: %d my_col = %d my_row = %d\n",rank,my_col,my_row );
  

  if(!readpgm(pgm))
  {
    if( rank==0 )
    pprintf( "An error occured while reading the pgm file\n" );
    MPI_Finalize();
    return 1;
  }
  //printf("Rank= %d field_width= %d field_height= %d local_width= %d local_height= %d\n",rank,field_width, field_height, local_width,local_height );
  
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
   
    int total_initial;
    MPI_Allreduce( &j, &total_initial, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    if( rank==0 )
      printf( "%d Checkerboard Implementation : total initial buggies\n", total_initial );
  int flag;
  for(steps = 0; steps<total_steps; steps++){
    flag = flag +1;
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

       
      //printf("Rank : %d coordinates[0] = %d and coordinates[1] = %d\n",rank,coordinates[0],coordinates[1] );
     /***************** MPI communication for checkerboard implementation ************************/
      /************** row communication starts here **********************************************************/
      if(coordinates[0] < q-1){
      //printf("Rank %d is inside MPI_Send\n",rank );
      MPI_Send(&field_a[(field_height-2)*field_width],field_width,MPI_INT,rank+q,0,MPI_COMM_WORLD);
      }
      if(coordinates[0] > 0){
        //printf("Rank %d is inside MPI_Recv\n",rank );
        MPI_Recv(&field_a[0], field_width, MPI_INT, rank-q, 0, MPI_COMM_WORLD, &status);
        for(int m=0; m<field_width*field_height; m++){
        //printf("Inside MPI_Recv: Rank: %d field_a[%d] = %d\n", rank,m,field_a[m]);
        }
      }
      for( int x=0; x<field_width * field_height; x++ ){
        //printf("Rank: %d field_a[%d] = %d\n",rank,x,field_a[x] );
      }
      if(coordinates[0] >0){
        //printf("Rank %d is inside MPI_Send\n",rank );
        MPI_Send(&field_a[1 * field_width],field_width, MPI_INT, rank-q, 1 ,MPI_COMM_WORLD);
      }
      if(coordinates[0] < q-1){
        //printf("Rank %d is inside MPI_Recv\n",rank );
        MPI_Recv(&field_a[(field_height-1) * field_width], field_width, MPI_INT, rank+q, 1, MPI_COMM_WORLD, &status);
        for(int m=0; m<field_width*field_height; m++){
        //printf("Inside MPI_Recv: Rank: %d field_a[%d] = %d\n", rank,m,field_a[m]);
        }
      }
      /************************* row communication ends here *************************************/
      /************************* column communication starts here ********************************/
      
      if(coordinates[1] < q-1){
        //printf("Rank %d is inside MPI_Send\n",rank );
        for(int z=0; z<field_height; z++)
          MPI_Send(&field_a[(z* field_width) + field_width-2],1,MPI_INT,rank+1,z,MPI_COMM_WORLD);//sedond last column to rank + 1 
          //MPI_Send(&z,1,MPI_INT,rank+1,z,MPI_COMM_WORLD);
        }
      if(coordinates[1] > 0){
        //printf("Rank %d is inside MPI_Recv\n",rank );
        
        for(int c=0; c<field_height; c++)
          MPI_Recv(&field_a[c*field_width],1,MPI_INT,rank-1,c,MPI_COMM_WORLD,&status); //receive second last column of rank - 1 into frist column
        /*for(int m=0; m<field_height; m++){
          //printf("Inside MPI_Recv: Rank: %d dummy_field_a[%d] = %d\n", rank,m,dummy_field_a[m]);
        }*/
        for(int m=0; m<field_width*field_height; m++){
        //printf("Inside MPI_Recv: Rank: %d field_a[%d] = %d\n", rank,m,field_a[m]);
        }
      }
      if(coordinates[1] > 0){
        //printf("Rank %d is inside MPI_Send\n",rank );
        for(int z=0; z<field_height; z++)
        MPI_Send(&field_a[(z*field_width) + 1],1,MPI_INT,rank-1,z,MPI_COMM_WORLD);//send first column to rank -1
      }
      if(coordinates[1] < q-1){
        //printf("Rank %d is inside MPI_Recv\n",rank );
        

        for(int t=0; t<field_height; t++)
        MPI_Recv(&field_a[t*field_width + field_width-1],1,MPI_INT,rank+1,t,MPI_COMM_WORLD,&status);// receive first column of rank +1 into second last column
        /*for(int m=0; m<field_height; m++){
          printf("Inside MPI_Recv: Rank: %d dummy_field_a1[%d] = %d\n", rank,m,dummy_field_a1[m]);
        }*/
        for(int m=0; m<field_width*field_height; m++){
        //printf("Inside MPI_Recv: Rank: %d field_a[%d] = %d\n", rank,m,field_a[m]);
        }
      }
      /*for(int m=0; m<field_width*field_height; m++){
        printf("Inside MPI_Recv: Rank: %d field_a[%d] = %d\n", rank,m,field_a[m]);
        }*/
     /************************ column communication ends here **********************************/
      

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
          //neighbors[ll] = neighbors[ll] + field_a[ll_N]; 

          int lx_S = (x);
          int ly_S = (y+1);
          int ll_S = (ly_S * field_width + lx_S );
         // printf("Rank : %d ll_S = %d field_a[ll_S] = %d\n",rank,ll_S,field_a[ll_S] );
          //neighbors[ll] = neighbors[ll] + field_a[ll_S];


          int lx_E = (x+1);
          int ly_E = (y);
          int ll_E = (ly_E * field_width + lx_E );
          //printf("Rank : %d ll_E = %d field_a[ll_E] = %d\n",rank,ll_E ,field_a[ll_E]);
          //neighbors[ll] = neighbors[ll] + field_a[ll_E];
          
          int lx_W = (x-1);
          int ly_W = (y);
          int ll_W = (ly_W * field_width + lx_W );
          //printf("Rank : %d ll_W = %d field_a[ll_W] = %d\n",rank,ll_W ,field_a[ll_W]);
          //neighbors[ll] = neighbors[ll] + field_a[ll_W];

          int lx_NW = (x-1);
          int ly_NW = (y-1);
          int ll_NW = (ly_NW * field_width + lx_NW );
          //printf("Rank : %d ll_NW = %d field_a[ll_NW] = %d\n",rank,ll_NW,field_a[ll_NW] );
          //neighbors[ll] = neighbors[ll] + field_a[ll_NW];

          int lx_NE = (x+1);
          int ly_NE = (y-1);
          int ll_NE = (ly_NE * field_width + lx_NE );
          //printf("Rank : %d ll_NE = %d field_a[ll_NE] = %d\n",rank,ll_NE ,field_a[ll_NE]);
          //neighbors[ll] = neighbors[ll] + field_a[ll_NE];


          int lx_SW = (x-1);
          int ly_SW = (y+1);
          int ll_SW = (ly_SW * field_width + lx_SW );
          //printf("Rank : %d ll_SW = %d field_a[ll_SW] = %d\n",rank,ll_SW,field_a[ll_SW] );
          //neighbors[ll] = neighbors[ll] + field_a[ll_SW];

          
          int lx_SE = (x+1);
          int ly_SE = (y+1);
          int ll_SE = (ly_SE * field_width + lx_SE );
          //printf("Rank : %d ll_SE = %d field_a[ll_SE] = %d\n",rank,ll_SE,field_a[ll_SW] );
          //neighbors[ll] = neighbors[ll] + field_a[ll_SE];

          
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
   // pprintf( "%i local buggies\n", i );
    int total;
    MPI_Allreduce( &i, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    if( rank==0 && flag == count)
      printf( "%d Checkerboard Implementation total buggies and count = %d\n", total,count );
    }
   // Free the fields
  if( field_a != NULL ) free( field_a );
  if( field_b != NULL ) free( field_b );

  // Finalize MPI and terminate
  if( rank==0 )
    pprintf( "Checkerboard Implementation : Terminating normally\n" );
  }

  MPI_Barrier(MPI_COMM_WORLD); 
  end = MPI_Wtime();
  local_elapsed = end - start;
  MPI_Reduce(&local_elapsed, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(rank == 0)
    printf("Time Elapsed = %f\n", elapsed_time);
  MPI_Finalize();
  return 0;
} 

