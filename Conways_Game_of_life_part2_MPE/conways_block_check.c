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
* 5 -> Either -Set <set of iterations eg. 1 2 3 4 5 and so on but should be less than max no of iterations>or -Range <Two numbers from to to eg. 1 to 100>
* Note: for this implementation 2nd Argument has to be greater than 4th argument
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include <mpe.h>
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
#define SIZE 10000


int main(int argc, char* argv[]) {
  // Initialize MPI
  double start,end,elapsed_time,local_elapsed;
  char *pgm;
  int count;
  MPI_Init(&argc, &argv);
  MPI_Pcontrol( 0 );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  int total_steps;
  int event1a, event1b, event2a, event2b,event3a, event3b, event4a, event4b,event5a,event5b,event6a,event6b;//events for MPE
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
  //command line arguments for set and Range
  int array[SIZE];
  int start_cl,end_cl;

  if(strcmp(argv[5],"-Set") == 0){
      //printf("Specify Set Of Iterations\n");
      int temp = 0;
      for(int i=6;i<argc;i++){
        array[temp] = atoi(argv[i]);
        //printf("argv[%d] = %s and array[%d] = %d\n",i,argv[i],temp,array[temp] );
        temp ++;
      }
      start_cl = 0;
      end_cl = argc - 6;
      //printf("start_cl %d and end_cl %d\n",start_cl,end_cl );
    }
    else if(strcmp(argv[5],"-Range") == 0){
      //printf("Specify Range Of Iterations\n");
      start_cl = atoi(argv[6]);
      end_cl = atoi(argv[7]);
      for(int i=start_cl;i<end_cl;i++){
        array[i] = i;
        //printf("array[%d] = %d\n",i,array[i] );  
      }
      //printf("Start_cl =%d to end_cl = %d\n",start_cl,end_cl );
    }// end of command line arguments Set and Range
MPE_Init_log();//Initialize MPE log
 /*
 *   *         user should NOT assign eventIDs directly in MPE_Describe_state()
 *     *                 Get the eventIDs for user-defined STATES(rectangles) from
 *       *                         MPE_Log_get_state_eventIDs() instead of the deprecated function
 *         *                                 MPE_Log_get_event_number().
 *           *                                     */
  MPE_Log_get_state_eventIDs( &event1a, &event1b );
  MPE_Log_get_state_eventIDs( &event2a, &event2b );
  MPE_Log_get_state_eventIDs( &event3a, &event3b );
  MPE_Log_get_state_eventIDs( &event4a, &event4b );
  MPE_Log_get_state_eventIDs( &event5a, &event5b );
  MPE_Log_get_state_eventIDs( &event6a, &event6b );
  if (rank == 0 ) {
        MPE_Describe_state( event1a, event1b, "Preprocessing", "red" );
        MPE_Describe_state( event2a, event2b, "MPI Communication", "orange" );
        MPE_Describe_state( event3a, event3b, "Nighbours calculations and updating field_a", "blue" );
        MPE_Describe_state( event4a, event4b, "MPI-IO", "green" );
        MPE_Describe_state( event5a, event5b, "creating temp vectors", "white" );
        MPE_Describe_state( event6a, event6b, "Global Communication", "cyan" );
  }

  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Pcontrol( 1 );
  /*
 *   *     MPE_Start_log();
 *     *         */


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
  MPE_Log_event( event1a, 0, NULL );//preprocessing start
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
    MPE_Log_event( event6a, 0, NULL );//global comm
    MPI_Allreduce( &j, &total_initial, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    MPE_Log_event( event6b, 0, NULL );
    if( rank==0 )
      printf( "%d Block Implementation : total initial buggies\n", total_initial );
    MPE_Log_event( event1b, 0, NULL );//preprocessing done
    int flag = 0;
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
      MPE_Log_event( event2a, 0, NULL );//MPI communication starts
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
      MPE_Log_event( event2b, 0, NULL );//MPI communication ends
        int x,y;
        for(y=0; y<field_height; y++ ){
          for(x=0; x<field_width; x++ ){
            //printf("Rank:%d field_a[%d]=%d\t", rank, (y*field_width + x),field_a[y*field_width + x]);
          }
            //printf("\n");
        }
     MPE_Log_event( event3a, 0, NULL );//neighbours calculations start 
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
   MPE_Log_event( event3b, 0, NULL );//neighbours calculations end and updating field_a ends
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
    MPE_Log_event( event6a, 0, NULL );
    MPI_Allreduce( &i, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    MPE_Log_event( event6b, 0, NULL );
    if( rank==0 && flag == count)
      printf( "%d Block implementation : total buggies and count = %d\n", total,flag );
    
    // create a temp field that takes local_heightxlocal_width
    MPE_Log_event( event5a, 0, NULL );
    int *temp_field_a = (int *)malloc( local_width * local_height * sizeof(int)); 
    int temp = 0;
    for(int j=1; j<local_height+1; j++){
      for(int k=1;k<local_width+1; k++){
        temp_field_a[temp] = field_a[j * field_width + k];
        //printf("steps : %d Rank : %d temp_field_a[%d] = %d\n",steps,rank,temp,temp_field_a[temp] );
        temp ++;  
      }
    }// end of for loop for temp field
    //create vector for darray
    int *anim_2d = (int*)malloc(local_height * local_width * sizeof(int)) ;
    int count = 0;
    for(int y=1;y<=local_height;y++){
            for(int x=1;x<=local_width;x++){
              anim_2d[count] =  field_a[y*field_width +x]+48;
              //printf("Rank %d and anim_2d[%d] = %d\n",rank,count,anim_2d[count] );
              count ++;            
            }
          }

    int height_main = local_height * np;//calculating total grid height
    int width_main = local_width * 1;//calculating total grid width
    int *anim = (int *)malloc( width_main * height_main * sizeof(int));
    MPE_Log_event( event5b, 0, NULL );
    //collect all local fields
    MPE_Log_event( event6a, 0, NULL ); 
    MPI_Gather(temp_field_a,local_height*local_width,MPI_INT,anim,local_height*local_width,MPI_INT,0,MPI_COMM_WORLD); //Gather local matrices
    MPE_Log_event( event6b, 0, NULL );
    if(rank == 0){
      for(int k=0;k< width_main*height_main; k++){
        //printf("steps : %d anim[%d] = %d\n",steps,k,anim[k] );
      }
    }
      //create .pbm file for animation block implementation
      //printf("start_cl %d and end_cl %d\n",start_cl,end_cl );
      //check if command line argument matches with current iteration
      //if matches create PGM file
      for(int k=start_cl; k<end_cl;k++){
      //printf("inside for loop\n");
      //printf("array === %d\n",array[k]);
    
      if(steps == array[k]){
      MPE_Log_event( event4a, 0, NULL );//MPI-IO
      char filename[1000];
      FILE *pgmfile;
      sprintf(filename, "test%d", steps);
      pgmfile = fopen(filename, "wb");
      if (pgmfile == NULL) {
          perror("cannot open file to write");
          exit(EXIT_FAILURE);
      }
      if(rank ==0){
      fprintf(pgmfile, "P5 \n");
      fprintf(pgmfile, "%d %d\n", width_main, height_main);
      fprintf(pgmfile, "%d \n", 1);
      }
      
      MPI_File  fh;
      MPI_Datatype filetype;
      int gsizes[2];
      int distribs[2];
      int dargs[2];
      int psizes[2];

      gsizes[0] = height_main; // no of rows in global array
      gsizes[1] = width_main; // no of columns in global array
      //printf("height_main %d and width_main %d\n",height_main,width_main );
      distribs[0] = MPI_DISTRIBUTE_BLOCK;
      distribs[1] = MPI_DISTRIBUTE_BLOCK;

      dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;

      psizes[0] = np;
      psizes[1] = 1;


      MPI_Type_create_darray(np,rank,2,gsizes,distribs,dargs,psizes,MPI_ORDER_C,MPI_INT,&filetype);
      MPI_Type_commit(&filetype);
      MPI_File_open(MPI_COMM_WORLD,filename , MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh) ; 
      MPI_File_set_view(fh,11,MPI_INT,filetype,"native",MPI_INFO_NULL);
      int local_array_size = local_height * local_width;
      MPI_File_write_all(fh,anim_2d,local_array_size,MPI_INT,&status);
      MPI_File_close(&fh);
      MPE_Log_event( event4b, 0, NULL );
      }
      }

    MPE_Log_sync_clocks();
    }//end of steps loop
       
  // Free the fields
  if( field_a != NULL ) free( field_a );
  if( field_b != NULL ) free( field_b );
  // Finalize MPI and terminate
  if( rank==0 )
    pprintf( "Block Implementation Terminating normally\n" );
  }//end of block implementation
  
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
  
  MPE_Log_event( event1a, 0, NULL );//pre-processing
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
    MPE_Log_event( event6b, 0, NULL );
    MPI_Allreduce( &j, &total_initial, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    MPE_Log_event( event6b, 0, NULL );
    if( rank==0 )
      printf( "%d Checkerboard Implementation : total initial buggies\n", total_initial );
    MPE_Log_event( event1b, 0, NULL );//pre-processing done

  int flag = 0;
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

      MPE_Log_event( event2a, 0, NULL );
 
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
      
      MPI_Datatype col;
      MPI_Type_vector(field_height,1,field_width,MPI_INT,&col);
      MPI_Type_commit(&col);




      /************************* column communication starts here ********************************/
     
      if(coordinates[1] < q-1){
        
        MPI_Send(field_a + field_width-2,1,col,rank+1,0,comm);
        }
      if(coordinates[1] > 0){
        
        MPI_Send(field_a +1,1,col,rank-1,1,comm);
      }
      if(coordinates[1] > 0){
        
        MPI_Recv(field_a +0,1,col,rank-1,0,comm,&status);
      }
      if(coordinates[1] < q-1){
        
          MPI_Recv(field_a + field_width-1,1,col,rank+1,1,comm,&status);
      }
     /************************ column communication ends here **********************************/
      MPE_Log_event( event2b, 0, NULL );
      MPE_Log_event( event3a, 0, NULL );
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
    MPE_Log_event( event3b, 0, NULL );

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
    MPE_Log_event( event6a, 0, NULL );
    MPI_Allreduce( &i, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    MPE_Log_event( event6b, 0, NULL );
    if( rank==0 && flag == count)
      printf( "%d Checkerboard Implementation total buggies and count = %d\n", total,count );
    //code for animation
    MPE_Log_event( event5a, 0, NULL );
    int *temp_field_a = (int *)malloc( local_width * local_height * sizeof(int)); 
    int temp = 0;
    for(int j=1; j<local_height+1; j++){
      for(int k=1;k<local_width+1; k++){
        temp_field_a[temp] = field_a[j * field_width + k];
        //printf("steps : %d Rank : %d temp_field_a[%d] = %d\n",steps,rank,temp,temp_field_a[temp] );
        temp ++;  
      }
    }// end of for loop for temp field
    //create vector for darray
    int *anim_2d = (int*)malloc(local_height * local_width * sizeof(int)) ;
    int count = 0;
    for(int y=1;y<=local_height;y++){
            for(int x=1;x<=local_width;x++){
              anim_2d[count] =  field_a[y*field_width +x]+48;
              //printf("Rank %d and anim_2d[%d] = %d\n",rank,count,anim_2d[count] );
              count ++;            
            }
          }
    int height_main = local_height * q;//calculating total grid height
    int width_main = local_width * q;//calculating total grid width
    int *anim = (int *)malloc( width_main * height_main * sizeof(int));
    MPE_Log_event( event5b, 0, NULL );
    //collect all local fields
    MPE_Log_event( event6a, 0, NULL );
    MPI_Gather(temp_field_a,local_height*local_width,MPI_INT,anim,local_height*local_width,MPI_INT,0,MPI_COMM_WORLD); //Gather local matrices
    MPE_Log_event( event6b, 0, NULL );
    if(rank == 0){
      for(int k=0;k< width_main*height_main; k++){
        //printf("rank : %d steps : %d anim[%d] = %d\n",rank,steps,k,anim[k] );
      }
    }

    //create .pbm file for animation checkerboard implementation
    for(int k=start_cl; k<end_cl;k++){
      //printf("inside for loop\n");
      //printf("array === %d\n",array[k]);
    //check if command line argument matches with current iteration
    //if matches create PGM file
    if(steps == array[k]){
    MPE_Log_event( event4a, 0, NULL );
    char filename[1000];
    FILE *pgmfile;
    sprintf(filename, "test%d", steps);
    pgmfile = fopen(filename, "wb");
    if (pgmfile == NULL) {
        perror("cannot open file to write");
        exit(EXIT_FAILURE);
    }
    if(rank ==0){
    fprintf(pgmfile, "P1 \n");
    fprintf(pgmfile, "%d %d\n", width_main, height_main);
    fprintf(pgmfile, "%d \n", 1);
    }
      
    MPI_File  fh;
    MPI_Datatype filetype;
    int gsizes[2];
    int distribs[2];
    int dargs[2];
    int psizes[2];

    gsizes[0] = height_main; // no of rows in global array
    gsizes[1] = width_main; // no of columns in global array
    //printf("height_main %d and width_main %d\n",height_main,width_main );
    distribs[0] = MPI_DISTRIBUTE_BLOCK;
    distribs[1] = MPI_DISTRIBUTE_BLOCK;

    dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;

    psizes[0] = q;
    psizes[1] = q;


    MPI_Type_create_darray(np,rank,2,gsizes,distribs,dargs,psizes,MPI_ORDER_C,MPI_INT,&filetype);
    MPI_Type_commit(&filetype);
    MPI_File_open(MPI_COMM_WORLD,filename , MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh) ; 
    MPI_File_set_view(fh,11,MPI_INT,filetype,"native",MPI_INFO_NULL);
    int local_array_size = local_height * local_width;
    MPI_File_write_all(fh,anim_2d,local_array_size,MPI_INT,&status);
    MPI_File_close(&fh);
    MPE_Log_event( event4b, 0, NULL );
    }
    }
 MPE_Log_sync_clocks();
    }//end of steps loop for checkerboard
   // Free the fields
  if( field_a != NULL ) free( field_a );
  if( field_b != NULL ) free( field_b );

  // Finalize MPI and terminate
  if( rank==0 )
    pprintf( "Checkerboard Implementation : Terminating normally\n" );
  }//end of checkerboard implementation
  MPE_Finish_log( "Conways_Game_Of_Life" );
  MPI_Barrier(MPI_COMM_WORLD); 
  end = MPI_Wtime();
  local_elapsed = end - start;
  MPI_Reduce(&local_elapsed, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(rank == 0)
    printf("Time Elapsed = %f\n", elapsed_time);
  MPI_Finalize();
  return 0;
} 

