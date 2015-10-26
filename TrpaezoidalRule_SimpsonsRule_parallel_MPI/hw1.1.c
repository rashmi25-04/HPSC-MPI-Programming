#include "mpi.h"
#include <stdio.h>
#include <string.h>

#define BUFSIZE 1024    /* maximum message size in chars */


int  my_rank, size;
int  pro_len;
char pro_name [MPI_MAX_PROCESSOR_NAME ];

void ring_hello_world() {

  char       recv_msg[ BUFSIZE ]; /* MPI_Recv buffer */
  char       msg[] = "Hello World !!!";
  int        prev,next;      /* next     processor ordinal in ring */
  char       send_msg[ BUFSIZE  ]; /* MPI_Send buffer */
  int        bytes_received;     /* byte count of received data */
  MPI_Status status[ sizeof( MPI_Status ) ];

  prev = ( my_rank - 1 + size ) % size;
  next = ( my_rank + 1 )        % size;

  if( my_rank == 0 ) {   /* initiate send, wait for last process */

    printf( "First message is '%s'.\n", msg );
    MPI_Send( msg, strlen( msg ), MPI_CHAR, next, 0, MPI_COMM_WORLD );
    MPI_Recv( recv_msg, BUFSIZE, MPI_CHAR, prev, 0, MPI_COMM_WORLD, status );
    MPI_Get_count( status, MPI_CHAR, &bytes_received );
    recv_msg[ bytes_received ] = '\0';
    printf("bytes_received: %d\n",bytes_received);
    printf( "Final mesage is '%s'.\n", recv_msg );

  } else {      /* receive message, foreward to next process */

    MPI_Recv( recv_msg, BUFSIZE, MPI_CHAR, prev, 0, MPI_COMM_WORLD, status );
    MPI_Get_count( status, MPI_CHAR, &bytes_received );
    recv_msg[ bytes_received ] = '\0';
    printf("message received is %s and received from process  %d\n",recv_msg, my_rank);

    sprintf( send_msg, "%s (%d/%s)", recv_msg, my_rank, pro_name );
    strcpy( send_msg, recv_msg );
    MPI_Send( send_msg, strlen( send_msg ), MPI_CHAR, next, 0, MPI_COMM_WORLD );

  }
}
int main( int argc, char *argv[]) {
  
  double t1, t2; 
  t1 = MPI_Wtime(); 

    MPI_Init( &argc, &argv );
    
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    
    MPI_Get_processor_name( pro_name, &pro_len );
    printf( "Process : %d processor name : %s\n", my_rank, pro_name );
    ring_hello_world();
    MPI_Finalize();
    t2 = MPI_Wtime(); 
printf( "Elapsed time is %f\n", t2 - t1 ); 
    
}
