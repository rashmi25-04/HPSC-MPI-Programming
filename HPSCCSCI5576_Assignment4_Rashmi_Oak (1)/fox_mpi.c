/* Fox's Algorithm */
/* Rashmi Oak */
/* References : http://www.lac.inpe.br/~stephan/CAP-372/Fox_example.pdf 
	          : Textbook -  Parallel Programming with MPI by Peter S Pacheco
   				Chapter 7: Communicators & Topology

 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "mpi.h"

#define SIZE 4 /* dimension of the input matrix */

int A[SIZE][SIZE];
int B[SIZE][SIZE];

typedef struct {
	int p; /* number of processors */
	MPI_Comm comm; /* handle to global grid communicator */
	MPI_Comm row_comm; /* row communicator */
	MPI_Comm col_comm; /* column communicator */
	int q; /* dimension of the grid, = sqrt(p) */
	int my_row; /* row position of a processor in a grid */
	int my_col; /* column position of a procesor in a grid */
	int my_rank; /* rank within the grid */
}grid_info;

/* normal matrix multiplication stuff */

void mm_mult(int **a, int **b, int **c, int size)
{
	int i,j,k;
       
	int **flag = (int**) malloc(size*sizeof(int*));
	for(i=0;i<size;i++)
		*(flag+i)=(int*) malloc(size*sizeof(int));

	for(i=0;i<size;i++)
	{
			for(j=0;j<size;j++)
			{
				flag[i][j]=0;
				for(k=0;k<size;k++){
					flag[i][j]=flag[i][j]+ (a[i][k] * b[k][j]);
				}
			}
	}
	
	for(i=0;i<size;i++)
		for(j=0;j<size;j++)
			c[i][j]+=flag[i][j];
	
}

void SetupGrid(grid_info *grid)
{
	int old_rank;
	int dimensions[2];
	int wrap_around[2];
	int coordinates[2];
	int free_coords[2];
	
	/* get the overall information before overlaying cart_grid */

	MPI_Comm_size(MPI_COMM_WORLD,&(grid->p));
	MPI_Comm_rank(MPI_COMM_WORLD,&old_rank);
	
	/* Assumption: p is a perfect square */
	grid->q=(int)sqrt((double)grid->p);
	/* set the dimensions */
	dimensions[0]=dimensions[1]=grid->q;
	
	/* we want a torus on the second dimension, so set it appropriately */

	wrap_around[0]=0;
	wrap_around[1]=1;
	
	MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,wrap_around,1,&(grid->comm));
	/* since we have set reorder to true, this might have changed the ranks */
	MPI_Comm_rank(grid->comm,&(grid->my_rank));
	/* get the cartesian coordinates for the current process */
	MPI_Cart_coords(grid->comm,grid->my_rank,2,coordinates);
	/* set the coordinate values for the current coordinate */
	grid->my_row=coordinates[0];
	grid->my_col=coordinates[1];

        /* create row communicators */
	free_coords[0]=0;
	free_coords[1]=1; /* row is gonna vary */
	MPI_Cart_sub(grid->comm,free_coords,&(grid->row_comm));
	
        /* create column communicators */
	free_coords[0]=1;
	free_coords[1]=0; /* row is gonna vary */
	MPI_Cart_sub(grid->comm,free_coords,&(grid->col_comm));
	
}


void frombuff(int *buff,int **a,int buffsize, int row, int col){
	
  	if(buffsize!=row*col)
	{
		//printf("transfer_data_from_buf: buffer size does not match matrix size!\n");
		exit(1);
	}
	int count=0, i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++){
			a[i][j]=buff[count];
			count++;
		}
	}
}

void tobuff(int *buff,int **a,int buffsize, int row, int col){
	
  	if(buffsize!=row*col)
	{
		//printf("transfer_data_to_buf: buffer size does not match matrix size!");
		exit(1);
	}
	int count=0, i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++){
			buff[count]=a[i][j];
			count++;
		}
	}
}




void Fox(int n,grid_info *grid,int **a, int **b, int **c)
{
	int **tempa;
	int *buffer; /* buffer for Bcast & send_recv */
	
	int stage;
	int root;
	int new_size; /* = n/q */
	int source;
	int dest;
	int i;
	MPI_Status status;
	
	new_size=n/grid->q;
	
	/* Initialize tempa */
	tempa=(int**) malloc(new_size*sizeof(int*));
	for(i=0;i<new_size;i++)
		*(tempa+i)=(int*) malloc(new_size*sizeof(int));
	/* initialize buffer */
	buffer=(int*)malloc(new_size*new_size*sizeof(int));

        /* we are gonna shift the elements of matrix b upwards with the column fixed */
	source = (grid->my_row+1) % grid->q; /* pick the emmediately lower element */
	dest= (grid->my_row+grid->q-1) % grid->q; /* move current element to immediately upper row */
	
	
	for(stage=0;stage<grid->q;stage++)
	{
		root=(grid->my_col+stage)%grid->q;
		if(root==grid->my_col)
		{
			tobuff(buffer,a,new_size*new_size, new_size,new_size);
			MPI_Bcast(buffer, new_size*new_size, MPI_INT,root,grid->row_comm);
			frombuff(buffer,a,new_size*new_size, new_size,new_size);
		
			mm_mult(a,b,c,new_size);
		}else
		{
			tobuff(buffer,tempa,new_size*new_size, new_size,new_size);
			MPI_Bcast(buffer,new_size*new_size, MPI_INT,root,grid->row_comm);
			frombuff(buffer,tempa,new_size*new_size, new_size,new_size);
			
			mm_mult(tempa,b,c, new_size);
		}
		tobuff(buffer,b,new_size*new_size, new_size,new_size);
		MPI_Sendrecv_replace(buffer,new_size*new_size, MPI_INT,dest,0,source,0,grid->col_comm,&status);
		frombuff(buffer,b,new_size*new_size, new_size,new_size);
	}


}




int main(int argc, char *argv[])
{	
	
	int i,j,size_g;
	int **localA;
	int **localB;
	int **localC;
	int my_rank;
	 double start,end,elapsed_time,local_elapsed;

	 
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

	grid_info grid;
	/*initialize Grid */

	SetupGrid(&grid);
     /* Initialize matrix A & B */
	static int count = 0;
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			A[i][j]=count++;
			B[i][j]=count++;
		}
	}
	
	
	
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<SIZE;j++)
		{
			//printf("A[%d][%d] = %d\n",i,j,A[i][j]);
			//printf("B[%d][%d] = %d\n",i,j,B[i][j]);
		}
	}


	size_g=SIZE/grid.q;
	/* allocate space for the three matrices */		

	
	localA=(int**) malloc(size_g*sizeof(int*));
	localB=(int**) malloc(size_g*sizeof(int*));
	localC=(int**) malloc(size_g*sizeof(int*));
	
	for(i=0;i<size_g;i++)
	{
		*(localA+i)=(int*) malloc(size_g*sizeof(int));
		*(localB+i)=(int*) malloc(size_g*sizeof(int));
		*(localC+i)=(int*) malloc(size_g*sizeof(int));
	}


/* Compute local matrices - Ideally the master should do this & pass it onto all the slaves */
/* At the same time initialize localC to all zeros */

	int row=grid.my_row*size_g;
	int col=grid.my_col*size_g;

	for(i=row;i<row+size_g;i++)
	{
		for(j=col;j<col+size_g;j++)
		{
		     localA[i-(row)][j-(col)]=A[i][j];
			 localB[i-(row)][j-(col)]=B[i][j];
			 localC[i-(row)][j-(col)]=0;
		}
	}
	Fox(SIZE,&grid,localA, localB, localC);

/* print results */
	//printf("rank=%d, row=%d col=%d\n",grid.my_rank,grid.my_row,grid.my_col);
	for(i=0;i<size_g;i++)
	{
		for(j=0;j<size_g;j++)
		{
			//printf("localC[%d][%d]=%d ", i,j,localC[i][j]);
			printf("%d ", localC[i][j]);
		}
		printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD); 
    end = MPI_Wtime();

    local_elapsed = end - start;
    
    MPI_Reduce(&local_elapsed, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(my_rank == 0)
      printf("Time Elapsed = %f\n", elapsed_time);
	MPI_Finalize ();
	exit(0);

}		