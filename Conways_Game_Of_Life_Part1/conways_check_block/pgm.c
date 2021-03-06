// System includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

// User includes
#include "globals.h"
#include "pprintf.h"

typedef enum { false, true } bool; // Provide C++ style 'bool' type in C

bool readpgm( char *filename )
{
  // Read a PGM file into the local task
  //
  // Input: char *filename, name of file to read
  // Returns: True if file read successfully, False otherwise
  //
  // Preconditions:
  //  * global variables nrows, ncols, my_row, my_col must be set
  //
  // Side effects:
  //  * sets global variables local_width, local_height to local game size
  //  * sets global variables field_width, field_height to local field size
  //  * allocates global variables field_a and field_b

  pp_set_banner( "pgm:readpgm" );

  // Open the file
  if( rank==0 )
    pprintf( "Opening file %s\n", filename );
  FILE *fp = fopen( filename, "r" );
  if( !fp )
  {
    pprintf( "Error: The file '%s' could not be opened.\n", filename );
    return false;
  }

  // Read the PGM header, which looks like this:
  //  |P5        magic version number
  //  |900 900       width height
  //  |255         depth
  char header[10];
  int width, height, depth;
  int rv = fscanf( fp, "%6s\n%i %i\n%i\n", header, &width, &height, &depth );
  if( rv != 4 )
  {
    if(rank==0) 
      pprintf( "Error: The file '%s' did not have a valid PGM header\n", 
        filename );
    return false;
  }
  if( rank==0 )
    pprintf( "%s: %s %i %i %i\n", filename, header, width, height, depth );

  // Make sure the header is valid
  if( strcmp( header, "P5") )
  {
    if(rank==0) 
      pprintf( "Error: PGM file is not a valid P5 pixmap.\n" );
    return false;
  }
  if( depth != 255) //for life.pgm 255
  {
    if(rank==0) 
      pprintf( "Error: PGM file has depth=%i, require depth=255 \n", 
        depth );
    return false;
  }

  // Make sure that the width and height are divisible by the number of
  // processors in x and y directions

  if( width % ncols )
  {
    if( rank==0 )
      pprintf( "Error: %i pixel width cannot be divided into %i cols\n", 
          width, ncols );
    return false;
  }
  if( height % nrows )
  {
    if( rank==0 )
      pprintf( "Error: %i pixel height cannot be divided into %i rows\n",
          height, nrows );
    return false;
  }

  // Divide the total image among the local processors
  local_width = width / ncols;
  local_height = height / nrows;

  // Find out where my starting range it
  int start_x = local_width * my_col;
  int start_y = local_height * my_row;

  pprintf( "Hosting data for x:%03i-%03i y:%03i-%03i\n", 
      start_x, start_x + local_width,
      start_y, start_y + local_height );

  // Create the array!
  field_width = local_width + 2;
  field_height = local_height + 2;
  field_a = (int *)malloc( field_width * field_height * sizeof(int));
  field_b = (int *)malloc( field_width * field_height * sizeof(int));
  //neighbors = (int *)malloc( field_width * field_height * sizeof(int));
  // Read the data from the file. Save the local data to the local array.
  int b, ll, lx, ly;
  for( int y=0; y<height; y++ )
  {
    for( int x=0; x<width; x++ )
    {
      // Read the next character
      b = fgetc( fp );
      if( b == EOF )
      {
        pprintf( "Error: Encountered EOF at [%i,%i]\n", y,x );
        return false;
      }

      // From the PGM, black cells (b=0) are bugs, all other 
      // cells are background 
      if( b==0 )
      {
        b=1;
      }
      else
      {
        b=0;
      }
      /*if( b==48 || b==49 )
      {
        b=1;
        //printf("[%d][%d]b = %d\n",y,x,b );
        //old[y][x] = b;
      }
      else
      {
        b=0;
        //printf("[%d][%d]b = %d\n",y,x,b );
        //old[y][x] = b;
    
      }*/
        

      // If the character is local, then save it!
      if( x >= start_x && x < start_x + local_width &&
        y >= start_y && y < start_y + local_height )
      {
        // Calculate the local pixels (+1 for ghost row,col)
        lx = x - start_x + 1;
        ly = y - start_y + 1;
        ll = (ly * field_width + lx );
        field_a[ ll ] = b;
        field_b[ ll ] = b;
       
        //printf("Rank : %d field_a[%d] = %d\n",rank,ll,field_a[ll] );
      } // save local point

    } // for x
  } // for y

  fclose( fp );

  pp_reset_banner();
  return true;
}
