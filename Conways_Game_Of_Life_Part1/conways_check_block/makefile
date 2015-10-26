########################################################################### 
# Compiler
########################################################################### 

CC = mpicc
CCFLAGS = -g -Wall -std=c99
ifeq ($(DEBUG),on)
	CCFLAGS += -DDEBUG
endif

########################################################################### 
# Files
########################################################################### 

C_FILES = conways_block_check.c pgm.c pprintf.c
O_FILES = conways_block_check.o pgm.o pprintf.o

########################################################################### 
# Targets
########################################################################### 

all: conways_block_check serial

conways_block_check: $(O_FILES)
	$(CC) -o conways_block_check $(O_FILES) $(LDFLAGS)

.PHONY: clean
clean:		
	/bin/rm -f core $(O_FILES) conways_block_check serial 

conways_base.o: pgm.o pprintf.o
.c.o:		
	$(CC) $(CCFLAGS) -c -o $*.o $*.c

.cpp.o:		
	$(CC) $(CCFLAGS) -c -o $*.o $*.cpp

########################################################################### 
# Dependencies
########################################################################### 

# All of the object files depend on the globals, so rebuild everything if they
# change!
*.o: globals.h

# Nothing really depends on the pprintf prototypes, but just be safe
*.o: pprintf.h

# Conway depends on PGM utilities
conways_base.o: pgm.h pprintf.h

serial:
	gcc -o serial serial_gol.c
