CC = /opt/matlab/opt/app-sync/matlab/bin/mex
CFLAGS = -ggdb 
#-pedantic -Wall -Wextra  -Wno-unused-parameter 
#LIB =  

OBJS   = hadamard.o 
HFILES = ../hadamard.h 


all: fastwht.mexa64 

fastwht.mexa64: $(HFILES) $(OBJS) fastwht.cpp  Makefile
	$(CC)   $(OBJS) fastwht.cpp 
hadamard.o : ../hadamard.cpp  
	$(CC) CFLAGS='-ggdb' -c ../hadamard.cpp 

