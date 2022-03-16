#/*************************************************************************
#*
#*                     MAKEFILE FOR SPIN.C                                
#*                                                                         
#*************************************************************************/



#######################################################################
#                                                                     #
# 1. Specify C compiler and ANSI option:                              #
#                                                                     #      
####################################################################### 

#Linux
CC=gcc

#/*************************************************************************
#*                     SPECIFY GRID SIZE
#*************************************************************************/

#STANDARD
LOWSIZE=-DMDIV=65 -DSDIV=129

#HIGH
#SIZE=-DMDIV=101 -DSDIV=201

#VERY HIGH
HIGHSIZE=-DMDIV=151 -DSDIV=301
#HIGHSIZE=-DMDIV=101 -DSDIV=201
#HIGHSIZE=-DMDIV=65 -DSDIV=129

#VERY VERY HIGH
SIZE=-DMDIV=201 -DSDIV=401
#HIGHSIZE=-DMDIV=201 -DSDIV=401
#401 601

#/*************************************************************************
#*                     COMPILING FLAGS
#*************************************************************************/


# DEBUGGING OPTION
MY_OWN =-g3 -Wall 
#MY_OWN = -g3 -Wall -DDEBUG

#MY_OWN =-O3

#/*************************************************************************
#*                    SOURCE AND OBJECT MACROS
#*************************************************************************/


SOBJ=nsss.o findmodel.o equil.o equil_util.o nrutil.o stableorbit.o surface.o interpol.o
MOBJ=findmodel.o equil.o equil_util.o nrutil.o stableorbit.o surface.o interpol.o maxmass.o

#/*************************************************************************
#*                    MAIN COMPILING INSTRUCTIONS
#*************************************************************************/
nsss: $(SOBJ)
	$(CC) $(MY_OWN)  $(SIZE) -o nsss $(SOBJ) -lm

maxmass: $(MOBJ)
	$(CC) $(MY_OWN) $(SIZE) -o maxmass $(MOBJ) -lm 

nsss.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h stableorbit.h nsss.c makefile interpol.h
	$(CC) -c $(MY_OWN) $(CFLAGS) $(COPTFLAGS) $(SIZE)  nsss.c 

maxmass.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h maxmass.c makefile
	$(CC) -c $(MY_OWN) $(CFLAGS) $(COPTFLAGS) $(SIZE)  maxmass.c 

findmodel.o: consts.h struct.h nrutil.h equil.h equil_util.h surface.h findmodel.h findmodel.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(SIZE)   findmodel.c

equil.o:equil.h equil_util.h nrutil.h consts.h equil.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(SIZE)   equil.c

equil_util.o:equil_util.h nrutil.h consts.h equil_util.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(SIZE)   equil_util.c

nrutil.o:nrutil.h nrutil.c makefile
	$(CC) -c $(COPTFLAGS)   nrutil.c

surface.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h surface.h stableorbit.h surface.c makefile
	$(CC) -c $(COPTFLAGS) $(SIZE)  surface.c 

stableorbit.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h surface.h stableorbit.h stableorbit.c makefile
	$(CC) -c $(COPTFLAGS) $(SIZE)  stableorbit.c

interpol.o: consts.h struct.h nrutil.h equil.h equil_util.h surface.h findmodel.h interpol.h interpol.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(SIZE)   interpol.c


