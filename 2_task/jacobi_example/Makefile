# Copyright (c) 2016-2017 The University of Tennessee and The University
#                         of Tennessee Research Foundation.  All rights
#                         reserved.
# AUTHOR: George Bosilca
# 
CC=${HOME}/SOFT/OpenMPi/bin/mpicc
LD=${HOME}/SOFT/OpenMPi/bin/mpicc

MPIDIR=${HOME}/SOFT/OpenMPi
MPIINC=-I$(MPIDIR)/include
MPILIB=-lpthread -L$(MPIDIR)/lib -lmpi

CFLAGS=-g
LDFLAGS= $(MPILIB) -g

LINK=$(LD)

APPS=jacobi_noft jacobi_bckpt

all: $(APPS)

jacobi_noft: jacobi_cpu_noft.o main.o
	$(LINK) -o $@ $^ -lm

jacobi_bckpt: jacobi_cpu_bckpt.o main.o
	$(LINK) -o $@ $^ -lm

%.o: %.c header.h
	$(CC) -c $(CFLAGS) -o $@ $<

clean:
	rm -f *.o $(APPS) *~
