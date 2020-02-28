# Makefile for scaff10x 
CC= gcc
CFLAGS= -O2 -march=x86-64 -mtune=generic
LFLAGS= -lm -pthread -lz

SOURCES=scaffHiC_screen.c scaffHiC_output.c scaffHiC_fastq.c scaffHiC_rename.c scaffHiC_outbreak.c scaffHiC_superAGP.c scaffhic.c breakhic.c scaffHiC_PCRdup.c scaffHiC_PCRdup2.c scaffHiC_PCRdup3.c scaffHiC_RDplace.c scaffHiC_proPair.c scaffHiC_AGPbuild.c scaffHiC_pmatrix.c scaffHiC_insert.c scaffHiC_lengthdis.c scaffHiC_orient.c scaffHiC_scf2agp.c scaffHiC_break.c scaffHiC_breakClean.c scaffHiC_BKplace.c scaffHiC_BKplace2.c scaffHiC_CIplace.c scaffHiC_scf2scf.c scaffHiC_getPair.c scaffHiC_breakPair.c scaffHiC_breakPair2 scaffHiC_breakAGP.c scaffHiC_breakTags.c scaffHiC_reads.c scaffHiC_screen-cover.c scaffHiC_translo.c scaffHiC_translo-chr.c

OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,scaff-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS) 


all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o $(LFLAGS); cp $$exe scaff-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '


clean: 
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp scaff-bin/scaffhic .
	@cp scaff-bin/breakhic .

cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
