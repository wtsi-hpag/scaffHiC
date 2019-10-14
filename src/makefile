# Makefile for scaff10x 
CC= gcc
CFLAGS= -O2 -march=x86-64 -mtune=generic
LFLAGS= -lm -pthread -lz

SOURCES=scaffHiC.c scaffHiC_agp2AGP.c scaffHiC_AGPbuild.c scaffHiC_cover.c scaffHiC_ctg2AGP.c scaffHiC_fastq.c scaffHiC_insert.c scaffHiC_lengthdis.c scaffHiC_motif1.c scaffHiC_motif2.c scaffHiC_motifm.c scaffHiC_motreads.c scaffHiC_offset.c scaffHiC_orient.c scaffHiC_output.c scaffHiC_pairs.c scaffHiC_pairs2.c scaffHiC_pairs-chr.c scaffHiC_pairs-cov.c scaffHiC_pairs-map.c scaffHiC_PCRdup2.c scaffHiC_PCRdup3.c scaffHiC_PCRdup.c scaffHiC_pmatrix.c scaffHiC_proPair.c scaffHiC_RDplace.c scaffHiC_reads.c scaffHiC_rename.c scaffHiC_scaffsort.c scaffHiC_scf2agp.c scaffHiC_scf2scf.c scaffHiC_screen.c scaffHiC_screen-cover.c scaffHiC_translo.c scaffHiC_translo-chr.c

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
	@cp scaff-bin/scaffHiC .

cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
