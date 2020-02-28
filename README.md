# scaffhic v1.1
Pipeline for genome scaffolding by modelling distributions of HiC pairs. Also breakhic is added to the pipeline.

Pipeline steps:
        
    Scaffolding with scaffhic:
      1 HiC pairs are aligned against the target assembly using BWA mem with a tag -5SPM for HiC reads
      2 PCR duplications and other artefacts are screened and removed from further data processing
      3 Based on read distributions, Contig Distance Index (CDI) is used to quantify the likelihood of contig joins
      4 Order and orientation are determined by pair numbers on contigs
      5 Relation matrix and partner matrix are used to reduce join errors
      6 6 iterations are carried out for genome scaffolding
      7 Contigs are arranged into a scaffolding structure
      8 Distributions of HiC insert mapping sizes are plotted as a png image file  
      9 HiC map is produced based on the scaffolded assembly

### Download and Compile:

    $ git clone  https://github.com/wtsi-hpag/scaffhic.git 
    $ cd scaffhic
    $ ./install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		


#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) is downloaded and compiled by scaffhic
The compression code pigz (https://github.com/madler/pigz) is downloaded and compiled by scaffhic
The HiC map builder PretextMap (https://github.com/wtsi-hpag/PretextMap) is downloaded by scaffhic
The HiC map viewer PretextView (https://github.com/wtsi-hpag/PretextView) is downloaded by scaffhic


### Run the pipelines

#### Run scaffhic for genome scaffolding:
           $ /full/path/to/scaffhic/src/scaffhic -nodes <nodes> -score <score> \
	   	 -depth <depth_of_matrix> -length <contig_length> \
                 -fq1 Input_read_1 -fq2 Input_read_2 draft-assembly.fasta genome-hic.fasta \
           
	       Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             depth:        number of columns of Partner Matrix [ default = 50 ]
             score:        score of Contig Distance Index (CDI) [ default CDI = 200 ]
             length:       minimum contig length considering for scaffolding [ default = 500000 ]


#### Run scaffhic for HiC map:
           $ /full/path/to/scaffhic/src/scaffhic -nodes <nodes> -score <score> \
	   	 -depth <depth_of_matrix> -length <contig_length> \
                 -map genome-hic.map -plot genome-hic-length.png -file 0 \ 
                 -fq1 Input_read_1 -fq2 Input_read_2 genome-hic.fasta test_scaffold_file \
           
	       Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             depth:        number of columns of Partner Matrix [ default = 50 ]
             score:        score of Contig Distance Index (CDI) [ default CDI = 200 ]
             length:       minimum contig length considering for scaffolding [ default = 500000 ]


#### View HiC map:
           PretextView is a desktop application for viewing pretext contact maps
           The package has different linux and Mac versions (linux version included in scaffhic)
           
	   $ /full/path/to/scaffhic/src/seqbin-bin/PretextView genome-hic.map \
	    
           More information, go to https://github.com/wtsi-hpag/PretextView


#### Run scaffhic on new parameter settings with alignment file align.dat:
           $ /full/path/to/scaffhic/src/scaffhic -nodes <nodes> -score <score> \
	   	 -depth <depth_of_matrix> -length <contig_length> \
                 -data /full/path/to/tmp_rununik_27152/align.dat  \
                 draft-assembly.fasta genome-hic.fasta \
          
           This saves alignment time and can get the results quickly. 


#### Example if you have an assembly and HiC read files 
##### genome_assembly.fasta  GM12878-HiC_1.fastq.gz GM12878-HiC_2.fastq.gz 
           $ /full/path/to/scaffhic/src/scaffhic -nodes 30 -score 200 -depth 50 -length 500000 \
                 -fq1 GM12878-HiC_1.fastq.gz -fq2 GM12878-HiC_2.fastq.gz genome_assembly.fasta genome-scaffhic.fasta \
          
           $ /full/path/to/scaffhic/src/scaffhic -nodes 30 -score 200 -depth 50 -length 500000 \
                 -map genome-hic.map -plot genome-hic-length.png -file 0 \
                 -fq1 GM12878-HiC_1.fastq.gz -fq2 GM12878-HiC_2.fastq.gz genome-scaffhic.fasta genome-hic2.fasta \
          
           After the above two steps, You now have  
              1. Scaffolded assembly file: genome-scaffhic.fasta; 
              2. HiC map: genome-hic.map;  
              3. HiC length distribution image: genome-hic-length.png 
              4. Assembly file produced in the visualization step:  genome-hic2.fasta 
                 This file is not used, but it might be slightly better than genome-scaffhic.fasta 
                 The visualization results are based on genome-scaffhic.fasta 
                 
          To view it on a desktop linux or Mac laptop use PretextView https://github.com/wtsi-hpag/PretextView  \
          $ ./PretextView genome-hic.map

#### Run breakhic to identify assembly breakpoints 
           $ /full/path/to/scaffhic/src/breakhic -nodes <nodes> -grid <grid size> \
                 -fq1 Input_read_1 -fq2 Input_read_2 draft-assembly.fasta break-hic.fasta \
           
               Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             grid:         grid size to search breakpoints, the smaller, the more CPU time [ default = 100 ]

#### Run breakhic on new parameter settings with alignment file align.dat:
           $ /full/path/to/scaffhic/src/breakhic -nodes <nodes> -grid <grid size> \
                 -data /full/path/to/tmp_rununik_27152/align.dat  \
                 draft-assembly.fasta break-hic.fasta \
          
           This saves alignment time and can get the results quickly. 

#### Example if you have an assembly and HiC read files 
##### genome_assembly.fasta  GM12878-HiC_1.fastq.gz GM12878-HiC_2.fastq.gz 
           $ /full/path/to/scaffhic/src/breakhic -nodes 30 -grid 100 \
                 -fq1 GM12878-HiC_1.fastq.gz -fq2 GM12878-HiC_2.fastq.gz genome_assembly.fasta genome-break.fasta \

           After the above step, You now have
              1. Assembly file with broken contigs: genome-break.fasta
              2. Breakpoint information file: genome-break.fasta.brk

Break1: 70 547700 217464711 3156 5 0      \
Break1: 70 854600 217464711 3447 5 0      \
Break1: 70 1155600 217464711 2736 6 0     \
Break2: 70 119259859 217464711 3431 5 200 \
Break2: 117 1340007 72126970 389 38 200   \
Break2: 117 4228188 72126970 1157 12 200  \
Break2: 117 28301466 72126970 556 26 200  \
Break1: 118 46043000 76603601 1388 9 0    \
Break2: 118 60007606 76603601 179 60 200  \
Break2: 118 69709038 76603601 116 60 200  \

2 - Contig/scaffold index;     \
3 - Breakpoint offset;         \
4 - Contig/scaffold length;    \
5 - Average HiC coverage;      \
6 - Break likelihood value QV; \
7 - Breakpoint nature: contig break (0); scaffold break or break at a gap (200) \
 

