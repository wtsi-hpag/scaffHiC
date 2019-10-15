# scaffHiC v1.0
Pipeline for genome scaffolding by modelling distributions of HiC pairs.

Pipeline steps:
        
    Scaffolding with scaffHiC:
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

    $ git clone  https://github.com/wtsi-hpag/scaffHiC.git 
    $ cd scaffHiC
    $ ./install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		


#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) is downloaded and compiled by scaffHiC.

### Run the pipelines

#### Run scaffHiC for genome scaffolding:
           $ /full/path/to/scaffHiC/src/scaffHiC -nodes <nodes> -score <score> \
	   	 -depth <depth_of_matrix> -length <contig_length> \
                 draft-assembly.fasta Input_read_1 Input_read_2 genome-hic.fasta \
           
	       Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             depth:        number of columns of Partner Matrix [ default = 50 ]
             score:        score of Contig Distance Index (CDI) [ default CDI = 200 ]
             length:       minimum contig length considering for scaffolding [ default = 500000 ]


#### Run scaffHiC for HiC map:
           $ /full/path/to/scaffHiC/src/scaffHiC -nodes <nodes> -score <score> \
	   	 -depth <depth_of_matrix> -length <contig_length> \
                 -map genome-hic.map -plot genome-hic-length.png -file 0 \ 
                 genome-hic.fasta Input_read_1 Input_read_2 test_scaffold_file \
           
	       Parameters:
             nodes:        number of CPUs requested  [ default = 30 ]
             depth:        number of columns of Partner Matrix [ default = 50 ]
             score:        score of Contig Distance Index (CDI) [ default CDI = 200 ]
             length:       minimum contig length considering for scaffolding [ default = 500000 ]


#### View HiC map:
           PretextView is a desktop application for viewing pretext contact maps and it has linux and Mac versions.
           
	   $ /full/path/to/scaffHiC/src/seqbin-bin/PretextView genome-hic.map \
	    
           More information, go to https://github.com/wtsi-hpag/PretextView
        
#### Run scaffHiC on new parameter settings with alignment file align.dat:
           $ /full/path/to/scaffHiC/src/scaffHiC -nodes <nodes> -score <score> \
	   	 -depth <depth_of_matrix> -length <contig_length> \
                 -data /full/path/to/tmp_rununik_27152/align.dat  \
                 draft-assembly.fasta genome-hic.fasta \
          
           This saves alignment time and can get the results quickly. 

