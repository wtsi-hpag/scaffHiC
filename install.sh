#!/bin/bash


projdir=`pwd`

bindir=$projdir/src/scaff-bin/
mkdir -p $bindir
mkdir -p $projdir/src/log/

errs=0

##### Download and install BWA ######

echo "Downloading and installing BWA"
if [[ ! -s $bindir/bwa ]]; then

    if [[ ! -d $projdir/src/bwa ]]; then
	cd $projdir/src/
	git clone https://github.com/lh3/bwa.git &> $projdir/src/log/bwa_cloning.log
    fi

    if [[ ! -s $projdir/src/bwa/bwa ]]; then
	cd $projdir/src/bwa
	make &> $projdir/src/log/bwa_installation.log
    fi

    cp bwa $bindir
fi

if  [[ ! -s $bindir/bwa ]]; then
    echo " !! Error: bwa not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if bwa was downloaded properly:" $projdir/src/log/bwa_cloning.log 
    echo "   Check if the bwa was compiled properly:" $projdir/src/log/bwa_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/bwa/bwa $bindir/bwa 
    
    errs=$(($errs+1))
else
    echo " BWA succesfully installed!"
    rm -rf $projdir/src/bwa/
fi


### Download and install PretextMap  ######

echo "Downloading and installing PretextMap"
if [[ ! -s $bindir/PretextMap ]]; then

    if [[ ! -d $projdir/src/PretextMap ]]; then
	cd $projdir/src/
        wget -r -np -nd https://github.com/wtsi-hpag/PretextMap/releases/download/0.0.41/PretextMap_Linux-x86-64.zip &> $projdir/src/log/Pretext_wetmap.log
        unzip PretextMap_Linux-x86-64.zip &> $projdir/src/log/Pretext_unzip.log
        rm -f PretextMap_Linux-x86-64.zip 
    fi

    chmod 755 PretextMap
    cp PretextMap $bindir
fi

if  [[ ! -s $bindir/PretextMap ]]; then
    echo " !! Error: PretextMap not installed properly!"; 
    echo "   Please check the log files:" 

    # Cleaning up
    cd $projdir/src
    rm -rf $bindir/PretextMap 
    
    errs=$(($errs+1))
else
    echo " PretextMap succesfully installed!"
    rm -rf $projdir/src/PretextMap/
fi

##### Download and install PretextView ######

echo "Downloading and installing PretextView"
if [[ ! -s $bindir/PretextView ]]; then

    if [[ ! -d $projdir/src/PretextView ]]; then
	cd $projdir/src/
        wget -r -np -nd https://github.com/wtsi-hpag/PretextView/releases/download/0.0.41/PretextView-Linux-x86-64.zip &> $projdir/src/log/PretextView_wget.log
        unzip PretextView-Linux-x86-64.zip &> $projdir/src/log/PretextView_unzip.log
        rm -f PretextView-Linux-x86-64.zip 
    fi

    chmod 755 PretextView
    cp PretextView $bindir
fi

if  [[ ! -s $bindir/PretextView ]]; then
    echo " !! Error: PretextView not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if PretextView was downloaded properly:" $projdir/src/log/PretextView_wget.log 

    # Cleaning up
    cd $projdir/src
    rm -rf $bindir/PretextView 
    
    errs=$(($errs+1))
else
    rm -rf $projdir/src/PretextView/
fi

###### Compile scaffhic sources ######

echo; echo "Compiling scaffHiC sources"

srcs=( scaffHiC_screen scaffHiC_output scaffHiC_fastq scaffHiC_rename scaffHiC_outbreak scaffHiC_superAGP scaffhic breakhic scaffHiC_PCRdup scaffHiC_PCRdup2 scaffHiC_PCRdup3 scaffHiC_RDplace scaffHiC_proPair scaffHiC_AGPbuild scaffHiC_pmatrix scaffHiC_insert scaffHiC_lengthdis scaffHiC_orient scaffHiC_scf2agp scaffHiC_break scaffHiC_breakClean scaffHiC_BKplace scaffHiC_BKplace2 scaffHiC_CIplace scaffHiC_scf2scf scaffHiC_getPair scaffHiC_breakPair scaffHiC_breakPair2 scaffHiC_breakAGP scaffHiC_breakTags scaffHiC_reads scaffHiC_screen-cover scaffHiC_translo scaffHiC_translo-chr)

cd $projdir/src
make &> $projdir/src/log/sources_compilation.log

echo; echo "Checking installation:"
for src in "${srcs[@]}"; do
    if [[ ! -s $bindir/$src ]]; then 
        echo " !! Error: executable $src missing in $bindir"
	echo "    Please check for errors the log file:" $projdir/src/log/sources_*	
        errs=$(($errs+1))
    fi
done

if [  $errs -gt 0 ]; then echo; echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi




