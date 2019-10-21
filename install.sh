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


##### Download and install pigz ######

echo "Downloading and installing pigz"
if [[ ! -s $bindir/pigz ]]; then

    if [[ ! -d $projdir/src/pigz ]]; then
	cd $projdir/src/
        wget -r -np -nd https://zlib.net/pigz/pigz-2.4.tar.gz &> $projdir/src/log/pigz_wget.log
        tar -xvzf pigz-2.4.tar.gz &> $projdir/src/log/pigz_untar.log
        rm -f pigz-2.4.tar.gz
    fi

    if [[ ! -s $projdir/src/pigz/pigz ]]; then
	cd $projdir/src/pigz-2.4
	make &> $projdir/src/log/pigz_installation.log
    fi

    cp pigz $bindir
fi

if  [[ ! -s $bindir/pigz ]]; then
    echo " !! Error: pigz not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if bwa was downloaded properly:" $projdir/src/log/pigz_cloning.log 
    echo "   Check if the bwa was compiled properly:" $projdir/src/log/pigz_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/pigz/pigz $bindir/pigz 
    
    errs=$(($errs+1))
else
    echo " pigz succesfully installed!"
    rm -rf $projdir/src/pigz/
fi

### Download and install PretextMap  ######

echo "Downloading and installing PretextMap/PretextView"
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
    echo "   Check if PretextMap was downloaded properly:" $projdir/src/log/pigz_cloning.log 

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
    echo " pigz succesfully installed!"
    rm -rf $projdir/src/PretextView/
fi

##### Download and install pigz ######

echo "Downloading and installing pigz"
if [[ ! -s $bindir/pigz ]]; then

    if [[ ! -d $projdir/src/pigz ]]; then
	cd $projdir/src/
        wget -r -np -nd https://zlib.net/pigz/pigz-2.4.tar.gz &> $projdir/src/log/pigz_wget.log
        tar -xvzf pigz-2.4.tar.gz &> $projdir/src/log/pigz_untar.log
        rm -f pigz-2.4.tar.gz
    fi

    if [[ ! -s $projdir/src/pigz/pigz ]]; then
	cd $projdir/src/pigz-2.4
	make &> $projdir/src/log/pigz_installation.log
    fi

    cp pigz $bindir
fi

if  [[ ! -s $bindir/pigz ]]; then
    echo " !! Error: pigz not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if bwa was downloaded properly:" $projdir/src/log/pigz_cloning.log 
    echo "   Check if the bwa was compiled properly:" $projdir/src/log/pigz_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/pigz/pigz $bindir/pigz 
    
    errs=$(($errs+1))
else
    echo " pigz succesfully installed!"
    rm -rf $projdir/src/pigz/
fi

###### Compile scaffHiC sources ######

echo; echo "Compiling scaffHiC sources"

srcs=( scaffHiC scaffHiC_agp2AGP scaffHiC_AGPbuild scaffHiC_cover scaffHiC_ctg2AGP scaffHiC_fastq scaffHiC_insert scaffHiC_lengthdis scaffHiC_motif1 scaffHiC_motif2 scaffHiC_motifm scaffHiC_motreads scaffHiC_offset scaffHiC_orient scaffHiC_output scaffHiC_pairs scaffHiC_pairs2 scaffHiC_pairs-chr scaffHiC_pairs-cov scaffHiC_pairs-map scaffHiC_PCRdup2 scaffHiC_PCRdup3 scaffHiC_PCRdup scaffHiC_pmatrix scaffHiC_proPair scaffHiC_RDplace scaffHiC_reads scaffHiC_rename scaffHiC_scaffsort scaffHiC_scf2agp scaffHiC_scf2scf scaffHiC_screen scaffHiC_screen-cover scaffHiC_translo scaffHiC_translo-chr)

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




