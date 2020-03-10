#!/bin/bash
function cleanFile
{
	FILE=$1
	shift
	cat $FILE | awk '{print $2 "\t" $3 }' | egrep -v "^0	"
}

#cleanFile human-hic-len.freq > bE378K21-screen.dat1.cleaned
#cleanFile hummingbird-hic.freq > bE378K21-screen.dat2.cleaned
#cleanFile fMasArm1-hic.freq > bE378K21-screen.dat3.cleaned
#cleanFile tdevil-hic.freq > bE378K21-screen.dat4.cleaned
#cleanFile human-mp15-hic.freq > bE378K21-screen.dat4.cleaned 

function plotcmd
{
        printf "set logscale x\n"
        printf "set logscale y\n"
	printf "set terminal svg\n"
        printf "set style line 1 lt 1 lw 3 pt 3 linecolor rgb \"black\"\n"
        printf "set style line 2 lt 1 lw 3 pt 3 linecolor rgb \"green\"\n"
        printf "set style line 3 lt 1 lw 3 pt 3 linecolor rgb \"blue\"\n"
        printf "set style line 4 lt 1 lw 3 pt 3 linecolor rgb \"red\"\n"
        printf "set style line 5 lt 1 lw 3 pt 3 linecolor rgb \"cyan\"\n"
	printf "set xlabel \"Paired Ends Mapping Length\"\n"
	printf "set ylabel \"Frequency / 100\"\n"
        printf "plot [ 1 to 100000000 ] [ 0.0001 to 50.0 ] \"human-hic-len.freq\" title \"Human-HiC\" with lines ls 1,\"hummingbird-hic.freq\" title \"Hummingbird-HiC\" with lines ls 2,\"fMasArm1-hic.freq\" title \"Fish MasArm-HiC\" with lines ls 3,\"tdevil-hic.freq\" title \"Tasmanian Devil-HiC\" with lines ls 4,\"sample-hic.freq\" title \"Test sample-HiC\" with lines ls 5"
}

plotcmd | gnuplot > hicplot.svg
inkscape -z --export-text-to-path --export-pdf hicplot.pdf hicplot.svg
gs -r600 -dNOPAUSE -dBATCH -sDEVICE=png256 -sOutputFile=hicplot.png hicplot.pdf

#rm -f bE378K21-raw.dat.cleaned bE378K21-screen.dat.cleaned data.svg
