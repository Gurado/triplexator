#!/bin/sh

DEMOS=`dirname $0`

[ ! -d ${DEMOS}/examples ] && echo "creating folder ${DEMOS}/examples" && mkdir ${DEMOS}/examples
echo "Searching for matching TFO-TTS pairs using ${DEMOS}/P00374.fasta and ${DEMOS}/P00374.fasta as input"
echo "ATTENTION: parameter setting spans large searchspace - use with care when running on a large dataset!"
echo "Command:"
echo "./bin/triplexator -l 15 -e 20 -c 2 -fr off -g 20 -m R -fm 0 -of 1 -od demos/examples -o P00374.tpx -ss demos/P00374.fasta  -ds demos/P00374.fasta "
./bin/triplexator -l 15 -e 20 -c 2 -fr off -g 20 -m R -fm 0 -of 1 -od demos/examples -o P00374.tpx -ss demos/P00374.fasta  -ds demos/P00374.fasta 
echo "...finished"
echo "Results are in ${DEMOS}/examples/P00374.tpx"
