#!/bin/sh

DEMOS=`dirname $0`

[ ! -d ${DEMOS}/examples ] && echo "creating folder ${DEMOS}/examples" && mkdir ${DEMOS}/examples
echo "Searching for putative triplex-forming oligonucleotides in ${DEMOS}/P00374.fasta"
echo "./bin/triplexator -l 15 -e 15 -m R,M -fr on -mrl 7 -mrp 1 -of 0 -po -od demos/examples -o P00374.tfo -ss demos/P00374.fasta "
echo "Command:"
./bin/triplexator -l 15 -e 15 -m R,M -fr on -mrl 7 -mrp 1 -of 0 -po -od demos/examples -o P00374.tfo -ss demos/P00374.fasta 
echo "...finished"
echo "Results are in ${DEMOS}/examples/P00374.tfo"