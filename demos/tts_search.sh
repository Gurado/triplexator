#!/bin/sh

DEMOS=`dirname $0`

[ ! -d ${DEMOS}/examples ] && echo "creating folder ${DEMOS}/examples" && mkdir ${DEMOS}/examples
echo "Searching for putative triplex target sites in ${DEMOS}/P00374.fasta"
echo "./bin/triplexator -l 15 -g 50 -e 10 -fr on -mrl 7 -mrp 3 -of 0 -od ${DEMOS}/examples -o P00374.tts -ds ${DEMOS}/P00374.fasta"
echo "Command:"
./bin/triplexator -l 15 -g 50 -e 10 -fr on -mrl 7 -mrp 3 -of 0 -od ${DEMOS}/examples -o P00374.tts -ds ${DEMOS}/P00374.fasta 
echo "...finished"
echo "Results are in ${DEMOS}/examples/P00374.tts"