#!/bin/sh -e 

DIR=`dirname $0`

# -----------------------------------------------------
# Location of programs used by triplex inspector
# -----------------------------------------------------
# triplex-inspector
INSPECTOR=${DIR}/../
# triplexator executable
TRIPLEXATOR=${DIR}/../../bin/triplexator
# cicos executable
CIRCOS=/opt/circos-0.62-1/bin/circos
# BEDtools folder containing the executables
BEDTOOLS=/opt/bedtools-2.17.0/bin/

# -----------------------------------------------------
# Location of the data
# -----------------------------------------------------
DATA=${DIR}/data
# annotation file
ANNOTATION=${DATA}/M_tuberculosis_H37Rv_2_exons.gtf
# genome in fasta format
GENOME=${DATA}/M_tuberculosis_H37Rv.fasta
# chromosome sizes
GENOMESIZE=${DATA}/M_tuberculosis_H37Rv.genome
# region of interest
TARGETREGION=${DATA}/nucleotide_biosynthesis.bed

# -----------------------------------------------------
# Location where the results are expected
# -----------------------------------------------------
RESULTS=${DIR}/results

# -----------------------------------------------------
# Parameters
# -----------------------------------------------------
# setting for primary targets
TTSOPTIONS="--lower-length-bound 12 --consecutive-errors 1 --error-rate 10 --runtime-mode 2 --output-format 0 --filter-repeats off --filtering-mode 0"
# setting for off-targets (relaxed constraints)
TPXOPTIONS="--lower-length-bound 10 --consecutive-errors 1 --error-rate 10 --runtime-mode 1 --output-format 0 --filter-repeats off --filtering-mode 0"

# -----------------------------------------------------
# test-related preparations
# -----------------------------------------------------
# extract test data
[ ! -d ${DIR}/data ] && tar xf ${DIR}/data.tar.gz

# -----------------------------------------------------
# start triplex inspector
# -----------------------------------------------------
${INSPECTOR}/run_inspector.sh \
	-g ${GENOME} \
	-s ${GENOMESIZE} \
	-1 ${TARGETREGION} \
	-a ${ANNOTATION} \
	-B ${BEDTOOLS} \
	-T ${TRIPLEXATOR} \
	-l "${TTSOPTIONS}" \
	-o "${TPXOPTIONS}" \
	-v \
	${RESULTS}
	
	
# -----------------------------------------------------
# finished
# -----------------------------------------------------
echo "Hooray! Everything went through smoothly!"
