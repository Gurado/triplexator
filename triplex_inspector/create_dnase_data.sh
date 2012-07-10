#!/bin/sh -e

USAGEMSG="usage: $(basename $0) -c chromosomes -o outputdirectory -G -P -m mode -k [ -f fseq -p picard -b bedtools -s samtools ] dnasedirectory

Prepare DNase hyper- or hypo-sensitivity data for the TFO design process from already aligned tagfiles (bam) 
This involves:
	1) duplicate removal
	2) tag shink (shink to 5' cutsite)
	3) replica merge

Requirements (e.g. in PATH environment):
	Java
	BEDtools >= v2.12.0
	Picard >= 1.48
	Samtools >= 0.1.12a

* dnase dnasedirector- The full path to directory containing folder of dnase bam files of the following structure:
	-dnase
		|- tissue/condition 1
		|	|- Replica 1.bam
		|	|- Replica 2.bam
		|- tissue/condition 2 
		|	|- ...
		...

* -c chromosomes - The full path to a file holding the chromosome sizes in tab-separated format, e.g. 'chr1    249250621'
* -p picard - Path to the folder containing the picard toolset (with final /; ommit when in path)
* -b bedtools - Path to the folder containing the bedtools binaries (with final /; ommit when in path)
* -s samtools - Path to samtools binary (with final /; ommit when in path)
* -f fseq - Path to fseq binary (with final /; ommit when in path)
* -G bedGraph - create bedGraph (using fseq)
* -P narrowPeak - create narrow peak format (using fseq)
* -k keep intermediate files - if enabled does not remove intermediate files (default off)
* -o outputdirectory - where you expect to find the files (default is the input folder)
* -m mode - either HYPER or HYPO for DNase hypersensitivty or hyposensitivity, hyper-mode will collapse the tags to the 5' nucleotide while hypo will keep the original tag size

Resulting collapsed, unique-mapping tags (shrinked to the 5' nucleotide in case of hypersensitivity) are placed in the dnase folder with
a bam file and the other requested formats per folder, e.g.:
	-dnase
	    |- tissue/condition_1.bam
	    |- tissue/condition_1.bedGraph
		|- tissue/condition 1
		|	|- Replica 1.bam
		|	|- Replica 2.bam
		|	|- ...
	    |- tissue/condition_2.bam
	    |- tissue/condition_2.bedGraph
		|- tissue/condition 2 
		|	|- Replica 1.bam
		|	|- ...
		|..


e.g. ./scripts/create_dnase_data.sh -k -c /Users/Shared/BEDTools/genomes/mouse.mm9.genome -p /Users/Shared/picard-tools-1.48/picard-tools-1.48/ -b /Users/Shared/BEDTools/bin/ -o /Users/<username>/TFO_design/data/dnase_mm9_8wk/ /Users/<username>/TFO_design/data/dnase_mm9_8wk/

Copyright Fabian Buske 2012
Contact: f.buske@uq.edu.au
"

[ $# -lt 4 ] && echo "$USAGEMSG" >&2 && exit 1

CHROMOSOMES="NONE"
DNASE="NONE"
PICARD=""
BEDTOOLS=""
SAMTOOLS=""
FSEQ=""
KEEPFILES="FALSE"
OUTPUT="NONE"
MODE="HYPER"
BEDGRAPH="FALSE"
NARROWPEAK="FALSE"

while getopts "c:p:b:s:f:o:m:kGP" opt
do
    case $opt in
	(c) CHROMOSOMES=$OPTARG;;
	(p) PICARD=$OPTARG;;
	(b) BEDTOOLS=$OPTARG;;
	(s) SAMTOOLS=$OPTARG;;
	(f) FSEQ=$OPTARG;;
	(o) OUTPUT=$OPTARG;;
	(m) MODE=$OPTARG;;
	(k) KEEPFILES="TRUE";;
	(G) BEDGRAPH="TRUE";;
	(P) NARROWPEAK="TRUE";;
	(*) echo "$0: error - unrecognized option $1" >&2; exit 1;;
	esac
done
shift $(($OPTIND - 1))
DNASE=$1

DIR=`dirname $0`

[ "$DNASE" = "NONE" ] && echo "[ERROR] DNASE folder missing" >&2 && exit 1
[ "$MODE" != "HYPER" ] && [ "$MODE" != "HYPO" ] && echo "[ERROR] MODE not known" >&2 && exit 1
[ "$CHROMOSOMES" = "NONE" ] && echo "[ERROR] CHROMOSOMES parameter missing" >&2 && exit 1
[ "$OUTPUT" == "NONE" ] && echo "[WARN] OUTPUTDIRECTORY not specified, using input folder instead" >&2 && OUTPUT=${DNASE}/

[ "$BEDTOOLS" != "" ] && echo "BEDTOOLS: ${BEDTOOLS}"
[ "$PICARD" != "" ] && echo "PICARD: ${PICARD}"
[ "$SAMTOOLS" != "" ] && echo "SAMTOOLS: ${SAMTOOLS}"
[ "$FSEQ" != "" ] && echo "FSEQ: ${FSEQ}"
[ "$KEEPFILES" == "FALSE" ] && echo "Keep intermediate files: No" || echo "Keep intermediate files: Yes"

echo "DNAse: ${DNASE}"
echo "CHROMOSOMES: ${CHROMOSOMES}"
echo "Output to: ${OUTPUT}"
echo "MODE: $MODE"
echo "create BEDgraph : ${BEDGRAPH}"
echo "create narrowPeak : ${NARROWPEAK}"

[ ! -d ${OUTPUT} ] && echo "create output directory ${OUTPUT}" && mkdir -p ${OUTPUT}

echo "processing DNase data (BAM files)"
for d in `ls -d ${DNASE}/*/` ; do
	dn=`basename $d`
	echo $dn
	if [ ! -f ${OUTPUT}/${dn}.bam ]; then 
	
		for BAM in `ls ${d}/*.bam`; do
			BAMNAME="${BAM%.*}" 
			echo  ${BAMNAME}
			echo "sort w/r to coordinates"
			[ ! -f ${BAMNAME}_sorted ] && java -Dsnappy.disable=true -Xmx2g -jar ${PICARD}SortSam.jar CREATE_INDEX=true SORT_ORDER=coordinate INPUT=${BAMNAME}.bam OUTPUT=${BAMNAME}_sorted
			echo "remove duplicates"
			[ ! -f ${BAMNAME}_noduplicate ] && java -Dsnappy.disable=true -Xmx2g -jar ${PICARD}MarkDuplicates.jar REMOVE_DUPLICATES=true INPUT=${BAMNAME}_sorted OUTPUT=${BAMNAME}_noduplicate METRICS_FILE=${BAMNAME}_noduplicate.met
			
			if [ "$MODE" == "HYPER" ]; then
				echo "shrink tag to 5' nucleotide"
				[ ! -f ${BAMNAME}_processed.bed ] && cat ${BAMNAME}_noduplicate | ${BEDTOOLS}bamToBed | awk -F\\t  '{OFS = "\t"; $4=(FNR); if ($6 =="+") $3=$2+1; else {$2=$3-1;}; print $0}' > ${BAMNAME}_processed.bed
			elif [ "$MODE" == "HYPO" ]; then
				echo "process bam to bed"
				[ ! -f ${BAMNAME}_processed.bed ] && cat ${BAMNAME}_noduplicate | ${BEDTOOLS}bamToBed | awk -F\\t  '{OFS = "\t"; $4=(FNR); print $0}' > ${BAMNAME}_processed.bed			
			fi
		done
		echo "merge all replicas"
		[ ! -f ${OUTPUT}/${dn}.bam ] && (cat ${d}/*_processed.bed | ${BEDTOOLS}bedToBam -g ${CHROMOSOMES} > ${OUTPUT}/${dn}_processed.bam.tmp && ${SAMTOOLS}samtools sort ${OUTPUT}/${dn}_processed.bam.tmp ${OUTPUT}/${dn} && rm ${OUTPUT}/${dn}_processed.bam.tmp ) 

	fi

	if [ $BEDGRAPH = "TRUE" ] || [ $NARROWPEAK = "TRUE" ]; then 
		[ -d ${OUTPUT}/_tmp_bed ] && rm -r -f ${OUTPUT}/_tmp_bed
		mkdir -p ${OUTPUT}/_tmp_bed
		${BEDTOOLS}bamToBed -i ${OUTPUT}/${dn}.bam > ${OUTPUT}/_tmp_bed/${dn}.bed
	fi

	if [ $BEDGRAPH = "TRUE" ] && [ ! -f  ${OUTPUT}/${dn}.bedGraph ]; then	
		echo "create bedgraph"
		mkdir -p ${OUTPUT}/_tmp_wig
		${FSEQ}fseq -d ${OUTPUT}/_tmp_bed -l 600 -f 0 -v -of wig -o ${OUTPUT}/_tmp_wig
		cat ${OUTPUT}/_tmp_wig/*.wig | python ${DIR}/scripts/wig2bedGraph.py > ${OUTPUT}/${dn}.bedGraph
		rm -r -f $OUTPUT/_tmp_wig
	fi

	if [ $NARROWPEAK = "TRUE" ] && [ ! -f ${OUTPUT}/${dn}.npf ] ; then
		echo "create narrow peak file"
		[ -d ${OUTPUT}/_tmp_npf ] && rm -r -f ${OUTPUT}/_tmp_npf
		mkdir -p ${OUTPUT}/_tmp_npf
		${FSEQ}fseq -d ${OUTPUT}/_tmp_bed -l 600 -f 0 -v -of npf -o ${OUTPUT}/_tmp_npf
		cat ${OUTPUT}/_tmp_npf/*.npf > ${OUTPUT}/${dn}.npf
		rm -r -f ${OUTPUT}/_tmp_npf
	fi

	if [ $BEDGRAPH = "TRUE" ] || [ $NARROWPEAK = "TRUE" ]; then
        	rm -r -f ${OUTPUT}/_tmp_bed
	fi
	
	if [ $KEEPFILES = "FALSE" ]; then
		echo "clean up"
		for BAM in `ls ${d}/*.bam`; do
			[ -f ${BAMNAME}_processed.bed ] && rm ${BAMNAME}_sorted* ${BAMNAME}_noduplicate* ${BAMNAME}_processed.bed 
		done
	fi
done

echo "finished creating DNase data"
