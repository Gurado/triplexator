#!/bin/sh -e

USAGEMSG="usage: $(basename $0) -l \"localConstraints\" -o \"offTargetConstraints\" -c chromatinData -f chromatinFormat -g genome -s genomeSize -1 lociOfInterest.bed -2 alternativeLoci.bed -T TRIPLEXATOR -B BEDtools -C Circos -v <outputdirectory>

Starts the triplex target analysis workflow on the given dataset. 

Requirements (in PATH environment or specified):
	Python 2.7 + Biopython
	Triplexator v1.2+
	BEDtools v2.12.0
	Samtools and pysam (bam) or bx-python (bigwig) need to be installed when chromatin data (e.g. DNase I hypersensitivity) are taken into account
	Circos v0.62+


* outputdirectory - where you expect to find the result files 
* -l localConstraints - a quoted string containing any options to pass to Triplexator for assessing the loci of interest for putative target sites.
* -o offTargetConstraints - a quoted string containing any options to pass to Triplexator for assessing off-targets (e.g. relaxing size and mismatch constraint). If not provided the same parameters as before will be used.
* -1 lociOfInterest.bed - full path to a 6-field bed file specifying the loci of interest.
* -2 alternativeLoci.bed - full path to a bed file specifying all loci to be considered in the off-target analysis. If not specified the whole genome is used.
* -c chromatinData - The full path to the file containing the chromatin organization data (e.g. DNase I cutsites).
* -f chromatinFormat - the file format of the chromatin data (bam or bigwig) .
* -g genome - The full path to a fasta file containing the genomic sequence data.
* -s genomeSize - The full path to a file the chromosome sizes in tab-separated format, e.g. 'chr1    249250621'.
* -P python - path to python executable.
* -T triplexator - path to triplexator executable.
* -B BEDtools - path to BEDtools executables with terminal directory delimiter (i.e. "/"), if no path is given BEDtool binaries are expected to be found in the PATH variable.
* -C Circos - path to Circos executable, if no path is given circos binaries are expected to be found in the PATH variable.
* -a annotation - The full path to a gtf file specifying (gene) annotation data
* -x - skips commands for which there already exists an result  (e.g. from a previous, disrupted run). This is to reduce runtime, but can result in undefined behavior.
* -v - print progress information (verbose).

"

DIR=`dirname $0`
VERSION="0.1.0"

[ $# -lt 4 ] && echo "$USAGEMSG" >&2 && exit 1

TTSOPTIONS="--lower-length-bound 17 --consecutive-errors 1 --error-rate 8 --runtime-mode 2 --output-format 0 --min-guanine 30 --filter-repeats off --filtering-mode 0"
TPXOPTIONS="NONE"
OUTPUT="NONE"
GENOME="NONE"
GENOMESIZE="NONE"
LOC="NONE"
LOISEQ="NONE"
LOCSEQ="NONE"
CHROMATIN="NONE"
CHROMATINFORMAT="NONE"
TRIPLEXATOR="triplexator"
PYTHON="python"
BEDTOOLS=""
CIRCOS="circos"
ANNOTATION="NONE"
MAXIMALSITES=50 		# threshold for the number of primary target regions to be tested for off-targets 
SKIPIFEXISTS="FALSE"
VERBOSE="--quiet"
JSON="primary_targets.json"
FLANKS=0 # length of flanking sequence shown

while getopts "1:2:g:s:c:f:l:o:P:T:a:B:C:xv" opt;
do
	case ${opt} in
	1) LOI="$OPTARG";;
	2) LOC="$OPTARG";;
	g) GENOME="$OPTARG";;
	s) GENOMESIZE="$OPTARG";;
	c) CHROMATIN="$OPTARG";;
    f) CHROMATINFORMAT="$OPTARG";;
	l) TTSOPTIONS="$OPTARG";;
	o) TPXOPTIONS="$OPTARG";;
	P) PYTHON="$OPTARG";;
	T) TRIPLEXATOR="$OPTARG";;
	B) BEDTOOLS="$OPTARG";;
	C) CIRCOS="$OPTARG";;
	a) ANNOTATION="$OPTARG";;
	x) SKIPIFEXISTS="TRUE";;
	v) VERBOSE="--verbose";;
	\?) print >&2 "$0: error - unrecognized option $1" 
		exit 1;;
	esac
done

shift $(($OPTIND-1))
OUTPUT=$1

if [ ! -d ${OUTPUT} ]; then  
	mkdir -p ${OUTPUT} 
else
	echo "[WARN] output directory already exists, content will be overwritten" 
fi

#---------------------------
# init logfile
#---------------------------
LOGFILE=${OUTPUT}/log.txt
cat /dev/null > ${LOGFILE}

echo "Output directory: $@" >> ${LOGFILE}

LOISHORT=${LOI##*/}
LOISHORT=${LOISHORT%.*}
LOCSHORT=${LOC##*/}
LOISHORT=${LOISHORT%.*}
GENOMESHORT=${GENOME##*/}
GENOMESHORT=${GENOMESHORT%.*}
CHROMATINSHORT=${CHROMATIN##*/}
CHROMATINSHORT=${CHROMATINSHORT%.*}

#---------------------------
# check pre-conditions
#---------------------------
echo "-------     " >> ${LOGFILE}
if [ "${TRIPLEXATOR}" != "triplexator" ]; then
	[ ! -x ${TRIPLEXATOR} ]  && echo "[ERROR] Triplexator executables not found (parameter -T)" >> ${LOGFILE} && exit 1
else
	[ "which ${TRIPLEXATOR}" = "" ] && echo "[ERROR] Triplexator executables not found (parameter -T)" >> ${LOGFILE} && exit 1
fi
if [ "${PYTHON}" != "python" ]; then
	[ ! -x ${PYTHON} ]  && echo "[ERROR] Python executables not found (parameter -P)" >> ${LOGFILE} && exit 1
else
	[ "which ${PYTHON}" = "" ] && echo "[ERROR] Python executables not found (parameter -P)" >> ${LOGFILE} && exit 1
fi
if [ "${BEDTOOLS}" != "" ]; then
	[ ! -x ${BEDTOOLS}fastaFromBed ]  && echo "[ERROR] BEDtools executables not found (parameter -B)" >> ${LOGFILE} && exit 1
else
	[ "which fastaFromBed" = "" ] && echo "[ERROR] BEDtools executable not in path (parameter -B)" >> ${LOGFILE} && exit 1
fi 
if [ "${CIRCOS}" != "circos" ]; then
	[ ! -x ${CIRCOS} ]  && echo "[NOTICE] Circos not found in path (parameter -C), no circos plots will be generated" >> ${LOGFILE}
else
	[ "which ${CIRCOS}" = "" ] && echo "[ERROR] Circos executable not found (parameter -C)" >> ${LOGFILE} && exit 1
fi
[ ! -f ${GENOME} ] && echo "[ERROR] genome file does not exist (parameter -g)" >> ${LOGFILE} && exit 1
[ ! -f ${LOI} ]  && echo "[ERROR] input file 1 (loci of interest) does not exist (parameter -1)" >> ${LOGFILE}  && exit 1
[ "${ANNOTATION}" != "NONE" ] && [ ! -f ${ANNOTATION} ] && echo "[ERROR] specified annotation file does not exist (parameter -a)" >> ${LOGFILE} && exit 1
[ "${LOC}" = "NONE" ] && echo "[NOTICE] input file 2 (alternative loci to consider) does not exist (parameter -2) using complete genome instead" >> ${LOGFILE} 
[ "${TPXOPTIONS}" = "NONE" ] && TPXOPTIONS=${TTSOPTIONS} && echo "[NOTICE] Using localConstraints as collateralConstraints as well: TPXOPTIONS=TTSOPTIONS" >> ${LOGFILE}

# add specific options the off-target scan crucial for this pipeline (screening with a pyrimidine TFO with errors corresponding to the TFO position)
TPXOPTIONS="${TPXOPTIONS} --triplex-motifs Y --error-reference 2"

[ "${CHROMATIN}" != "NONE" ] && [ "${CHROMATINFORMAT}" == "NONE" ] &&  echo "[ERROR] format of chromatin data not provided (parameter -f)" >> ${LOGFILE} && exit 1

[ "${CHROMATINFORMAT}" != "NONE" ] && [ "${CHROMATINFORMAT}" != "bam" ] && [ "${CHROMATINFORMAT}" != "bigwig" ] && echo "[ERROR] format of chromatin data not known (parameter -f)" >> ${LOGFILE} && exit 1


#---------------------------
# output parameter and version to logfile
#---------------------------

echo "------------------" >> ${LOGFILE}
echo "INSPECTOR VERSION:"${VERSION} >> ${LOGFILE}
echo "EXE DIR:          "${DIR} >> ${LOGFILE}
echo "LOISHORT:         "${LOISHORT} >> ${LOGFILE}
echo "LOCSHORT:         "${LOCSHORT} >> ${LOGFILE}
echo "GENOMESHORT:      "${GENOMESHORT} >> ${LOGFILE}
echo "GENOME:           "${GENOME} >> ${LOGFILE}
echo "GENOMESIZE:       "${GENOMESIZE} >> ${LOGFILE}
echo "TRIPLEXATOR:      "${TRIPLEXATOR} >> ${LOGFILE}
echo "PYTHON:           "${PYTHON} >> ${LOGFILE}
echo "TTSOPTIONS:       "${TTSOPTIONS} >> ${LOGFILE}
echo "TPXOPTIONS:       "${TPXOPTIONS} >> ${LOGFILE}
echo "CHROMATIN:        "${CHROMATIN} >> ${LOGFILE}
echo "CHROMATINFORMAT:  "${CHROMATINFORMAT} >> ${LOGFILE}
echo "OUTPUT:           "${OUTPUT} >> ${LOGFILE}
echo "ANNOTATION:       "${ANNOTATION} >> ${LOGFILE}
echo "BEDTOOLS:         "${BEDTOOLS} >> ${LOGFILE}
echo "CIRCOS:           "${CIRCOS} >> ${LOGFILE}
echo "SKIPIFEXISTS:     "${SKIPIFEXISTS} >> ${LOGFILE}
echo "VERBOSE:          "${VERBOSE} >> ${LOGFILE}
echo "------------------" >> ${LOGFILE}


#---------------------------
# start the workflow
#---------------------------
echo "-------     " >> ${LOGFILE}
echo "*** "`date +%H:%M:%S`" inspector started" >> ${LOGFILE} 

#---------------------------
# prepare sequence data
#---------------------------
echo "-------     " >> ${LOGFILE}
echo "*** "`date +%H:%M:%S`" get sequence data" >> ${LOGFILE}
echo "1) loci of interest" >> ${LOGFILE}

if	[ ! -f ${OUTPUT}/${LOISHORT}.fasta ] || [ ! ${SKIPIFEXISTS} = "TRUE" ]; then	
	echo "...getting sequence data for loci of interest" >> ${LOGFILE} 
	${BEDTOOLS}fastaFromBed -name -s -fi ${GENOME} -bed ${LOI} -fo ${OUTPUT}/${LOISHORT}.fasta
fi
LOISEQ=${OUTPUT}/${LOISHORT}.fasta
echo "LOISEQ:     "${LOISEQ} >> ${LOGFILE}

echo "2) regions to consider for off-targets" >> ${LOGFILE}
if  [ "${LOC}" = "NONE" ]; then
	echo "...consider whole genome for alternative regions" >> ${LOGFILE}
	LOCSEQ=${GENOME}
	cat ${GENOMESIZE} | awk '{print $1"\t"0"\t"$2"\t"$1"\t0\t+"}' > ${OUTPUT}/${GENOMESHORT}.bed
	LOC=${OUTPUT}/${GENOMESHORT}.bed
else
	if	[ -f ${LOC} ]  && ([ ! -f ${OUTPUT}/${LOCSHORT}.fasta ] || [ ! ${SKIPIFEXISTS} = "TRUE" ]); then
		LOCSEQ=${OUTPUT}/${LOCSHORT}.fasta 
		echo "...getting sequence data for alternative regions" >> ${LOGFILE} 
		${BEDTOOLS}fastaFromBed -name -fi ${GENOME} -bed ${LOC} -fo ${OUTPUT}/${LOCSHORT}.fasta 
		LOCSEQ=${OUTPUT}/${LOCSHORT}.fasta
	fi
fi
echo "LOCSEQ:     "${LOCSEQ} >> ${LOGFILE}
echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}

#---------------------------
# prepare annotation data (split into individual files grouped by column 2 entries)
#---------------------------
if [ "$ANNOTATION" != "NONE" ]; then
	echo "-------     " >> ${LOGFILE}
	echo "*** "`date +%H:%M:%S`" prepare annotation data" >> ${LOGFILE} 
	echo "...convert annotation into BED files" >> ${LOGFILE} 
	[ ! -d ${OUTPUT}/annotation ] && mkdir ${OUTPUT}/annotation
	for ANNO in `cat ${ANNOTATION} | awk -F\\t '{print $2}' | sort -u`; do
		if [ ! -f ${OUTPUT}/annotation/${ANNO}.bed ] || [ ! ${SKIPIFEXISTS} = "TRUE" ]; then
			echo "... processing annotation" ${ANNO} >> ${LOGFILE} 
			cat ${ANNOTATION} | awk -F\\t -v d1="${ANNO}" '$2==d1 {print $1"\t"$4"\t"$5"\t"$3"\t"$6"\t"$7}' > ${OUTPUT}/annotation/${ANNO}.bed
		fi
	done
	echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}
fi

#---------------------------
# screen region of interest for putative primary targets
#---------------------------
[ ! -d ${OUTPUT}/tts ] && mkdir ${OUTPUT}/tts
if [ ! -f ${OUTPUT}/tts${LOISHORT}.TTS ] || [ ! ${SKIPIFEXISTS} = "TRUE" ]; then
	echo "-------     " >> ${LOGFILE}
	echo "*** "`date +%H:%M:%S`" search primary targets (TTSs)" >> ${LOGFILE}
	${TRIPLEXATOR} ${TTSOPTIONS} --merge-features -od ${OUTPUT}/tts -o ${LOISHORT}.TTS -ds ${LOISEQ}
	echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}
else
	echo "...skipping primary target region detection due to pre-existing results" >> ${LOGFILE}
fi
# have we found any region that could be a putative primary target at all?
FOUNDTTS=`cat ${OUTPUT}/tts/${LOISHORT}.TTS | awk '{if (NR!=1) {print $0}}' | awk '/./{n++}; END {print n+0}'` 
test ${FOUNDTTS} == 0 && \
	( echo "[WARNING] No primary targets found. \nEither increase the loci of interest or relax the constraints for finding primary targets" >> ${LOGFILE} ) \
	&& exit 1

# are there more targets that we want to process at a time?
test ${FOUNDTTS} -gt ${MAXIMALSITES} && \
	( echo "[WARNING] Number of potential targets exceeds threshold MAXIMALSITES \n Either increase the MAXIMALSITES threshold or reduce the number of primary targets by applying more stringent constraints (parameter localConstraints) " >> ${LOGFILE} ) \
	&& exit 1

#---------------------------
# screen for off-targets
#---------------------------

[ ! -d ${OUTPUT}/tfo ] && mkdir ${OUTPUT}/tfo
echo "-------     " >> ${LOGFILE}
echo "*** "`date +%H:%M:%S`" convert primary targets to genomic locations" >> ${LOGFILE} 
cat ${OUTPUT}/tts/${LOISHORT}.TTS \
	| awk -F\\t '{if (NR!=1) {printf  "%s\t%s\t%s\tt_%02d\t%s\t%s\t%s\t%s\n", $1, $2, $3, (NR-1), $4, $5, $10, $7}}' \
	| sort -k1,1 -k2,3n -k6,6 -u \
	| ${PYTHON} ${DIR}/scripts/convertTTStoGenomicLoci.py --bed=6 ${VERBOSE} ${LOI} \
	| sort -k1,1 -k2,3n -k6,6 \
	> ${OUTPUT}/tts/${LOISHORT}.TTSpool 

echo "...design TFOs against primary targets" >> ${LOGFILE} 
cat ${OUTPUT}/tts/${LOISHORT}.TTSpool | awk -F\\t '{print ">"$4"\n"$8}' > ${OUTPUT}/tts/${LOISHORT}.TTSpool.fa
cat ${OUTPUT}/tts/${LOISHORT}.TTSpool.fa | ${PYTHON} ${DIR}/scripts/translateTTStoTFO.py ${VERBOSE} --is-fasta 1> ${OUTPUT}/tfo/${LOISHORT}.TFO 2>> ${LOGFILE}

echo "...find alternative targets for all TFOs" >> ${LOGFILE} 
[ ! -d ${OUTPUT}/tpx ] && mkdir ${OUTPUT}/tpx
TFO=${OUTPUT}/tfo/${LOISHORT}.TFO
TFOSHORT=${TFO##*/}
TFOSHORT=${TFOSHORT%.*}
echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}

if [ ! -f ${OUTPUT}/tpx/${TFOSHORT}.TPX ] || [ ! ${SKIPIFEXISTS} = "TRUE" ]; then
	echo "-------     " >> ${LOGFILE}
	echo "*** "`date +%H:%M:%S`" processing ${TFOSHORT}" >> ${LOGFILE} 
	${TRIPLEXATOR} ${TPXOPTIONS} -od ${OUTPUT}/tpx -o ${TFOSHORT}.TPX -ss ${TFO} -ds ${LOCSEQ} 

	# remove overlap with primary target
	echo "...removing primary target entries from off-target list" >> ${LOGFILE} 
	${PYTHON} ${DIR}/scripts/_removeOverlapWithPrimaryTarget.py ${VERBOSE} ${OUTPUT}/tts/${LOISHORT}.TTSpool ${OUTPUT}/tpx/${TFOSHORT}.TPX 1> ${OUTPUT}/tpx/${TFOSHORT}.TPXcleaned 2>> ${LOGFILE}
	mv ${OUTPUT}/tpx/${TFOSHORT}.TPXcleaned ${OUTPUT}/tpx/${TFOSHORT}.TPX

	echo "...annotating ${TFOSHORT}" >> ${LOGFILE}
	
	if [ "${ANNOTATION}" != "NONE" ]; then
		# discard previous annotation in case this is run repeatedly
		awk -F\\t ' {for (i=1; i<=12; i++) 
			printf("%s\t", $i) 
			printf("\n")}
		' ${OUTPUT}/tpx/${TFOSHORT}.TPX > ${OUTPUT}/tpx/${TFOSHORT}.TPXreset
		mv ${OUTPUT}/tpx/${TFOSHORT}.TPXreset ${OUTPUT}/tpx/${TFOSHORT}.TPX
		
		# get number of off-targets and extract their coordinates
		cat ${OUTPUT}/tpx/${TFOSHORT}.TPX | awk -F\\t '{if (NR!=1) {print  $4"\t"$5"\t"$6"\t"NR"\t"$7"\t"$11}}' > ${OUTPUT}/tpx/${TFOSHORT}.bed
		FOUNDOFFS=`cat ${OUTPUT}/tpx/${TFOSHORT}.bed | awk '/./{n++}; END {print n+0}'` 
		echo "...Number of off-targets found: ${FOUNDOFFS}" >> ${LOGFILE}
	
		if [ ! "${FOUNDOFFS}" = "0" ]; then
			# intersect with each annotation file
			for ANNO in `cat ${ANNOTATION} | awk -F\\t '{print $2}' | sort -u`; do
				echo "...Intersecting with ${ANNO}" >> ${LOGFILE}
				${BEDTOOLS}intersectBed -wo -a ${OUTPUT}/tpx/${TFOSHORT}.bed -b ${OUTPUT}/annotation/${ANNO}.bed | sort -k4,4n > ${OUTPUT}/tpx/${TFOSHORT}.${ANNO}
				# append annotation to tpx file (additional columns, multiple entries are comma separated)
				${PYTHON} ${DIR}/scripts/_annotationTPX.py ${OUTPUT}/tpx/${TFOSHORT}.TPX ${OUTPUT}/tpx/${TFOSHORT}.${ANNO} ${ANNO} > ${OUTPUT}/tpx/${TFOSHORT}.TPXanno
				mv ${OUTPUT}/tpx/${TFOSHORT}.TPXanno ${OUTPUT}/tpx/${TFOSHORT}.TPX
				rm ${OUTPUT}/tpx/${TFOSHORT}.${ANNO}
			done	
		fi	
	fi
	echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}
else
	echo "...skipping off-target detection due to pre-existing results" >> ${LOGFILE}
fi

#---------------------------
# find optimal primary targets
#---------------------------

[ ! -d ${OUTPUT}/json ] && mkdir -p ${OUTPUT}/json
if [ ! -f ${OUTPUT}/json/${JSON} ] || [ ! ${SKIPIFEXISTS} = "TRUE" ]; then
	echo "-------     " >> ${LOGFILE}
	echo "*** "`date +%H:%M:%S`" find optimal primary targets" >> ${LOGFILE}
	echo "...get all eligible primary targets" >> ${LOGFILE} 
	cat ${OUTPUT}/tts/${LOISHORT}.TTSpool | awk -F\\t '{print  $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ${OUTPUT}/tts/${LOISHORT}.TTSpool.bed
	${BEDTOOLS}fastaFromBed -name -s -fi ${GENOME} -bed ${OUTPUT}/tts/${LOISHORT}.TTSpool.bed -fo ${OUTPUT}/tts/${LOISHORT}.submatches
	if [ ! -f ${OUTPUT}/tts/${LOISHORT}.submatches.TTS ] || [ ! ${SKIPIFEXISTS} = "TRUE" ]; then
		echo "...collect all primary targets"  >> ${LOGFILE} 
		${TRIPLEXATOR} ${TTSOPTIONS} --all-matches -od ${OUTPUT}/tts -o ${LOISHORT}.submatches.TTS -ds ${OUTPUT}/tts/${LOISHORT}.submatches
	fi
	echo "...intersect primary targets with previously obtained off-targets"  >> ${LOGFILE} 
	${PYTHON} ${DIR}/scripts/_submatches.py --verbose ${OUTPUT}/tts/${LOISHORT}.TTSpool ${OUTPUT}/tts/${LOISHORT}.submatches.TTS ${OUTPUT}/tpx/${LOISHORT}.TPX ${OUTPUT}/tpx/${LOISHORT}.TPX.log > ${OUTPUT}/json/${JSON} 2>> ${LOGFILE}
	echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}
else
	echo "...skipping primary target detection due to pre-existing results" >> ${LOGFILE}
fi

#---------------------------
# generate JSON files
#---------------------------

if [ ! -f ${OUTPUT}/json/primary_target_regions.json ] || [ ! ${SKIPIFEXISTS} = "TRUE" ]; then
	echo "-------     " >> ${LOGFILE}
	echo "*** "`date +%H:%M:%S`" produce JSON files" >> ${LOGFILE}
	if [ "${CHROMATIN}" != "NONE" ]; then
		[ "$CHROMATINFORMAT" == "bam" ] && echo "...create bam index " >> ${LOGFILE} && samtools index ${CHROMATIN}
		${PYTHON} ${DIR}/scripts/collectOutput.py ${VERBOSE} --chromatin ${CHROMATIN} --chromatin-format ${CHROMATINFORMAT} --output-dir ${OUTPUT}/json/ --flanks ${FLANKS} --fasta-file ${LOISEQ} ${LOI} ${LOC} ${OUTPUT}/tts/${LOISHORT}.TTSpool ${OUTPUT}/tpx/${LOISHORT}.TPX 2>> ${LOGFILE}
	else
		${PYTHON} ${DIR}/scripts/collectOutput.py ${VERBOSE} --output-dir ${OUTPUT}/json/ --flanks ${FLANKS} --fasta-file ${LOISEQ} ${LOI} ${LOC} ${OUTPUT}/tts/${LOISHORT}.TTSpool ${OUTPUT}/tpx/${LOISHORT}.TPX  2>> ${LOGFILE}
	fi
	echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}
fi

#---------------------------
# generate html output
#---------------------------

echo "-------     " >> ${LOGFILE}
echo "*** "`date +%H:%M:%S`" generate HTML file" >> ${LOGFILE}
[ ! -d ${OUTPUT}/includes/ ] && mkdir ${OUTPUT}/includes/
cp -r ${DIR}/includes/* ${OUTPUT}/includes/
${PYTHON} ${DIR}/scripts/_make_html.py ${VERBOSE} --with-chromatin ${CHROMATIN} --output ${OUTPUT}/ ${LOI} ${OUTPUT}/tts/${LOISHORT}.TTSpool ${OUTPUT}/tpx/${LOISHORT}.TPX ${JSON} 2>> ${LOGFILE}
cp ${DIR}/FAQs.txt ${OUTPUT}/FAQs.txt
echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}

#---------------------------
# generate circos plots (optional)
#---------------------------

if hash ${CIRCOS} 2>&- ; then
	echo "-------     " >> ${LOGFILE}
	echo "*** "`date +%H:%M:%S`" create circos plots" >> ${LOGFILE}
	mkdir -p ${OUTPUT}/circos
	if [ "${CHROMATIN}" != "NONE" ]; then
		${PYTHON} ${DIR}/scripts/_make_circos_data.py ${VERBOSE} --chromatin ${CHROMATIN} --chromatin-format ${CHROMATINFORMAT} --output-dir ${OUTPUT}/circos/ --fasta-file ${LOISEQ} ${LOI} ${LOC} ${OUTPUT}/tts/${LOISHORT}.TTSpool ${OUTPUT}/tpx/${LOISHORT}.TPX ${OUTPUT}/json/primary_targets.json ${GENOMESIZE} 2>> ${LOGFILE}
	else
		${PYTHON} ${DIR}/scripts/_make_circos_data.py ${VERBOSE} --output-dir ${OUTPUT}/circos/ --fasta-file ${LOISEQ} ${LOI} ${LOC} ${OUTPUT}/tts/${LOISHORT}.TTSpool ${OUTPUT}/tpx/${LOISHORT}.TPX ${OUTPUT}/json/primary_targets.json ${GENOMESIZE} 2>> ${LOGFILE}
	fi
	
	ORIGINDIR=`pwd` # save current directory
	cat /dev/null > ${OUTPUT}/circos/log.txt
	# make one circos plot for each target
	for FOLDER in `ls -d ${OUTPUT}/circos/t*/`; do
		echo ${FOLDER} >> ${LOGFILE}
		cd ${FOLDER} && ${CIRCOS} -conf ./circos.conf >> ../log.txt && cd ${ORIGINDIR}
	done

	echo `pwd`

	# augment html	
	echo "<html><head><meta content='text/html; charset=UTF-8' http-equiv='content-type'>
	<script type="text/javascript" src="includes/js/jquery.min.js"></script>
	<script type="text/javascript" src="includes/js/circos.js"></script>
	<style type='text/css'>.transformed {-webkit-transform: scale(0.01, 0.01);-webkit-transform-origin: top left;-moz-transform: scale(0.01, 0.01);-moz-transform-origin: top left;-ms-transform: rscale(0.01, 0.01);-ms-transform-origin: top left;transform: scale(0.01, 0.01);}</style>
	</head><body>" > ${OUTPUT}/inspector_circos.html
	
	mkdir -p ${OUTPUT}/img
	mv ${OUTPUT}/circos/*.png ${OUTPUT}/img/
	
	# adjust html		
	for HTML in `ls -f ${OUTPUT}/circos/t*.html`; do
		TID=${HTML##*/}
		TID=${TID%.*}
				
		echo "Augmenting HTML for ${TID}" >> ${LOGFILE}
		echo "<img src='img/${TID}.png' class='transformed' usemap='#${TID}_map'>" >> ${OUTPUT}/inspector_circos.html
		cat ${HTML} >> ${OUTPUT}/inspector_circos.html
		rm ${HTML}
	done
	echo '</body></html>' >> ${OUTPUT}/inspector_circos.html
	echo "*** "`date +%H:%M:%S`" done     " >> ${LOGFILE}
fi

#---------------------------
# finish
#---------------------------
echo "-------     " >> ${LOGFILE}
echo "*** "`date +%H:%M:%S`" inspector finished" 	>> ${LOGFILE}

