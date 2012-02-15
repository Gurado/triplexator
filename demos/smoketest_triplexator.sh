#!/bin/sh

USAGEMSG="usage: $(basename $0) triplexator

with "triplexator" indicating the location of the executable
"

[ $# -lt 1 ] && echo "$USAGEMSG" >&2 && exit 1
TRIPLEXATOR="NONE"

shift $(($OPTIND - 1))
TRIPLEXATOR=$1

[ "$TRIPLEXATOR" = "NONE" ] && echo "ERROR: Specify path to Triplexator executable " >&2 && exit 1
[ ! -f ${TRIPLEXATOR} ] && echo "ERROR: Triplexator executable does not exist" >&2 && exit 1


DEMOS=`dirname $0`

FAILED=0
PASSED=0

[ ! -d ${DEMOS}/reference ] && echo "ERROR: Could not find ${DEMOS}/reference data " >&2 && exit 1

[ ! -d ${DEMOS}/tests ] && mkdir ${DEMOS}/tests


echo "______________ START TESTING ... TFOs _______________"

$TRIPLEXATOR -o test_default.tfo -od ${DEMOS}/tests -of 0 --pretty-output -ss ${DEMOS}/single-stranded.fasta
if  [ -f ${DEMOS}/tests/test_default.tfo ] && [ $(diff ${DEMOS}/reference/test_default.tfo ${DEMOS}/tests/test_default.tfo | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: default settings TFO.........................OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: default settings TFO.........................FAILED"
fi

$TRIPLEXATOR --lower-length-bound 14 --filtering-mode 0 --error-rate 9 -o test_minimum_size.tfo -od ${DEMOS}/tests -of 0 --pretty-output -ss ${DEMOS}/single-stranded.fasta 
if [ -f ${DEMOS}/tests/test_minimum_size.tfo ] && [ $(diff ${DEMOS}/reference/test_minimum_size.tfo ${DEMOS}/tests/test_minimum_size.tfo | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: minimum size TFO.............................OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: minimum size TFO.............................FAILED"
fi

$TRIPLEXATOR --lower-length-bound 14 --error-rate 10 --filtering-mode 0 --max-guanine 50 --min-guanine 50 -o test_guanine50.tfo -od ${DEMOS}/tests -of 0 --pretty-output -ss ${DEMOS}/single-stranded.fasta
if [ -f ${DEMOS}/tests/test_guanine50.tfo ] && [ $(diff ${DEMOS}/reference/test_guanine50.tfo ${DEMOS}/tests/test_guanine50.tfo | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: 50% guanine rate settings TFO................OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: 50% guanine rate settings TFO................FAILED"
fi


$TRIPLEXATOR --lower-length-bound 14 --error-rate 0  --triplex-motifs R -o test_perfect_R.tfo -od ${DEMOS}/tests -of 0 --pretty-output -ss ${DEMOS}/single-stranded.fasta
if [ -f ${DEMOS}/tests/test_perfect_R.tfo ] && [ $(diff ${DEMOS}/reference/test_perfect_R.tfo ${DEMOS}/tests/test_perfect_R.tfo | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: perfect TFO using purine motif...............OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: perfect TFO using purine motif...............FAILED"
fi

$TRIPLEXATOR --lower-length-bound 14 --error-rate 0  --triplex-motifs Y -o test_perfect_Y.tfo -od ${DEMOS}/tests -of 0 --pretty-output -ss ${DEMOS}/single-stranded.fasta 
if [ -f ${DEMOS}/tests/test_perfect_Y.tfo ] && [ $(diff ${DEMOS}/reference/test_perfect_Y.tfo ${DEMOS}/tests/test_perfect_Y.tfo | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: perfect TFO using pyrimidine motif...........OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: perfect TFO using pyrimidine motif...........FAILED"
fi

$TRIPLEXATOR --lower-length-bound 14 --error-rate 0  --triplex-motifs M -o test_perfect_M.tfo -od ${DEMOS}/tests -of 0 --pretty-output -ss ${DEMOS}/single-stranded.fasta 
if [ -f ${DEMOS}/tests/test_perfect_M.tfo ] && [ $(diff ${DEMOS}/reference/test_perfect_M.tfo ${DEMOS}/tests/test_perfect_M.tfo | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: perfect TFO using purine-pyrimidine motif....OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: perfect TFO using purine-pyrimidine motif....FAILED"
fi

$TRIPLEXATOR  --error-rate 10 --duplicate-locations --detect-duplicates 1 --duplicate-cutoff 2 -o test_duplicates_filtering_permissive.tfo -od ${DEMOS}/tests -of 0 --pretty-output -ss ${DEMOS}/single-stranded.fasta 
if [ -f ${DEMOS}/tests/test_duplicates_filtering_permissive.tfo ] && [ $(diff ${DEMOS}/reference/test_duplicates_filtering_permissive.tfo ${DEMOS}/tests/test_duplicates_filtering_permissive.tfo | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: TFO permissive duplicate filtering...........OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: TFO permissive duplicate filtering...........FAILED"
fi

$TRIPLEXATOR  --error-rate 10 --duplicate-locations --detect-duplicates 2 --duplicate-cutoff 2 -o test_duplicates_filtering_strict.tfo -od ${DEMOS}/tests -of 0 --pretty-output -ss ${DEMOS}/single-stranded.fasta 
if [ -f ${DEMOS}/tests/test_duplicates_filtering_strict.tfo ] && [ $(diff ${DEMOS}/reference/test_duplicates_filtering_strict.tfo ${DEMOS}/tests/test_duplicates_filtering_strict.tfo | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: TFO strict duplicate filtering...............OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: TFO strict duplicate filtering...............FAILED"
fi

echo "______________ START TESTING ... TTSs _______________"

$TRIPLEXATOR -o test_default.tts -od ${DEMOS}/tests -of 0 --pretty-output -ds ${DEMOS}/double-stranded.fasta
if [ -f ${DEMOS}/tests/test_default.tts ] && [ $(diff ${DEMOS}/reference/test_default.tts ${DEMOS}/tests/test_default.tts | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: default settings TTS.........................OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: default settings TTS.........................FAILED"
fi

$TRIPLEXATOR --lower-length-bound 14 --filtering-mode 0 --error-rate 9 -o test_minimum_size.tts -od ${DEMOS}/tests -of 0 --pretty-output -ds ${DEMOS}/double-stranded.fasta 
if [ -f ${DEMOS}/tests/test_minimum_size.tts ] && [ $(diff ${DEMOS}/reference/test_minimum_size.tts ${DEMOS}/tests/test_minimum_size.tts | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: minimum size TTS.............................OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: minimum size TTS.............................FAILED"
fi

$TRIPLEXATOR --error-rate 10 --detect-duplicates 1 --duplicate-cutoff 0 -o test_duplicates_filtering.tts -od ${DEMOS}/tests -of 0 --pretty-output -ds ${DEMOS}/double-stranded.fasta
if [ -f ${DEMOS}/tests/test_duplicates_filtering.tts ] && [ $(diff ${DEMOS}/reference/test_duplicates_filtering.tts ${DEMOS}/tests/test_duplicates_filtering.tts | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: TTS duplicate filtering......................OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: TTS duplicate filtering......................FAILED"
fi

$TRIPLEXATOR  --lower-length-bound 14 --filtering-mode 0 --error-rate 20 --filter-repeats on --minimum-repeat-length 4 --maximum-repeat-period 1 -o test_low_complexity_region_filtering.tts -od ${DEMOS}/tests -of 0 --pretty-output -ds ${DEMOS}/double-stranded.fasta
if [ -f ${DEMOS}/tests/test_low_complexity_region_filtering.tts ] && [ $(diff ${DEMOS}/reference/test_low_complexity_region_filtering.tts ${DEMOS}/tests/test_low_complexity_region_filtering.tts | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: TTS low complexity region filtering..........OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: TTS low complexity region filtering..........FAILED"
fi

echo "______________ START TESTING ... TFO/TTS matching ___"

$TRIPLEXATOR -o test_default.triplex -od ${DEMOS}/tests -of 1 --pretty-output -ss ${DEMOS}/single-stranded.fasta -ds ${DEMOS}/double-stranded.fasta  
if [ -f ${DEMOS}/tests/test_default.triplex ] && [ $(diff ${DEMOS}/reference/test_default.triplex ${DEMOS}/tests/test_default.triplex | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: default settings triplex.....................OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: default settings triplex.....................FAILED"
fi

$TRIPLEXATOR --error-rate 10 --lower-length-bound 14 --filtering-mode 0 --error-rate 9 -o test_minimum_size.triplex -od ${DEMOS}/tests -of 1 --pretty-output -ss ${DEMOS}/single-stranded.fasta -ds ${DEMOS}/double-stranded.fasta 
if [ -f ${DEMOS}/tests/test_minimum_size.triplex ] && [ $(diff ${DEMOS}/reference/test_minimum_size.triplex ${DEMOS}/tests/test_minimum_size.triplex | wc -l) -eq 0 ]
then
	PASSED=`expr ${PASSED} + 1`
	echo "Test: minimum size triplex.........................OK"
else
	FAILED=`expr ${FAILED} + 1`
	echo "Test: minimum size triplex.........................FAILED"
fi



echo "============== FINISHED TESTING ====================="
echo "Total number of tests passed: ${PASSED}"
echo "Total number of tests failed: ${FAILED}"
echo ""
echo "Run triplexator from : ${TRIPLEXATOR}"

