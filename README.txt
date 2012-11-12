*** Triplexator - Finding Nucleic Acid Triple Helices ***

---------------------------------------------------------------------------
Table of Contents
---------------------------------------------------------------------------
  1.   Overview
  2.   Installation
  3.   Usage
  4.   Output Format
  5.   Examples
  6.   Contact

---------------------------------------------------------------------------
1. Overview
---------------------------------------------------------------------------

Triplexator is a tool for detecting nucleic acid triple helices and triplex
features in nucleotide sequences using the canonical triplex-formation
rules.

---------------------------------------------------------------------------
2. Installation 
---------------------------------------------------------------------------
---------------------------------------------------------------------------
2.1 Installation - binaries
---------------------------------------------------------------------------

Triplexator binaries are available for some plattforms from 
  http://code.google.com/p/triplexator/downloads/list

  1) untar:
     tar xf triplexator.<plattform>.tar.gz
  2) change directory:
     cd triplexator
  3) run triplexator:
     ./bin/triplexator --help
     This should output a brief usage message.

---------------------------------------------------------------------------
2.2 Installation - from source
---------------------------------------------------------------------------

Triplexator sources can be obtained from googlecode using git and build 
using cmake:

    1) obtain triplexator:
    >git clone https://code.google.com/p/triplexator/ triplexator
    
    2) change directory:
      >cd triplexator
    
    3) create directory and change into it:
      >mkdir -p build/Release && cd build/Release
    
    4) run cmake and make:
      >cmake ../.. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -G "Unix Makefiles" && make
    
    5) change directory:
      >cd ../..
    
    6) try the binary:
      >./bin/triplexator --help
    
    7) run the smoketest:
      >./demos/smoketest_triplexator.sh ./bin/triplexator 

On success, an executable file triplexator was build and a brief usage 
description has been dumped.

---------------------------------------------------------------------------
3. Usage
---------------------------------------------------------------------------

To get a short usage description of Triplexator, you can execute 
triplexator -h or  triplexator --help.

Usage: triplexator [OPTION]... -ss <SINGLE-STRANDED FILE> -ds <DUPLEX FILE>

Triplexator expects the names of one or two DNA/RNA (multi-) Fasta files and 
runs in different operative modes depending on which data are input. 


  [ -ss <FILE> ],  [ --single-strand-file <FILE> ]

  File in FASTA format that is searched for triplex-forming capability 
  (e.g. RNA or DNA). If only this file is supplied, Triplexator will search 
  and output putative Triplex Forming Oligonucleotides (TFOs) only.

  [ -ds <FILE> ],  [ --duplex-file <FILE> ]

  File in FASTA format that is searched for triplex-forming capability (e.g. DNA)
  If only these files are supplied <span id="sc">Triplexator</span> will search 
  and output Triplex Target Sites (TTSs) only.


---------------------------------------------------------------------------
3.1. Main Options
---------------------------------------------------------------------------

  [ -l NUM ],  [ --lower-length-bound NUM ]
  
  Specifies the minimum length of a TFO, TTS or triplex (TTS-TFO pair)
  
  [ -u NUM ],  [ --upper-length-bound NUM ]
  
  Specifies the maximum length of a TFO, TTS or triplex (TTS-TFO pair), 
  -1 = unrestricted  (default -1)
  
  [ -m ],  [ --triplex-motifs MOTIF1,MOTIF2,... ]
  
  Specifies the motifs from the canonical triplex-formation rules to be 
  used when searching for TFOs in the third strand:
  R - the purine motif that permit guanines (G) and adenines (A).
  Y - the pyrimidine motif that permit cytosines (C) and thymines (T).
  M - the mixed motif, purine-pyrimdine, that permit guanines (G) and 
      thymines (T).
  P - parallel binding, i.e. motifs facilitating Hoogsten bonds;
      this covers the pyrimidine motif and the purine-pyrimidine motif
	  in parallel configuration.
  A - anti-parallel binding, i.e. motifs facilitating reverse Hoogsten 
      bonds; this covers the purine motif and the purine-pyrimidine motif
	  in anti-parallel configuration.
	  
  By default all motifs are used.

  [ -mpmg NUM ],  [ --mixed-parallel-max-guanine NUM ]
  
  Specifies the maximum guanine proportion (in %) in a mixed-motif triplexes 
  (GT) to consider this feature for parallel binding (Hoogsteen bonds) 
  (default 100).
  
  As GT-TFOs can bind in either orientation this parameter can be used 
  to specify at which guanine content a GT-TFO should not be able to 
  bind in parallel orientation
  (because anti-parallel binding will always dominate due to high G content).
  
  Works in conjunction with --mixed-antiparallel-min-guanine but 
  parameters are keep separate as there may be a smooth transition
  between the binding modes.
  
  [ -mamg NUM ], [ --mixed-antiparallel-min-guanine NUM ]
  
  Specifies the minimum guanine proportion (in %) in a mixed-motif triplexes 
  (GT) to consider this feature for anti-parallel binding (reverse Hoogsteen 
  bonds) (default 0).

  As GT-TFOs can bind in either orientation this parameter can be used 
  to specify at which guanine content a GT-TFO should not be able to 
  bind in anti-parallel orientation.
  (because parallel binding will always dominate due to the low G content).
  
  Works in conjunction with --mixed-parallel-max-guanine but 
  parameters are keep separate as there may be a smooth transition
  between the binding modes.
  
  [ -e NUM ],  [ --error-rate NUM ]

  Set the maximal error-rate in % tolerated (default 5). 
  Triplexator searches for matches with an error-rate percent 
  of at most NUM. A match of a feature R with E errors has 
  error-rate of 100*(E/|R|), whereby |R| is the feature length. 
  In other words, a feature is allowed to have not more than 
  |R|*ceil(NUM)/100 errors.

  [ -E NUM ],  [ --maximal-error NUM ]
  
  Set the maximal overall error tolerated, disable with -1  (default -1).
  The maximal overall error is a hard threshold that can be used in 
  conjunction with the error-rate (see above).
  For example, in a scenario using an error-rate of 10% and a 
  maximal error of 3, the error-rate will be the limiting factor up to
  features of length 30, after which the maximal error takes over.
  
  [ -c NUM ],  [ --consecutive-errors NUM ]
  
  Sets the tolerated number of consecutive errors with respect to the
  canonical triplex rules as such were found to greatly destabilize 
  triplexes in vitro. The maximum permitted number is 3.
  
  [ -g NUM ],  [ --min-guanine NUM ]
  
  Set the minimum guanine conten in triplex features. NUM must be a value
  between 0 and 100 (default is 10). 
  The minimum guanine rate controls the ratio of guanines required in the 
  any triplex target site. For triplex-forming oligonucleotide this  
  constraint will be applied to their respective target. 
 
  [ -G NUM ],  [ --max-guanine NUM ]
  
  Set the maximum guanine conten in triplex features. NUM must be a value
  between 0 and 100 (default is 100). 
  The maximum guanine rate controls the ratio of guanines required in the 
  any triplex target site. For triplex-forming oligonucleotide this  
  constraint will be applied to their respective target. 
 
  [ -b NUM ],  [ --minimum-block-run NUM ]
  
  Sets the number of consecutive matches required in a feature 
  discarding any feature that violates this constrait (default 1).
  The rational behind this parameter is that a seed of consecutive matching
  positions of a given length is required to initiate triplex formation.
  Given the observation that central errors are more disruptive than 
  errors at the flanks of a triplex, this parameter will be especially
  effective for short features.
  
  Example: 
  feature 			valid -b	discarded with -b
  AGGAGAGtGAGAAAGA    <= 8	         >= 9
  AGGAGAGGAGAAAtGA    <= 13	         >= 14	
  
  [ -a ],   [--all-matches ]
  
  Flag indicates that all qualifying sub-matches should be processed and
  reported in addition to the longest match.
  Careful! This can result in hugh output files when searching for 
  TFO-TTS pairs (i.e. providing single-stranded and doubles-stranded input).
  		

  [ -mf NUM],  [ --merge-features NUM ]

   merge overlapping features into a cluster and report the spanning region
   Only supported for TFO and TTS detection, respectively. 
   For TFO-TTS pairs (triplexes) features are merged in the TFO and TTS
   detection phase on default.
   Any merge is performed before duplicate detection (-dd).
  		
  [ -dd NUM],  [ --detect-duplicates NUM ]

  Indiates whether and how duplicates should be detected (default 0). 
  Choices are:
  
  0 = off         do not detect any duplicates
  1 = permissive  detect duplicates in feature space, 
                  e.g. AGGGAcGAGGA != AGGGAtGAGGA
  2 = strict      detect duplicates in target space, 
                  e.g. AGGGAcGAGGA == AGGGAtGAGGA == AGGGAYGAGGA
				
  Detection of duplicates requires all input sequence to be present in 
  memory at the same time, which will increase memory consumption 
  particularly when whole genomes are under investigation. 
  
  It is further advised to enable filtering of repeat and low complexity
  regions to minimize the workload during duplicate detection. 

  [ -ssd [on|off] ], [--same-sequence-duplicates [on|off] ]
  
  Whether to count a feature copy in the same sequence as duplicates 
  or not. (default off)
 
  [ -v ],  [ --verbose ]
  
  Verbose. Print extra information and running times.

  [ -vv ],  [ --vverbose ]

  Very verbose. Like -v, but also print filtering statistics like true and
  false positives (TP/FP).

  [ -V ],  [ --version ]
  
  Print version information.

  [ -h ],  [ --help ]

  Print a brief usage summary.

---------------------------------------------------------------------------
3.2. Output Format Options
---------------------------------------------------------------------------

  [ -z ],  [ --zip ]
  
  Compress output with gzip on-the-fly.
  Requires gzip and boost libraries during compilation.
  
  [ -dl ],  [ --duplicate-locations ]
  
  If enabled, the locations of duplicates are reported for individual
  triplex features. Requires the setting of --duplicates-cutoff to > 0.

  [ -ns [on|off] ], [ --normalized-score [on|off] ]
  
  Whether to compute the triplex potenial normalized over the sequences 
  (default off).

  [ -o FILE ],  [ --output FILE ]

  Change the output filename to FILE. By default, this is the input file
  name extended by the suffix ".TFO", ".TTS" or ".TRIPLEX" depending on
  the operative mode.

  [ -od FILEDIR ],  [ --output-directory FILEDIR ]

  Specifies the output directory where the result files will be written.
  By default the current directory is used.
  
  [ -of NUM ],  [ --output-format NUM ]

  Select the output format the matches should be stored in. See section 4.
 
  [ -po ],  [ --pretty-output ]

  Pretty output indicates matches with capital letters and deviations from 
  the triplex-formation rules by small letters.
  
  [ -er NUM ],  [ --error-reference NUM ]
  
  Sets the reference to which the error should correspond (default 0)
  0 = the Watson strand of the target (TTS)
  1 = the purine strand of the target (TTS)
  2 = the third strand (TFO)
  
---------------------------------------------------------------------------
3.3. Filtration Options
---------------------------------------------------------------------------

  [ -fm NUM],  [ --filtering-mode NUM ]
  
  Method to quickly discard non-hits (default 1). 
  
  0 = brute-force approach      use no filtering, go the extra mile
  1 = q-gram filtering          filter hits using qgrams
  
  G-gram filtering will use more memory but can improve runtime.
  The greedy approach, however, will catch up on runtime
  when the q-grams get very small, due to an high error-rate,
  small minimum triplex length or disabled repeat filtering.

  The qgram weight is calculated as followed:
  min(14.0,floor((qgramThreshold -1 -minLength)/-(ceil(errorRate*minLength)+1)))

  [ -t NUM ],  [ --qgram-threshold NUM ]
  
  Minimal number of q-grams required per potential hit (default 2).
  A higher threshold means more stringent filtering therefore requiring
  fewer validations but also leads to shorter qgrams, which increases the
  number of lookups.
                                             
  [ -fr ],  [ --filter-repeats NUM ]
  
  Activates the filtering of low complexity regions and repeats in the 
  sequence data. This option can greatly decrease the memory consumption 
  and runtime of Triplexator. However, many repeat regions comply to 
  triplex-formation rules. Hence this option is deactivated by default.
  
  [ -mrl NUM ] [ --minimum-repeat-length NUM ]
  
  Only considered with -r. Specifies the minimum length of a repeat 
  
  [ -mrp NUM ],  [ --maximum-repeat-period NUM ]
  
  Only considered with -r. Maximum period that defines a repeat or low 
  complexity region. 
  
  [ -dc NUM ],  [ --duplicates-cutoff NUM ]
  
  Feature is disregarded if it occurs more often than specified with
  this cutoff. Disable filtering by setting cutoff to -1. (default -1)

---------------------------------------------------------------------------
3.4. Performance Options
---------------------------------------------------------------------------

  Performance Options are only available if Triplexator has been compiled
  with OpenMP enabled.

  [ -rm NUM ],  [ --runtime-mode NUM ]
  
  The computational bottle-neck of triplexator is the matching of TFOs with 
  their putative targets when detecting the triplexes. Therefore, 
  Triplexator can leverage multi-processor architectures on bases of OpenMP.
  
  Depending on the dataset and the computational resources different 
  runtime modes can be chosen from. See below for additional information.
  
  0 = Serial (default)
  1 = Parallelize TTSs
  2 = Parallelize duplexes
  
  In case of memory capacity issues it can be helpful to divide the 
  single-strand sequence file into several smaller chunks and to execute
  Triplexator on each of them. Is can also be helpful to use this approach
  and distribute the smaller task over a cluster.
  
  [ -p NUM ],  [ --processors NUM ]
  
  Number of processors used when executed in parallel mode.
  Specify -1 to detect automatically. (default -1)
  
---------------------------------------------------------------------------
3.4.1 Serial
---------------------------------------------------------------------------

  The default option performs the triplex-matching serially. This option
  is a good tradeof if memory is a constraint or Triplexator is run on a 
  one processor achitecture.
  
---------------------------------------------------------------------------
3.4.2 Parallelize triplex target sites (TTSs)
---------------------------------------------------------------------------

  In case the duplex sequences are rather long, i.e. chromosomes, this 
  is the appropriate runtime-mode. Parallelize triplex matches evaluates  
  one duplex sequence at a time but parallelize the matching of all its  
  putative triplex target sites when searching for suitable partners in the 
  single-strand sequence set.
 
---------------------------------------------------------------------------
3.4.3 Parallelize duplexes
---------------------------------------------------------------------------

  This is the appropriate runmode-option in case many rather small duplex 
  sequences are searched for their triplex potential. Parallelize duplexes
  reads all duplex sequences into memory and performs the triplex search 
  in parallel trading runtime for memory consumption.  
 
---------------------------------------------------------------------------
4. Output Formats
---------------------------------------------------------------------------

  [ -of NUM ],  [ --output-format NUM ]

  Triplexator supports currently 3 different output formats:

  0 = Tab-separated Format + Summary Format
  1 = Triplexator Format + Summary Format
  2 = Summary Format only
  
  All output formats are sensitive to the operative mode that Triplexator
  runs in, i.e. the results for the search of TFOs, TTSs and triplexes.
  
  In addition a log file will be generated each time Triplexator is run.
 
---------------------------------------------------------------------------
4.1. Tab-separated Format
---------------------------------------------------------------------------

  The tabulator separated file contains a header line indicated with an "#"
  followed by the meaning of each column (header) as indicated below by a 
  separating pipe symbol "|". Each line corresponds to one idividual entry 
  of a triplex feature (TFO/TTS) or triplex match.

  Searching putative triplex-forming oligonucleotides:
    Sequence-ID|Start|End|Score|Motif|Error-rate|Errors|Guanine-rate|
    			...Duplicates|TFO
	
  Searching putative triplex target sites:
	Duplex-ID|Start|End|Score|Strand|Error-rate|Errors|Guanine-rate|
				...Duplicates|TTS
	
  Searching putative triplexes:
	Sequence-ID|TFO start|TFO end|Duplex-ID|TTS start|TTS end|Score|Error-rate|
	            ...Motif|Strand|Orientation|Guanine-rate
	
  with the following meaning:
  - Sequence-ID     : id of the single-stranded sequence providing the TFO
  - Duplex-ID       : id of the double-stranded sequence providing the TTS
  - Start [TFO/TTS] : start index of the feature
  - End [TFO/TTS]   : end index of the feature
  - Score           : score of the feature/triplex (number of matches)
  - Motif           : Binding motif of the canonical triplex rules used  
                      to find this TFO (R|Y|M)
  - Error-rate      : rate of deviations from the canonical rules for this 
                      feature/triplex
  - Errors          : Deviation from the triplex rulesets are encoded with
                      respect to the participant that contains the 
                      deviation, i.e "d3" means that in the duplex 
                      nucleotide 3 does not match the rules (starting at 0). 
                      "o5" means that in the third strand oligonucleotide 
                      the 5th position deviates and "b3" means that both 
                      participants contain an error, while "t3" means that 
                      the participants are conform to their individual rules
                      but don't match. 
                      Coordinates depend on the setting of the parameter 
                      --error-reference.
  - Duplicates      : Number of times this triplex feature or respective
                      target does occur in the corresponding sequence set
					  NOT supported when searching for triplexes.
  - TFO/TTS         : sequence of the TFO/TTS
  - Strand          : strand of the duplex providing the poly-purine tract
  - Orientation     : parallel or anti-parallel binding of the TFO to the 
                      TTS [P|A]
  - Guanine-rate    : proportion of guanines w/r/t the target able to 
                      participate in triplex formation                  
---------------------------------------------------------------------------
4.2. Triplexator Format
---------------------------------------------------------------------------

  This format is complies to fasta format when searching either for TFOs or
  for TTSs. When searching TFO-TTS pairs a visual representation of the
  alignment is given as well.

---------------------------------------------------------------------------
4.3. Summary Format
---------------------------------------------------------------------------

  The summary file will be generated each time Triplexator is started. 
  The summary file is a tab-separated file that aggregates triplex feature 
  and triplex hits for the corresponding sequence(s).
  
  The absolute number of features/triplexes is given for all motifs 
  individually (abs), while the relative measure (triplex potential) is 
  adjusted for the sequence and feature length of the sequence(s).
  
  Note, to save disk space only sequences that have at least one triplex
  feature are output.
  
  Searching triplex-forming oligonucleotides:
    Sequence-ID|TFOs (abs)|TFOs (rel)|GA (abs)|GA (rel)|TC (abs)|TC (rel)|
	            ...GT (abs)|GT (rel)
	
  Searching triplex target sites:
	Duplex-ID|TTSs (abs)|TTSs (rel)
	
  Searching triplexes:
	Duplex-ID|Sequence-ID|Total (abs)|Total (rel)|GA (abs)|GA (rel)|TC (abs)|
	            ...TC (rel)|GT (abs)|GT (rel)
  
  with the following meaning:
  - Sequence-ID     : id of the single-stranded sequence providing the TFO
  - Duplex-ID       : id of the double-stranded sequence providing the TTS
  - TFOs (abs)      : absolute number of maximal TFOs found in the sequence
                      over all triplex motifs   
  - TFOs (rel)      : length-adjusted triplex potential wrt. TFO features
                      over all triplex motifs 
  - TTSs (abs)      : absolute number of maximal TTSs found in the sequence
  - TTSs (rel)      : length-adjusted triplex potential wrt. TTS features
  - Total (abs)     : absolute number of maximal TFO-TTS pairs found in the 
                      two sequences
  - Total (rel)     : length-adjusted triplex potential wrt. TFO-TTS pairs
                      found in the two sequences
  - GA/TC/GT (abs)  : absolute number of maximal TFOs or TFO-TTS pairs 
                      (depending on the context)
  - GA/TC/GT (rel)  : length-adjusted triplex potential or TFOs or TFO-TTS 
                      pairs wrt. the specified motif   
                      (depending on the context)

---------------------------------------------------------------------------
5. Examples
---------------------------------------------------------------------------

---------------------------------------------------------------------------
5.1. Identify TFOs in single-strand sequences
---------------------------------------------------------------------------

  We want to find all putative triplex-forming oligonucleotide in a set of 
  <transcripts> subject to the following specifics:
  
	- at least 20 bps in length "-l 20"
	- having at most 15% errors in the motif "-e 15"
	- we may be only interested in TFOs that form triplexes of the 
	  purine-pyrimidine motif "-m M"
	- we want to remove low complexity regions of length >= 7 and period <=1 
	  (e.g. for polyA filtering) "-fr on -mrl 7 -mrp 1" 
	- output all sites "-of 0"
	- output to the file names <transcripts>.TFO "-o <transcripts>.TFO"
	- place the results in a specific forder "-od <folder>"

  >triplexator -l 20 -e 15 -m M -fr on -mrl 7 -mrp 1 -of 0 -od <folder> 
		-o <transcripts>.TFO -ss <transcripts>.fasta 
  
---------------------------------------------------------------------------
5.2. Identify high quality putative TTSs in a genome 
---------------------------------------------------------------------------

  We want to find all putative target sites in <genome>, which comply to 
  the following specifics:
  
	- at least 15 bps in length "-l 15"
	- containing at least 50% guanines "-g 50"
	- having at most 10% pyrimidine interruptions "-e 10"
	- filtered for low complexity regions of length >= 7 and period <=3 
	  "-fr on -mrl 7 -mrp 3"
	- at most 5 duplicates in the whole genome "-dc 5"
	- output all sites "-of 0"
	- output to the file names <genome>.TTS "-o <genome>.TTS" 
	
  >triplexator -l 15 -g 50 -e 10 -fr on -mrl 7 -mrp 3 -dc 5 -of 0 
        -o <genome>.TTS -ds <genome>.fasta 

---------------------------------------------------------------------------
5.3. Identify TFO-TTS pairs in single-strand and duplex sequences
---------------------------------------------------------------------------

  We want to find all putative triplexes that can form between a set of 
  <transcripts> and <promoters>, which comply to the following specifics:
  
	- at least 20 bps in length "-l 20"
	- having at most 5% mismatches and errors "-e 5"
	- filtered for low complexity regions of length >= 7 and period <=3 
	  "-fr on -mrl 7 -mrp 3"
	- output the alignments "-of 1"
	- we like to look at the alignments so make them pretty "-po "
	- output to the file names <transcripts_promoters>.TRIPLEX 
	  "-o <transcripts_promoters>.TRIPLEX" 
	- we don't have much time but lots of memory so run in parallel, 
	  promoters are fairly short, so parallelize on duplexes "-rm 3"
	- but don't use all my processors I still have to work, I'll give you 3 
	  "-p 3"

  >triplexator -l 20 -e 5 -fr on -mrl 7 -mrp 3 -dc 5 -of 1 
		-o <transcripts_promoters>.TRIPLEX -po -rm 3 -p 3
        -ss <transcripts>.fasta -ds <promoters>.fasta 

---------------------------------------------------------------------------
6. Contact
---------------------------------------------------------------------------

For questions or comments, contact:
  Fabian Buske <fbuske@uq.edu.au>
