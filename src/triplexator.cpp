// ==========================================================================
//                                triplexator
// ==========================================================================
// Copyright (c) 2011,2012, Fabian Buske, UQ
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Fabian Buske or the University of Queensland nor 
//       the names of its contributors may be used to endorse or promote products 
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Fabian Buske <fbuske@uq.edu.au>
// ==========================================================================

#define SEQAN_PROFILE					// enable time measuring
//#define TRIPLEX_DEBUG					// print verification regions

//#ifndef SEQAN_ENABLE_PARALLELISM
//#define SEQAN_ENABLE_PARALLELISM 1		// disable parallelism on default
//#endif

#include <seqan/platform.h>
#ifdef PLATFORM_WINDOWS
	#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
	#define SEQAN_DEFAULT_TMPDIR "./"
#endif

#if SEQAN_ENABLE_PARALLELISM
#include <seqan/parallel.h>
#include <seqan/parallel/parallel_macros.h>
#endif  // #if SEQAN_ENABLE_PARALLELISM

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include "triplexator.h"
#include "triplex.h"

#include <iostream>
#include <sstream>

using namespace std;


namespace SEQAN_NAMESPACE_MAIN
{

	void _setupCommandLineParser(CommandLineParser & parser, Options & options){
		::std::string rev = "$Revision: 12251 $";
		addVersionLine(parser, "Version 1.3.2 (30/03/2012) SeqAn Revision: " + rev.substr(11, 4) + "");
		append(options.version, "Version 1.3.2 (30/03/2012) SeqAn Revision: " + rev.substr(11, 4) + "");
		
		addTitleLine(parser, "***********************************************************************************");
		addTitleLine(parser, "*** Triplexator - Finding nucleic acid triple helices with approximate matching ***");
		addTitleLine(parser, "***              (c) Copyright 2011, 2012, 2013 by Fabian Buske                 ***");
		addTitleLine(parser, "***                 Comments, Bugs, Feedback: f.buske@uq.edu.au                 ***");
		addTitleLine(parser, "***********************************************************************************");
		addUsageLine(parser, "[OPTIONS] -ss <FASTA FILE> -ds <FASTA FILE>");
		addSection(parser, "Input:");
		addHelpLine(parser, "Triplexator will run in different modes depending on the input. Providing:");
		addHelpLine(parser, "1) only the third strand (-ss) - search for putative triplex-forming oligonucleotides (TFO)");
		addHelpLine(parser, "2) only the duplex (-ds) - search for putative triplex target sites (TTS)");
		addHelpLine(parser, "3) both - search for triplexes (matching TFO-TTS pairs)");
		addHelpLine(parser, "");
		addOption(parser, addArgumentText(CommandLineOption("ss",  "single-strand-file",    "File in FASTA format that is searched for TFOs (e.g. RNA or DNA)", OptionType::String), "<FILE>"));
		addOption(parser, addArgumentText(CommandLineOption("ds", "duplex-file", 			"File in FASTA format that is searched for TTSs (e.g. DNA)", OptionType::String), "<FILE>"));
		addSection(parser, "Main Options:");
		addOption(parser, CommandLineOption("l",  "lower-length-bound",						"minimum triplex feature length required", OptionType::Int| OptionType::Label, options.minLength));
		addOption(parser, CommandLineOption("L",  "upper-length-bound",						"maximum triplex feature length permitted, -1 = unrestricted ", OptionType::Int | OptionType::Label, options.maxLength ));
		addOption(parser, CommandLineOption("e",  "error-rate",								"set the maximal error-rate in % tolerated", OptionType::Double | OptionType::Label, (100.0 * options.errorRate)));
		addOption(parser, CommandLineOption("E",  "maximal-error",							"set the maximal overall error tolerated, disable with -1 ", OptionType::Int | OptionType::Label, options.maximalError));
		addOption(parser, CommandLineOption("c",  "consecutive-errors",						"maximum number of consecutive errors", OptionType::Int | OptionType::Label, options.maxInterruptions));
		addOption(parser, CommandLineOption("g",  "min-guanine",							"set the minimal guanine proportion required in %", OptionType::Double | OptionType::Label, 100.0 * options.minGuanineRate));
		addOption(parser, CommandLineOption("G",  "max-guanine",							"set the maximal guanine proportion allowed in %", OptionType::Double | OptionType::Label, 100.0 * options.maxGuanineRate));
		addOption(parser, addArgumentText(CommandLineOption("m", "triplex-motifs", 			"Triplex motifs allowed [R,Y,M,P,A] (default all)", OptionType::String), "MOTIF1,MOTIF2,..."));
		addHelpLine(parser, "R = purine (G and A)");
		addHelpLine(parser, "Y = pyrimidine (C and T)");
		addHelpLine(parser, "M = mixed (G and T)");
		addHelpLine(parser, "P = parallel (Hoogsteen bonds)");
		addHelpLine(parser, "A = anti-parallel (reverse Hoogsteen bonds)");
		addOption(parser, CommandLineOption("mpmg",  "mixed-parallel-max-guanine",			"maximum guanine content (%) in mixed-motif (GT) to consider parallel binding", OptionType::Double | OptionType::Label, (100.0 * options.mixed_parallel_max_guanine)));
		addOption(parser, CommandLineOption("mamg",  "mixed-antiparallel-min-guanine",		"minimum guanine content (%) in mixed-motif (GT) to consider anti-parallel binding", OptionType::Double | OptionType::Label, (100.0 * options.mixed_antiparallel_min_guanine)));
		
		addOption(parser, CommandLineOption("b",  "minimum-block-run",						"required number of consecutive matches, disable with -1", OptionType::Int | OptionType::Label, options.minBlockRun));
		addOption(parser, CommandLineOption("a",  "all-matches",							"process and report all sub-matches in addition to the longest match", OptionType::Boolean));
		addHelpLine(parser, "Careful! This can result in hugh output files when searching for TFO-TTS pairs.");
		addOption(parser, CommandLineOption("dd", "detect-duplicates",						"indicates whether and how duplicates should be detected", OptionType::Int | OptionType::Label, options.detectDuplicates));
		addHelpLine(parser, "0 = off         do not detect duplicates");
		addHelpLine(parser, "1 = permissive  detect duplicates in sequence space, e.g. AGGGAcGAGGA != AGGGAtGAGGA");	
		addHelpLine(parser, "2 = strict      detect duplicates in target space, e.g. AGGGAcGAGGA == AGGGAtGAGGA == AGGGAnGAGGA");
		addOption(parser, addArgumentText(CommandLineOption("ssd", "same-sequence-duplicates",	"whether to count a feature copy in the same sequence as duplicates or not.", OptionType::String | OptionType::Label, (options.sameSequenceDuplicates?"on":"off")), "[on|off]"));
		addOption(parser, CommandLineOption("v",  "verbose",			"verbose mode", OptionType::Boolean));
		addOption(parser, CommandLineOption("vv", "vverbose",			"very verbose mode", OptionType::Boolean));
		addSection(parser, "Filtration Options:");
		addOption(parser, CommandLineOption("fm", "filtering-mode",		"filtering mode - method to quickly discard non-hits", OptionType::Int | OptionType::Label, options.filterMode));
		addHelpLine(parser, "0 = brute-force approach      use no filtering, go the extra mile");
		addHelpLine(parser, "1 = q-gram filtering          filter hits using qgrams (benefical for features > 20 nt)");
		addOption(parser, CommandLineOption("t", "qgram-threshold",		"number of q-grams (must be > 0)", OptionType::Int | OptionType::Label, options.qgramThreshold));
		addHelpLine(parser, "A higher threshold means more stringent filtering therefore requiring fewer validations but also leads to shorter qgrams, which increases the number of lookups.");
		addOption(parser, addArgumentText(CommandLineOption("fr",  "filter-repeats",         "if enabled, disregards repeat and low-complex regions ", OptionType::String | OptionType::Label, (options.filterRepeats?"on":"off")), "[on|off]"));
		addOption(parser, CommandLineOption("mrl",  "minimum-repeat-length","minimum length requirement for low-complex regions to be filtered", OptionType::Int | OptionType::Label, options.minRepeatLength));
		addOption(parser, CommandLineOption("mrp",  "maximum-repeat-period","maximum repeat period for low-complex regions to be filtered", OptionType::Int | OptionType::Label, options.maxRepeatPeriod));
		addOption(parser, CommandLineOption("dc",   "duplicate-cutoff",		"disregard feature if it occurs more often than this cutoff, disable with -1.", OptionType::Int | OptionType::Label, options.duplicatesCutoff));
		addSection(parser, "Output Options:");
#ifdef BOOST
		addOption(parser, CommandLineOption("z", "zip",					"compress output with gzip (requires gzip & boost)", OptionType::Boolean));
#endif
		addOption(parser, CommandLineOption("mf", "merge-features","merge overlapping features into a cluster and report the spanning region", OptionType::Boolean));
		addHelpLine(parser, "Supported for TFO and TTS detection only. Merge is performed before duplicate detection.");
		addOption(parser, CommandLineOption("dl", "duplicate-locations","Report the location of duplicates", OptionType::Boolean));
		addHelpLine(parser, "Only works when duplicate cutoff is set to greater than 0.");
		addOption(parser, addArgumentText(CommandLineOption("o", "output",	"output filename (default standard out)", OptionType::String), "FILE"));
		addOption(parser, addArgumentText(CommandLineOption("od", "output-directory", 	"output will be written to this location", OptionType::String), "FILEDIR"));
		addOption(parser, CommandLineOption("of", "output-format",     "set output format", OptionType::Int | OptionType::Label, options.outputFormat));
		addHelpLine(parser, "0 = Tab-separated");
		addHelpLine(parser, "1 = Triplexator format (contains sequence/alignment)");
		addHelpLine(parser, "2 = Summary only");
		addOption(parser, CommandLineOption("po", "pretty-output",		"indicate matching/mismatching characters with upper/lower case", OptionType::Boolean));
		addOption(parser, CommandLineOption("er", "error-reference",	"reference to which the error should correspond", OptionType::Int | OptionType::Label, options.errorReference));
		addHelpLine(parser, "0 = the Watson strand of the target");
		addHelpLine(parser, "1 = the purine strand of the target");
		addHelpLine(parser, "2 = the third strand");

#if SEQAN_ENABLE_PARALLELISM
		addSection(parser, "Performance Options:");
		addOption(parser, CommandLineOption("rm", "runtime-mode",		"execution mode - parallel modes require OpenMP support during compilation", OptionType::Int | OptionType::Label, options.runtimeMode));
		addHelpLine(parser, "0 = Serial               process in serial (most memory efficient)");
		addHelpLine(parser, "1 = Parallelize TTSs     process targets per duplex in parallel (for long duplex sequences)");	
		addHelpLine(parser, "2 = Parallelize duplex   process duplexes in parallel (for short duplex sequences)");
		addHelpLine(parser, "Note: potential runtime speedup is at the cost of higher memory usage. ");
		addHelpLine(parser, "Option 2 is prone to consume lots of memory especially with heaps of duplex sequence.");
		addOption(parser, CommandLineOption("p", "processors",			"number of processors used in parallel mode. -1 = detect automatically.", OptionType::Int | OptionType::Label, options.processors));
#endif
		requiredArguments(parser, 0);
	}
	
	int _parseCommandLineAndCheck(Options & options,
								 CommandLineParser & parser,
								 int argc,
								 char const ** argv)
	{
		bool stop = !parse(parser, argc, argv);
		if (stop)
			return 1;
		if (isSetLong(parser, "help")) {
			options.showHelp = true;
			return 0;
		} else {
			options.showHelp = false;
		}
		
		if (isSetLong(parser, "version")) {
			options.showVersion = true;
			return 0;
		} else {
			options.showVersion = false;
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// Extract options
		getOptionValueLong(parser, "error-rate", options.errorRate);
		getOptionValueLong(parser, "maximal-error", options.maximalError);
		getOptionValueLong(parser, "min-guanine", options.minGuanineRate);
		getOptionValueLong(parser, "max-guanine", options.maxGuanineRate);
		if (isSetLong(parser, "consecutive-errors")){
			getOptionValueLong(parser, "consecutive-errors", options.maxInterruptions);
		}
		
		getOptionValueLong(parser, "mixed-parallel-max-guanine", options.mixed_parallel_max_guanine);
		getOptionValueLong(parser, "mixed-antiparallel-min-guanine", options.mixed_antiparallel_min_guanine);
		
		
		
		::std::string tmpVal;
		getOptionValueLong(parser, "filter-repeats", tmpVal);
		if (tmpVal == "off"){
			options.filterRepeats = false;
		} else if (tmpVal == "on"){
			options.filterRepeats = true;
		} else {
			cerr << "Unknown specification for the option filter repeats." << ::std::endl;
			stop = true;
		}
		
		if (isSetLong(parser, "minimum-repeat-length")){
			getOptionValueLong(parser, "minimum-repeat-length", options.minRepeatLength);
		}
		if (isSetLong(parser, "maximum-repeat-period")){
			getOptionValueLong(parser, "maximum-repeat-period", options.maxRepeatPeriod);
		}
		
		getOptionValueLong(parser, "output", options.output);
		getOptionValueLong(parser, "output-directory", options.outputFolder);
		
		if (empty(options.outputFolder)){
			options.outputFolder = "./";
		} else {
			// make sure the last character is a forward slash
			if (options.outputFolder[length(options.outputFolder)-1]!= '/'){
				append(options.outputFolder, '/');
			}
		}
		getOptionValueLong(parser, "output-format", options.outputFormat);
		if (isSetLong(parser, "lower-length-bound")){
			getOptionValueLong(parser, "lower-length-bound", options.minLength);
		}
		if (isSetLong(parser, "upper-length-bound")){
			getOptionValueLong(parser, "upper-length-bound", options.maxLength);
		}
		if(options.maxLength >= options.minLength){
			options.applyMaximumLengthConstraint = true;
		}
		
#if SEQAN_ENABLE_PARALLELISM
		getOptionValueLong(parser, "processors", options.processors);
		if (options.processors < 1){
			options.processors = omp_get_max_threads();
		}
		getOptionValueLong(parser, "runtime-mode", options.runtimeMode);
		if (options.processors == 1){
            options.runtimeMode = RUN_SERIAL;
        } else if(options.runtimeMode == RUN_SERIAL){
			options.processors = 1;
		} else {
			options.processors = min(options.processors, omp_get_max_threads());
			omp_set_num_threads(options.processors);
		}
#else
		options.processors = 1;
        options.runtimeMode = RUN_SERIAL;
#endif 
		
		getOptionValueLong(parser, "filtering-mode", options.filterMode);
		getOptionValueLong(parser, "error-reference", options.errorReference);
		getOptionValueLong(parser, "qgram-threshold", options.qgramThreshold);
		
		if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit
		if (isSetLong(parser, "verbose")) options._debugLevel = max(options._debugLevel, 1);
		if (isSetLong(parser, "vverbose")) options._debugLevel = max(options._debugLevel, 3);
		if (isSetLong(parser, "pretty-output")) options.prettyString = true;	
		
		getOptionValueLong(parser, "single-strand-file", tmpVal);
		if (tmpVal.length()>0){
			appendValue(options.tfoFileNames, tmpVal, Generous());
			options.tfoFileSupplied = true;
		}
		
		getOptionValueLong(parser, "duplex-file", tmpVal);
		if (tmpVal.length()>0){
			appendValue(options.duplexFileNames, tmpVal, Generous());
			options.ttsFileSupplied = true;
		}
		
		//	getOptionValueLong(parser, "duplex-file", tmpVal);
		//	unsigned int beg = 0;
		//	for(unsigned int i = 0; i<tmpVal.length(); ++i) {
		//		if (tmpVal[i] == ',') {
		//			appendValue(options.duplexFileNames, tmpVal.substr(beg, i - beg));
		//			beg = i + 1;
		//			options.ttsFileSupplied = true;
		//		}
		//	}
		//	if (beg != tmpVal.length()){
		//		appendValue(options.duplexFileNames, tmpVal.substr(beg, tmpVal.length() - beg));
		//		options.ttsFileSupplied = true;
		//	}
		
		if (options.ttsFileSupplied && options.tfoFileSupplied)
			options.runmode=TRIPLEX_TRIPLEX_SEARCH;
		else if (options.ttsFileSupplied && !options.tfoFileSupplied)
			options.runmode=TRIPLEX_TTS_SEARCH;
		else if (!options.ttsFileSupplied && options.tfoFileSupplied)
			options.runmode=TRIPLEX_TFO_SEARCH;
		else
			options.runmode = 0;
		
#ifdef BOOST
		getOptionValueLong(parser, "zip", options.compressOutput);
#endif
		if (isSetLong(parser, "duplicate-cutoff")){
			getOptionValueLong(parser, "duplicate-cutoff", options.duplicatesCutoff);
		}
		if (isSetLong(parser, "all-matches")) options.allMatches = true;	
		getOptionValueLong(parser, "minimum-block-run", options.minBlockRun);
		getOptionValueLong(parser, "detect-duplicates", options.detectDuplicates);
		getOptionValueLong(parser, "duplicate-locations", options.reportDuplicateLocations);
		if (isSetLong(parser, "merge-features")) options.mergeFeatures = true;
		
		getOptionValueLong(parser, "same-sequence-duplicates", tmpVal);
		if (tmpVal == "off"){
			options.sameSequenceDuplicates = false;
		} else if (tmpVal == "on"){
			options.sameSequenceDuplicates = true;
		} else {
			::std::cerr << "Unknown specification for the option same-sequence-duplicates." << ::std::endl;
			stop = true;
		}
		
		
		// get triplex motif if set
		getOptionValueLong(parser, "triplex-motifs", tmpVal);
		unsigned beg = 0;
		if (tmpVal.length()>0){
			options.motifGA = false;
			options.motifTC = false;
			options.motifGT_p = false;
			options.motifGT_a = false;
		}
		for(unsigned int i = 0; i<tmpVal.length(); ++i) {
			if (tmpVal[i] == ',') {
				string motif =  tmpVal.substr(beg, i - beg);
				beg = i + 1;
				if (motif == "M"){
					options.motifGT_p = true;
					options.motifGT_a = true;
				} else if (motif == "R"){
					options.motifGA = true;
				} else if (motif == "Y"){
					options.motifTC = true;
				} else if (motif == "P") {
					options.motifTC = true;
					options.motifGT_p = true;
				} else if (motif == "A") {
					options.motifGA = true;
					options.motifGT_a = true;
				}
			}
		}
		
		if (beg != tmpVal.length()){
			string motif =  tmpVal.substr(beg, tmpVal.length() - beg);
			if (motif == "M"){
				options.motifGT_p = true;
				options.motifGT_a = true;
			} else if (motif == "R"){
				options.motifGA = true;
			} else if (motif == "Y"){
				options.motifTC = true;
			} else if (motif == "P") {
				options.motifTC = true;
				options.motifGT_p = true;
			} else if (motif == "A") {
				options.motifGA = true;
				options.motifGT_a = true;
			}
		}
		//////////////////////////////////////////////////////////////////////////////
		// Check options
		if (options.runmode == 0){
			::std::cerr << "At least one type of input files has to be supplied." << endl << "Run 'triplexator --help' for more information." << ::std::endl;
			options.showHelp = true;
			return 0;
		}
		if ((options.errorRate > 20 || options.errorRate < 0) && (stop = true))
			::std::cerr << "Error-rate must be a value between 0 and 20" << ::std::endl;
		if ((options.minGuanineRate < 0 || options.minGuanineRate > 100) && (stop = true))
			::std::cerr << "Minimum guanine proportion in the triplex target site must be a value between 0 and 100" << ::std::endl;
		if ((options.maxGuanineRate < 0 || options.maxGuanineRate > 100) && (stop = true))
			::std::cerr << "Maximum guanine proportion in the triplex target site must be a value between 0 and 100" << ::std::endl;
		if ((options.minGuanineRate > options.maxGuanineRate) && (stop = true))
			::std::cerr << "Maximum guanine proportion cannot be smaller than minimum guanine proportion" << ::std::endl;
		if ((options.mixed_antiparallel_min_guanine < 0 || options.mixed_antiparallel_min_guanine > 100) && (stop = true))
			::std::cerr << "Min guanine proportion antiparallel mixed motif TFOs must be a value between 0 and 100" << ::std::endl;
		if ((options.mixed_parallel_max_guanine < 0 || options.mixed_parallel_max_guanine > 100) && (stop = true))
			::std::cerr << "Max guanine proportion parallel mixed motif TFOs must be a value between 0 and 100" << ::std::endl;
		if ((options.minLength < 10) && (stop = true))
			::std::cerr << "Minimum triplex length should be greater or equal than 10. " << options.minLength << ::std::endl;
		if ((options.maxLength > 1000) && (stop = true))
			::std::cerr << "Maximum triplex length needs to be smaller or equal than 1000. " << options.maxLength << ::std::endl;
		if ((options.maxInterruptions > 3) && (stop = true))
			::std::cerr << "Maximum consecutive interruptions needs to be smaller or equal than 3." << options.maxInterruptions << ::std::endl;
		if ((options.outputFormat > 2) && (stop = true))
			::std::cerr << "Invalid output format option." << ::std::endl;
		if (! (options.runtimeMode==RUN_SERIAL || options.runtimeMode==RUN_PARALLEL_DUPLEX || options.runtimeMode==RUN_PARALLEL_TRIPLEX || options.runtimeMode==RUN_PARALLEL_STRANDS) && (stop = true))
			::std::cerr << "Runtime mode not known" << ::std::endl;
		if (options.duplicatesCutoff >= 0 && options.detectDuplicates == DETECT_DUPLICATES_OFF && (stop = true))
			::std::cerr << "Duplicate filtering with specified cutoff requires duplicate detection mode to be enabled" << ::std::endl;
		if (! (options.filterMode==BRUTE_FORCE || options.filterMode==FILTERING_GRAMS) && (stop = true))
			::std::cerr << "Filtering mode not known" << ::std::endl;
		if (! (options.errorReference==WATSON_STAND || options.errorReference==PURINE_STRAND || options.errorReference==THIRD_STRAND) && (stop = true))
			::std::cerr << "Error reference not known" << ::std::endl;
		if ((options.errorReference==WATSON_STAND || options.errorReference==PURINE_STRAND) && options.runmode==TRIPLEX_TFO_SEARCH)
			::std::cerr << "Note: reference defaulted to thrid strand for TFO search" << ::std::endl;
		if (options.errorReference==THIRD_STRAND && options.runmode==TRIPLEX_TTS_SEARCH)
			::std::cerr << "Note: reference defaulted to Watson strand for TTS search" << ::std::endl;
		if (options.qgramThreshold <= 0 && (stop = true))
			::std::cerr << "qgram theshhold needs to be positive, otherwise filtering is void" << ::std::endl;
		
		options.errorRate = options.errorRate / 100.0;
		options.minGuanineRate = options.minGuanineRate / 100.0;
		options.maxGuanineRate = options.maxGuanineRate / 100.0;
		options.mixed_parallel_max_guanine = options.mixed_parallel_max_guanine / 100.0;
		options.mixed_antiparallel_min_guanine = options.mixed_antiparallel_min_guanine / 100.0;
		options.minGuanine = static_cast<unsigned>(ceil(options.minLength * options.minGuanineRate));
		options.tolError = static_cast<unsigned>(floor(options.errorRate*options.minLength));
		
		if (options.errorRate == 0 || options.maximalError == 0)
			options.maxInterruptions=0;
		
		if (options.applyMaximumLengthConstraint){
			if (options.maximalError < 0){
				options.maximalError = static_cast<int>(floor(options.errorRate * options.maxLength));
			} else {
				options.maximalError = min(options.maximalError, static_cast<int>(floor(options.errorRate * options.maxLength)));
			}
		}
		
		//	optimizing shape/q-gram for threshold >= 2 
		if (options.filterMode == FILTERING_GRAMS && options.runmode==TRIPLEX_TRIPLEX_SEARCH){
			int qgram = _calculateShape(options);
			if (qgram <= 4 && (stop = true)){
				::std::cerr << "Error-rate, minimum length and qgram-threshold settings do not allow for efficient filtering with q-grams of weight >= 5 (currently " << qgram << ")." << ::std::endl;
				::std::cerr << "Consider disabling filtering-mode (brute-force approach)" << ::std::endl;
			}
		}
		if ((options.minBlockRun > options.minLength - 2*options.tolError) && (stop = true)) {
			::std::cerr << "Block match too large given minimum length constraint and error rate." << ::std::endl;
		} 
		
		if (stop)
		{
			::std::cerr << "Exiting ..." << ::std::endl;
			return TRIPLEX_INVALID_OPTIONS;
		} else {
			return 0;	
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Find triplexes in many duplex sequences (import from Fasta)
	void _populateLogFile(int argc, const char *argv[], Options	&options)
	{
		options.logFileHandle << "*************************************************************" << ::std::endl;
		options.logFileHandle << "*** Triplexator - Finding nucleic acid triple helices     ***" << ::std::endl;
		options.logFileHandle << "***         (c) Copyright 2011 by Fabian Buske            ***" << ::std::endl;
		options.logFileHandle << "***     Comments, Bugs, Feedback: f.buske@uq.edu.au       ***" << ::std::endl;
		options.logFileHandle << "*************************************************************" << ::std::endl;
		options.logFileHandle << "*** " << options.version << ::std::endl;
		
		options.logFileHandle << "*** COMMAND:" << ::std::endl;
		options.logFileHandle << ">";
		for (int i=0;i<argc;++i){
			options.logFileHandle << argv[i] << " ";
		}
		options.logFileHandle << ::std::endl;
		options.logFileHandle << "*************************************************************" << ::std::endl;
		options.logFileHandle << "*** INPUT:" << ::std::endl;
		options.logFileHandle << "- single-stranded file supplied : " << (options.tfoFileSupplied?"Yes":"No") << ::std::endl;
		options.logFileHandle << "- duplex file supplied : " << (options.ttsFileSupplied?"Yes":"No") << ::std::endl;
		options.logFileHandle << "-> ";
		switch (options.runmode) {
			case TRIPLEX_TTS_SEARCH:
				options.logFileHandle  << "search putative triplex target sites" << ::std::endl;
				break;
			case TRIPLEX_TFO_SEARCH:
				options.logFileHandle  << "search putative triplex-forming oligonucleotides" << ::std::endl;
				break;
			case TRIPLEX_TRIPLEX_SEARCH:
				options.logFileHandle  << "search putative triplexes (matching triplex-forming oligonucleotides and target sites)" << ::std::endl;
				break;
			default:
				break;
		}		
		options.logFileHandle << "*************************************************************" << ::std::endl;
		options.logFileHandle << "*** Output Options:" << ::std::endl;
		options.logFileHandle << "- output directory : " << options.outputFolder << ::std::endl;
		options.logFileHandle << "- output file : " << options.output << ::std::endl;
		options.logFileHandle << "- output format : ";
		switch (options.outputFormat) {
			case FORMAT_BED:
				options.logFileHandle << FORMAT_BED << " = Triplex" << ::std::endl;
				break;
			case FORMAT_TRIPLEX:
				options.logFileHandle << FORMAT_TRIPLEX << " = extended Triplex (+Alignment)" << ::std::endl;
				break;
			case FORMAT_SUMMARY:
				options.logFileHandle << FORMAT_SUMMARY << " = Summary (tsv)" << ::std::endl;
				break;
			default:
				break;
		}
		options.logFileHandle << "- merge features : " << (options.runmode == TRIPLEX_TRIPLEX_SEARCH || options.mergeFeatures?"Yes":"No") << ::std::endl;
		options.logFileHandle << "- report duplicate locations : " << (options.reportDuplicateLocations?"Yes":"No") << ::std::endl;
	#ifdef BOOST
		options.logFileHandle << "- compress output : " << (options.compressOutput?"Yes":"No") << ::std::endl;
	#endif	
		options.logFileHandle << "- error reference : " ;
		switch (options.errorReference) {
			case WATSON_STAND:
				options.logFileHandle << WATSON_STAND << " = Watson strand (TTS)" << ::std::endl;
				break;
			case PURINE_STRAND:
				options.logFileHandle << PURINE_STRAND << " = purine strand (TTS)" << ::std::endl;
				break;
			case THIRD_STRAND:
				options.logFileHandle << THIRD_STRAND << " = third strand (TFO)" << ::std::endl;
				break;
			default:
				break;
		}		
		options.logFileHandle << "*************************************************************" << ::std::endl;
		options.logFileHandle << "*** Main Options:" << ::std::endl;
	//	options.logFileHandle << "- consider forward strand in duplex : " << (options.forward?"Yes":"No") << ::std::endl;
	//	options.logFileHandle << "- consider reverse strand in duplex : " << (options.reverse?"Yes":"No") << ::std::endl;
		options.logFileHandle << "- maximum error-rate : " << (options.errorRate*100) << "%" << ::std::endl;
		if (options.maximalError>=0)
			options.logFileHandle << "- maximum total error : " << (options.maximalError) << ::std::endl;	
		else
			options.logFileHandle << "- maximum total error : " << "not specified" << ::std::endl;	
		
		options.logFileHandle << "- minimum guanine content with respect to the target : " << (options.minGuanineRate*100) << "%" << ::std::endl;
		options.logFileHandle << "- maximum guanine content with respect to the target : " << (options.maxGuanineRate*100) << "%" << ::std::endl;
		
		options.logFileHandle << "- minimum length : " << options.minLength << " nucleotides" << ::std::endl;
		if (!options.applyMaximumLengthConstraint)
			options.logFileHandle << "- maximum length : omitted" << ::std::endl;
		else 
			options.logFileHandle << "- maximum length : " << options.maxLength << " nucleotides" << ::std::endl;


		options.logFileHandle << "- maximum number of tolerated consecutive pyrimidine interruptions in a target: " << options.maxInterruptions << ::std::endl;
		if (options.runmode == TRIPLEX_TRIPLEX_SEARCH || options.runmode == TRIPLEX_TFO_SEARCH){
			options.logFileHandle << "- include GT-motif : " << (options.motifGT_a || options.motifGT_p?"Yes":"No") << ::std::endl;
			options.logFileHandle << "- include GA-motif : " << (options.motifGA?"Yes":"No") << ::std::endl;
			options.logFileHandle << "- include TC-motif : " << (options.motifTC?"Yes":"No") << ::std::endl;
		}
		
		if (options.motifGT_a || options.motifGT_p){
			options.logFileHandle << "- minimum guanine-percentage in anti-parallel mixed motif TFOs : " << (options.mixed_antiparallel_min_guanine*100) << "%" << ::std::endl;
			options.logFileHandle << "- maximum guanine-percentage in parallel mixed motif TFOs : " << (options.mixed_parallel_max_guanine*100) << "%" << ::std::endl;
		}
		
		options.logFileHandle << "- number of consecutive matches required in a feature : " << (options.minBlockRun) << "" << ::std::endl;
		
		
		if (!options.allMatches)
			options.logFileHandle << "- longest match only : yes" << ::std::endl;
		else 
			options.logFileHandle << "- longest match only : no ( report all matches )" << ::std::endl;
		
		options.logFileHandle << "- detect duplicates : ";
		switch (options.detectDuplicates) {
			case DETECT_DUPLICATES_OFF:
				options.logFileHandle << DETECT_DUPLICATES_OFF << " = off" << ::std::endl;
				break;
			case DETECT_DUPLICATES_PERMISSIVE:
				options.logFileHandle << DETECT_DUPLICATES_PERMISSIVE << " = permissive" << ::std::endl;
				break;
			case DETECT_DUPLICATES_STRICT:
				options.logFileHandle << DETECT_DUPLICATES_STRICT << " = strict" << ::std::endl;
				break;
			default:
				break;
		}
		options.logFileHandle << "- same sequence duplicates : " << (options.sameSequenceDuplicates?"on":"off") << ::std::endl;	
		
		
		options.logFileHandle << "*************************************************************" << ::std::endl;
		options.logFileHandle << "*** Filtration Options :" << ::std::endl;

		options.logFileHandle << "- filter repeats : " << (options.filterRepeats?"Yes":"No") << ::std::endl;
		if (options.filterRepeats){
			options.logFileHandle << "- minimum repeat length : " << options.minRepeatLength << ::std::endl;
			options.logFileHandle << "- maximum repeat period : " << options.maxRepeatPeriod << ::std::endl;
		}
		options.logFileHandle << "- duplicate cutoff : " << options.duplicatesCutoff << ::std::endl;
		if (options.runmode == TRIPLEX_TRIPLEX_SEARCH){
			if (options.filterMode == FILTERING_GRAMS){
				options.logFileHandle << "- filtering : qgrams" << ::std::endl;
				options.logFileHandle << "- weight : " << length(options.shape) << ::std::endl;
				options.logFileHandle << "- min. threshold specified: " << options.qgramThreshold << ::std::endl;
				int minSeedsThreshold = static_cast<int>(options.minLength+1-(min(static_cast<int>(ceil(options.errorRate*options.minLength)), options.maximalError)+1)*length(options.shape));
				options.logFileHandle << "- min. threshold actual: " << minSeedsThreshold << ::std::endl;			
			} else {
				options.logFileHandle << "- filtering : none - brute force" << ::std::endl;
			}
		}
		options.logFileHandle << "*************************************************************" << ::std::endl;
		options.logFileHandle << "*** Runtime mode:" << ::std::endl;
		options.logFileHandle << "- OpenMP support : ";
	#if SEQAN_ENABLE_PARALLELISM	
		options.logFileHandle << "Yes" << ::std::endl;
	#endif 
	#ifndef SEQAN_ENABLE_PARALLELISM	
		options.logFileHandle << "No" << ::std::endl;
	#endif 	
		
		options.logFileHandle << "- runtime mode : ";
		switch (options.runtimeMode) {
			case RUN_SERIAL:
				options.logFileHandle << RUN_SERIAL << " = serial" << ::std::endl;
				break;
	#if SEQAN_ENABLE_PARALLELISM			
			case RUN_PARALLEL_TRIPLEX:
				options.logFileHandle << RUN_PARALLEL_TRIPLEX << " = parallel (target sites) - " << options.processors << " processors" << ::std::endl;
				break;
			case RUN_PARALLEL_STRANDS:
				options.logFileHandle << RUN_PARALLEL_STRANDS << " = parallel (strands) - " << options.processors << " processors" << ::std::endl;
				break;
			case RUN_PARALLEL_DUPLEX:
				options.logFileHandle << RUN_PARALLEL_DUPLEX << " = parallel (duplex sequences) - " << options.processors << " processors" << ::std::endl;
				break;
	#endif
			default:
				break;
		}
		options.logFileHandle << "*************************************************************" << ::std::endl;
		options.logFileHandle << "*** Log messages:" << ::std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Find triplexes in many duplex sequences (import from Fasta)
	template <
	typename TMotifSet,
	typename TFile,
	typename TShape>
	int _findTriplex(TMotifSet						&tfoMotifSet,
					 StringSet<CharString> const	&tfoNames,
					 TFile							&outputfile,
					 Options						&options,
					 TShape const					&shape)
	{
		typedef Index<TMotifSet, IndexQGram<TShape, OpenAddressing> >				TQGramIndex;
		typedef Pattern<TQGramIndex, QGramsLookup< TShape, Standard_QGramsLookup > > TPattern;
		
		typedef __int64															TId;
		typedef Gardener<TId, GardenerUngapped>									TGardener;
		
		unsigned errorCode = TRIPLEX_NORMAL_PROGAM_EXIT;
		
		SEQAN_PROTIMESTART(find_time);
		options.logFileHandle << _getTimeStamp() << " * Started searching for triplexes" << ::std::endl;
		
		TId duplexSeqNo = 0;
		// open duplex file
		options.logFileHandle << _getTimeStamp() << " * Processing " << options.duplexFileNames[0] << ::std::endl;
	#if SEQAN_ENABLE_PARALLELISM	
		// run in parallel if requested
		if (options.runtimeMode==RUN_PARALLEL_DUPLEX){
			if (options.filterMode == FILTERING_GRAMS){
				// create index
				if (options._debugLevel >= 1)
					options.logFileHandle << _getTimeStamp() <<  " - Started creating q-gram index for all TFOs" << ::std::endl;
				TQGramIndex index_qgram(tfoMotifSet);
				resize(indexShape(index_qgram), weight(shape));
				// create pattern	
				TPattern pattern(index_qgram,shape);
				options.timeFindTriplexes = 0;
				// create index
				if (options._debugLevel >= 1)
					options.logFileHandle << _getTimeStamp() <<  " - Finised creating q-gram index for all TFOs" << ::std::endl;
				
				errorCode = startTriplexSearchParallelDuplex(tfoMotifSet, tfoNames, pattern, outputfile, duplexSeqNo, options, TGardener());
			} else {
                TQGramIndex pattern;
				errorCode = startTriplexSearchParallelDuplex(tfoMotifSet, tfoNames, pattern, outputfile, duplexSeqNo, options, BruteForce());
			}
		} else {
		// otherwise go for serial processing
	#endif	
			if (options.filterMode == FILTERING_GRAMS){
				// create index
				if (options._debugLevel >= 1)
					options.logFileHandle << _getTimeStamp() <<  " - Started creating q-gram index for all TFOs" << ::std::endl;
				TQGramIndex index_qgram(tfoMotifSet);
				resize(indexShape(index_qgram), weight(shape));
				// create pattern	
				TPattern pattern(index_qgram,shape);
				options.timeFindTriplexes = 0;
				// create index
				if (options._debugLevel >= 1)
					options.logFileHandle << _getTimeStamp() <<  " - Finised creating q-gram index for all TFOs" << ::std::endl;
				
				errorCode = startTriplexSearchSerial(tfoMotifSet, tfoNames, pattern, outputfile, duplexSeqNo, options, TGardener());
			} else {
                TQGramIndex pattern;
				errorCode = startTriplexSearchSerial(tfoMotifSet, tfoNames, pattern, outputfile, duplexSeqNo, options, BruteForce());
			}	
	#if SEQAN_ENABLE_PARALLELISM	
		}
	#endif
		
		if (errorCode == TRIPLEX_NORMAL_PROGAM_EXIT){
			options.logFileHandle << _getTimeStamp() << " * Finished processing " << options.duplexFileNames[0] << ::std::endl; 
			options.timeFindTriplexes += SEQAN_PROTIMEDIFF(find_time);	
			options.logFileHandle << _getTimeStamp() << " * Finished searching for triplexes  within " << ::std::setprecision(3) << options.timeFindTriplexes << " seconds (summed over all cpus)" << ::std::endl;
		}
		return errorCode;
	}

	////////////////////////////////////////////////////////////////////////////////
	//// Main triplex mapper function
	template <typename TOligoSet, typename TMotifSet>
	int mapTriplexes(Options &options)
	{
		typedef typename Iterator<TOligoSet, Standard>::Type 	TOligoIter;
		typedef typename Iterator<TMotifSet, Standard>::Type 	TIterMotifSet;

		TOligoSet				oligoSequences;
		StringSet<CharString>	oligoNames;		// tfo names, taken from the Fasta file
		StringSet<CharString> 	duplexNames;	// tts names, taken from the Fasta file
		String< Pair<CharString, unsigned> >  ttsnoToFileMap;
		typedef Repeat<unsigned, unsigned>						TRepeat;
		typedef String<TRepeat>									TRepeatString;
		typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
		
		// circumvent numerical obstacles
	//	options.errorRate += 0.0000001;

		//////////////////////////////////////////////////////////////////////////////
		// Step 1: verify duplex file
		options.logFileHandle << _getTimeStamp() << " * Started checking duplex file" << ::std::endl;
		// try opening each duplex file once before running the whole searching procedure
		int filecount = 0;
		int numTTSFiles = 0;
		while(filecount < numTTSFiles){
			::std::ifstream file;
			file.open(toCString(options.duplexFileNames[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
			if (!file.is_open())
				return TRIPLEX_READFILE_FAILED;
			file.close();
			++filecount;
		}
		options.logFileHandle << _getTimeStamp() << " * Finished checking duplex file" << ::std::endl;
		
		//////////////////////////////////////////////////////////////////////////////
		// Step 2: read in TFO files

		options.logFileHandle << _getTimeStamp() << " * Started reading single-stranded file:" << options.tfoFileNames[0] << ::std::endl;
		
		if (!_loadOligos(oligoSequences, oligoNames, toCString(options.tfoFileNames[0]), options)) {
			options.logFileHandle << "ERROR: Failed to load single-stranded sequence set" << ::std::endl;
			cerr << "Failed to load sequence set containing the triplex forming oligonucleotides" << endl;
			return TRIPLEX_TFOREAD_FAILED;
		}

		options.logFileHandle << _getTimeStamp() << " * Finished reading single-stranded file (" << length(oligoSequences) << " sequences read)" << ::std::endl;
		
		//////////////////////////////////////////////////////////////////////////////
		// Step 3:  pre-process all sequences with the requested TFO motifs

		SEQAN_PROTIMESTART(find_time);
		options.logFileHandle << _getTimeStamp() << " * Started detecting triplex-forming oligonucleotides in single-stranded sequences" << ::std::endl;
		
		Shape<Triplex, SimpleShape > ungappedShape;
		
	#ifdef TRIPLEX_DEBUG
		::std::cout << options.shape << ::std::endl;
	#endif
		
		if (!stringToShape(ungappedShape, options.shape)){
			return TRIPLEX_SHAPE_FAILED;
		}

		bool reduceSet = true; // merge overlapping features
		
		unsigned oligoSeqNo = 0;
		TMotifSet tfoMotifSet;
		for (TOligoIter it = begin(oligoSequences); it != end(oligoSequences); ++it){
			
			// find low complexity regions and mask sequences if requested
			if (options.filterRepeats){
				TRepeatString data_repeats;
				findRepeats(data_repeats, *it, options.minRepeatLength, options.maxRepeatPeriod);
				for (TRepeatIterator rbeg = begin(data_repeats); rbeg != end(data_repeats); ++rbeg){
					TRepeat repeat = *rbeg;
					CharString replacement = string(repeat.endPosition-repeat.beginPosition, 'N' );
					replace(*it, repeat.beginPosition, repeat.endPosition, replacement);
				}
			}
			
			// process TC motif
			if (options.motifTC) {
				processTCMotif(tfoMotifSet, *it, oligoSeqNo, reduceSet, options);
			}
			// process GA motif
			if (options.motifGA) {
				processGAMotif(tfoMotifSet, *it, oligoSeqNo, reduceSet, options);
			}
			// process GT motif
			if (options.motifGT_p && options.motifGT_a) {
				processGTMotif(tfoMotifSet, *it, oligoSeqNo, TRIPLEX_ORIENTATION_BOTH, reduceSet, options);
			} else if (options.motifGT_p){
				processGTMotif(tfoMotifSet, *it, oligoSeqNo, TRIPLEX_ORIENTATION_PARALLEL, reduceSet, options);
			} else if (options.motifGT_a){
				processGTMotif(tfoMotifSet, *it, oligoSeqNo, TRIPLEX_ORIENTATION_ANTIPARALLEL, reduceSet, options);
			}
			++oligoSeqNo;
		}
		
		// any business with duplicates?
		if (options.detectDuplicates != DETECT_DUPLICATES_OFF){

			// detect duplicates if requested
			if (options.detectDuplicates == DETECT_DUPLICATES_STRICT)
				_countDuplicatesStrict(tfoMotifSet, options, TFO() );
			else
				_countDuplicatesPermissive(tfoMotifSet, options, TFO() );
			
			if (options.duplicatesCutoff >= 0){
				unsigned removed = _filterDuplicatesWithCutoff(tfoMotifSet, options);
				options.logFileHandle << _getTimeStamp() << " * Duplicate filtering removed " << ::std::setprecision(3)  << removed << " entries (" << length(tfoMotifSet) << " remain)" << ::std::endl;
			}
		}
		
		options.timeFindTfos += SEQAN_PROTIMEDIFF(find_time);	

		options.logFileHandle << _getTimeStamp() << " * Finished detecting TFOs within " << ::std::setprecision(3)  << options.timeFindTfos << " seconds (" << length(tfoMotifSet) << " TFOs detected)" << ::std::endl;
		
	#ifdef TRIPLEX_DEBUG
		TIterMotifSet itr = begin(tfoMotifSet,Standard());
		TIterMotifSet itrEnd = end(tfoMotifSet,Standard());
		::std::cout << "printing all tfo segments" << ::std::endl;
		while (itr != itrEnd){
			::std::cout << "tfo: " << tfoString(*itr) << " type: " << (*itr).motif << " length: "<< length(*itr) <<  " position: "<< beginPosition(*itr) << " " << (*itr).parallel << ::std::endl;
			++itr;
		}
	#endif

		//////////////////////////////////////////////////////////////////////////////
		// Step 4: prepare output file & scan duplex sequences with tfo patterns
		// 
	
		unsigned errorCode = TRIPLEX_NORMAL_PROGAM_EXIT;
	#ifdef BOOST
		if (options.compressOutput){
			io::filtering_ostream filterstream;
			filterstream.push(io::gzip_compressor());
			::std::ofstream filehandle;
			if (!empty(options.output) && options.outputFormat!=2){
				openOutputFile(filterstream, options);
			} else {
				filterstream.push(::std::cout);	
			}
			printTriplexHeader(filterstream, options);
			errorCode = _findTriplex(tfoMotifSet, oligoNames, filterstream, options, ungappedShape);
			closeOutputFile(filterstream, options);
		} else {
	#endif
			::std::ofstream filehandle;
			if (!empty(options.output) && options.outputFormat!=2){
				openOutputFile(filehandle, options);
				printTriplexHeader(filehandle, options);
				errorCode = _findTriplex(tfoMotifSet, oligoNames, filehandle, options, ungappedShape);
				closeOutputFile(filehandle, options);
			} else {
				printTriplexHeader(::std::cout, options);
				errorCode = _findTriplex(tfoMotifSet, oligoNames, ::std::cout, options, ungappedShape);
			}
	#ifdef BOOST
		}
	#endif
		return errorCode;
	}

	template <
	typename TDuplexName, 
	typename TInput, 
	typename TSeqNo, 
	typename TMap, 
	typename TString, 
	typename TOutput
	>
	inline void _investigateTTSconsecutively(TDuplexName filename, 
											 TInput &file, 
											 TSeqNo &seqNo, 
											 TMap &ttsnoToFileMap,
											 TString &duplexNames, 
											 TOutput &outputhandle, 
											 Options &options)
	{
		typedef Repeat<unsigned, unsigned>						TRepeat;
		typedef String<TRepeat>									TRepeatString;
		typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
		typedef ::std::vector<unsigned>							THitList;
		typedef TriplexPotential<TId>							TPotential;
		typedef typename ::std::list<TPotential>				TPotentials;
		
		options.logFileHandle << _getTimeStamp() << " * Starting search for targets in serial mode " << ::std::endl;
		
		TDuplex	duplexString;
		CharString	id;
		unsigned seqNoWithinFile = 0;
		THitList ttsMaxMatchList;
		TPotentials potentials;
		
		for(; !_streamEOF(file); ++seqNo,++seqNoWithinFile){
			resize(ttsMaxMatchList, seqNoWithinFile+1);
			TPotential potential(seqNoWithinFile);
			
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Started reading next duplex sequence " << ::std::endl;

			TTargetSet ttsSet;
			readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
			appendValue(duplexNames, id, Generous());
			
			read(file, duplexString, Fasta());			// read Fasta sequence
			ttsnoToFileMap.insert(::std::make_pair<unsigned,::std::pair< ::std::string,unsigned> >(seqNo,::std::make_pair< ::std::string,unsigned>(filename,seqNoWithinFile)));
			
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Finished reading next duplex sequence" << ::std::endl;

			bool reduceSet = false || options.mergeFeatures; //merge overlapping features on request
			
	#ifdef TRIPLEX_DEBUG
			::std::cout << ">" << id << ::std::endl << duplexString << ::std::endl;
	#endif
			
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Started processing duplex " << id << ::std::endl;

			// find low complexity re_countDuplicatesgions and mask sequences if requested
			if (options.filterRepeats){
				if (options._debugLevel > 1 )
					options.logFileHandle << _getTimeStamp() << "   ... Started low complexity filtering of duplex sequence" << ::std::endl;

				TRepeatString	data_repeats;
				findRepeats(data_repeats, duplexString, options.minRepeatLength, options.maxRepeatPeriod);
				for (TRepeatIterator rbeg = begin(data_repeats); rbeg != end(data_repeats); ++rbeg){
					TRepeat repeat = *rbeg;
					CharString replacement = string(repeat.endPosition-repeat.beginPosition, 'N' );
					replace(duplexString, repeat.beginPosition, repeat.endPosition, replacement);
				}
				
				if (options._debugLevel > 1 )
					options.logFileHandle << _getTimeStamp() << "   ... Finished low complexity filtering of duplex sequence" << ::std::endl;

			}
			
			if (options.forward){
				if (options._debugLevel > 1 )
					options.logFileHandle << _getTimeStamp() << "   ... Started processing the Watson strand" << ::std::endl;

				unsigned ttsMaxMatches = processDuplex(ttsSet, duplexString, seqNoWithinFile, true, reduceSet, options);
				addCount(potential, ttsMaxMatches, '+');
				
				if (options._debugLevel > 1 )
					options.logFileHandle << _getTimeStamp() << "   ... Finished processing the Watson strand" << ::std::endl;

			}
			
			if (options.reverse){
				
				if (options._debugLevel > 1 )
					options.logFileHandle << _getTimeStamp() << "   ... Started processing the Crick strand" << ::std::endl;

				unsigned ttsMaxMatches = processDuplex(ttsSet, duplexString, seqNoWithinFile, false, reduceSet, options);
				addCount(potential, ttsMaxMatches, '-');

				if (options._debugLevel > 1 )
					options.logFileHandle << _getTimeStamp() << "   ... Finished processing the Crick strand" << ::std::endl;
				
			}
			
	#ifdef TRIPLEX_DEBUG
			typedef typename Iterator<TTargetSet>::Type  TIterMotifSet;
			::std::cout << "printing all tts segments" << ::std::endl;
			for (TIterMotifSet itr=begin(ttsSet); itr != end(ttsSet);++itr){
				::std::cout << "tts: " << ttsString(*itr) << " type: " << (*itr).motif << " length: "<< length(*itr) <<  " position: "<< beginPosition(*itr) << " " << ::std::endl;
			}
	#endif
			
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Started outputing results " << ::std::endl;

			
			// get norm for sequence
			setNorm(potential, length(duplexString), options);
			// add potential to list
			appendValue(potentials, potential);
			
			dumpTtsMatches(outputhandle, ttsSet, duplexNames, options);
			
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Finished outputing results " << ::std::endl;

			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Finished processing duplex " << id << ::std::endl;

		}
		dumpSummary(potentials, duplexNames, options, TTS());
		file.close();	
	}
	

	template <
	typename TDuplexName, 
	typename TInput, 
	typename TSeqNo, 
	typename TMap, 
	typename TString, 
	typename TOutput
	>
	inline void _investigateTTSsimultaneous(TDuplexName filename, 
											TInput &file, 
											TSeqNo &seqNo, 
											TMap &ttsnoToFileMap,
											TString &duplexNames, 
											TOutput &outputhandle, 
											Options &options)
	{
		typedef Repeat<unsigned, unsigned>						TRepeat;
		typedef String<TRepeat>									TRepeatString;
		typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
		typedef typename Iterator<TTriplexSet, Standard>::Type 	TIter;
		typedef typename Iterator<TTargetSet, Standard>::Type 	TTIter;	
		typedef TriplexPotential<TId>							TPotential;
		typedef typename ::std::vector<TPotential>				TPotentials;
		
		options.logFileHandle << _getTimeStamp() << " * Starting search for targets in parallel mode" << ::std::endl;
		
		TTriplexSet duplexSet;
		CharString	id;
		unsigned seqNoWithinFile = 0;

		for(; !_streamEOF(file); ++seqNo,++seqNoWithinFile){
			TDuplex	duplexString;
			readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
			appendValue(duplexNames, id, Generous());
			read(file, duplexString, Fasta());		// read Fasta sequence
			ttsnoToFileMap.insert(::std::make_pair<unsigned,::std::pair< ::std::string,unsigned> >(seqNo,::std::make_pair< ::std::string,unsigned>(filename,seqNoWithinFile)));
			appendValue(duplexSet, duplexString);	
			
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Finished reading sequence " << id << ".\n";
		}
		options.logFileHandle << _getTimeStamp() << " * Finished reading " << length(duplexSet) << " duplex sequences " << ::std::endl;
		
		::std::vector<TTargetSet> tmpTtsSets;
		resize(tmpTtsSets, length(duplexSet));
		TPotentials potentials;
		resize(potentials, length(duplexSet));
		
		options.logFileHandle << _getTimeStamp() << " * Detecting triplex targets in all duplexes (" << options.processors << " threads)" << ::std::endl;
		
		if (options._debugLevel >= 1 )
			options.logFileHandle << _getTimeStamp() << "   ... Started TTS detection " << ::std::endl;

		bool reduceSet = false || options.mergeFeatures; //merge overlapping features on request
		
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel)
		{
			SEQAN_PRAGMA_IF_PARALLEL(omp for schedule(dynamic) )
			for (int duplexSeqNo=0;duplexSeqNo<(int)length(duplexSet);++duplexSeqNo) {
				
				// find low complexity regions and mask sequences if requested
				if (options.filterRepeats){
					TRepeatString	data_repeats;
					findRepeats(data_repeats, value(duplexSet, duplexSeqNo), options.minRepeatLength, options.maxRepeatPeriod);
					for (TRepeatIterator rbeg = begin(data_repeats); rbeg != end(data_repeats); ++rbeg){
						TRepeat repeat = *rbeg;
						CharString replacement = string(repeat.endPosition-repeat.beginPosition, 'N' );
						replace(value(duplexSet, duplexSeqNo), repeat.beginPosition, repeat.endPosition, replacement);
					}
					if (options._debugLevel > 1 )
						options.logFileHandle << _getTimeStamp() << "   ... Finished filtering sequence " << duplexSeqNo << " for low complexity regions" << ::std::endl;

				}
				
				TPotential potential(duplexSeqNo);
				if (options.forward){
					unsigned totalNumberOfMatches = processDuplex(tmpTtsSets[duplexSeqNo], value(duplexSet, duplexSeqNo), duplexSeqNo, true, reduceSet, options);
					addCount(potential, totalNumberOfMatches, '+');
					if (options._debugLevel > 1 )
						options.logFileHandle << _getTimeStamp() << "   ... Finished processing Watson strand of sequence " << duplexSeqNo << ::std::endl;			
				}
				
				if (options.reverse){
					unsigned totalNumberOfMatches = processDuplex(tmpTtsSets[duplexSeqNo], value(duplexSet, duplexSeqNo), duplexSeqNo, false, reduceSet, options);
					addCount(potential, totalNumberOfMatches, '-');
					if (options._debugLevel > 1 )
						options.logFileHandle << _getTimeStamp() << "   ... Finished processing Crick strand of sequence " << duplexSeqNo << ::std::endl;
				}
				// get norm for sequence
				setNorm(potential, length(value(duplexSet, duplexSeqNo)), options);
				// add potential to list
				potentials[duplexSeqNo] = potential;
			}
		}
		if (options._debugLevel >= 1 )
			options.logFileHandle << _getTimeStamp() << "   ... Finished TTS detection " << ::std::endl;


		if (options._debugLevel >= 1 )
			options.logFileHandle << _getTimeStamp() << "   ... Started accumulating TTSs in common datastructure " << ::std::endl;

		TTargetSet ttsSet;
		for (int duplexSeqNo=0;duplexSeqNo<(int)length(duplexSet);++duplexSeqNo){
			// copy values
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Add  " << length(tmpTtsSets[duplexSeqNo]) << " to " << length(ttsSet) ;
			
			for (TTIter itb = begin(tmpTtsSets[duplexSeqNo]); itb != end(tmpTtsSets[duplexSeqNo]); ++itb){
				appendValue(ttsSet, *itb);
			}
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << " to finally  " << length(ttsSet) << ::std::endl;
				
		}

#ifdef TRIPLEX_DEBUG
		typedef typename Iterator<TTargetSet>::Type  TIterMotifSet;
		::std::cout << "printing all tts segments" << ::std::endl;
		for (TIterMotifSet itr=begin(ttsSet); itr != end(ttsSet);++itr){
			::std::cout << "tts: " << ttsString(*itr) << " type: " << (*itr).motif << " length: "<< length(*itr) <<  " position: "<< beginPosition(*itr) << " " << ::std::endl;
		}
#endif
		
		if (options._debugLevel >= 1 )
			options.logFileHandle << _getTimeStamp() << "   ... Finished accumulating TTSs in common datastructure " << ::std::endl;
		
		// any business with duplicate detection?
		if (options.detectDuplicates != DETECT_DUPLICATES_OFF){
			if (options._debugLevel >= 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Started duplicate detection " << ::std::endl;

			if (options.detectDuplicates == DETECT_DUPLICATES_STRICT)
				_countDuplicatesStrict(ttsSet, options, TTS() );
			else
				_countDuplicatesPermissive(ttsSet, options, TTS() );
			
			// remove duplicates with respect to cutoff
			if (options.duplicatesCutoff >= 0){
				unsigned removed = _filterDuplicatesWithCutoff(ttsSet, options);
				options.logFileHandle << _getTimeStamp() << " * Duplicate filtering removed " << ::std::setprecision(3)  << removed << " entries (" << length(ttsSet) << " remain)" << ::std::endl;
			}
			if (options._debugLevel >= 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Finished duplicate detection " << ::std::endl;

		}

		options.logFileHandle << _getTimeStamp() << " * Finished detecting targets "  << ::std::endl;

		dumpTtsMatches(outputhandle, ttsSet, duplexNames, options);	
		dumpSummary(potentials, duplexNames, options, TTS());	
		
		file.close();	
	}

	template <
	typename TDuplexName, 
	typename TInput, 
	typename TSeqNo, 
	typename TMap, 
	typename TString, 
	typename TOutput
	>
	inline void _investigateTTS(TDuplexName filename, 
								TInput &file, 
								TSeqNo &seqNo, 
								TMap &ttsnoToFileMap,
								TString &duplexNames, 
								TOutput &outputhandle, 
								Options &options)
	{
		if (options.detectDuplicates == DETECT_DUPLICATES_OFF && options.runtimeMode != RUN_PARALLEL_DUPLEX ){
			_investigateTTSconsecutively(filename, file, seqNo, ttsnoToFileMap, duplexNames, outputhandle, options);
		} else {
			_investigateTTSsimultaneous(filename, file, seqNo, ttsnoToFileMap, duplexNames, outputhandle, options);
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Main TTS mapper function
	template <typename TTargetSet>
	int investigateTTS(Options &options)
	{
		MultiFasta				duplexSet;
		StringSet<CharString>	duplexNames;		// duplex names, taken from the Fasta file
		
		map<unsigned,pair< string,unsigned> > ttsnoToFileMap;
		
		// circumvent numerical obstacles
	//	options.errorRate += 0.0000001;
		
		//////////////////////////////////////////////////////////////////////////////
		// Step 1: search for TTS
		
		options.logFileHandle << _getTimeStamp() << " * Started searching for triplex target sites " << ::std::endl;
		unsigned seqNo = 0;
		
		// open duplex file
		::std::ifstream file;
		options.logFileHandle << _getTimeStamp() << " * Processing " << options.duplexFileNames[0] << ::std::endl;
		file.open(toCString(options.duplexFileNames[0]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return TRIPLEX_DUPLEXREAD_FAILED;
		
		// remove the directory prefix of current duplex file
		::std::string ttsFile(toCString(options.duplexFileNames[0]));
		size_t lastPos = ttsFile.find_last_of('/') + 1;
		if (lastPos == ttsFile.npos) lastPos = ttsFile.find_last_of('\\') + 1;
		if (lastPos == ttsFile.npos) lastPos = 0;
		::std::string duplexFileName = ttsFile.substr(lastPos);
		
		// iterate over duplex sequences
		SEQAN_PROTIMESTART(find_time);
		
		// create output file
	#ifdef BOOST
		if (options.compressOutput){
			io::filtering_ostream filterstream;
			filterstream.push(io::gzip_compressor());
			::std::ofstream filehandle;
			if (!empty(options.output) && options.outputFormat!=2){
				openOutputFile(filterstream, options);
			} else {
				filterstream.push(::std::cout);	
			}
			printTTSHeader(filterstream, options);
			_investigateTTS(duplexFileName, file, seqNo, ttsnoToFileMap, duplexNames, filterstream, options);
			closeOutputFile(filterstream, options);
		} else {
	#endif
			::std::ofstream filehandle;
			if (!empty(options.output) && options.outputFormat!=2){
				openOutputFile(filehandle, options);
				printTTSHeader(filehandle, options);
				_investigateTTS(duplexFileName, file, seqNo, ttsnoToFileMap, duplexNames, filehandle, options);
				closeOutputFile(filehandle, options);
			} else {
				printTTSHeader(::std::cout, options);
				_investigateTTS(duplexFileName, file, seqNo, ttsnoToFileMap, duplexNames, ::std::cout, options);
			}
	#ifdef BOOST
		}
	#endif
				
		CharString sfName;
		_getShortFilename(sfName, options.duplexFileNames[0]);

		options.timeFindTtss += SEQAN_PROTIMEDIFF(find_time);

		options.logFileHandle << _getTimeStamp() << " * Finished processing " << options.duplexFileNames[0] << ::std::endl;
			
		options.logFileHandle << _getTimeStamp() << " * Finished searching for triplex target sites within " << ::std::setprecision(3)  << options.timeFindTtss << " seconds." << ::std::endl;
		return TRIPLEX_NORMAL_PROGAM_EXIT;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Main TFO mapper function
	template <typename TOligoSet, typename TMotifSet>
	int investigateTFO(Options &options
	){
		typedef typename Iterator<TOligoSet, Standard>::Type 	TOligoIter;
		typedef typename Iterator<TMotifSet, Standard>::Type	Titer;
		typedef Repeat<unsigned, unsigned>						TRepeat;
		typedef String<TRepeat>									TRepeatString;
		typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
		typedef TriplexPotential<TId>							TPotential;
		typedef typename ::std::list<TPotential>				TPotentials;
		
		TOligoSet				oligoSequences;
		StringSet<CharString>	oligoNames;				// tfs names, taken from the Fasta file
		
		// circumvent numerical obstacles
	//	options.errorRate += 0.0000001;
		
		//////////////////////////////////////////////////////////////////////////////
		// Step 1: read in TFO files
		
		options.logFileHandle << _getTimeStamp() << " * Started reading single-stranded file:" << options.tfoFileNames[0] << ::std::endl;
		
		if (!_loadOligos(oligoSequences, oligoNames, toCString(options.tfoFileNames[0]), options)) {
			options.logFileHandle << "ERROR: Failed to load single-stranded sequence set" << ::std::endl;
			cerr << "Failed to load single-stranded sequence set" << endl;
			return TRIPLEX_TFOREAD_FAILED;
		}
		
		options.logFileHandle << _getTimeStamp() << " * Finished reading single-stranded file (" << length(oligoSequences) << " sequences read)" << ::std::endl;
		
		
		//////////////////////////////////////////////////////////////////////////////
		// Step 2:  process all sequences with the requested TFO motifs
		
		SEQAN_PROTIMESTART(find_time);
		options.logFileHandle << _getTimeStamp() << " * Started detecting triplex-forming oligonucleotides (TFOs) in single-stranded sequences" << ::std::endl;
		
		unsigned oligoSeqNo = 0;
		TMotifSet tfoMotifSet;
		
		bool reduceSet = false || options.mergeFeatures; //merge overlapping features on request
		TPotentials potentials;
		
		for (TOligoIter it = begin(oligoSequences, Standard()); it != end(oligoSequences, Standard()); ++ it){
			// find low complexity regions and mask sequences if requested
			if (options.filterRepeats){
				TRepeatString	data_repeats;
				findRepeats(data_repeats, *it, options.minRepeatLength, options.maxRepeatPeriod);
				for (TRepeatIterator rbeg = begin(data_repeats); rbeg != end(data_repeats); ++rbeg){
					TRepeat repeat = *rbeg;
					CharString replacement = string(repeat.endPosition-repeat.beginPosition, 'N' );
					replace(*it, repeat.beginPosition, repeat.endPosition, replacement);
				}
			}
			
			TPotential potential(oligoSeqNo);
			// process TC motif
			if (options.motifTC) {
				unsigned totalNumberOfMatches = processTCMotif(tfoMotifSet, *it, oligoSeqNo, reduceSet, options);
				addCount(potential, totalNumberOfMatches, 'Y');
			}
			// process GA motif
			if (options.motifGA) {
				unsigned totalNumberOfMatches = processGAMotif(tfoMotifSet, *it, oligoSeqNo, reduceSet, options);
				addCount(potential, totalNumberOfMatches, 'R');
			}
			// process GT motif
			if (options.motifGT_p || options.motifGT_a) {
				// for TFO search the parallel GT motifs need to be recorded only
				unsigned totalNumberOfMatches = processGTMotif(tfoMotifSet, *it, oligoSeqNo, TRIPLEX_ORIENTATION_PARALLEL, reduceSet, options);
				addCount(potential, totalNumberOfMatches, 'M');
			}
			// get norm for sequence
			setNorm(potential, length(*it), options);
			// add potential to list
			appendValue(potentials, potential);
			// increase oligo id counter
			++oligoSeqNo;
		}
		
		// any business with duplicates?
		if (options.detectDuplicates != DETECT_DUPLICATES_OFF){
			
			// detect duplicates if requested
			if (options.detectDuplicates == DETECT_DUPLICATES_STRICT)
				_countDuplicatesStrict(tfoMotifSet, options, TFO() );
			else
				_countDuplicatesPermissive(tfoMotifSet, options, TFO() );
			
			if (options.duplicatesCutoff >= 0){
				unsigned removed = _filterDuplicatesWithCutoff(tfoMotifSet, options);
				options.logFileHandle << _getTimeStamp() << " * Duplicate filtering removed " << ::std::setprecision(3)  << removed << " entries (" << length(tfoMotifSet) << " remain)" << ::std::endl;
			}
		}
		
		options.timeFindTfos += SEQAN_PROTIMEDIFF(find_time);	

		options.logFileHandle << _getTimeStamp() << " * Finished detecting TFOs within " << ::std::setprecision(3)  << options.timeFindTfos << " seconds (" << length(tfoMotifSet) << " TFOs detected)" << ::std::endl;
		
		//////////////////////////////////////////////////////////////////////////////
		// Step 4: output
		
		options.logFileHandle << _getTimeStamp() << " * Started printing results " << ::std::endl;
		// create output file
		
		
	#ifdef BOOST
		if (options.compressOutput){
			io::filtering_ostream filterstream;
			filterstream.push(io::gzip_compressor());
			::std::ofstream filehandle;
			if (!empty(options.output) && options.outputFormat!=2){
				openOutputFile(filterstream, options);
			} else {
				filterstream.push(::std::cout);	
			}
			printTFOHeader(filterstream, options);
			dumpTfoMatches(filterstream, tfoMotifSet, oligoNames, options);
			dumpSummary(potentials, oligoNames, options, TFO());
			closeOutputFile(filterstream, options);
		} else {
	#endif
			::std::ofstream filehandle;
			if (!empty(options.output) && options.outputFormat!=2){
				openOutputFile(filehandle, options);
				printTFOHeader(filehandle, options);
				dumpTfoMatches(filehandle, tfoMotifSet, oligoNames, options);
				dumpSummary(potentials, oligoNames, options, TFO());
				closeOutputFile(filehandle, options);
			} else {
				printTFOHeader(::std::cout, options);
				dumpTfoMatches(::std::cout, tfoMotifSet, oligoNames, options);
				dumpSummary(potentials, oligoNames, options, TFO());
			}
	#ifdef BOOST
		}
	#endif
		
		options.logFileHandle << _getTimeStamp() << " * Finished outputing results " << ::std::endl;
		return TRIPLEX_NORMAL_PROGAM_EXIT;
	}

	// starting the main program
	int _mainWithOptions(int argc, char const ** argv, Options &options)
	{
		SEQAN_PROTIMESTART(runtime);
		
		openLogFile(options);
		openSummaryFile(options);
		
		_populateLogFile(argc, argv, options);
		
		int result = TRIPLEX_NORMAL_PROGAM_EXIT;
		if (options.runmode == TRIPLEX_TTS_SEARCH){ // investigate TTS only
			result = investigateTTS<TTargetSet>(options);
		} else if (options.runmode == TRIPLEX_TFO_SEARCH){ // investigate TFO only
			result = investigateTFO<TTriplexSet, TMotifSet>(options);
		} else if (options.runmode == TRIPLEX_TRIPLEX_SEARCH){ // map TFO and TTSs
			result = mapTriplexes<TTriplexSet, TMotifSet>(options);
		} else {
			cerr << "Exiting ... invalid runmode" << endl;
			options.logFileHandle << "ERROR: Exit due to invalid options " << ::std::endl;	
			closeLogFile(options);
			return TRIPLEX_INVALID_OPTIONS;
		}
		
		if (result != TRIPLEX_NORMAL_PROGAM_EXIT){
			cerr << "Exiting ... invalid result" << endl;
			options.logFileHandle << "ERROR: Exit with errors " << ::std::endl;	
		} else {
			options.logFileHandle << _getTimeStamp() << " * Exit without errors " << ::std::endl;
		}
		
		closeSummaryFile(options);
		
		options.logFileHandle << _getTimeStamp() << " * Finished program within " <<  ::std::setprecision(3)  << SEQAN_PROTIMEDIFF(runtime) << " seconds" << ::std::endl;
		closeLogFile(options);
		
		return result;
	}
} // namespace SEQAN_NAMESPACE_MAIN 

//////////////////////////////////////////////////////////////////////////////
// Command line parsing and parameter choosing
// Program entry point
int main(int argc, char const ** argv)
{
    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    _setupCommandLineParser(parser, options);
    
    // Then, parser the command line and handle the cases where help display
    // is requested or errornoeous parameters were given.
    int ret = _parseCommandLineAndCheck(options, parser, argc, argv);
    if (ret != 0)
        return ret;
        
    if (options.showHelp || options.showVersion)
        return 0;
    
    // Finally, launch the program.
    ret = _mainWithOptions(argc, argv, options);
    return ret;
}



