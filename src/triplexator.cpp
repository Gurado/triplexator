// ==========================================================================
//                                triplexator
// ==========================================================================
// Copyright (c) 2011, Fabian Buske, UQ
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

#ifndef SEQAN_ENABLE_PARALLELISM
#define SEQAN_ENABLE_PARALLELISM 1		// disable parallelism on default
#endif

#include <seqan/platform.h>
#ifdef PLATFORM_WINDOWS
	#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
	#define SEQAN_DEFAULT_TMPDIR "./"
#endif

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/parallel.h> 
#include <seqan/parallel/parallel_macros.h> 
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include "triplexator.h"
#include "triplex_alphabet.h"
#include "gardener.h"
#include "output.h"

#include <iostream>
#include <sstream>

using namespace std;
using namespace seqan;


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
	options.logFileHandle << "- report duplicate locations : " << (options.reportDuplicateLocations?"Yes":"No") << ::std::endl;
#ifdef BOOST
	options.logFileHandle << "- compress output : " << (options.compressOutput?"Yes":"No") << ::std::endl;
#endif	
	options.logFileHandle << "*************************************************************" << ::std::endl;
	options.logFileHandle << "*** Main Options:" << ::std::endl;
//	options.logFileHandle << "- consider forward strand in duplex : " << (options.forward?"Yes":"No") << ::std::endl;
//	options.logFileHandle << "- consider reverse strand in duplex : " << (options.reverse?"Yes":"No") << ::std::endl;
	options.logFileHandle << "- maximum error-rate : " << (options.errorRate*100) << "%" << ::std::endl;
	options.logFileHandle << "- maximum total error : " << (options.maximalError) << ::std::endl;	
	options.logFileHandle << "- minimum guanine-ratio with respect to the target : " << (options.guanineRate*100) << "%" << ::std::endl;
	options.logFileHandle << "- relax guanine-ratio at flanks : " << (options.relaxGuanineRate?"Yes":"No") << ::std::endl;
	
	options.logFileHandle << "- minimum length : " << options.minLength << " nucleotides" << ::std::endl;
//	if (options.maxLength<options.minLength)
//		options.logFileHandle << "- maximum length : omitted" << ::std::endl;
//	else 
//		options.logFileHandle << "- maximum length : " << options.maxLength << " nucleotides" << ::std::endl;
	

	options.logFileHandle << "- maximum number of tolerated consecutive pyrimidine interruptions in a target: " << options.maxInterruptions << ::std::endl;
	if (options.runmode == TRIPLEX_TRIPLEX_SEARCH || options.runmode == TRIPLEX_TFO_SEARCH){
		options.logFileHandle << "- include GT-motif : " << (options.motifGT_a || options.motifGT_p?"Yes":"No") << ::std::endl;
		options.logFileHandle << "- include GA-motif : " << (options.motifGA?"Yes":"No") << ::std::endl;
		options.logFileHandle << "- include TC-motif : " << (options.motifTC?"Yes":"No") << ::std::endl;
	}
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
	options.logFileHandle << "- minimum repeat length : " << options.minRepeatLength << ::std::endl;
	options.logFileHandle << "- maximum repeat period : " << options.maxRepeatPeriod << ::std::endl;
	options.logFileHandle << "- duplicate cutoff : " << options.duplicatesCutoff << ::std::endl;
	if (options.runmode == TRIPLEX_TRIPLEX_SEARCH){
		if (options.filterMode == FILTERING_GRAMS){
			options.logFileHandle << "- filtering : qgrams" << ::std::endl;
			options.logFileHandle << "- weight(qgram) : " << length(options.shape) << ::std::endl;
			options.logFileHandle << "- threshold(qgram) : " << options.qgramThreshold << ::std::endl;
		} else {
			options.logFileHandle << "- filtering : none - greedy algorithm" << ::std::endl;
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
int _findTriplex(TMotifSet				&tfoMotifSet,
				StringSet<CharString>	&tfoNames,
				TFile					&outputfile,
				Options					&options,
				TShape const			&shape)
{
	typedef Index<TMotifSet, IndexQGram<TShape, OpenAddressing> >			TQGramIndex;
	typedef __int64															TId;
	typedef Gardener<TId, GardenerUngapped>									TGardener;
	
	SEQAN_PROTIMESTART(find_time);
	options.logFileHandle << _getTimeStamp() << " * Started searching for triplexes" << ::std::endl;
	
	// create index
	if (options._debugLevel >= 1)
		options.logFileHandle << _getTimeStamp() <<  " - Started creating q-gram index for all TFOs" << ::std::endl;
	
	TQGramIndex index_qgram(tfoMotifSet);
	resize(indexShape(index_qgram), weight(shape));
	// create pattern	
	Pattern<TQGramIndex, QGramsLookup< TShape, Standard_QGramsLookup > > pattern(index_qgram,shape);
	options.timeFindTriplexes = 0;
	
	// create index
	if (options._debugLevel >= 1)
		options.logFileHandle << _getTimeStamp() <<  " - Finised creating q-gram index for all TFOs" << ::std::endl;
	
	TId duplexSeqNo = 0;
	// open duplex file
	options.logFileHandle << _getTimeStamp() << " * Processing " << options.duplexFileNames[0] << ::std::endl;
#if SEQAN_ENABLE_PARALLELISM	
	// run in parallel if requested
	if (options.runtimeMode==RUN_PARALLEL_DUPLEX){
		if (options.filterMode == FILTERING_GRAMS){
			startTriplexSearchParallelDuplex(tfoMotifSet, tfoNames, pattern, outputfile, duplexSeqNo, options, TGardener());
		} else {
			startTriplexSearchParallelDuplex(tfoMotifSet, tfoNames, pattern, outputfile, duplexSeqNo, options, BruteForce());			
		}
	} else {
	// otherwise go for serial processing
#endif	
		if (options.filterMode == FILTERING_GRAMS){
			startTriplexSearchSerial(tfoMotifSet, tfoNames, pattern, outputfile, duplexSeqNo, options, TGardener());
		} else {
			startTriplexSearchSerial(tfoMotifSet, tfoNames, pattern, outputfile, duplexSeqNo, options, BruteForce());
		}	
#if SEQAN_ENABLE_PARALLELISM	
	}
#endif
	options.logFileHandle << _getTimeStamp() << " * Finished processing " << options.duplexFileNames[0] << ::std::endl; 

	options.timeFindTriplexes += SEQAN_PROTIMEDIFF(find_time);	

	options.logFileHandle << _getTimeStamp() << " * Finished searching for triplexes  within " << ::std::setprecision(3) << options.timeFindTriplexes << " seconds (summed over all cpus)" << ::std::endl;
	
	return 0;
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
	typedef std::list<TRepeat>								TRepeatString; //@TODO workaround for memory leak in seqan string
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
			TRepeatString	data_repeats;
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
		_findTriplex(tfoMotifSet, oligoNames, filterstream, options, ungappedShape);
		closeOutputFile(filterstream, options);
	} else {
#endif
		::std::ofstream filehandle;
		if (!empty(options.output) && options.outputFormat!=2){
			openOutputFile(filehandle, options);
			printTriplexHeader(filehandle, options);
			_findTriplex(tfoMotifSet, oligoNames, filehandle, options, ungappedShape);
			closeOutputFile(filehandle, options);
		} else {
			printTriplexHeader(::std::cout, options);
			_findTriplex(tfoMotifSet, oligoNames, ::std::cout, options, ungappedShape);
		}
#ifdef BOOST
	}
#endif
	
	return TRIPLEX_NORMAL_PROGAM_EXIT;
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
	typedef std::list<TRepeat>								TRepeatString; //@TODO workaround for memory leak in seqan string
	typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;

	options.logFileHandle << _getTimeStamp() << " * Starting search for targets in serial mode " << ::std::endl;
	
	TDuplex	duplexString;
	CharString	id;
	unsigned seqNoWithinFile = 0;
	for(; !_streamEOF(file); ++seqNo,++seqNoWithinFile){
		
		if (options._debugLevel > 1 )
			options.logFileHandle << _getTimeStamp() << "   ... Started reading next duplex sequence " << ::std::endl;

		TTargetSet ttsSet;
		readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
		appendValue(duplexNames, id, Generous());
		
		read(file, duplexString, Fasta());			// read Fasta sequence
		ttsnoToFileMap.insert(::std::make_pair<unsigned,::std::pair< ::std::string,unsigned> >(seqNo,::std::make_pair< ::std::string,unsigned>(filename,seqNoWithinFile)));
		
		if (options._debugLevel > 1 )
			options.logFileHandle << _getTimeStamp() << "   ... Finished reading next duplex sequence" << ::std::endl;

		bool reduceSet = false; // done merge overlapping features
		
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

			processDuplex(ttsSet, duplexString, seqNoWithinFile, true, reduceSet, options);
			
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Finished processing the Watson strand" << ::std::endl;

		}
		
		if (options.reverse){
			
			if (options._debugLevel > 1 )
				options.logFileHandle << _getTimeStamp() << "   ... Started processing the Crick strand" << ::std::endl;

			processDuplex(ttsSet, duplexString, seqNoWithinFile, false, reduceSet, options);

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

		dumpTtsMatches(outputhandle, ttsSet, duplexNames, seqNoWithinFile, options);				
		
		if (options._debugLevel > 1 )
			options.logFileHandle << _getTimeStamp() << "   ... Finished outputing results " << ::std::endl;

		if (options._debugLevel > 1 )
			options.logFileHandle << _getTimeStamp() << "   ... Finished processing duplex " << id << ::std::endl;

	}
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
	typedef std::list<TRepeat>								TRepeatString; //@TODO workaround for memory leak in seqan string
	typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
	typedef typename Iterator<TTriplexSet, Standard>::Type 	TIter;
	typedef typename Iterator<TTargetSet, Standard>::Type 	TTIter;	
	
	
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
	
	options.logFileHandle << _getTimeStamp() << " * Detecting triplex targets in all duplexes (" << options.processors << " threads)" << ::std::endl;
	
	if (options._debugLevel >= 1 )
		options.logFileHandle << _getTimeStamp() << "   ... Started TTS detection " << ::std::endl;

	bool reduceSet = false; // don't merge overlapping features
	
	SEQAN_OMP_PRAGMA(omp parallel)
	{
		SEQAN_OMP_PRAGMA(omp for schedule(dynamic) )
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
			
			if (options.forward){
				processDuplex(tmpTtsSets[duplexSeqNo], value(duplexSet, duplexSeqNo), duplexSeqNo, true, reduceSet, options);
				if (options._debugLevel > 1 )
					options.logFileHandle << _getTimeStamp() << "   ... Finished processing Watson strand of sequence " << duplexSeqNo << ::std::endl;			
			}
			
			if (options.reverse){
				processDuplex(tmpTtsSets[duplexSeqNo], value(duplexSet, duplexSeqNo), duplexSeqNo, false, reduceSet, options);
				if (options._debugLevel > 1 )
					options.logFileHandle << _getTimeStamp() << "   ... Finished processing Crick strand of sequence " << duplexSeqNo << ::std::endl;
			}
			
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
	
	// any business with duplicates?
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
		
	dumpTtsMatches(outputhandle, ttsSet, duplexNames, seqNoWithinFile, options);	
	
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
	if (options.detectDuplicates == DETECT_DUPLICATES_OFF){
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
	typedef std::list<TRepeat>								TRepeatString; //@TODO workaround for memory leak in seqan string
	typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
	
	
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
	
	bool reduceSet = false; // done merge overlapping features
	
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
		
		// process TC motif
		if (options.motifTC) {
			processTCMotif(tfoMotifSet, *it, oligoSeqNo, reduceSet, options);
		}
		// process GA motif
		if (options.motifGA) {
			processGAMotif(tfoMotifSet, *it, oligoSeqNo, reduceSet, options);
		}
		// process GT motif
		if (options.motifGT_p || options.motifGT_a) {
			// for TFO search the parallel GT motifs need to be recorded only
			processGTMotif(tfoMotifSet, *it, oligoSeqNo, TRIPLEX_ORIENTATION_PARALLEL, reduceSet, options);
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
		closeOutputFile(filterstream, options);
	} else {
#endif
		::std::ofstream filehandle;
		if (!empty(options.output) && options.outputFormat!=2){
			openOutputFile(filehandle, options);
			printTFOHeader(filehandle, options);
			dumpTfoMatches(filehandle, tfoMotifSet, oligoNames, options);
			closeOutputFile(filehandle, options);
		} else {
			printTFOHeader(::std::cout, options);
			dumpTfoMatches(::std::cout, tfoMotifSet, oligoNames, options);
		}
#ifdef BOOST
	}
#endif
	
	options.logFileHandle << _getTimeStamp() << " * Finished outputing results " << ::std::endl;
	return TRIPLEX_NORMAL_PROGAM_EXIT;
}

// starting the main program
int mainWithOptions(int argc, char const ** argv, Options &options)
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

//////////////////////////////////////////////////////////////////////////////
// Command line parsing and parameter choosing
// Program entry point
int main(int argc, char const ** argv)
{
    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parser the command line and handle the cases where help display
    // is requested or errornoeous parameters were given.
    int ret = parseCommandLineAndCheck(options, parser, argc, argv);
    if (ret != 0)
        return ret;
        
    if (options.showHelp || options.showVersion)
        return 0;
    
    // Finally, launch the program.
    ret = mainWithOptions(argc, argv, options);
    return ret;
}



