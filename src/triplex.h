// ==========================================================================
//                                triplexator
// ==========================================================================
// Copyright (c) 2011, Fabian Buske, UQTRIPLEXATOR_TRIPLEX_H_
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

#ifndef SANDBOX_FBUSKE_APPS_TRIPLEXATOR_TRIPLEX_H_
#define SANDBOX_FBUSKE_APPS_TRIPLEXATOR_TRIPLEX_H_

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/parallel.h> 
#include <seqan/index.h>
#include <seqan/modifier/modifier_view.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.
#include <seqan/find.h>
#include <seqan/score.h>

#include "triplex_alphabet.h"
#include "triplex_pattern.h"
#include "gardener.h"

#ifndef SEQAN_PRAGMA_IF_PARALLEL
#if SEQAN_ENABLE_PARALLELISM
#define STRINGIFY(a) #a
#define SEQAN_PRAGMA_IF_PARALLEL(code) \
_Pragma(STRINGIFY(code))
#else // SEQAN_ENABLE_PARALLELISM
#define SEQAN_PRAGMA_IF_PARALLEL(code)
#endif // SEQAN_ENABLE_PARALLELISM
#endif // SEQAN_PRAGMA_IF_PARALLEL

using namespace std;
namespace SEQAN_NAMESPACE_MAIN
{


// ============================================================================
// Tags, Classes, Enums
// ============================================================================
	
	enum MODE_OF_OPERATION
	{
		TRIPLEX_TTS_SEARCH		= 1,
		TRIPLEX_TFO_SEARCH		= 2,
		TRIPLEX_TRIPLEX_SEARCH	= 3
	};
	
	enum RUNTIME_MODE
	{
		RUN_SERIAL				= 0,
		RUN_PARALLEL_TRIPLEX	= 1,
		RUN_PARALLEL_DUPLEX  	= 2,
		RUN_PARALLEL_STRANDS	= 3

	};
	
	enum FILTER_MODE
	{
		BRUTE_FORCE				= 0,
		FILTERING_GRAMS			= 1
	};
	
	enum DETECT_DUPLICATES
	{
		DETECT_DUPLICATES_OFF			= 0,
		DETECT_DUPLICATES_PERMISSIVE	= 1,
		DETECT_DUPLICATES_STRICT		= 2
	};
	
	enum ORIENTATION
	{
		TRIPLEX_ORIENTATION_PARALLEL	 	= 1,
		TRIPLEX_ORIENTATION_ANTIPARALLEL  	= -1,
		TRIPLEX_ORIENTATION_BOTH			= 0
	};
	
	enum OUTPUT_FORMAT
	{
		FORMAT_BED	 	= 0,
		FORMAT_TRIPLEX  = 1,
		FORMAT_SUMMARY	= 2
	};

	enum TRIPLEX_ERROR
	{
		TRIPLEX_NORMAL_PROGAM_EXIT  =  0,
		TRIPLEX_INVALID_OPTIONS 	= -1,
		TRIPLEX_TFOREAD_FAILED  	= -2,
		TRIPLEX_READFILE_FAILED 	= -3,
		TRIPLEX_INVALID_SHAPE   	= -4,
		TRIPLEX_DUPLEXREAD_FAILED 	= -5,
		TRIPLEX_SHAPE_FAILED    	= -6,
		TRIPLEX_OUTPUTFILE_FAILED  	= -7
	};
	
	enum ERROR_REFERENCE
	{
		WATSON_STAND				= 0,
		PURINE_STRAND				= 1,
		THIRD_STRAND				= 2
	};
	
	
	struct _TTS;
	typedef Tag<_TTS> TTS; // tag for a triplex target site in the duplex
	
	struct _TFO;
	typedef Tag<_TFO> TFO; // tag for a triplex forming oligonucleotide
	
	struct _MIXEDMOTIF;
	typedef Tag<_MIXEDMOTIF> MIXEDMOTIF; // tag for a triplex forming oligonucleotide of the GT Motif

	struct _PURINEMOTIF;
	typedef Tag<_PURINEMOTIF> PURINEMOTIF; // tag for a triplex forming oligonucleotide of the GA Motif

	struct _PYRIMIDINEMOTIF;
	typedef Tag<_PYRIMIDINEMOTIF> PYRIMIDINEMOTIF; // tag for a triplex forming oligonucleotide of the TC Motif
	
	struct _BruteForce;
	typedef Tag<_BruteForce> BruteForce; // tag for a brute force approach

	
	// definition of a queue block entry
	struct QueueBlock
	{
		unsigned int	blockSize;		// size of entry
		char			blockClass;		// class of entry
	};

	// definition of a triplex match
	template <typename _TGPos, typename TSize, typename TScore>
	struct TriplexMatch
	{
		typedef typename MakeSigned_<_TGPos>::Type TGPos;
		
		TSize		tfoSeqNo;		// tfo seqNo
		TGPos		oBegin;			// begin position of the match in the third strand (oligo)
		TGPos		oEnd;			// end position of the match in the third strand (oligo)
		unsigned	ttsSeqNo;		// tts seqNo
		TGPos		dBegin;			// begin position of the match in the duplex
		TGPos		dEnd;			// end position of the match in the duplex
		TScore		mScore;			// matching score
		bool		parallel;		// false = anti-parallel, true = parallel
		char		motif;			// triplex motif [R,Y,M]
		char		strand;			// '+'..Watson, '-'..Crick strand
		
		TriplexMatch() {
			tfoSeqNo = -1;
			oBegin = -1;
			oEnd = -1;
			ttsSeqNo = 0;
			dBegin = -1;
			dEnd = -1;
			mScore = 0;
			parallel = false;
			motif = ' ';
			strand = ' ';
			SEQAN_CHECKPOINT
		}
		
		TriplexMatch(TSize _tfoSeqNo, TGPos _oBegin, TGPos _oEnd, TSize _ttsSeqNo, TGPos _dBegin, TGPos _dEnd, TScore _mScore, bool _parallel, char _motif, char _strand):
			tfoSeqNo(_tfoSeqNo),
			oBegin(_oBegin),
			oEnd(_oEnd),
			ttsSeqNo(_ttsSeqNo),
			dBegin(_dBegin),
			dEnd(_dEnd),
			mScore(_mScore),
			parallel(_parallel),
			motif(_motif),
			strand(_strand)
		{
			SEQAN_CHECKPOINT
		}
		
		
		~TriplexMatch(){
			SEQAN_CHECKPOINT
		}
		
		template <typename TSource>
		inline TriplexMatch &
		operator = (TSource const & source) {
			SEQAN_CHECKPOINT
			assign(*this, source);
			return *this;
		}
		
		void print(){
			::std::cerr << "Match:" << tfoSeqNo << " ("<< oBegin << ":" << oEnd << ")\t";
			::std::cerr << ttsSeqNo << " ("<< dBegin << ":" << dEnd << ")\tD:"<< (dBegin-oBegin)<< "\t";
			::std::cerr << motif << " " << strand << " " << (parallel?"P":"A") << " " << mScore << ::std::endl;
		}
	};	
	
	// ... to sort matches according to dBegin
	template <typename TTriplexMatch>
	struct LessRNoGBeginPos : public ::std::binary_function < TTriplexMatch, TTriplexMatch, bool >
	{
		inline bool operator() (TTriplexMatch const &a, TTriplexMatch const &b) const
		{
			// tts position and orientation
			if (a.ttsSeqNo < b.ttsSeqNo) return true;
			if (a.ttsSeqNo > b.ttsSeqNo) return false;
			// tfo seq number
			if (a.tfoSeqNo < b.tfoSeqNo) return true;
			if (a.tfoSeqNo > b.tfoSeqNo) return false;
			// tts begin position
			if (a.dBegin < b.dBegin) return true;
			if (a.dBegin > b.dBegin) return false;
			// tfo begin position
			if (a.oBegin < b.oBegin) return true;
			if (a.oBegin > b.oBegin) return false;
			// tfo end position
			if (a.oEnd < b.oEnd) return true;
			if (a.oEnd > b.oEnd) return false;
			// triplex motif
			if (a.motif < b.motif) return true;
			if (a.motif > b.motif) return false;
			// tts strand
			if (a.strand < b.strand) return true;
			if (a.strand > b.strand) return false;
			// triplex orientation
			if (a.parallel < b.parallel) return true;
			if (a.parallel > b.parallel) return false;
			// score
			return a.mScore > b.mScore;
		}
	};
	
	// definition of sequence location
	template <typename TSeq, typename TPos>
	struct SeqPos
	{
		TSeq		seqnr;		// sequence number
		TPos		position;	// offset position
		
		bool operator==(const SeqPos<TSeq, TPos>& b) const;
		bool operator!=(const SeqPos<TSeq, TPos>& b) const;
		bool operator<(const SeqPos<TSeq, TPos>& b) const;
		bool operator>(const SeqPos<TSeq, TPos>& b) const;
		
		
		SeqPos(TSeq _seqnr, TPos _pos):
		seqnr(_seqnr),
		position(_pos)
		{
			SEQAN_CHECKPOINT
		}
		
		
		SeqPos(){
			SEQAN_CHECKPOINT
		}
		
		template <typename TSource>
		inline SeqPos &
		operator = (TSource const & source) {
			SEQAN_CHECKPOINT
			assign(*this, source);
			return *this;
		}
		
		void print(){
			::std::cerr << "SeqPos:" << seqnr << ":"<< position << ::std::endl;
		}
	};	
	
	//____________________________________________________________________________
	
	template <typename TSeq, typename TPos>
	bool SeqPos<TSeq, TPos>::operator==(const SeqPos<TSeq, TPos>& b) const {
		if (seqnr != b.seqnr) return false;
		if (position != b.position) return false;
		return true;
	}
	
	//____________________________________________________________________________
	
	template <typename TSeq, typename TPos>	
	bool SeqPos<TSeq, TPos>::operator!=(const SeqPos<TSeq, TPos>& b) const {
		return !(*this == b);
	}
	
	//____________________________________________________________________________
	
	template <typename TSeq, typename TPos>
	bool SeqPos<TSeq, TPos>::operator<(const SeqPos<TSeq, TPos>& b) const {
		if (seqnr < b.seqnr) return true;
		if (seqnr > b.seqnr) return false;
		if (position < b.position) return true;
		return false;
	}
	
	//____________________________________________________________________________
	
	template <typename TSeq, typename TPos>
	bool SeqPos<TSeq, TPos>::operator>(const SeqPos<TSeq, TPos>& b) const {
		return b<*this;
	}
	
	//____________________________________________________________________________

	template <typename TSeq, typename TPos>
	inline TSeq getSequenceNo(SeqPos<TSeq, TPos> & me){
		return me.seqnr;
	}
	
	
	template <typename TSeq, typename TPos>
	inline TSeq getSequenceNo(SeqPos<TSeq, TPos> const & me){
		return me.seqnr;
	}
	
	//____________________________________________________________________________

	template <typename TSeq, typename TPos>
	inline TPos getPosition(SeqPos<TSeq, TPos> & me){
		return me.position;
	}
	
	template <typename TSeq, typename TPos>
	inline TPos getPosition(SeqPos<TSeq, TPos> const & me){
		return me.position;
	}
	
	//____________________________________________________________________________
	
	
	typedef Dna5String										Ttts;
	typedef StringSet<Ttts>									TttsSet;
	typedef StringSet<TriplexString>						TOligoMotifSet;
	typedef TriplexString									TTriplex;
	typedef TriplexString									TDuplex;
	typedef StringSet<ModStringTriplex<TDuplex, TDuplex> > 	TTargetSet;
	typedef	DnaRYString										TMaskRY;
	typedef	DnaKMString										TMaskKM;
	typedef StringSet<TMaskRY>								TMaskRYSet;
	typedef StringSet<TMaskKM>								TMaskKMSet;
	typedef StringSet<TTriplex>								TTriplexSet;
	typedef __int64											TId;
	typedef TriplexMatch<Difference<TriplexString>::Type, TId, double>		TMatch;		// a single match
	typedef ::std::set<TMatch, LessRNoGBeginPos<TMatch> > TMatchSet;	// set of matches
	typedef Graph<Automaton<Triplex, Triplex> >				TGraph;		// automaton graph
	typedef VertexDescriptor<TGraph>::Type					TVertexDescriptor;
	typedef ModStringTriplex<TTriplex, TTriplex>			TOligoMotif;
	typedef StringSet<TOligoMotif>							TMotifSet;
	typedef __int64											TPos;
	
	
	template <typename TTriplexMatch>
	struct BestScore : public ::std::binary_function < TTriplexMatch, TTriplexMatch, bool >
	{
		inline bool operator() (TTriplexMatch const &a, TTriplexMatch const &b) const
		{
			// tfo read number
			if (a.tfoSeqNo < b.tfoSeqNo) return true;
			if (a.tfoSeqNo > b.tfoSeqNo) return false;
			
			// quality
			return a.mScore > b.mScore;
		}
	};

	template<typename TPos, typename TSize, typename TScore>
	TPos length(TriplexMatch<TPos, TSize, TScore> m){
		return min(m.oEnd-m.oBegin, m.dEnd-m.dBegin);
	}
	
	template<typename TPos, typename TSize, typename TScore>
	struct Position<TriplexMatch<TPos, TSize, TScore> >
	{
		typedef TPos Type;
	};
	
	template<typename TPos, typename TSize, typename TScore>
	struct Value<TriplexMatch<TPos, TSize, TScore> >
	{
		typedef TScore Type;
	};
	
	template<typename TPos, typename TSize, typename TScore>
	struct Size<TriplexMatch<TPos, TSize, TScore> >
	{
		typedef TSize Type;
	};
	
	template <typename TTriplexSet, typename TShape, typename TSpec>
	struct Cargo< Index<TTriplexSet, IndexQGram<TShape, TSpec> > > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};
	
	template <>
	struct Size<TriplexString>
	{
		typedef unsigned Type;
	};
	
	template <typename TShape>
	struct Size< Index<TTriplexSet, IndexQGram<TShape> > >
	{
		typedef unsigned Type;
	};
	
	
	struct Options
	{
		// main options
		bool 		showHelp;
		bool 		showVersion;
   		unsigned 	hashsize;			// mod used for hashing
		unsigned	runmode;			// mode in which Triplexator is run
		unsigned	filterMode;			// what kind of filtering should be used ... or whether brut force the the mode of choice
		// 0..undetermined
		// 1..TRIPLEX_TTS_SEARCH
		// 2..TRIPLEX_TFO_SEARCH
		// 3..TRIPLEX_TRIPLEX_SEARCH
		bool		ttsFileSupplied;  	// indicates that at least one file as triplex target has been specified
		bool		tfoFileSupplied;  	// indicates that one file as triplex source has been specified
		bool		forward;			// compute forward oriented read matches
		bool		reverse;			// compute reverse oriented read matches
		double		errorRate;			// Criteria 1 threshold (max)
		int			maximalError;		// maximal absolute error tolerated
		double		guanineRate;		// Criteria 2 threshold (min)
		bool		motifTC;			// use triplex TC motifs
		bool		motifGA;			// use triplex AG motifs
		bool		motifGT_p;			// use triplex GT motifs (parallel configuration)
		bool		motifGT_a;			// use triplex GT motifs (anti-parallel configuration)
		int			qgramThreshold;		// the threshold used to calculate the weight of the qgram
		bool		bothTFOStrands;		// search both strands of the sequence for TFOs
		unsigned int minGuanine;		// minimum number of guanines required
		bool		filterRepeats;		// filter repeats 
		unsigned	minRepeatLength;	// minimum length of low complex region to be filtered out
		unsigned	maxRepeatPeriod;	// maximum repeat period defining a low complexity region
		int			duplicatesCutoff;	// threshold above which a feature will not be reported
		unsigned	minBlockRun;		// minimum number of consecutive matches (block) required for a feature
		unsigned	detectDuplicates;	// whether and how to detect duplicates 
		bool		reportDuplicateLocations; // whether to report the locations of duplicates
		bool		sameSequenceDuplicates; // whether to count a feature copy in the same sequence as duplicate or not
		CharString	output;				// name of result file
		CharString	outputFolder;		// name of result folder
		int			_debugLevel;		// level of verbosity
		bool		printVersion;		// print version number
		bool		allMatches;			// output all matches rahter than the longest one only
		CharString	version;			// version info
		int			runtimeMode;		// parallel runtime mode
		int			processors;         // number of processors in parallel mode (threads)
		__int64		minLength;			// minimum length of a triplex (default 14)
		__int64		maxLength;			// maximum length of a triplex (default 50)
		unsigned	maxInterruptions;	// maximum consecutive interruptions
		unsigned	tolError;			// tolerated numbers of error for the minimum length
		
		double		mixed_parallel_max_guanine;		// maximum guanine ratio to consider parallel binding of the GT motif
		double		mixed_antiparallel_min_guanine; // minimum guanine ratio to consider anti-parallel binding of the GT motif
		
		// log options
		CharString	logFileName;
		::std::ofstream logFileHandle;
		CharString	summaryFileName;
		::std::ofstream summaryFileHandle;
		
		// output format options
		bool		computeTpot;
		bool		prettyString;		// indicate matching/mismatching characters with upper/lower case if true
		unsigned	outputFormat;		// 0..FASTA format
		unsigned	errorReference;		// Reference the error is reported to (0=Watson strand of TTS, 1=purine strand of TTS, 2 TFO)
		bool		mergeFeatures;		// combine overlapping features into a feature cluster
		// 1..Triplex format
		// 2..summary only
		const char	*runID;				// runID needed for gff output	
#ifdef BOOST
		bool		compressOutput;
#endif
		
		// filtration parameters
		CharString shape;			// shape (e.g. 111111)
		
		// statistics
		double		timeLoadFiles;		// time for loading input files
		double		timeFindTriplexes;	// time for matching tfos and ttss
		double		timeFindTfos;		// time for finding tfos
		double		timeFindTtss;		// time for finding ttss
		double		timeDumpResults;	// time for dumping the results
		
		// flags
		bool applyMaximumLengthConstraint;
		
		// data
		StringSet<CharString>	tfoFileNames;
		StringSet<CharString>	duplexFileNames;
		
		// 
		TGraph triplexParser;
		
		Options()
		{
			// Set defaults.
			showHelp = true;
			hashsize = 127;
			ttsFileSupplied = false;
			tfoFileSupplied = false;
			runmode = TRIPLEX_TRIPLEX_SEARCH;
			forward = true;
			reverse = true;
			errorRate = 0.05;
			maximalError = -1;
			guanineRate = 0.1;
			motifTC = true;
			motifGA = true;
			motifGT_p = true;
			motifGT_a = true;
			bothTFOStrands = false;
			duplicatesCutoff = -1;
			minBlockRun = 1;
			detectDuplicates = DETECT_DUPLICATES_OFF;
			reportDuplicateLocations = false;
			sameSequenceDuplicates = true;
			qgramThreshold = 2;
			allMatches = false;
			output = "";
			outputFolder = "";
			_debugLevel = 0;
			printVersion = false;
			runtimeMode = RUN_SERIAL;
			filterMode = FILTERING_GRAMS;
			processors= -1;
			minLength = 16;
			maxLength = -1;
			tolError = 0;
			maxInterruptions = 1;
			filterRepeats = 1;
			minRepeatLength = 10;
			maxRepeatPeriod = 4;
			
			computeTpot = false;
			prettyString = false;
			outputFormat = 0;
			errorReference = WATSON_STAND;
			mergeFeatures = false;
			runID = "s";
			mixed_parallel_max_guanine     = 1.;
			mixed_antiparallel_min_guanine = 0.;
			
#ifdef BOOST
			compressOutput = false;
#endif
			version = "";
			shape = "11111";
			logFileName = "triplex_search.log";
			summaryFileName = "triplex_search.summary";
			
			timeLoadFiles = 0.0;
			timeFindTriplexes = 0.0;
			timeFindTfos = 0.0;
			timeFindTtss = 0.0;
			timeDumpResults = 0.0;
			
			applyMaximumLengthConstraint = false;
		}
	};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

	//////////////////////////////////////////////////////////////////////////////
	// produce a timestamp for the log file
	inline int _calculateShape(Options &options){	
		// number of errors when using error-rate
		int errors = static_cast<int>(ceil(options.errorRate*options.minLength));
		// check if this number is capped by a specified maximum
		if (options.maximalError >=0)
			errors  = min(errors, options.maximalError);
		int qgramWeight =static_cast<int>(min(14.0,floor((options.qgramThreshold-1-options.minLength)/-(errors+1))));
		resize(options.shape,qgramWeight);
		for (int i=0;i<qgramWeight;++i){
			options.shape[i]='1';
		}		
#ifdef TRIPLEX_DEBUG
		::std::cerr << "qgram-threshold:" << options.qgramThreshold << " weight(qgram):" << qgramWeight << ::std::endl;
#endif
		return qgramWeight;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// produce a timestamp for the log file
	CharString _getTimeStamp(){
		time_t rawtime;
		struct tm * timeinfo = new tm;
		char buffer [80];
		time(&rawtime);
		localtime_r(&rawtime, timeinfo);
		strftime (buffer,80,"[%x %X]",timeinfo);
		delete timeinfo;
		return buffer;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Search a sequence compatible with the TC-motif and add these for parallel binding
	template <typename TOligoMotifSet, typename TString, typename TId>
	inline void processTCMotif(TOligoMotifSet	&motifSet,
							   TString			&sequence,
							   TId const		&tfoSeqNo,
							   bool const		reduceSet,
							   Options const	&options
							   ){
		typedef typename Value<TOligoMotifSet>::Type			TTfoMotif;
		typedef typename Iterator<TString>::Type 				TIter;
		typedef typename Infix<TString>::Type					TSegment;
		typedef String<TSegment>								TSegString;
		typedef typename Iterator<TSegString, Standard>::Type	TSegStringIter;
		
		// parse TFOs for valid substrings with respect to maximum number of consecutive interruptions
		TString valid   = "TCY";   // the valid characters
		TString invalid = "GARN";  // the interrupting characters
		
		// create parser
		TGraph parser;
		_makeParser(parser, valid, invalid, options);
		
		// split tfo sequence into valid parts
		TSegString seqString;	// target segment container
		_parse(seqString,parser, sequence, options);
		
		// convert tfo sequences into matching tts to allow pattern search
		for (TSegStringIter it = begin(seqString, Standard()); it != end(seqString, Standard()); ++it){
			TTfoMotif tfomotif(*it, true, tfoSeqNo, true, 'Y');
			_filterWithGuanineAndErrorRate(motifSet, tfomotif, 'G', 'N', reduceSet, TRIPLEX_ORIENTATION_PARALLEL, options, PYRIMIDINEMOTIF());
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Search a sequence compatible with the GA-motif and add these for anti-parallel binding
	template <typename TOligoMotifSet, typename TString, typename TId>
	inline void processGAMotif(TOligoMotifSet	&motifSet,
							   TString			&sequence,
							   TId const		&tfoSeqNo,
							   bool const		reduceSet,
							   Options const	&options
							   ){
		typedef typename Value<TOligoMotifSet>::Type			TTfoMotif;
		typedef typename Iterator<TString>::Type				TIter;
		typedef typename Infix<TString>::Type					TSegment;
		typedef String<TSegment>								TSegString;
		typedef typename Iterator<TSegString, Standard>::Type	TSegStringIter;
		
		// parse TFOs for valid substrings with respect to maximum number of consecutive interruptions
		TString valid   = "GAR";  // the valid characters
		TString invalid = "TCYN"; // the interrupting characters
		
		// create parser
		TGraph parser;
		_makeParser(parser, valid, invalid, options);
		
		// split tfo sequence into valid parts
		TSegString seqString;	// target segment container
		_parse(seqString,parser, sequence, options);
		
		// convert tfo sequences into matching tts to allow pattern search
		for (TSegStringIter it = begin(seqString, Standard()); it != end(seqString, Standard()); ++it){
			TTfoMotif tfomotif(*it, false, tfoSeqNo, true, 'R');
			_filterWithGuanineAndErrorRate(motifSet, tfomotif, 'G', 'N', reduceSet, TRIPLEX_ORIENTATION_ANTIPARALLEL, options, PURINEMOTIF());
		}
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Search a sequence compatible with the GT-motif and add these for 
	// the requested binding orientation(s)
	template <typename TOligoMotifSet, typename TString, typename TId>
	inline void processGTMotif(TOligoMotifSet		&motifSet,
							   TString				&sequence,
							   TId const			&tfoSeqNo,
							   ORIENTATION const	orientation, // the orientation to be considered (>0 parallel, <0 antiparallel, 0=both
							   bool const			reduceSet,
							   Options const		&options
							   ){		
		typedef typename Value<TOligoMotifSet>::Type			TTfoMotif;
		typedef typename Iterator<TString>::Type				TIter;
		typedef typename Infix<TString>::Type					TSegment;
		typedef ModifiedString<TSegment, ModView< FunctorRYFilter > >  	TFilter;
		typedef String<TSegment>								TSegString;
		typedef typename Iterator<TSegString, Standard>::Type	TSegStringIter;
		
		// parse TFOs for valid substrings with respect to maximum number of consecutive interruptions
		TString valid   = "GTK";  // the valid characters
		TString invalid = "CAMN"; // the interrupting characters
		
		// create parser
		TGraph parser;
		_makeParser(parser, valid, invalid, options);
		
		// split tfo sequence into valid parts
		TSegString seqString;	// target segment container
		_parse(seqString,parser, sequence, options);
		
		// convert tfo sequences into matching tts to allow pattern search
		for (TSegStringIter it = begin(seqString, Standard()); it != end(seqString, Standard()); ++it){
#ifdef TRIPLEX_DEBUG
			::std::cerr << "processing:" << *it << ::std::endl;
#endif
			if ((orientation == TRIPLEX_ORIENTATION_BOTH || orientation == TRIPLEX_ORIENTATION_PARALLEL) && options.mixed_parallel_max_guanine >= options.guanineRate){
				TTfoMotif tfomotif(*it, true, tfoSeqNo, true, 'M');
				_filterWithGuanineAndErrorRate(motifSet, tfomotif, 'G', 'N', reduceSet, orientation, options, MIXEDMOTIF());
			}
			if (orientation == TRIPLEX_ORIENTATION_BOTH || orientation == TRIPLEX_ORIENTATION_ANTIPARALLEL){
				TTfoMotif tfomotif_rev(*it, false, tfoSeqNo, true, 'M');
				_filterWithGuanineAndErrorRate(motifSet, tfomotif_rev, 'G', 'N', reduceSet, orientation, options, MIXEDMOTIF());
			}
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Search a sequence for a putative triplex target
	template <typename TDuplexMotifSet, typename TString, typename TId>
	inline void processDuplex(TDuplexMotifSet	&ttsSet,
							  TString			&duplex,
							  TId const			&seqNo,
							  bool const		plusstrand,
							  bool const		reduceSet,
							  Options			&options
							  ){
		typedef typename Value<TDuplexMotifSet>::Type				TTtsMotif;
		typedef typename Iterator<TString>::Type 					TIter;
		typedef typename Infix<TString>::Type						TSegment;
		typedef String<TSegment>									TSegString;
		typedef typename Iterator<TSegString, Standard>::Type		TSegStringIter;
		
		// parse duplex for valid substrings with respect to maximum number of consecutive interruptions
		TString valid;		// the valid characters
		TString invalid;	// the interrupting characters
		if (plusstrand){
			valid = "GAR";
			invalid  = "TCYN";
		} else {
			valid = "TCY";
			invalid  = "GARN";
		}
		// create parser
		TGraph parser;		
		_makeParser(parser, valid, invalid, options);
			
		// split duplex into valid parts
		TSegString seqString;	// target segment container
		_parse(seqString, parser, duplex, options);
		
		// process one segment at a time
		for (TSegStringIter it = begin(seqString, Standard()); it != end(seqString, Standard()); ++it){
#ifdef TRIPLEX_DEBUG
			::std::cerr << "pTTS:" << *it << ::std::endl;
#endif
			if (plusstrand){
				TTtsMotif ttsfilter(*it, true, seqNo, false, '+');
				_filterWithGuanineAndErrorRate(ttsSet, ttsfilter, 'G', 'Y', reduceSet, TRIPLEX_ORIENTATION_BOTH, options, TTS());
			} else {
				TTtsMotif ttsfilter(*it, true, seqNo, false, '-');
				_filterWithGuanineAndErrorRate(ttsSet, ttsfilter, 'G', 'Y', reduceSet, TRIPLEX_ORIENTATION_BOTH, options, TTS());
			}
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Search for triplexes given a set of duplexes and a set of TFOs
	template<
	typename TId, 
	typename TGardenerSpec,
	typename TPattern,
	typename TQuery
	>
	inline void _detectTriplex(Gardener< TId, TGardenerSpec>	&gardener,
							   TPattern	const					&pattern,
							   TQuery							&ttsSet,
							   Options const					&options
							   ){
		
		// adjust errorRate if maximalError is set and caps the errorRate setting wrt the minimum length constraint
		double eR = options.errorRate;
		if (options.maximalError >= 0){
			eR = min(options.errorRate, max(double(options.maximalError)/options.minLength, 0.0));
		}
#if SEQAN_ENABLE_PARALLELISM	
		if (options.runtimeMode==RUN_PARALLEL_TRIPLEX){
			plant(gardener, pattern, ttsSet, eR, options.minLength, options.maxInterruptions+1, MULTIPLE_WORKER() );
		} else {
#endif
			plant(gardener, pattern, ttsSet, eR, options.minLength, options.maxInterruptions+1, SINGLE_WORKER() );
#if SEQAN_ENABLE_PARALLELISM
		}
#endif
	}

	
	
#if SEQAN_ENABLE_PARALLELISM
	
	//////////////////////////////////////////////////////////////////////////////
	// copy the matches across from the source to the sink
	template<
	typename TMatches
	>
	inline void _saveMatches(TMatches &match_sink,
							 TMatches & match_source
							 ){
		
		typedef typename Iterator<TMatches, Standard>::Type	TIter;
		typedef typename Value<TMatches>::Type				TMatch;
		
		for (TIter it = begin(match_source, Standard()); it != end(match_source, Standard()); ++it){
			appendValue(match_sink, *it);
		}
		
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Search for a triplex given a target string and a set of TFOs on both 
	// strands of the duplex in parallel, which requires about twice as much memory	
	template<
	typename TMatches,
	typename TId, 
	typename TPattern,
	typename TDuplex,
	typename TGardenerSpec
	>
	void _detectTriplexParallelStrands(TMatches			&matches,
									   TPattern const	&pattern,
									   TDuplex			&duplexString,
									   TId const		&duplexId,
									   Options			&options,
									   Gardener<TId, TGardenerSpec>
									   ){	
		typedef StringSet<ModStringTriplex<TDuplex, TDuplex> > 	TDuplexModSet;
		typedef Gardener< TId, TGardenerSpec>					TGardener;
		
		
		TGardener gardener_forward;
		TGardener gardener_reverse;
		TDuplexModSet ttsSet_forward;
		TDuplexModSet ttsSet_reverse;
		TMatches matches_forward;
		TMatches matches_reverse;
		
		bool reduceSet = true; // merge overlapping features
		
		omp_set_num_threads(2);
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel sections)
		//#pragma omp parallel sections // starts new team
		{	
			SEQAN_PRAGMA_IF_PARALLEL(omp section)
			//#pragma omp section
			{
				// prefilter for putative TTSs
				processDuplex(ttsSet_forward, duplexString, duplexId, true, reduceSet, options);
				if (length(ttsSet_forward)>0){
					_detectTriplex(gardener_forward, pattern, ttsSet_forward, options);
					_filterAndStore(matches_forward, gardener_forward, pattern, ttsSet_forward, duplexId, true, options);
					
				}
				
			}
			SEQAN_PRAGMA_IF_PARALLEL(omp section)
			//#pragma omp section
			{
				// prefilter for putative TTSs
				processDuplex(ttsSet_reverse, duplexString, duplexId, false, reduceSet, options);
				if (length(gardener_reverse)>0){
					_detectTriplex(gardener_reverse, pattern, ttsSet_reverse, options);
					_filterAndStore(matches_reverse, gardener_reverse, pattern, ttsSet_reverse, duplexId, false, options);
				}			
			}
			
		}  /* end of sections */
		
		_saveMatches(matches, matches_forward);
		_saveMatches(matches, matches_reverse);
		
		eraseAll(gardener_forward);
		eraseAll(gardener_reverse);
		
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Search for a triplex given a target string and a set of TFOs on both 
	// strands of the duplex in parallel, which requires about twice as much memory	
	template<
	typename TMatches,
	typename TId, 
	typename TPattern,
	typename TDuplex
	>
	void _detectTriplexParallelStrands(TMatches			&matches,
									   TPattern			&tfoSet,
									   TDuplex			&duplexString,
									   TId const		&duplexId,
									   Options			&options,
									   BruteForce
									   ){	
		typedef StringSet<ModStringTriplex<TDuplex, TDuplex> > 	TDuplexModSet;
		
		TDuplexModSet ttsSet_forward;
		TDuplexModSet ttsSet_reverse;
		TMatches matches_forward;
		TMatches matches_reverse;
		
		bool reduceSet = true; // merge overlapping features
		
		omp_set_num_threads(2);
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel sections)	
		//#pragma omp parallel sections // starts new team
		{	
			SEQAN_PRAGMA_IF_PARALLEL(omp section)
			//#pragma omp section
			{
				// prefilter for putative TTSs
				processDuplex(ttsSet_forward, duplexString, duplexId, true, reduceSet, options);
				_detectTriplexBruteForce(matches_forward, tfoSet, ttsSet_forward, duplexId, options);					
			}
			SEQAN_PRAGMA_IF_PARALLEL(omp section)
			//#pragma omp section
			{
				// prefilter for putative TTSs
				processDuplex(ttsSet_reverse, duplexString, duplexId, false, reduceSet, options);
				_detectTriplexBruteForce(matches_reverse, tfoSet, ttsSet_reverse, duplexId, options);					
			}
			
		}  /* end of sections */
		
		_saveMatches(matches, matches_forward);
		_saveMatches(matches, matches_reverse);
		
	}
#endif
	
	//////////////////////////////////////////////////////////////////////////////
	// Search for a triplex given a target string and a set of TFOs
	template<
	typename TMatches,
	typename TId, 
	typename TPattern,
	typename TDuplex,
	typename TGardenerSpec
	>
	void _detectTriplex(TMatches		&matches,
						TPattern const	&pattern,
						TDuplex			&duplexString,
						TId const		&duplexId,
						Options			&options,
						Gardener<TId, TGardenerSpec>
						){	
		typedef StringSet<ModStringTriplex<TDuplex, TDuplex> > 	TDuplexModSet;
		typedef Gardener<TId, TGardenerSpec>					TGardener;
		
		bool reduceSet = true; // merge overlapping features
		
		if (options.forward){
			TGardener gardener_forward;
			TDuplexModSet ttsSet_forward;
			// prefilter for putative TTSs
			processDuplex(ttsSet_forward, duplexString, duplexId, true, reduceSet, options);
#ifdef TRIPLEX_DEBUG
			typedef typename Iterator<TDuplexModSet>::Type  TIterMotifSet;
			::std::cerr << "printing all tts segments (forward)" << ::std::endl;
			for (TIterMotifSet itr=begin(ttsSet_forward); itr != end(ttsSet_forward);++itr){
				::std::cerr << "tts: " << ttsString(*itr) << " type: " << (*itr).motif << " length: "<< length(*itr) <<  " position: "<< beginPosition(*itr) << " " << ::std::endl;
			}
#endif
			
			if (length(ttsSet_forward)>0){
				_detectTriplex(gardener_forward, pattern, ttsSet_forward, options);
				_filterAndStore(matches, gardener_forward, pattern, ttsSet_forward, duplexId, true, options);
			}
			eraseAll(gardener_forward);
		}
		
		if (options.reverse) {
			TGardener gardener_reverse;
			TDuplexModSet ttsSet_reverse;
			// prefilter for putative TTSs
			processDuplex(ttsSet_reverse, duplexString, duplexId, false, reduceSet, options);
#ifdef TRIPLEX_DEBUG
			typedef typename Iterator<TDuplexModSet>::Type  TIterMotifSet;
			::std::cerr << "printing all tts segments (reverse)" << ::std::endl;
			for (TIterMotifSet itr=begin(ttsSet_reverse); itr != end(ttsSet_reverse);++itr){
				::std::cerr << "tts: " << ttsString(*itr) << " type: " << (*itr).motif << " length: "<< length(*itr) <<  " position: "<< beginPosition(*itr) << " " << ::std::endl;
			}
#endif
			if (length(ttsSet_reverse)>0){
				_detectTriplex(gardener_reverse, pattern, ttsSet_reverse, options);
				_filterAndStore(matches, gardener_reverse, pattern, ttsSet_reverse, duplexId, false, options);
			}
			eraseAll(gardener_reverse);
		}
	}

	
	//////////////////////////////////////////////////////////////////////////////
	// Search for a triplex given a target string and a set of TFOs
	template<
	typename TMatches,
	typename TId, 
	typename TPatterns,
	typename TDuplex
	>
	void _detectTriplex(TMatches		&matches,
						TPatterns		&tfoSet,
						TDuplex			&duplexString,
						TId const		&duplexId,
						Options			&options,
						BruteForce
								  ){	
		typedef ModStringTriplex<TDuplex, TDuplex>	TTts;
		typedef StringSet<TTts>						TTtsSet;
		typedef typename Iterator<TTtsSet>::Type	TTtsIter;
		
		bool reduceSet = true; // merge overlapping features
		TTtsSet ttsSet;
		// prefilter for putative TTSs
		if (options.forward) {
			processDuplex(ttsSet, duplexString, duplexId, true, reduceSet, options);
		}
		if (options.reverse) {
			processDuplex(ttsSet, duplexString, duplexId, false, reduceSet, options);
		}
#ifdef TRIPLEX_DEBUG
		::std::cerr << "printing all tts segments" << ::std::endl;
		for (TTtsIter itr=begin(ttsSet); itr != end(ttsSet);++itr){
			::std::cerr << "tts: " << ttsString(*itr) << " type: " << (*itr).motif << " length: "<< length(*itr) <<  " position: "<< beginPosition(*itr) << " " << ::std::endl;
		}
#endif
		
#if SEQAN_ENABLE_PARALLELISM	
		if (options.runtimeMode==RUN_PARALLEL_TRIPLEX){
			
			SEQAN_PRAGMA_IF_PARALLEL(omp parallel)
			{
				SEQAN_PRAGMA_IF_PARALLEL(omp for schedule(dynamic) )
				for (int tts=0; tts<length(ttsSet); ++tts){
					
					TTtsSet tmp_ttsSet;
					appendValue(tmp_ttsSet, ttsSet[tts]);
					TMatches tmp_matches;
					_detectTriplexBruteForce(tmp_matches, tfoSet, tmp_ttsSet, duplexId, options);
					
					if (length(tmp_matches)>0){
						SEQAN_PRAGMA_IF_PARALLEL(omp critical(addMatches) )
						_saveMatches(matches, tmp_matches);
					}
				}
			}
		} else		
#endif
		_detectTriplexBruteForce(matches, tfoSet, ttsSet, duplexId, options);
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Search for a triplex given a target string and a set of TFOs
	template<
	typename TMatches,
	typename TId, 
	typename TPatterns,
	typename TTtsSet
	>
	void _detectTriplexBruteForce(TMatches			&matches,
								  TPatterns			&tfoSet,
								  TTtsSet			&ttsSet,
								  TId const			&duplexId,
								  Options			&options
								  ){	
		typedef typename Value<TTtsSet>::Type		TTts;
		typedef typename Host<TTts>::Type			TDuplex;
		typedef typename Iterator<TTtsSet>::Type	TTtsIter;
		typedef typename Iterator<TPatterns>::Type	TPattIter;
		typedef typename Value<TPatterns>::Type		TPattern;
		typedef typename Host<TPattern>::Type		TTfo;
		typedef typename Position<TPatterns>::Type	TPos;
		typedef typename Infix<TDuplex>::Type		TSegment;
		typedef String<TSegment>					TSegString;
		typedef typename Iterator<TSegString>::Type	TSegStringIter;
		
		int minScore = options.minLength- static_cast<int>(ceil(options.errorRate * options.minLength));
		if (length(ttsSet)>0){
			// iterate over all TFO candidates
			int tfoNo = 0;
			for (TPattIter itO = begin(tfoSet);itO != end(tfoSet); ++itO, ++tfoNo){
#ifdef TRIPLEX_DEBUG
				TPattern pat = *itO;
				::std::cerr << "tfo candidate: " << pat << ::std::endl << "--- candidate: " << tfoString(*itO) << " " << getMotif(*itO) <<  " " << isParallel(*itO) << :: std::endl;
#endif
				TTfo tfoCandidate = ttsString(*itO);
				Finder<TTfo > finder(tfoCandidate);
				// iterate over all TTS candidates
				for (TTtsIter itD = begin(ttsSet); itD != end(ttsSet); ++itD){
					TDuplex ttsCandidate = ttsString(*itD);
#ifdef TRIPLEX_DEBUG
					::std::cerr << "tts candidate: " << ttsCandidate << ::std::endl;
#endif
					// iterate all suitable (>= minLength) diagonals
					for (int diag = -(length(ttsCandidate)-options.minLength); diag <= length(tfoCandidate)-options.minLength; ++diag){
						int offsetTts = 0;
						int offsetTfo = 0;
						if (diag < 0)
							offsetTts = -diag;
						else if (diag > 0)
							offsetTfo = diag;
						
						int lenS = min(length(ttsCandidate)-offsetTts, length(tfoCandidate)-offsetTfo);
						TDuplex triplex(infix(ttsCandidate, offsetTts, offsetTts+lenS));
						TTfo tfo(infix(tfoCandidate, offsetTfo, offsetTfo+lenS));
#ifdef TRIPLEX_DEBUG
						::std::cerr << "tts: " << triplex << ::std::endl;
						::std::cerr << "tfo: " << tfo << ::std::endl;
#endif	
						int m = 0; //matches between the sequences
						for (int k=0; k<lenS; ++k){
							if (triplex[k]!=tfo[k]){
								triplex[k]='N';
							} else {
								++ m;
							}
						}
#ifdef TRIPLEX_DEBUG	
						::std::cerr << "tts: " << triplex <<  " " << (m < minScore?"discarded":"to verify") << ::std::endl;
#endif					
						// check for minimum requirement of matching positions
						if (m >= minScore){
							// run through TTS parser
							TTtsSet triplexSet;
							// create parser once
							if (empty(options.triplexParser)){
								TDuplex valid("GAR");		// the valid characters
								TDuplex invalid("TCYN");	// the interrupting characters			
								_makeParser(options.triplexParser, valid, invalid, options);
							}
							
							// split duplex into valid parts
							TSegString seqString;	// target segment container
							_parse(seqString, options.triplexParser, triplex, options);
							
							// process one segment at a time
							for (TSegStringIter it = begin(seqString, Standard()); it != end(seqString, Standard()); ++it){
#ifdef TRIPLEX_DEBUG
								::std::cerr << "pTriplex:" << *it << ::std::endl;
#endif
								TTts ttsfilter(*it, true, getSequenceNo(*itD), false, '+');
								bool reduceSet = false; // don't merge overlapping triplexes	
								_filterWithGuanineAndErrorRate(triplexSet, ttsfilter, 'G', 'Y', reduceSet, TRIPLEX_ORIENTATION_BOTH, options, TTS());
							}
							
							// save all matches
							TPos tfoStart;
							TPos tfoEnd;
							TPos ttsStart;
							TPos ttsEnd;
							char strand;
							for (TTtsIter itr=begin(triplexSet); itr!=end(triplexSet); ++itr){
								// compute score = matching positions, which can be found in the triplex string (number of N's within the interval found)
								int score=0;
								for (unsigned int i=beginPosition(*itr); i<endPosition(*itr); ++i){
									if (triplex[i]!='N')
										score += 1;
								}
								
								// calculate tts positions according to strand in the duplex
								if (getMotif(*itD)=='+'){
									ttsStart = offsetTts + beginPosition(*itD)+ beginPosition(*itr);
									ttsEnd = ttsStart + length(*itr);
									strand = '+';			
								} else {
									ttsEnd = endPosition(*itD) - (offsetTts + beginPosition(*itr));
									ttsStart = ttsEnd - length(*itr);
									strand = '-';
								}
								
								// calculate tfo positions according to binding orientation
								if (isParallel(*itO)){
									tfoStart = offsetTfo + beginPosition(*itO)+beginPosition(*itr);
									tfoEnd = tfoStart + length(*itr);
								} else {
									tfoEnd = endPosition(*itO) - (offsetTfo + beginPosition(*itr));
									tfoStart = tfoEnd - length(*itr);
									
								}
								
								// save the corresponding triplex match 
								TMatch match(tfoNo,
											 tfoStart,
											 tfoEnd,
											 duplexId,
											 ttsStart,
											 ttsEnd,
											 score,
											 isParallel(*itO),
											 getMotif(*itO),
											 strand
											 );
								appendValue(matches, match);
								
#ifdef TRIPLEX_DEBUG
								::std::cerr << "tts: " << ttsString(*itr) << " length: "<< length(*itr) <<  " position: "<< beginPosition(*itr) << "-" <<endPosition(*itr) << " " << infix(triplex, beginPosition(*itr), endPosition(*itr)) << ::std::endl;
								TSegment ftfo = infix(host(*itO), tfoStart, tfoEnd);
								TSegment ftts = infix(host(*itD), ttsStart, ttsEnd);
								::std::cerr << "ttsC: " << tfoCandidate << ::std::endl << "tfoC: "<< tfoCandidate << ::std::endl;
								::std::cerr << "tts: " << ftts << ::std::endl << "tfo: "<< ftfo << ::std::endl;
#endif
							}
							
						}
						
					}
				}
				
			}
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Outputs the block queue 
	void _printQueue(std::vector<QueueBlock> queue){
		for (unsigned i=0;i<length(queue);++i){
			::std::cout << queue[i].blockSize << ":" << queue[i].blockClass << ::std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Count the number of times a particular string is present in the set of TTSs
	template <typename TSpec>
	inline void _countDuplicatesPermissive(StringSet<ModStringTriplex<TSpec, TSpec> > &strings, 
										  Options const &options, 
										  TTS const & 
	){
		typedef ModStringTriplex<TSpec, TSpec>		TString;
		typedef typename Host<TString>::Type		THost;
		typedef typename Value<TString>::Type		THostValue;		
		typedef StringSet<TString>					TStringSet;
		typedef StringSet<TriplexString>			TText;
		typedef typename Iterator<TStringSet>::Type	TIter;
		typedef typename Iterator<TText>::Type		TSIter;
		typedef Index<TText>						TIndex;
		typedef ModifiedString<ModifiedString<THost, ModView< FunctorComplement<THostValue> > >, ModReverse> TModRevCompl;
		typedef typename Position<TString>::Type	TPos;
		typedef typename Id<TString>::Type			TId;
		typedef SeqPos<TId, TPos>					TKey;
		typedef Map<TKey, Skiplist<> >				TMap;
		typedef typename Iterator<TMap>::Type		TMapIter;
		
		// create index structure based on underlying sequence
		TText sset;
		for (TIter it=begin(strings); it != end(strings); ++it){
			if (getMotif(*it)=='+')
				appendValue(sset, getSegment(*it));
			else
				appendValue(sset, TModRevCompl(getSegment(*it)));
		}
		
		if (options._debugLevel > 1)
			::std::cerr << "Finished preparing index structure for abundancy counter: " << length(sset) << " sequences\n";
		
		TIndex index_esa(sset);
		if (options._debugLevel > 1)
			::std::cerr << "Finished creating index structure for abundancy counter\n";
		
		
		if (options._debugLevel > 1)
			::std::cerr << "Counting duplicates in " << length(strings) <<  " sequences\n";
		
		
		for (TIter it=begin(strings); it != end(strings); ++it){
			TMap lociMap;
			int occ = 0;
			if (options._debugLevel > 2)
				::std::cerr << "Detecting duplicates of " << outputString(*it) << "\t" << ttsString(*it) << ":";
			
			Finder<TIndex> finder_esa(index_esa);
			TKey identity(getSequenceNo(*it), beginPosition(*it));
			
			TriplexString pattern;
			if (getMotif(*it)=='+')
				pattern = getSegment(*it);
			else
				pattern = TModRevCompl(getSegment(*it));
			
			while(find(finder_esa, pattern) ){
				TString ts = value(strings, position(finder_esa).i1);
				
				// skip same sequence duplicates if requested
				if (!options.sameSequenceDuplicates && getSequenceNo(*it) == getSequenceNo(ts))
					continue;
					
				if (getMotif(ts)=='+'){
					TKey key(getSequenceNo(ts), beginPosition(ts)+position(finder_esa).i2);
					if (key != identity && !hasKey(lociMap, key)){
						++occ;
						add(lociMap, key);
					} 
#ifdef TRIPLEX_DEBUG	
					else {
						::std::cerr << "Found duplicate " << getSequenceNo(ts) << " " << beginPosition(ts)+position(finder_esa).i2 << ::std::endl;
					}
					::std::cerr << ttsString(*it) << " " << position(finder_esa) << " " << getSequenceNo(ts) << " " << getMotif(ts) << ":" << beginPosition(ts)+position(finder_esa).i2 << " " << infix(host(ts), beginPosition(ts)+position(finder_esa).i2, beginPosition(ts)+position(finder_esa).i2+length(*it)) << ::std::endl;
#endif
				} else {
					TKey key(getSequenceNo(ts), endPosition(ts)-position(finder_esa).i2-length(*it));
					if (key != identity && !hasKey(lociMap, key)){
						++occ;
						add(lociMap, key);
					} 
#ifdef TRIPLEX_DEBUG	
					else {
						::std::cerr << "Found duplicate " << getSequenceNo(ts) << " " << beginPosition(ts)+position(finder_esa).i2 << ::std::endl;
					}
					::std::cerr << ttsString(*it) << " " << position(finder_esa) << " " << getSequenceNo(ts) << " " << getMotif(ts) << ":" << endPosition(ts)-position(finder_esa).i2 << " " << infix(host(ts), endPosition(ts)-position(finder_esa).i2-length(*it), endPosition(ts)-position(finder_esa).i2) << ::std::endl;
#endif
				}
				if (options._debugLevel > 2 && occ>0 && occ % 10000 == 0)
					::std::cerr << occ << ":";
				
			}
			duplicates(*it, occ);
			
			// record duplicate locations
			if (options.reportDuplicateLocations && options.duplicatesCutoff>=0 && occ <= options.duplicatesCutoff){
				for (TMapIter mit = begin(lociMap); mit != end(lociMap); ++mit){
					addDuplicate(*it, getSequenceNo(*mit), getPosition(*mit));
				}
			}
			
			if (options._debugLevel > 2)
				::std::cerr << occ << " duplicates found" << ::std::endl;
			
		}
		
		if (options._debugLevel > 1)
			::std::cerr << "Finished counting duplicates for all features\n";

	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Count the number of times a particular string is present in the set of TFOs
	template <typename TSpec>
	inline void _countDuplicatesPermissive(StringSet<ModStringTriplex<TSpec, TSpec> > &strings, 
										  Options const &options, 
										  TFO const & 
	){
		typedef ModStringTriplex<TSpec, TSpec>		TString;
		typedef typename Host<TString>::Type		THost;
		typedef typename Value<TString>::Type		THostValue;		
		typedef StringSet<TString>					TStringSet;
		typedef StringSet<TriplexString>			TText;
		typedef typename Iterator<TStringSet>::Type	TIter;
		typedef typename Iterator<TText>::Type		TSIter;
		typedef Index<TText>						TIndex;
		typedef ModifiedString<ModifiedString<THost, ModView< FunctorComplement<THostValue> > >, ModReverse> TModRevCompl;
		typedef typename Position<TString>::Type	TPos;
		typedef typename Id<TString>::Type			TId;
		typedef SeqPos<TId, TPos>					TKey;
		typedef Map<TKey, Skiplist<> >				TMap;
		typedef typename Iterator<TMap>::Type		TMapIter;
		
		// create index structure based on underlying sequence
		TText sset;
		for (TIter it=begin(strings); it != end(strings); ++it){
			appendValue(sset, tfoString(*it));
		}
		
		if (options._debugLevel > 1)
			::std::cerr << "Finished preparing index structure for abundancy counter: " << length(sset) << " sequences\n";
		
		TIndex index_esa(sset);
		if (options._debugLevel > 1)
			::std::cerr << "Finished creating index structure for abundancy counter\n";
		
		
		if (options._debugLevel > 1)
			::std::cerr << "Counting duplicates in " << length(strings) <<  " sequences\n";
		
		
		for (TIter it=begin(strings); it != end(strings); ++it){
			TMap lociMap;
			int occ = 0;
			if (options._debugLevel > 2)
				::std::cerr << "Detecting duplicates of " << outputString(*it) << "\t" << ttsString(*it) << ":";
			
			Finder<TIndex> finder_esa(index_esa);
			TKey identity(getSequenceNo(*it), beginPosition(*it) );
			
			while(find(finder_esa, tfoString(*it)) ){
				TString ts = value(strings, position(finder_esa).i1);
				
				// skip same sequence duplicates if requested
				if (!options.sameSequenceDuplicates && getSequenceNo(*it) == getSequenceNo(ts))
					continue;
					
				TKey key(getSequenceNo(ts), beginPosition(ts)+position(finder_esa).i2);
				if (key != identity && !hasKey(lociMap, key)){
					++occ;
					add(lociMap, key);
				} 
				// stop iterating if duplicatesCutoff is set and already exceeded
				if (options.duplicatesCutoff>=0 && occ > options.duplicatesCutoff)
					break;
				
			}
			duplicates(*it, occ);
			
			// record duplicate locations
			if (options.reportDuplicateLocations && options.duplicatesCutoff>=0 && occ <= options.duplicatesCutoff){
				for (TMapIter mit = begin(lociMap); mit != end(lociMap); ++mit){
					addDuplicate(*it, getSequenceNo(*mit), getPosition(*mit));
				}
			}
			
			if (options._debugLevel > 2)
				::std::cerr << occ << " duplicates found" << ::std::endl;
			
		}
		
		if (options._debugLevel > 1)
			::std::cerr << "Finished counting duplicates for all features\n";
		
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Count the number of times a particular string is present in the set wrt. target space.
	// Duplicates are detected with respect to the starting position of the corresponding) target site.
	template <typename TSpec, typename TTag>
	inline void _countDuplicatesStrict(StringSet<ModStringTriplex<TSpec, TSpec> > &strings, 
									  Options &options, 
									  TTag const & 
	){
		typedef ModStringTriplex<TSpec, TSpec>		TString;
		typedef StringSet<TString>					TStringSet;
		typedef StringSet<TriplexString>			TText;
		typedef typename Iterator<TStringSet>::Type	TIter;
		typedef typename Iterator<TText>::Type		TSIter;
		typedef Index<TText>						TIndex;
		typedef typename Position<TString>::Type	TPos;
		typedef typename Id<TString>::Type			TId;
		typedef SeqPos<TId, TPos>					TKey;
		typedef Map<TKey, Skiplist<> >				TMap;
		typedef typename Iterator<TMap>::Type		TMapIter;
		
		// create index structure in target space
		if (options._debugLevel >= 1)
			options.logFileHandle << _getTimeStamp() << "   ... Started preparing index structure for duplicate detection on " << length(strings) << " features" << ::std::endl;
		
		
		TText sset;
		for (TIter it=begin(strings); it != end(strings); ++it){
			// a feature has a duplicate, if another string in the ttsString space
			// has the same sequence
			appendValue(sset, ttsString(*it));
		}
		TIndex index_esa(sset);
		
		if (options._debugLevel >= 1){
			options.logFileHandle << _getTimeStamp() << "   ... Finished preparing index structure for duplicate detection on " << length(strings) << " features" << ::std::endl;
			options.logFileHandle << _getTimeStamp() << "   ... Started counting duplicates" << ::std::endl;
		}
		
		for (TIter it=begin(strings); it != end(strings); ++it){
			TMap lociMap;
			int occ = 0;
#ifdef TRIPLEX_DEBUG	
				::std::cerr << "Detecting duplicates of " << outputString(*it) << "\t" << ttsString(*it) << ":";
#endif			
			Finder<TIndex> finder_esa(index_esa);
			TKey identity(getSequenceNo(*it), beginPosition(*it) );
			
			while(find(finder_esa, ttsString(*it)) ){
				TString ts = value(strings, position(finder_esa).i1);

				// skip same sequence duplicates if requested
				if (!options.sameSequenceDuplicates && getSequenceNo(*it) == getSequenceNo(ts))
					continue;
				
				if ( isTFO(ts) || getMotif(ts)=='+'){
					TKey key(getSequenceNo(ts), beginPosition(ts)+position(finder_esa).i2);
					// check for identity or whether location has been reported already
					if (key != identity && !hasKey(lociMap, key)){
						++occ;
						add(lociMap, key);
					} 
#ifdef TRIPLEX_DEBUG	
					else {
						::std::cerr << "Found duplicate " << getSequenceNo(ts) << " " << beginPosition(ts)+position(finder_esa).i2 << ::std::endl;
					}
					::std::cerr << ttsString(*it) << " " << position(finder_esa) << " " << getSequenceNo(ts) << " " << getMotif(ts) << ":" << beginPosition(ts)+position(finder_esa).i2 << " " << infix(host(ts), beginPosition(ts)+position(finder_esa).i2, beginPosition(ts)+position(finder_esa).i2+length(*it)) << ::std::endl;
#endif
				} else { 
					TKey key(getSequenceNo(ts), endPosition(ts)-position(finder_esa).i2-length(*it));
					if (key != identity && !hasKey(lociMap, key)){
						++occ;
						add(lociMap, key);
					} 
#ifdef TRIPLEX_DEBUG
					else {
						::std::cerr << "Found duplicate " << getSequenceNo(ts) << " " << beginPosition(ts)+position(finder_esa).i2 << ::std::endl;
					}
					::std::cerr << ttsString(*it) << " " << position(finder_esa) << " " << getSequenceNo(ts) << " " << getMotif(ts) << ":" << endPosition(ts)-position(finder_esa).i2-length(*it) << " " << infix(host(ts), endPosition(ts)-position(finder_esa).i2-length(*it), endPosition(ts)-position(finder_esa).i2)  << ::std::endl;
#endif					
				}				
				// stop iterating if duplicatesCutoff is set and already exceeded
				if (options.duplicatesCutoff>=0 && occ > options.duplicatesCutoff)
					break;
				
			}
			duplicates(*it, occ);
			
			// record duplicate locations
			if (options.reportDuplicateLocations && options.duplicatesCutoff>=0 && occ <= options.duplicatesCutoff){
				for (TMapIter mit = begin(lociMap); mit != end(lociMap); ++mit){
					addDuplicate(*it, getSequenceNo(*mit), getPosition(*mit));
				}
			}
#ifdef TRIPLEX_DEBUG
				::std::cerr << " duplicates found for " << getSequenceNo(*it) << " : " << occ << ::std::endl;
#endif
			
		}
		
		if (options._debugLevel >= 1)
			options.logFileHandle << _getTimeStamp() << "   ... Finished counting duplicates" << ::std::endl;
		
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Filter a string with respect to the duplicate cutoff 
	template <typename TMotifSet>
	inline unsigned _filterDuplicatesWithCutoff(TMotifSet		&motifSet,
												Options const	&options
												)			
	{
		typedef typename Iterator<TMotifSet, Standard>::Type	TMotifSetIter;
		
		TMotifSet tmpMotifSet;
		
		for (TMotifSetIter it=begin(motifSet, Standard()); it != end(motifSet, Standard());++it){
			if (duplicates(*it) <= options.duplicatesCutoff){
				appendValue(tmpMotifSet, *it);
			}
		}
		unsigned removed = length(motifSet) - length(tmpMotifSet);
		
		// return if nothing has been removed
		if (removed == 0u){
			return removed;
		}
		
		// otherwise reset motifset with valid entries only
		clear(motifSet);
		for (TMotifSetIter it=begin(tmpMotifSet, Standard()); it != end(tmpMotifSet, Standard());++it){
			appendValue(motifSet, *it);
		}
		
		return removed;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Increase left handsite iterator 
	template< typename TNum, typename TArray, typename TPos>
	inline void _increaseLeft(TArray	&encoded_seq,
							  TPos		&pos,
							  TNum		&filter_chars,
							  TNum		&interrupt_chars,
							  TNum		&nonfilter_char
							  ){
		// increase leftmost pointer
		filter_chars -= encoded_seq[0][pos];
		interrupt_chars -= encoded_seq[1][pos];
		nonfilter_char -= encoded_seq[2][pos];
		++pos;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Increase right handsite iterator
	template< typename TNum, typename TArray, typename TPos>
	inline void _increaseRight(TArray	&encoded_seq,
							   TPos		&pos,
							   TNum		&filter_chars,
							   TNum		&interrupt_chars,
							   TNum		&nonfilter_char
							   ){
		// increase right pointer
		filter_chars += encoded_seq[0][pos];
		interrupt_chars += encoded_seq[1][pos];
		nonfilter_char += encoded_seq[2][pos];
		++pos;
	} 
	
	//////////////////////////////////////////////////////////////////////////////
	// Return whether pposition corresponds to interrupting character
	template<typename TArray, typename TPos>
	inline bool _isInterruptingChar(TArray	&encoded_seq,
									TPos	pos
									){
		return encoded_seq[1][pos];
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Encodes a sequence in a 2D array of dim [3][length(sequence)] where
	// where the first row contains a one at positions containing the filter char
	// the second row contains a 1 where a interrupt char 
	template <
	typename TString,
	typename TChar
	>
	inline bool** _encodeSeq(TString		&sequence,
							 TChar const	&filter_char,
							 TChar const	&interrupting_char,
							 bool**			blockruns, 
							 Options const	&options
							){
		typedef typename Iterator<TString>::Type	TIter;
		
		int rows = 3;
		int cols = length(sequence);

		// declaration
		bool** encoded_seq;
		
		// allocation
		encoded_seq = new bool*[rows];
		for(int i = 0; i < rows; i++){
			encoded_seq[i] = new bool[cols];
		}
		
		// initialization
		for(int i = 0; i < rows; ++i){
			for(int j = 0; j < cols; ++j){
				encoded_seq[i][j] = false;
			}
		}
		
		unsigned counter = 0;
		unsigned runCounter = 0;
		for (TIter it = begin(sequence); it != end(sequence); ++it, ++counter){
			if (*it == filter_char){
				encoded_seq[0][counter] = true;
			} else if (*it == interrupting_char){
				encoded_seq[1][counter] = true;
				encoded_seq[2][counter] = true;
				if (counter - runCounter >= options.minBlockRun){
					for (unsigned t1=0; t1 <= counter - options.minBlockRun; ++t1){
						for (unsigned t2=max(t1,runCounter) + options.minBlockRun; t2 <= length(sequence); ++t2){
							blockruns[t1][t2] = true;
						}
					}
				}
				runCounter = counter+1;
			} else {
				encoded_seq[2][counter] = true;
			}
		}
		// final segment
		if (counter - runCounter >= options.minBlockRun){
			for (unsigned t1=0; t1 <= counter - options.minBlockRun ; ++t1){
				for (unsigned t2=max(t1,runCounter) + options.minBlockRun; t2 <= length(sequence); ++t2){
					blockruns[t1][t2] = true;
				}
			}

		}
		return encoded_seq;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Filter a string with the requested guanine rate AND the error rate
	template <typename TMotifSet, typename TTag>
	inline int _filterWithGuanineAndErrorRate(TMotifSet							&patternString,
											  typename Value<TMotifSet>::Type	&pattern,
											  char								filter_char,
											  char								interrupting_char,
											  bool								reduceSet,
											  ORIENTATION const					orientation,
											  Options const						&options,
											  TTag const & 
											  ){
		typedef typename Value<TMotifSet>::Type			TPattern;
		typedef typename Host<TPattern>::Type			TString;
		typedef typename Iterator<TString>::Type		TStringIter;
		typedef ModifiedString<TString, ModView< FunctorRYFilter > >  	TFilter;
		typedef typename Iterator<TFilter>::Type		TFilterIter;
		typedef short									TEncode;
		
		TMotifSet tmp_pattern_set;
		TMotifSet* ptr_pattern_set;
		if (reduceSet){
			ptr_pattern_set = &tmp_pattern_set;
		} else {
			ptr_pattern_set = &patternString;
		}

		bool** encoded_seq;
		bool** blockruns; // flag indicating for each seq interval whether there exists a valid blockRun inbetween
		// allocation
		blockruns = new bool*[length(pattern)-options.minBlockRun+1];
		for(unsigned i = 0; i <= length(pattern)-options.minBlockRun; i++){
			blockruns[i] = new bool[length(pattern)+1];
		}
		// initialization
		for (unsigned t1=0; t1 <= length(pattern)-options.minBlockRun; ++t1){
			for (unsigned t2=0; t2 <= length(pattern); ++t2){
				blockruns[t1][t2] = false;
			}
		}

		// if no guanine rate restriction, then collapse filter and tolerated chars
		if (options.guanineRate <= 0.0){
			TFilter filter_seq(pattern);
			if (filter_char == 'G')
				filter_char='R';
			else
				filter_char='Y';
			encoded_seq = _encodeSeq(filter_seq, filter_char, interrupting_char, blockruns, options);
		} else {
			encoded_seq = _encodeSeq(pattern, filter_char, interrupting_char, blockruns, options);
		}
	
#ifdef TRIPLEX_DEBUG
		::std::cerr << pattern << ::std::endl;
		for (int r=0; r<3;++r){
			for (unsigned i=0; i<length(pattern) ; ++i){
				::std::cerr << encoded_seq[r][i] ;
			}
			::std::cerr << ::std::endl;
		}
		::std::cerr << ::std::endl;
		for (unsigned t1=0; t1 <= length(pattern)-options.minBlockRun; ++t1){
			for (unsigned t2=0; t2 <= length(pattern); ++t2){
				::std::cerr << blockruns[t1][t2];
			}
			::std::cerr  << ::std::endl;
		}
		
		::std::cerr << "Any hit: " << blockruns[0][length(pattern)] << " " << length(pattern) << ::std::endl;
#endif		
		
		double max_error = floor(length(pattern)*options.errorRate);
		if (options.maximalError >= 0)
			max_error = min(max_error, double(options.maximalError));
		double max_tolerated = (floor(length(pattern)*(1.0-options.guanineRate)));
		
		
		double cnt_filter_chars = 0.0;
		double cnt_interrupt_chars = 0.0;
		double cnt_nonfilter_chars = 0.0;
		bool is_match = false;
		unsigned max_length = length(pattern);
		if (options.maxLength >= options.minLength)
			max_length = options.maxLength;
		
		double tmp_error = 0.0;
		unsigned tmp_start = 0;
		unsigned tmp_end = 0;
		unsigned covered_end= 0;
		
		unsigned itLeft = 0;
		unsigned itRight = 0;
		
		// there must be another valid blockrun
		while (blockruns[itLeft][length(pattern)] && itLeft+options.minLength <= length(pattern)){
			// obey minimum length 
			while (itRight-itLeft < options.minLength && itRight < length(pattern) ){
				while (itRight-itLeft < options.minLength && itRight < length(pattern) ){
					_increaseRight(encoded_seq, itRight, cnt_filter_chars, cnt_interrupt_chars, cnt_nonfilter_chars);
				}
				
				// obey maximum error
				while (cnt_interrupt_chars > max_error)
					_increaseLeft(encoded_seq, itLeft, cnt_filter_chars, cnt_interrupt_chars, cnt_nonfilter_chars);
				
				// obey maximum filterchars 
				while (cnt_nonfilter_chars > max_tolerated)
					_increaseLeft(encoded_seq, itLeft, cnt_filter_chars, cnt_interrupt_chars, cnt_nonfilter_chars);
				
				// hit cannot start with interrupting char
				while (itLeft < length(pattern) && _isInterruptingChar(encoded_seq, itLeft))
					_increaseLeft(encoded_seq, itLeft, cnt_filter_chars, cnt_interrupt_chars, cnt_nonfilter_chars);

				if (itRight < itLeft){
					itRight = itLeft;
					cnt_filter_chars = 0.0;
					cnt_interrupt_chars = 0.0;
					cnt_nonfilter_chars = 0.0;	
				}
			}
			
			// cannot fulfil minimum length constraint
			if (itRight-itLeft < options.minLength){
				break;
				
			} else {
				is_match = false;
				// got a minimum length segment that does not violate maximum constraints
				// extend to the left as far as possible
				while (cnt_interrupt_chars <= max_error && cnt_nonfilter_chars <= max_tolerated && itRight-itLeft <= max_length){
					if (blockruns[itLeft][itRight] && !_isInterruptingChar(encoded_seq, itRight-1) 
						&& cnt_interrupt_chars/(itRight-itLeft) <= options.errorRate 
						&& cnt_filter_chars/(itRight-itLeft) >= options.guanineRate 
						&& _motifSpecificConstraints(double(itRight-itLeft), cnt_filter_chars, cnt_interrupt_chars, orientation, options, TTag()))
					{	
						is_match = true;
						tmp_start = itLeft;
						tmp_end = itRight;
						tmp_error = cnt_interrupt_chars;
						// add match straight away if all matches should be reported
						if (options.allMatches){
#ifdef TRIPLEX_DEBUG		
							::std::cerr << "add match:" << infix(pattern,tmp_start,tmp_end) << ::std::endl;
							double tmp_cnt_filter_chars = 0.0;
							double tmp_cnt_interrupt_chars = 0.0;
							double tmp_cnt_nonfilter_char = 0.0;
							unsigned it=tmp_start;
							while(it<tmp_end){
								_increaseRight(encoded_seq, it, tmp_cnt_filter_chars, tmp_cnt_interrupt_chars, tmp_cnt_nonfilter_char);
							}
							::std::cerr << "len : " << (tmp_end-tmp_start) << ::std::endl;
							::std::cerr << "err : " << (tmp_cnt_interrupt_chars/(tmp_end-tmp_start)) << ::std::endl;
							::std::cerr << "gua : " << (tmp_cnt_nonfilter_char/(tmp_end-tmp_start)) << ::std::endl;							
#endif	
							_addMatch(*ptr_pattern_set, pattern, tmp_start, tmp_end, tmp_error, TTag());
							covered_end = tmp_end;
							is_match = false; // prevent redundant addition
						}
					}
					
					if (itRight < length(pattern)){
						_increaseRight(encoded_seq, itRight, cnt_filter_chars, cnt_interrupt_chars, cnt_nonfilter_chars);
					} else {
						break;
					}
				}
				if (is_match && tmp_end > covered_end){
#ifdef TRIPLEX_DEBUG		
					::std::cerr << "add match:" << infix(pattern,tmp_start,tmp_end) << ::std::endl;
					double tmp_cnt_filter_chars = 0.0;
					double tmp_cnt_interrupt_chars = 0.0;
					double tmp_cnt_nonfilter_char = 0.0;
					unsigned it=tmp_start;
					while(it<tmp_end){
						_increaseRight(encoded_seq, it, tmp_cnt_filter_chars, tmp_cnt_interrupt_chars, tmp_cnt_nonfilter_char);
					}
					::std::cerr << "len : " << (tmp_end-tmp_start) << ::std::endl;
					::std::cerr << "err : " << (tmp_cnt_interrupt_chars/(tmp_end-tmp_start)) << ::std::endl;
					::std::cerr << "gua : " << (tmp_cnt_nonfilter_char/(tmp_end-tmp_start)) << ::std::endl;
#endif	
					_addMatch(*ptr_pattern_set, pattern, tmp_start, tmp_end, tmp_error, TTag());
					covered_end = tmp_end;
					// leave loop if there cannot be a subsequent hit
					if (tmp_end >= length(pattern))
						break;
				}
				is_match = false;
				// increase leftmost pointer & skip errors
				++itLeft;
				while (itLeft < length(pattern) && _isInterruptingChar(encoded_seq, itLeft))
					++itLeft;
				
				itRight = itLeft;
				cnt_filter_chars = 0.0;
				cnt_interrupt_chars = 0.0;
				cnt_nonfilter_chars = 0.0;
				while (!options.allMatches && itRight < tmp_end){
					_increaseRight(encoded_seq, itRight, cnt_filter_chars, cnt_interrupt_chars, cnt_nonfilter_chars);
				}					
			}
		}
		
		for (int r=0; r<3;++r)
			delete [] encoded_seq[r];
		delete [] encoded_seq;
		for (unsigned t1=0; t1<=length(pattern)-options.minBlockRun; ++t1)
			delete [] blockruns[t1];
		delete [] blockruns;
		
		// reduce motif set for triplex search
		if (reduceSet){
#ifdef TRIPLEX_DEBUG	
			::std::cerr << "# Elements before merging:" << length(tmp_pattern_set) << ::std::endl;
#endif			
			_reduceMotifSet(patternString, tmp_pattern_set);
		} 
#ifdef TRIPLEX_DEBUG
		::std::cerr << "# Elements final patterns:" << length(patternString) << ::std::endl;
#endif			
		return 0;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Processing motif specific constraints
	template <typename TSize, typename TTag>
	inline bool _motifSpecificConstraints(TSize const &length,
										  TSize const &filter,
										  TSize const &errors,
										  ORIENTATION const	&orientation,
										  Options const &options,
										  TTag
										  ){
#ifdef TRIPLEX_DEBUG
		::std::cerr << "length:" << length << " filter:" << filter << " errors:" << errors << " orientation:" << orientation << ::std::endl;
#endif
		return true;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Processing motif specific constraints
	template <typename TSize>
	inline bool _motifSpecificConstraints(TSize const &length,
										  TSize const &filter,
										  TSize const &errors,
										  ORIENTATION const	&orientation,
										  Options const &options,
										  MIXEDMOTIF
										  ){
#ifdef TRIPLEX_DEBUG
		::std::cerr << "length:" << length << " filter:" << filter << " errors:" << errors << " orientation:" << orientation << " qualifies" << ::std::endl;
#endif
		if (orientation != TRIPLEX_ORIENTATION_PARALLEL && filter/length >= options.mixed_antiparallel_min_guanine)
			return true;

		if (orientation != TRIPLEX_ORIENTATION_ANTIPARALLEL && filter/length <= options.mixed_parallel_max_guanine)
			return true;
		
#ifdef TRIPLEX_DEBUG
		::std::cerr << "NOT !" << ::std::endl;
#endif
		return false;
	}	
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Relay for adding new TFO match to the result set
	template <typename TMotifSet>
	inline void _addMatch(TMotifSet							&patternString,
						  typename Value<TMotifSet>::Type	&pattern,
						  unsigned							start,
						  unsigned							end,
						  unsigned							errors,
						  MIXEDMOTIF
						  ){
		_addMatch(patternString, pattern, start, end, errors, TFO());
	}

	//////////////////////////////////////////////////////////////////////////////
	// Relay for adding new TFO match to the result set
	template <typename TMotifSet>
	inline void _addMatch(TMotifSet							&patternString,
						  typename Value<TMotifSet>::Type	&pattern,
						  unsigned							start,
						  unsigned							end,
						  unsigned							errors,
						  PURINEMOTIF
						  ){
		_addMatch(patternString, pattern, start, end, errors, TFO());
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Relay for adding new TFO match to the result set
	template <typename TMotifSet>
	inline void _addMatch(TMotifSet							&patternString,
						  typename Value<TMotifSet>::Type	&pattern,
						  unsigned							start,
						  unsigned							end,
						  unsigned							errors,
						  PYRIMIDINEMOTIF
						  ){
		_addMatch(patternString, pattern, start, end, errors, TFO());
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Add a new TFO match to the result set
	template <typename TMotifSet>
	inline void _addMatch(TMotifSet							&patternString,
						  typename Value<TMotifSet>::Type	&pattern,
						  unsigned							start,
						  unsigned							end,
						  unsigned							errors,
						  TFO
						  ){
		typedef typename Value<TMotifSet>::Type	TPattern;
		typedef typename Position<TPattern>::Type	TPos;
		
		if (pattern.parallel){
			TPattern seqpattern(host(pattern),beginPosition(pattern)+start, beginPosition(pattern)+end, isParallel(pattern), getSequenceNo(pattern), isTFO(pattern), getMotif(pattern));
			setScore(seqpattern,end-start-errors);
			appendValue(patternString, seqpattern);
		} else{
			TPattern seqpattern(host(pattern),beginPosition(pattern)+(length(pattern)-end), beginPosition(pattern)+(length(pattern)-start), isParallel(pattern), getSequenceNo(pattern), isTFO(pattern), getMotif(pattern));
			setScore(seqpattern,end-start-errors);
			appendValue(patternString, seqpattern);
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Add a new TTS match to the result set
	template <typename TMotifSet>
	inline void _addMatch(TMotifSet						&patternString,
						  typename Value<TMotifSet>::Type	&pattern,
						  unsigned							start,
						  unsigned							end,
						  unsigned							errors,
						  TTS
						  ){
		typedef typename Value<TMotifSet>::Type	TPattern;
		typedef typename Position<TPattern>::Type	TPos;
		
		if (pattern.motif=='+'){
			TPattern seqpattern(host(pattern),beginPosition(pattern)+start, beginPosition(pattern)+end, isParallel(pattern), getSequenceNo(pattern), isTFO(pattern), getMotif(pattern));
			setScore(seqpattern,end-start-errors);
			appendValue(patternString, seqpattern);
		} else {
			TPattern seqpattern(host(pattern),beginPosition(pattern)+(length(pattern)-end), beginPosition(pattern)+(length(pattern)-start), isParallel(pattern), getSequenceNo(pattern), isTFO(pattern), getMotif(pattern));
			setScore(seqpattern,end-start-errors);
			appendValue(patternString, seqpattern);
		}
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	//  Merge motifs that overlap with respect to the underlying covered sequence
	//  Merging is performed based on IntervalTrees
	//	1.) all intervals are inserted into an interval tree
	//	2.) an undirected graph is constructed containing the intervals as nodes and
	//		edges if the intervals overlap (edge computation is done using the interval tree)
	//	3.) all connected components are calculated to detect overlapping intervals
	//	4.) the minimum and maximum index are determined as the new interval border
	//		for a set of overlapping intervals
	template <typename TMotifSet>
	void _reduceMotifSet(TMotifSet &reducedMotifSet,	// OUT
						 TMotifSet &motifSet)			// IN
	{
		typedef typename Value<TMotifSet>::Type					TOligoMotif;
		typedef typename Iterator<TMotifSet, Standard>::Type	TMotifSetIter;
		typedef typename Position<TOligoMotif>::Type			TIntervalValue;
		typedef IntervalTree<TIntervalValue,TVertexDescriptor>	TIntervalTree;
		
		typedef Graph<Undirected<void> >						TGraph;
		typedef VertexDescriptor<TGraph>::Type					TVertexDescriptor;
		typedef EdgeDescriptor<TGraph>::Type					TEdgeDescriptor;
		typedef typename Size<TGraph>::Type						TSize;
		typedef Iterator<TGraph, VertexIterator>::Type			TVertexIterator; 
		typedef IntervalAndCargo<TIntervalValue, TVertexDescriptor>	TInterval;
		typedef String<TVertexDescriptor>						TCounts;
		typedef typename Iterator<TCounts>::Type				TCountIter;
		typedef typename Iterator<String<unsigned int> >::Type	TCompIter;
		typedef Map<Pair<TVertexDescriptor, TOligoMotif*>, Skiplist<> > TPropMap;
		typedef String<TSize>									TComponent;
		typedef Map<Pair<TSize, TOligoMotif*>, Skiplist<> >		TCompMap;
		
		// if there is just one motif there is nothing to reduce
		if (length(motifSet) == 1){
			appendValue(reducedMotifSet, *begin(motifSet));
			return;
		}
		
		String<TInterval> intervals;
		resize(intervals, length(motifSet));
		int count = 0;
		
		// create graph from interval tree where intervals correspond to nodes,
		// edges correspond to direct overlaps and clusters of connected graphs
		// correspond to all direct or indirect connected graphs
		TGraph g;
		
		// fill the containers and create the nodes in the graph
		TPropMap vPropMap;
		
		for (TMotifSetIter it=begin(motifSet, Standard()); it != end(motifSet, Standard());++it,++count){
			TVertexDescriptor v = addVertex(g);
			insert(vPropMap, v, &(*it));
			// create interval and link vertex as cargo
			intervals[count].i1 = static_cast<TIntervalValue>(beginPosition(*it));
			intervals[count].i2 = static_cast<TIntervalValue>(endPosition(*it));
			intervals[count].cargo = v;
		}
		
		// create interval trees
		TIntervalTree itree(intervals);
		
		// add edges to the graph; edges correspond to overlaps of intervals 
		TVertexIterator itv(g); 
		while(!atEnd(itv)) { 
			TCounts itreeResult;
			TOligoMotif* motif = cargo(vPropMap, getValue(itv));
			findIntervals(itree, (TIntervalValue)beginPosition(*motif), (TIntervalValue)endPosition(*motif), itreeResult);
			for (TCountIter countIt=begin(itreeResult);countIt!=end(itreeResult);++countIt){
				// ignore trivial hit to same node
				if (getValue(itv) == *countIt)
					continue;
				// add edge between overlapping interval nodes
				addEdge(g,getValue(itv),*countIt);
			}
			++itv;
		}
		
		// compute connected components
		TComponent components;
		TSize num_compartments = connectedComponents(g, components);
		
		// iterate over all components and merge the intervals
		goBegin(itv); 
		TCompMap compMap;
		
		while(!atEnd(itv)) { 
			TVertexDescriptor v = getValue(itv);
			if (hasKey(compMap, getProperty(components, v))){
				// get stored motif
				TOligoMotif* oldmotif = cargo(compMap, getProperty(components, v));
				// extend motif boundaries if required
				merge(*oldmotif, *cargo(vPropMap, v));
			} else {
				// copy motif
				TOligoMotif* motif(cargo(vPropMap, v));
				// insert into compartment map
				insert(compMap, getProperty(components, v), motif);
			}
			++itv;
		}
		
		// add all new motifs into the reduced set
		for (TSize it=0; it < num_compartments;++it){
			appendValue(reducedMotifSet, *cargo(compMap, it));
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Parse a given string using an automaton adding the resulting segments to
	// the parse result object
	template <typename TString, typename TStringSegment, typename TGraph>
	inline int _parse(TStringSegment	&parseresult,	// the results of the parser
					  TGraph			&parser,		// parser object (automaton)
					  TString			&sequence,		// sequence to parse
					  Options const		&options)
	{
		typedef typename Value<TStringSegment>::Type			TSegment;
		typedef typename Iterator<TString, Standard>::Type		TIter;
		TIter it = begin(sequence,Standard());
		TIter run = begin(sequence,Standard());
		TIter itEnd = end(sequence,Standard());
		
		parseString(parser,getRoot(parser),run,itEnd);
		while (run != itEnd){
			// shift iterator back before the interruptions
			unsigned shift = min(options.maxInterruptions,static_cast<unsigned>(run-it));
			run -= shift;
			long size = max(run-run,run-it); // run-run = 0 for type compatility reasons
			
			if (size >= options.minLength){
				TSegment segment = infix(sequence, it, run);
				appendValue(parseresult, segment, Generous());
			}
			run += shift+1;
			it = run;
			if (run != itEnd)
				parseString(parser,getRoot(parser),run,itEnd);
		}
		// last entry
		unsigned size = max(itEnd-itEnd,itEnd-it);
		if (size >= options.minLength){
			TSegment segment = infix(sequence, it, itEnd);
			appendValue(parseresult, segment, Generous());
		}
		return 0;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Create a parser object using an automaton that allows segments with up to
	// options.maxInterruptions interruptions of invalid characters
	template <typename TGraph, typename TString>
	int _makeParser(TGraph			&parser, 	// the automaton to be created
					TString			&valids,	// the valid characters
					TString			&invalids,	// the invalid characters (interruptions)
					Options const	&options)
	{
		typedef typename Iterator<TString, Rooted>::Type		TIter;
		
		// create the root
		TVertexDescriptor root = addVertex(parser);
		assignRoot(parser, root);
		
		TIter itValid = begin(valids,Rooted());
		TIter itInvalid = begin(invalids,Rooted());
		TVertexDescriptor validChild = addVertex(parser);
		
		// add all valid edges create a child
		for (goBegin(itValid); !atEnd(itValid); goNext(itValid))
		{
			addEdge(parser, root, validChild , *itValid);
			addEdge(parser, validChild, validChild , *itValid);
		}
		
		TVertexDescriptor lastChild = validChild;
		TVertexDescriptor invalidChild;
		
		for (unsigned int i = 1; i<=options.maxInterruptions;++i){
			// add an interruption
			invalidChild = addVertex(parser);
			itInvalid = begin(invalids);
			// add all invalid edges
			for (goBegin(itInvalid); !atEnd(itInvalid); goNext(itInvalid)){
				addEdge(parser, lastChild, invalidChild , *itInvalid);
			}
			// add all valid edges
			itValid = begin(valids);
			for (goBegin(itValid); !atEnd(itValid); goNext(itValid)){
				addEdge(parser, invalidChild, validChild , *itValid);
			}
			lastChild = invalidChild;
		}
		
#ifdef TRIPLEX_DEBUG
		if (options._debugLevel > 2 ){
			::std::cerr << "Parser: vertices=" << numVertices(parser) << " edges=" << numEdges(parser) << ::std::endl;
//			write(::std::cerr, parser, DotDrawing());
		}
#endif
		
		return 0;
	}
		
	//////////////////////////////////////////////////////////////////////////////
	// Store all matches in a vector and convert the references accordingly
	template<
	typename TMatches,
	typename TId, 
	typename TGardenerSpec,
	typename TPattern,
	typename TStringSet
	>
	void _filterAndStore(TMatches						&matches, 
						 Gardener<TId, TGardenerSpec>	&gardener,
						 TPattern	const				&pattern,
						 TStringSet	const				&ttsSet, 
						 TId	const					&duplexId, 
						 bool	const					plusstrand,
						 Options						&options
						 ){
		typedef Gardener<TId, TGardenerSpec>				TGardener;
		typedef typename Value<TGardener>::Type				THitMap;
		typedef typename Value<THitMap>::Type				THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type		THitSetPointer;
		typedef typename Value<THitSetPointer>::Type		THitSet;
		typedef typename Value<THitSet>::Type				THit;
		typedef typename Iterator<THitSet>::Type			TIter;
		typedef typename Position<TMatch>::Type				TPos;
		typedef typename Value<TStringSet>::Type			TString;
		typedef typename Iterator<TStringSet>::Type			TStringIter;
		typedef typename Host<TString>::Type				THost;
		typedef typename Iterator<THost>::Type				THostIter;
		
		typedef typename Infix<THost>::Type					TSegment;
		typedef String<TSegment>							TSegString;
		typedef typename Iterator<TSegString>::Type			TSegStringIter;
		
		// check all queries for hits	
		for (TId queryid=0; queryid<(TId)length(ttsSet); ++queryid){
			for (TIter it = harvestBegin(gardener,queryid); it != harvestEnd(gardener, queryid); ++it){
				THit hit = *it;
				TPos tfoStart;
				TPos tfoEnd;
				TPos ttsStart;
				TPos ttsEnd;
				char strand;
				
				THost tfo = infix(ttsString(getSequenceByNo(hit.getNdlSeqNo(),needle(pattern))), hit.getNdlPos(), hit.getNdlPos() + hit.getHitLength());
				THost triplex(infix(ttsString(value(ttsSet,hit.getHstId())), hit.getHstkPos(), hit.getHstkPos()+hit.getHitLength()));
#ifdef TRIPLEX_DEBUG
				::std::cerr << "transform: " << triplex << :: std::endl;
#endif
				THostIter itTTS=begin(triplex);
				for (THostIter itTFO=begin(tfo); itTFO!=end(tfo); ++itTFO,++itTTS){
					if (*itTFO!=*itTTS){
						*itTTS = 'N';
					}	
				}
#ifdef TRIPLEX_DEBUG
				::std::cerr << "to       : " << triplex << :: std::endl;
#endif
				
				// run through TTS parser
 				TStringSet triplexSet;
				// create parser once
				if (empty(options.triplexParser)){
					THost valid("GAR");		// the valid characters
					THost invalid("TCYN");	// the interrupting characters			
					_makeParser(options.triplexParser, valid, invalid, options);
				}
				
				// split duplex into valid parts
				TSegString seqString;	// target segment container
				_parse(seqString, options.triplexParser, triplex, options);
				
				// process one segment at a time
				for (TSegStringIter it = begin(seqString, Standard()); it != end(seqString, Standard()); ++it){
#ifdef TRIPLEX_DEBUG
					::std::cerr << "pTriplex:" << *it << ::std::endl;
#endif
					TString ttsfilter(*it, true, hit.getHstId(), false, '+');
					bool reduceSet = false; // don't merge overlapping triplexes	
					_filterWithGuanineAndErrorRate(triplexSet, ttsfilter, 'G', 'Y', reduceSet, TRIPLEX_ORIENTATION_BOTH, options, TTS());
				}
				
				
				for (TStringIter itr=begin(triplexSet); itr!=end(triplexSet); ++itr){
					// compute score = matching positions, which can be found in the triplex string (number of N's within the interval found)
					int score=0;
					for (unsigned int i=beginPosition(*itr); i<endPosition(*itr); ++i){
						if (triplex[i]!='N')
							score += 1;
					}
					
					// calculate tts positions according to strand in the duplex
					if (plusstrand){
						ttsStart = hit.getHstkPos() + beginPosition(value(ttsSet,hit.getHstId()))+ beginPosition(*itr);
						ttsEnd = ttsStart + length(*itr);
						strand = '+';			
					} else {
						ttsEnd = endPosition(value(ttsSet,hit.getHstId())) - (hit.getHstkPos() + beginPosition(*itr));
						ttsStart = ttsEnd - length(*itr);
						strand = '-';
					}
					
					// calculate tfo positions according to binding orientation
					if (isParallel(getSequenceByNo(hit.getNdlSeqNo(),needle(pattern)))){
						tfoStart = hit.getNdlPos() + beginPosition(getSequenceByNo(hit.getNdlSeqNo(),needle(pattern)))+beginPosition(*itr);
						tfoEnd = tfoStart + length(*itr);
					} else {
						tfoEnd = endPosition(getSequenceByNo(hit.getNdlSeqNo(),needle(pattern))) - (hit.getNdlPos() + beginPosition(*itr));
						tfoStart = tfoEnd - length(*itr);

					}
					
					// save the corresponding triplex match 
					TMatch match(hit.getNdlSeqNo(),
								 tfoStart,
								 tfoEnd,
								 duplexId,
								 ttsStart,
								 ttsEnd,
								 score,
								 isParallel(getSequenceByNo(hit.getNdlSeqNo(),needle(pattern))),
								 getMotif(getSequenceByNo(hit.getNdlSeqNo(),needle(pattern))),
								 strand
								 );
					appendValue(matches, match);
					
#ifdef TRIPLEX_DEBUG
 					::std::cerr << "tts: " << ttsString(*itr) << " length: "<< length(*itr) <<  " position: "<< beginPosition(*itr) << "-" <<endPosition(*itr) << " " << infix(triplex, beginPosition(*itr), endPosition(*itr)) << ::std::endl;
					THost ftfo = infix(ttsString(getSequenceByNo(hit.getNdlSeqNo(),needle(pattern))), tfoStart, tfoEnd);
					THost ftts = infix(ttsString(value(ttsSet,hit.getHstId())), ttsStart, ttsEnd);
					::std::cerr << "tts: " << ftts << ::std::endl << "tfo: "<< ftfo << ::std::endl;
#endif
					
				}
			}
		}
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Take a sequence and generate a block queue (tolerated, matching, interrupting chars)
	// return the pointer to the queue after which the size constraint (minLength)
	// could never be fulfilled for a sequence starting behind that pointer
	// (essential datastructure for error-rate and guanine content filtering)	
	template <typename TString>
	inline int _makeBlockQueue(TString						&sequence,
							   vector<QueueBlock>			&queue,
							   Options const	&options
							   ){
		typedef typename Iterator<TString>::Type	TIter;
		TIter itStart = begin(sequence);
		TIter itEnd = end(sequence);
		unsigned queue_size_restricted = 0;
		char last = '\0';
		int count = 0;
		stack<int> blockSizes;
		
		if (itStart != itEnd){
			last = *itStart;
			count = 1;
			++itStart;
		}
		for (;itStart != itEnd;++itStart){
			if (last == (*itStart))
				count += 1;
			else{
				if (last != '>'){
					QueueBlock tb = { count, last };
					queue.push_back(tb);
					blockSizes.push(count);
				}
				last = *itStart;
				count = 1;
			}
		}
		// add last entry
		if (last != '>'){
			QueueBlock tb = { count, last };
			queue.push_back(tb);
			blockSizes.push(count);
		}
		
		queue_size_restricted = queue.size();
		unsigned tmp_sum=0;
		while (tmp_sum+blockSizes.top() < options.minLength){
			--queue_size_restricted;
			tmp_sum += blockSizes.top();
			blockSizes.pop();
		}
		return queue_size_restricted;
	}
	
	
	// extract the core file name
	template <typename TString>
	void _getShortFilename(TString &sfName, 	// short file name
						   TString const &fName 	// file name
						   ){
		// remove the directory prefix file
		::std::string _fName;
		assign(_fName, fName);
		size_t lastPos = _fName.find_last_of('/') + 1;
		if (lastPos == _fName.npos) lastPos = _fName.find_last_of('\\') + 1;
		if (lastPos == _fName.npos) lastPos = 0;
		sfName = _fName.substr(lastPos);
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Load multi-Fasta sequences
	template <typename TSequenceSet, typename TNameSet, typename TOptions>
	bool _loadOligos(TSequenceSet		&sequences,
				     TNameSet			&fastaIDs,
				     const char *		fileName,
				     TOptions const	&options)
	{
		
		MultiSeqFile multiFasta;
		if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
		
		AutoSeqFormat format;
		guessFormat(multiFasta.concat, format);
		split(multiFasta, format);
		
		unsigned seqCount = length(multiFasta);
		resize(sequences, seqCount, Exact());
		resize(fastaIDs, seqCount, Exact());
		
		TriplexString seq;
		for(unsigned i = 0; i < seqCount; ++i)
		{
			assignCroppedSeqId(fastaIDs[i], multiFasta[i], format);	// read Fasta id
			assignSeq(seq, multiFasta[i], format);				// read third strand
			assign(sequences[i], seq, Exact());
		}
		if (options._debugLevel > 1 )
			::std::cerr << "read " << length(sequences) << " sequences.\n";
		return (seqCount > 0);
	}
		
	//////////////////////////////////////////////////////////////////////////////
	// Find triplexes in many duplex sequences (import from Fasta) in parallel
	// by reading in all duplex sequences and storing the results on memory
	template <
	typename TMotifSet,
	typename TFile,
	typename TPattern,
	typename TId,
	typename TGardenerSpec
	>
	int inline startTriplexSearchSerial(TMotifSet					&tfoMotifSet,
										StringSet<CharString>		&tfoNames,
										TPattern const				&pattern,			 
										TFile						&outputfile,
										TId							duplexSeqNo,
										Options						&options,
										Gardener<TId, TGardenerSpec>
										){
		typedef TriplexString										TDuplex;
		typedef StringSet<ModStringTriplex<TDuplex, TDuplex> >		TDuplexModSet;
		typedef Gardener<TId, TGardenerSpec>						TGardener;
		typedef ::std::list<TMatch>									TMatches;
		typedef Repeat<unsigned, unsigned>						TRepeat;
		typedef std::list<TRepeat>								TRepeatString; //@TODO workaround for memory leak in seqan string
		typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
		
		// open duplex file
		::std::ifstream file;
		file.open(toCString(options.duplexFileNames[0]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return TRIPLEX_READFILE_FAILED;
		
		// remove the directory prefix of current duplex file
		::std::string duplexFile(toCString(options.duplexFileNames[0]));
		size_t lastPos = duplexFile.find_last_of('/') + 1;
		if (lastPos == duplexFile.npos) lastPos = duplexFile.find_last_of('\\') + 1;
		if (lastPos == duplexFile.npos) lastPos = 0;
		::std::string duplexName = duplexFile.substr(lastPos);
		TId duplexSeqNoWithinFile = 0;
		
		if (options._debugLevel >= 1)
			::std::cerr << "Starting on duplex file " << duplexName << ::std::endl;
		
		// iterate over duplex sequences
		for(; !_streamEOF(file); ++duplexSeqNo,++duplexSeqNoWithinFile){
			TMatches matches;
			TDuplex	duplexString;
			CharString duplexId;
			//readID(file, id, Fasta());			// read Fasta id
			readShortID(file, duplexId, Fasta());	// read Fasta id up to first whitespace
			if (options._debugLevel >= 2)
				::std::cerr << "Processing:\t" << duplexId << "\t(seq " << duplexSeqNoWithinFile << ")\r" << ::std::flush;
			
			read(file, duplexString, Fasta());			// read Fasta sequence

			// find low complexity regions and mask sequences if requested
			if (options.filterRepeats){
				TRepeatString	data_repeats;
				findRepeats(data_repeats, duplexString, options.minRepeatLength, options.maxRepeatPeriod);
				for (TRepeatIterator rbeg = begin(data_repeats); rbeg != end(data_repeats); ++rbeg){
					TRepeat repeat = *rbeg;
					CharString replacement = string(repeat.endPosition-repeat.beginPosition, 'N' );
					replace(duplexString, repeat.beginPosition, repeat.endPosition, replacement);
				}
			}
			
#if SEQAN_ENABLE_PARALLELISM	
			// run in parallel if requested and both strands are actually searched
			if (options.runtimeMode==RUN_PARALLEL_STRANDS && options.forward && options.reverse)
				_detectTriplexParallelStrands(matches, pattern, duplexString, duplexSeqNoWithinFile, options, TGardener());
				// otherwise go for serial processing
#endif	
			_detectTriplex(matches, pattern, duplexString, duplexSeqNoWithinFile, options, TGardener());
			
			// output all entries
			printTriplexEntry(matches, duplexId, duplexString, tfoMotifSet, tfoNames, outputfile, options);
			
			// clean up
			clear(matches);
		}
		file.close();
		return TRIPLEX_NORMAL_PROGAM_EXIT;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Find triplexes in many duplex sequences (import from Fasta) in parallel
	// by reading in all duplex sequences and storing the results on memory
	template <
	typename TMotifSet,
	typename TFile,
	typename TPattern,
	typename TId
	>
	int inline startTriplexSearchSerial(TMotifSet					&tfoMotifSet,
										StringSet<CharString>		&tfoNames,
										TPattern const				&pattern,			 
										TFile						&outputfile,
										TId							duplexSeqNo,
										Options						&options,
										BruteForce
										){
		typedef TriplexString										TDuplex;
		typedef StringSet<ModStringTriplex<TDuplex, TDuplex> >		TDuplexModSet;
		typedef ::std::list<TMatch>									TMatches;
		typedef Repeat<unsigned, unsigned>						TRepeat;
		typedef std::list<TRepeat>								TRepeatString; //@TODO workaround for memory leak in seqan string
		typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
		
		// open duplex file
		::std::ifstream file;
		file.open(toCString(options.duplexFileNames[0]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return TRIPLEX_READFILE_FAILED;
		
		// remove the directory prefix of current duplex file
		::std::string duplexFile(toCString(options.duplexFileNames[0]));
		size_t lastPos = duplexFile.find_last_of('/') + 1;
		if (lastPos == duplexFile.npos) lastPos = duplexFile.find_last_of('\\') + 1;
		if (lastPos == duplexFile.npos) lastPos = 0;
		::std::string duplexName = duplexFile.substr(lastPos);
		TId duplexSeqNoWithinFile = 0;
		
		if (options._debugLevel >= 1)
			::std::cerr << "Starting on duplex file " << duplexName << ::std::endl;
		
		// iterate over duplex sequences
		for(; !_streamEOF(file); ++duplexSeqNo,++duplexSeqNoWithinFile){
			TMatches matches;
			TDuplex	duplexString;
			CharString duplexId;
			readShortID(file, duplexId, Fasta());	// read Fasta id up to first whitespace
			if (options._debugLevel >= 2)
				::std::cerr << "Processing:\t" << duplexId << "\t(seq " << duplexSeqNoWithinFile << ")\r" << ::std::flush;
			
			read(file, duplexString, Fasta());			// read Fasta sequence
			
			// find low complexity regions and mask sequences if requested
			if (options.filterRepeats){
				// find low complexity regions and mask sequences if requested
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
			
#if SEQAN_ENABLE_PARALLELISM	
			// run in parallel if requested and both strands are actually searched
			if (options.runtimeMode==RUN_PARALLEL_STRANDS && options.forward && options.reverse)
				_detectTriplexParallelStrands(matches, tfoMotifSet, duplexString, duplexSeqNoWithinFile, options, BruteForce());
			// otherwise go for serial processing
#endif	
			_detectTriplex(matches, tfoMotifSet, duplexString, duplexSeqNoWithinFile, options, BruteForce());
			
			// output all entries
			printTriplexEntry(matches, duplexId, duplexString, tfoMotifSet, tfoNames, outputfile, options);
			
			// clean up
			clear(matches);
		}
		file.close();
		return TRIPLEX_NORMAL_PROGAM_EXIT;
	}
	
#if SEQAN_ENABLE_PARALLELISM	
	//////////////////////////////////////////////////////////////////////////////
	// Find triplexes in many duplex sequences (import from Fasta) in parallel
	// by reading in all duplex sequences and storing the results on memory
	template <
	typename TMotifSet,
	typename TFile,
	typename TPattern,
	typename TId,
	typename TTag
	>
	int inline startTriplexSearchParallelDuplex(TMotifSet					&tfoMotifSet,
												StringSet<CharString> const	&tfoNames,
												TPattern const				&pattern,			 
												TFile						&outputfile,
												TId							duplexSeqNo,
												Options		&options,
												TTag
												){
		typedef TriplexString									TDuplex;
		typedef StringSet<ModStringTriplex<TDuplex, TDuplex> >	TDuplexModSet;
		typedef ::std::list<TMatch>								TMatches;
		typedef Triple<TId, CharString, TDuplex>				TSeq;
		typedef ::std::vector<TSeq>								TDataContainer;
		typedef Repeat<unsigned, unsigned>						TRepeat;
		typedef std::list<TRepeat>								TRepeatString; //@TODO workaround for memory leak in seqan string
		typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
		
		TDataContainer data;
		
		// open duplex file and read all fasta files
		::std::ifstream file;
		file.open(toCString(options.duplexFileNames[0]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return TRIPLEX_READFILE_FAILED;
		
		// remove the directory prefix of current duplex file
		::std::string duplexFile(toCString(options.duplexFileNames[0]));
		size_t lastPos = duplexFile.find_last_of('/') + 1;
		if (lastPos == duplexFile.npos) lastPos = duplexFile.find_last_of('\\') + 1;
		if (lastPos == duplexFile.npos) lastPos = 0;
		::std::string duplexName = duplexFile.substr(lastPos);
		TId duplexSeqNoWithinFile = 0;
		
		if (options._debugLevel >= 1)
			::std::cerr << "Starting on duplex file " << duplexName << ::std::endl;
		
		// iterate over duplex sequences
		SEQAN_PROTIMESTART(find_time);
		if (options._debugLevel >= 1)
			::std::cerr << "Reading: duplex sequences" << ::std::endl;
		
		// read duplex data
		options.logFileHandle << _getTimeStamp() << " * Reading all sequences " << ::std::endl;
		for(; !_streamEOF(file); ++duplexSeqNo,++duplexSeqNoWithinFile){
			TDuplex	duplexString;
			CharString duplexId;
			readShortID(file, duplexId, Fasta());	// read Fasta id up to first whitespace
			read(file, duplexString, Fasta());		// read Fasta sequence
			
			// find low complexity regions and mask sequences if requested
			if (options.filterRepeats){
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
			
			TSeq seq(duplexSeqNoWithinFile, duplexId, duplexString);
			appendValue(data, seq);

		}
		
		// search
		if (options._debugLevel >= 1)
			::std::cerr << "Searching triplexes\r" << ::std::endl;
		
		// run parallel
		_invokeParallelSequenceProcessing(data, pattern, tfoMotifSet, tfoNames, outputfile, options, TTag());
		
		
		if (options._debugLevel >= 1)
			::std::cerr << "Outputting results\r" << ::std::endl;
		
		options.timeFindTriplexes += SEQAN_PROTIMEDIFF(find_time);
		file.close();
		return TRIPLEX_NORMAL_PROGAM_EXIT;
	}
#endif 
	
	
#if SEQAN_ENABLE_PARALLELISM	
	template <
	typename TDataContainer,
	typename TPattern,
	typename TMotifSet,
	typename TFile,
	typename TId,
	typename TGardenerSpec
	>
	inline void _invokeParallelSequenceProcessing(TDataContainer				&data,
												  TPattern const				&pattern,
												  TMotifSet						&tfoSet,
												  StringSet<CharString>	const	&tfoNames,
												  TFile							&outputfile,
												  Options						&options,
												  Gardener<TId, TGardenerSpec>
												  ){
		typedef typename Value<TDataContainer>::Type	TSeq;
		typedef typename Value<TSeq,1>::Type			TCounter;
		typedef typename Value<TSeq,2>::Type			TString;
		typedef typename Value<TSeq,3>::Type			TDuplex;
		typedef Gardener<TId, TGardenerSpec>			TGardener;
		typedef ::std::list<TMatch>						TMatches;
		
		// parallel section 
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel)
		{
			SEQAN_PRAGMA_IF_PARALLEL(omp for schedule(dynamic) )
			for (int i=0; i<(int)length(data);++i){
				TMatches matches;
				TSeq seq = data[i];
				TCounter duplexCounter = seq.i1;
				TString duplexId = seq.i2;
				TDuplex duplex = seq.i3;
				
				_detectTriplex(matches, pattern, duplex, duplexCounter, options, TGardener());

				if (length(matches)>0){
					SEQAN_PRAGMA_IF_PARALLEL(omp critical(printTriplexEntry) )
					printTriplexEntry(matches, duplexId, duplex, tfoSet, tfoNames, outputfile, options);
				}
			}
		}
	}
#endif
	
#if SEQAN_ENABLE_PARALLELISM	
	template <
	typename TDataContainer,
	typename TPattern,
	typename TMotifSet,
	typename TFile
	>
	inline void _invokeParallelSequenceProcessing(TDataContainer				&data,
												  TPattern const				&pattern,
												  TMotifSet						&tfoSet,
												  StringSet<CharString>	const	&tfoNames,
												  TFile							&outputfile,
												  Options						&options,
												  BruteForce
												  ){
		typedef typename Value<TDataContainer>::Type	TSeq;
		typedef typename Value<TSeq,1>::Type			TCounter;
		typedef typename Value<TSeq,2>::Type			TString;
		typedef typename Value<TSeq,3>::Type			TDuplex;
		typedef ::std::list<TMatch>						TMatches;
		
		// parallel section 
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel)
		{
			SEQAN_PRAGMA_IF_PARALLEL(omp for schedule(dynamic) )
			for (int i=0; i<(int)length(data);++i){
				TMatches matches;
				TSeq seq = data[i];
				TCounter duplexCounter = seq.i1;
				TString duplexId = seq.i2;
				TDuplex duplex = seq.i3;
				
				_detectTriplex(matches, tfoSet, duplex, duplexCounter, options, BruteForce());
				
				if (length(matches)>0){
					SEQAN_PRAGMA_IF_PARALLEL(omp critical(printTriplexEntry) )
					printTriplexEntry(matches, duplexId, duplex, tfoSet, tfoNames, outputfile, options);
				}
			}
		}
	}
#endif



}
#endif  // #ifndef SANDBOX_FBUSKE_APPS_TRIPLEXATOR_TRIPLEX_H_
