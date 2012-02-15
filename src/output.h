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

#ifndef SEQAN_HEADER_OUTPUT_FORMAT_H
#define SEQAN_HEADER_OUTPUT_FORMAT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>

#include "triplexator.h"
#include "triplex_alphabet.h"
#include "triplex_pattern.h"
#include <seqan/align.h>
#include <seqan/misc/priority_type_heap.h>


#ifdef BOOST
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
namespace io = boost::iostreams;
#endif

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// Less-operators ...	
	template <typename TTriplexMatch>
	struct PotentialComparator : public ::std::binary_function < TTriplexMatch, TTriplexMatch, bool >
	{
		inline bool operator() (TTriplexMatch const &a, TTriplexMatch const &b) const
		{
			// tts position and orientation
			if (a.ttsSeqNo < b.ttsSeqNo) return true;
			if (a.ttsSeqNo > b.ttsSeqNo) return false;
			// tfo seq number
			if (a.tfoNo < b.tfoNo) return true;
			if (a.tfoNo > b.tfoNo) return false;
			// diagonal
			if (a.dBegin-a.oBegin < b.dBegin-b.oBegin) return true;
			if (a.dBegin-a.oBegin > b.dBegin-b.oBegin) return false;
			// tts start position
			if (a.dBegin < b.dBegin) return true;
			if (a.dBegin > b.dBegin) return false;
			// tts end position
			if (a.dEnd > b.dEnd) return true;
			if (a.dEnd < b.dEnd) return false;
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
		
	template <typename TModString>
	struct ModTriplexPotentialComparator : public ::std::binary_function < TModString, TModString, bool >
	{
		inline bool operator() (TModString const &a, TModString const &b) const
		{
			// seqNo
			if (getSequenceNo(a) < getSequenceNo(b)) return true;
			if (getSequenceNo(a) > getSequenceNo(b)) return false;
			// beginPosition
			if (beginPosition(a) < beginPosition(b)) return true;
			if (beginPosition(a) > beginPosition(b)) return false;
			// endPosition
			if (endPosition(a) > endPosition(b)) return true;
			if (endPosition(a) < endPosition(b)) return false;
			// isParallel
			if (isParallel(a) < isParallel(b)) return true;
			if (isParallel(a) > isParallel(b)) return false;
			// getMotif
			if (getMotif(a) < getMotif(b)) return true;
			if (getMotif(a) > getMotif(b)) return false;
			// score
			return score(a) > score(b);
		}
	};

	
	//////////////////////////////////////////////////////////////////////////////
	// prepare output
	template <typename TFile>
	void openOutputFile(TFile 			&filehandle,		// file handle
						Options const	&options)
	{
		// create output file
		CharString fileName = options.outputFolder;
		CharString tmp = "tmp_";
		if (!empty(options.output) && options.outputFormat!=2){
			append(tmp, options.output);
			append(fileName,tmp);
			if (options._debugLevel >= 1)		
				::std::cerr << "open " << fileName << ::std::endl;
			
			filehandle.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::trunc);
			if (!filehandle.is_open()) {
				::std::cerr << "Failed to open temporary output file:" << fileName << ::std::endl;
				return;
			}
		}
	}

#ifdef BOOST
	void openOutputFile(io::filtering_ostream		&filterstream,
						Options const	&options)
	{
		// create output file
		CharString fileName = options.outputFolder;
		CharString tmp = "tmp_";
		if (!empty(options.output) && options.outputFormat!=2){
			append(tmp, options.output);
			append(fileName,tmp);
			append(fileName,".gz");
			
			if (options._debugLevel >= 1)		
				::std::cerr << "open " << fileName << ::std::endl;
			filterstream.push(io::file_sink(toCString(fileName), std::ios::binary));
		}
	}
#endif
	
	//////////////////////////////////////////////////////////////////////////////
	// finish output
	template <typename TFile>
	inline int finishOutputFile(TFile			&filehandle,		// file handle
						 Options const	&options)
	{
		// rename temporary file to final location
		CharString workFileName = options.outputFolder;
		CharString fileName = options.outputFolder;
		if (!empty(options.output) && options.outputFormat!=2){
			CharString tmp = "tmp_";
			append(tmp,options.output);
			append(workFileName,tmp);
			append(fileName,options.output);
			if (filehandle.is_open()){
				filehandle.close();
			}
			// remove existing file first
			if (options._debugLevel >= 1)		
				::std::cerr << "rename temorary file " << workFileName << " to " << fileName << ::std::endl;
			::std::remove(toCString(fileName));
			// rename tmporary file to final destination
			int result = ::std::rename( toCString(workFileName) , toCString(fileName) );
			if (result != 0){
				::std::cerr << "Failed to rename output file " << workFileName << " to " << fileName << ::std::endl;
				return 1;
			}
		} 
		return 0;
	}
	
#ifdef BOOST	
	inline int finishOutputFile(io::filtering_ostream	&filterstream,		
						 Options const			&options)
	{
		// rename temporary file to final location
		CharString workFileName = options.outputFolder;
		CharString fileName = options.outputFolder;
		if (!empty(options.output) && options.outputFormat!=2){
			CharString tmp = "tmp_";
			append(tmp,options.output);
			append(tmp,".gz");
			append(workFileName,tmp);
			append(fileName,options.output);
			append(fileName,".gz");
			
			filterstream.flush();
			
			// remove existing file first
			if (options._debugLevel >= 1)		
				::std::cerr << "rename temorary file " << workFileName << " to " << fileName << ::std::endl;
			::std::remove(toCString(fileName));
			// rename tmporary file to final destination
			int result = ::std::rename( toCString(workFileName) , toCString(fileName) );
			if (result != 0){
				::std::cerr << "Failed to rename output file " << workFileName << " to " << fileName << ::std::endl;
				return 1;
			}
		} 
		return 0;
	}
#endif
	
	//////////////////////////////////////////////////////////////////////////////
	// close and finish output file
	template <typename TFile>
	inline void closeOutputFile(TFile	&filehandle,		// file handle
								Options	&options)
	{
		if (!empty(options.output) && options.outputFormat!=2){
			// close output file
			if (filehandle.is_open()){
				filehandle.close();
			}
			finishOutputFile(filehandle, options);
		}
	}
	
#ifdef BOOST	
	inline void closeOutputFile(io::filtering_ostream	&filterstream,
								Options					&options)
	{
		if (!empty(options.output) && options.outputFormat!=2){
			// flush filterstream
			filterstream.flush();		
			finishOutputFile(filterstream, options);
		}
	}
#endif
	
	//////////////////////////////////////////////////////////////////////////////
	// prepare output
	void openLogFile(Options &options)
	{
		// create output file
		CharString fileName = options.outputFolder;
		CharString tmp = "tmp_";
		if (!empty(options.output)){
			append(tmp, options.output);
			append(tmp, ".log");
			append(fileName,tmp);
		} else {
			append(fileName,tmp);
			append(fileName,options.logFileName);
		}
		options.logFileHandle.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::trunc);
		if (!options.logFileHandle.is_open()) {
			::std::cerr << "Failed to create log file:" << fileName << ::std::endl;
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// finish log file
	inline int finishLogFile(Options &options)
	{
		// rename temporary file to final location
		CharString workFileName = options.outputFolder;
		CharString fileName = options.outputFolder;
		CharString tmp = "tmp_";
		if (!empty(options.output)){
			append(tmp, options.output);
			append(tmp, ".log");
			append(workFileName,tmp);
			append(fileName,options.output);
			append(fileName,".log");
		} else {
			append(workFileName,tmp);
			append(workFileName,options.logFileName);
			append(fileName,options.logFileName);

		}
		if (options.logFileHandle.is_open()){
			::std::cerr << "File still open. Renaming aborted: " << workFileName << ::std::endl;
			return 1;
		}
		// remove existing file first
		::std::remove(toCString(fileName));
		// rename tmporary file to final destination
		int result = ::std::rename( toCString(workFileName) , toCString(fileName) );
		if (result != 0){
			::std::cerr << "Failed to rename output file " << workFileName << " to " << fileName << ::std::endl;
			return 1;
		}
		return 0;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	// close logfile
	void closeLogFile(Options &options)
	{
		// close output file
		if (options.logFileHandle.is_open()){
			options.logFileHandle.close();
			finishLogFile(options);
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// prepare summary file
	void openSummaryFile(Options &options)
	{
		// create output file
		CharString fileName = options.outputFolder;
		CharString tmp = "tmp_";
		if (!empty(options.output)){
			append(tmp, options.output);
			append(tmp, ".summary");
			append(fileName,tmp);
		} else {
			append(fileName,tmp);
			append(fileName,options.summaryFileName);
		}
		options.summaryFileHandle.open(toCString(fileName), ::std::ios_base::out | ::std::ios_base::trunc);
		if (!options.summaryFileHandle.is_open()) {
			::std::cerr << "Failed to create temporary summary file:" << fileName << ::std::endl;
		}
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	// finish summary file
	inline int finishSummaryFile(Options &options)
	{
		// rename temporary file to final location
		CharString workFileName = options.outputFolder;
		CharString fileName = options.outputFolder;
		CharString tmp = "tmp_";
		if (!empty(options.output)){
			append(tmp, options.output);
			append(tmp, ".summary");
			append(workFileName,tmp);
			append(fileName,options.output);
			append(fileName,".summary");
		} else {
			append(workFileName,tmp);
			append(workFileName,options.summaryFileName);
			append(fileName,options.summaryFileName);
		}
		if (options.summaryFileHandle.is_open()){
			::std::cerr << "File still open. Renaming aborted: " << workFileName << ::std::endl;
			return 1;
		}
		// remove existing file first
		::std::remove(toCString(fileName));
		// rename tmporary file to final destination
		int result = ::std::rename( toCString(workFileName) , toCString(fileName) );
		if (result != 0){
			::std::cerr << "Failed to rename output file " << workFileName << " to " << fileName << ::std::endl;
			return 1;
		}
		return 0;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// close summary file
	void closeSummaryFile(Options &options)
	{
		// close output file
		if (options.summaryFileHandle.is_open()){
			options.summaryFileHandle.close();
			finishSummaryFile(options);
		}
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	// print header for triplex file
	template <typename TFile>
	void printTriplexHeader(TFile	&filehandle,		// file handle
							Options	&options)
	{
		char _sep_ = '\t';
		switch (options.outputFormat)
		{
			case 0:	// brief Triplex Format
				filehandle << "# Sequence-ID" << _sep_ << "TFO start" << _sep_ << "TFO end" << _sep_ << "Duplex-ID" << _sep_ << "TTS start" << _sep_ << "TTS end" << _sep_ << "Score" << _sep_ << "Error-rate" << _sep_ << "Errors" << _sep_ << "Motif" << _sep_ << "Strand" << _sep_ << "Orientation" << _sep_ << "Guanine-rate" << ::std::endl;
				break;
			case 1:	// brief Triplex Format
				filehandle << "# Sequence-ID" << _sep_ << "TFO start" << _sep_ << "TFO end" << _sep_ << "Duplex-ID" << _sep_ << "TTS start" << _sep_ << "TTS end" << _sep_ << "Score" << _sep_ << "Error-rate" << _sep_ << "Errors" << _sep_ << "Motif" << _sep_ << "Strand" << _sep_ << "Orientation" << _sep_ << "Guanine-rate" << ::std::endl;
				break;	
			default:
				break;
		}
		// summary file 
		options.summaryFileHandle << "# Duplex-ID" << _sep_ << "Sequence-ID" << _sep_ << "Total (abs)" << _sep_ << "Total (rel)" << _sep_ << "GA (abs)" << _sep_ << "GA (rel)" << _sep_ << "TC (abs)" << _sep_ << "TC (rel)" << _sep_ << "GT (abs)" << _sep_ << "GT (rel)" << ::std::endl;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// prepare tfo output
	template <typename TFile>
	void printTFOHeader(TFile	&filehandle,		// file handle
						Options	&options)
	{
		char _sep_ = '\t';
		switch (options.outputFormat)
		{
			case 0:	// TFO Format
				filehandle << "# Sequence-ID" << _sep_ << "Start" << _sep_ << "End" << _sep_ << "Score" << _sep_ << "Motif" << _sep_ << "Error-rate" << _sep_ << "Errors"  << _sep_ << "Guanine-rate" << _sep_ <<  "Duplicates"  << _sep_ <<  "TFO"  << _sep_ <<  "Duplicate locations" << ::std::endl;
				break;
			case 1:	// FASTA
				break;
			default:
				break;				
		}
		// summary file 
		options.summaryFileHandle << "# Sequence-ID" << _sep_ << "TFOs (abs)" << _sep_ << "TFOs (rel)" << _sep_ << "GA (abs)" << _sep_ << "GA (rel)" << _sep_ << "TC (abs)" << _sep_ << "TC (rel)" << _sep_ << "GT (abs)"  << _sep_ << "GT (rel)" << ::std::endl;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	// prepare tts output
	template <typename TFile>
	void printTTSHeader(TFile		&filehandle,		// file handle
						Options 	&options)
	{
		char _sep_ = '\t';
		switch (options.outputFormat)
		{
			case 0:	// TTS Format
				filehandle << "# Duplex-ID" << _sep_ << "Start" << _sep_ << "End" << _sep_ << "Score" << _sep_ << "Strand" << _sep_ << "Error-rate" << _sep_ << "Errors"  << _sep_ << "Guanine-rate" << _sep_ <<  "Duplicates"  << _sep_ << "TTS"  << _sep_ << "Duplicate locations"  << ::std::endl;
				break;
			case 1:	// FASTA
				break;
			default:
				break;				
		}
		// summary file 
		options.summaryFileHandle << "# Duplex-ID" << _sep_ << "TTSs (abs)" << _sep_ << "TTSs (rel)" << ::std::endl;
	}

	
	//////////////////////////////////////////////////////////////////////////////
	// Output triplex match alignments
	template <
	typename TMatch,
	typename TString,
	typename TMotifSet,
	typename TFile
	>
	void dumpAlignment(TMatch		&match,				// forward/reverse matches
					   TString		&duplex,			// duplex string
					   TMotifSet	&tfoSet,
					   TFile		&filehandle,
					   Options		&options)
	{	
		if (options.outputFormat != FORMAT_TRIPLEX)
			return;
		
		typedef typename Iterator<TString, Standard>::Type	TIter;
		typedef typename Value<TMotifSet>::Type				TMotif;
		TMotif tts_(duplex, match.dBegin, match.dEnd, match.parallel, match.ttsSeqNo, false, match.strand);		
		TMotif tfo_(host(value(tfoSet,match.tfoNo)), match.oBegin, match.oEnd, match.parallel, value(tfoSet,match.tfoNo).seqNo, true, match.motif);		

		CharString psTFO = prettyString(tfo_);
		CharString psTTS = prettyString(tts_);
		TString opp(psTTS);
		complement(opp);
		
		if (match.strand == '-'){
			reverse(opp);
			reverse(psTTS);
			filehandle << "     5'- " << opp << " -3'" << ::std::endl;
			filehandle << "TTS: 3'- " << psTTS << " -5'" << ::std::endl;
			TIter itTts = begin(tts_);
			TIter itTtsEnd = end(tts_);
			TIter itTfo = begin(tfo_);
			TIter itTfoEnd = end(tfo_);
			filehandle << "         ";
			while(itTtsEnd != itTts and itTfoEnd != itTfo){
				--itTtsEnd;
				--itTfoEnd;
				if (*itTtsEnd == *itTfoEnd)
					filehandle << "|";
				else
					filehandle << "*";
			}
			filehandle << ::std::endl;
			if (!value(tfoSet,match.tfoNo).parallel){
				filehandle << "TFO: 5'- " << psTFO << " -3'" << ::std::endl;
			} else {
				reverse(psTFO);
				filehandle << "TFO: 3'- " << psTFO << " -5'" << ::std::endl;
			}
		} else { // '+' strand
			if (!value(tfoSet,match.tfoNo).parallel){
				reverse(psTFO);
				filehandle << "TFO: 3'- " << psTFO << " -5'" << ::std::endl;
			} else {
				filehandle << "TFO: 5'- " << psTFO << " -3'" << ::std::endl;
			}
			TIter itTts = begin(tts_);
			TIter itTtsEnd = end(tts_);
			TIter itTfo = begin(tfo_);
			TIter itTfoEnd = end(tfo_);
			filehandle << "         ";
			while(itTts != itTtsEnd and itTfo != itTfoEnd){
				if (*itTts == *itTfo)
					filehandle << "|";
				else
					filehandle << "*";
				++itTts;
				++itTfo;
			}
			filehandle << ::std::endl;
			filehandle << "TTS: 5'- " << psTTS <<  " -3'" << ::std::endl;
			filehandle << "     3'- " << opp <<  " -5'" << ::std::endl;
		}
		filehandle << ::std::endl;
	}
		
	
	//////////////////////////////////////////////////////////////////////////////
	// get error string for triplex matches
	template <
	typename TMatch,
	typename TString,
	typename TMotifSet
	>
	CharString _errorString(TMatch		&match,				// forward/reverse matches
							TString		&duplex,			// duplex string
							TMotifSet 	&tfoSet,
							Options		&options)
	{	
		std::ostringstream errors;
		
		typedef typename Iterator<TString, Standard>::Type	TIter;
		typedef typename Value<TMotifSet>::Type				TMotif;
		TMotif tts_(duplex, match.dBegin, match.dEnd, match.parallel, match.ttsSeqNo, false, match.strand);		
		TMotif tfo_(host(value(tfoSet,match.tfoNo)), match.oBegin, match.oEnd, match.parallel, value(tfoSet,match.tfoNo).seqNo, true, match.motif);		
		
		CharString psTFO = prettyString(tfo_);
		CharString psTTS = prettyString(tts_);
				
		if (options.errorReference == WATSON_STAND){
			// requires consideration of pruine tract and parallel/anti-parallel triplex formation
			if (match.strand == '-'){
				reverse(psTTS);
				if (value(tfoSet,match.tfoNo).parallel) 
					reverse(psTFO);
//				::std::cerr << "\n" << psTTS << "\n" << psTFO << "\n" << tts_ << "\n" << tfo_ << "\n";
				TIter itTts = begin(tts_);
				TIter itTtsEnd = end(tts_);
				TIter itTfo = begin(tfo_);
				TIter itTfoEnd = end(tfo_);
				unsigned i = 0;
				while(itTtsEnd != itTts and itTfoEnd != itTfo){
					--itTtsEnd;
					--itTfoEnd;
					if (*itTtsEnd != *itTfoEnd){
						if (! isupper(value(psTTS,i)) && ! isupper(value(psTFO,i))) errors << 'b' << i;
						else if (! isupper(value(psTTS,i))) errors << 'd' << i;
						else if (! isupper(value(psTFO,i))) errors << 'o' << i;	
						else errors << 't' << i;
					}
					++i;
				}			
			} else { // '+' strand
				if (!value(tfoSet,match.tfoNo).parallel)
					reverse(psTFO);
//				::std::cerr << "\n" << psTTS << "\n" << psTFO << "\n" << tts_ << "\n" << tfo_ << "\n";
				TIter itTts = begin(tts_);
				TIter itTtsEnd = end(tts_);
				TIter itTfo = begin(tfo_);
				TIter itTfoEnd = end(tfo_);
				unsigned i = 0;
				while(itTts != itTtsEnd && itTfo != itTfoEnd){
					if (*itTts != *itTfo){
						if (! isupper(value(psTTS,i)) && ! isupper(value(psTFO,i))) errors << 'b' << i;
						else if (! isupper(value(psTTS,i))) errors << 'd' << i;
						else if (! isupper(value(psTFO,i))) errors << 'o' << i;		
						else errors << 't' << i;
					}
					++itTts;
					++itTfo;
					++i;
				}
			}
		} else if (options.errorReference == PURINE_STRAND){ 
			// only requires adjustment of anti-parallel binding third strands
			if (!value(tfoSet,match.tfoNo).parallel) 
				reverse(psTFO);
//			::std::cerr << "\n" << psTTS << "\n" << psTFO << "\n" << tts_ << "\n" << tfo_ << "\n";
			
			TIter itTts = begin(tts_);
			TIter itTtsEnd = end(tts_);
			TIter itTfo = begin(tfo_);
			TIter itTfoEnd = end(tfo_);
			unsigned i = 0;
			while(itTts != itTtsEnd && itTfo != itTfoEnd){
				if (*itTts != *itTfo){
					if (! isupper(value(psTTS,i)) && ! isupper(value(psTFO,i))) errors << 'b' << i;
					else if (! isupper(value(psTTS,i))) errors << 'd' << i;
					else if (! isupper(value(psTFO,i))) errors << 'o' << i;	
					else errors << 't' << i;
				}
				++itTts;
				++itTfo;
				++i;
			}
		} else if (options.errorReference == THIRD_STRAND){
			// requires consideration of parallel/anti-parallel binding
			if (value(tfoSet,match.tfoNo).parallel) {
//				::std::cerr << "\n" << psTTS << "\n" << psTFO << "\n" << tts_ << "\n" << tfo_ << "\n";
				TIter itTts = begin(tts_);
				TIter itTtsEnd = end(tts_);
				TIter itTfo = begin(tfo_);
				TIter itTfoEnd = end(tfo_);
				unsigned i = 0;
				while(itTts != itTtsEnd && itTfo != itTfoEnd){
					if (*itTts != *itTfo){
						if (! isupper(value(psTTS,i)) && ! isupper(value(psTFO,i))) errors << 'b' << i;
						else if (! isupper(value(psTTS,i))) errors << 'd' << i;
						else if (! isupper(value(psTFO,i))) errors << 'o' << i;		
						else errors << 't' << i;
					}
					++itTts;
					++itTfo;
					++i;
				}
			} else {
				reverse(psTTS);
//				::std::cerr << "\n" << psTTS << "\n" << psTFO << "\n" << tts_ << "\n" << tfo_ << "\n";
				TIter itTts = begin(tts_);
				TIter itTtsEnd = end(tts_);
				TIter itTfo = begin(tfo_);
				TIter itTfoEnd = end(tfo_);
				unsigned i = 0;
				while(itTtsEnd != itTts and itTfoEnd != itTfo){
					--itTtsEnd;
					--itTfoEnd;
					if (*itTtsEnd != *itTfoEnd){
						if (! isupper(value(psTTS,i)) && ! isupper(value(psTFO,i))) errors << 'b' << i;
						else if (! isupper(value(psTTS,i))) errors << 'd' << i;
						else if (! isupper(value(psTFO,i))) errors << 'o' << i;	
						else errors << 't' << i;
					}
					++i;
				}
			}
		}
//		::std::cerr << errors.str() << "\n";
		return errors.str();
	}
	
	
	//////////////////////////////////////////////////////////////////////////////
	// Output triplex matches
	template <
	typename TMatches,
	typename TPotentials,
	typename TString,
	typename TMotifSet,
	typename TFile
	>
	void printTriplexEntry(TMatches		&matches,			// forward/reverse matches
						   TPotentials	&tpots, 
						   CharString	&duplexId,			// tts names (read from Fasta file, currently unused)
						   TString		&duplex,			// list of tts names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
						   TMotifSet const				&tfoSet,	// set of tfos
						   StringSet<CharString> const	&tfoNames,	// tfo names (read from Fasta file, currently unused)
						   TFile		&filehandle,				// output file
						   Options		&options
						   ){
		typedef typename Iterator<TMatches, Standard>::Type		TIter;
		typedef typename Value<TMatches>::Type					TMatch;
		typedef ::std::list<TMatch>								TMatchList;
//		typedef Triple<unsigned, unsigned, unsigned>			TCounts;
//		typedef Pair<int,int>									TPair;
//		typedef ::std::set<TPair, LessRPair<TPair> >			TPairList;	// set of matches
//		typedef TPairList *										TPairListP;
//		typedef Triple<unsigned, TCounts, TPairListP>			TBipartiteGraphValue;
		typedef unsigned										TBipartiteGraphKey;
//		typedef Pair<TBipartiteGraphKey,TBipartiteGraphValue>	TBipartiteGraphPair; 
//		typedef Map<TBipartiteGraphPair, Skiplist<> >			TBipartiteGraph;
//		typedef typename Iterator<TBipartiteGraph>::Type		TBipartiteGraphIter;
		
		typedef typename Iterator<TPotentials>::Type			TPotIter;
		typedef typename Value<TPotentials>::Type				TPotValue;
		typedef typename Key<TPotValue>::Type					TPotKey;
		typedef typename Cargo<TPotValue>::Type					TPotCargo;
		
		char _sep_ = '\t';
		for (TIter it = begin(matches, Standard()); it != end(matches, Standard());++it){
			TMatch match = (*it);
			TBipartiteGraphKey seqNo = value(tfoSet,match.tfoNo).seqNo;
			switch (options.outputFormat)
			{
				case 0:	// brief Triplex Format
					filehandle << tfoNames[seqNo] << _sep_ << match.oBegin << _sep_ << match.oEnd << _sep_ ;
					filehandle << duplexId << _sep_ << match.dBegin << _sep_ << match.dEnd << _sep_ ;
					filehandle << match.mScore << _sep_ << ::std::setprecision(2) << (1.0-match.mScore/(match.dEnd-match.dBegin)) << _sep_ ;
					filehandle << _errorString(match, duplex, tfoSet, options) << _sep_ << match.motif << _sep_ << match.strand << _sep_ <<  (match.parallel?'P':'A') << _sep_ <<  (match.guanines/(match.dEnd-match.dBegin)) << ::std::endl;
					break;
				case 1:	// extended Triplex Format
					filehandle << tfoNames[seqNo] << _sep_ << match.oBegin << _sep_ << match.oEnd << _sep_ ;
					filehandle << duplexId << _sep_ << match.dBegin << _sep_ << match.dEnd << _sep_ ;
					filehandle << match.mScore << _sep_ << ::std::setprecision(2) << (1.0-match.mScore/(match.dEnd-match.dBegin)) << _sep_ ;
					filehandle << _errorString(match, duplex, tfoSet, options) << _sep_ << match.motif << _sep_ << match.strand << _sep_ <<  (match.parallel?'P':'A') << _sep_ <<  (match.guanines/(match.dEnd-match.dBegin)) << ::std::endl;
					dumpAlignment(match, duplex, tfoSet, filehandle, options);
					break;
				default:
					break;
			}
		}

		// output triplex potentials
		for(TPotIter it = begin(tpots); it != end(tpots); ++it){
			TPotValue tpotvalue = *it;
			TPotCargo tpotkey = key(tpotvalue);
			TPotCargo tpot = cargo(tpotvalue);
			options.summaryFileHandle << duplexId << _sep_ << tfoNames[getKey(tpot).i1] << _sep_ << getCounts(tpot) << _sep_ << ::std::setprecision(3) << (getCounts(tpot)/getNorm(tpot)) << _sep_;
			options.summaryFileHandle << getCount(tpot,'R') << _sep_ << ::std::setprecision(3) << (getCount(tpot,'R')/getNorm(tpot)) << _sep_;
			options.summaryFileHandle << getCount(tpot,'Y') << _sep_ << ::std::setprecision(3) << (getCount(tpot,'Y')/getNorm(tpot)) << _sep_;
			options.summaryFileHandle << getCount(tpot,'M') << _sep_ << ::std::setprecision(3) << (getCount(tpot,'M')/getNorm(tpot)) << _sep_ << ::std::endl;
		}	
		options.summaryFileHandle.flush();
	}
	
	// output single TTS hit
	template<typename TFile, typename TEntry, typename TDuplexNames>
	inline void printTtsEntry(TFile					&filehandle,
							  TEntry				&entry,
							  unsigned				&counter,
							  TDuplexNames const	&ttsIDs,	// duplex names (read from Fasta file, currently unused)
							  Options const			&options
							  ){
		
		char _sep_ = '\t';
		if (options.outputFormat == 0){ //  Triplex Format
			filehandle << ttsIDs[getSequenceNo(entry)] << _sep_ << beginPosition(entry) << _sep_ <<  endPosition(entry) << _sep_ ;
			filehandle << score(entry) << _sep_ << entry.motif << _sep_ << ::std::setprecision(2) << (1.0-score(entry)/(endPosition(entry)-beginPosition(entry))) << _sep_ << errorString(entry) << _sep_ << guanineRate(entry) << _sep_ << duplicates(entry) << _sep_ ;
			if (options.prettyString)
				filehandle << prettyString(entry) << _sep_ ;
			else
				filehandle << outputString(entry) << _sep_ ;
			
			if (!options.reportDuplicateLocations || duplicates(entry) < 1 || options.duplicatesCutoff <= duplicates(entry) ){
				filehandle << "-" ;
			} else {
				for (int d = 0; d<duplicates(entry); ++d){
					filehandle << ttsIDs[getDuplicateAt(entry, d).i1] << ":" << getDuplicateAt(entry, d).i2 << "-" << getDuplicateAt(entry, d).i2+length(entry) << ";";
				}
			}
			filehandle << ::std::endl;
			
		} else if (options.outputFormat == 1){ // FASTA
			filehandle << ">" ;
			filehandle << ttsIDs[getSequenceNo(entry)] << "_" << counter << _sep_;
			filehandle << beginPosition(entry) << "-" <<  endPosition(entry) << " " << entry.motif << _sep_  << score(entry) << _sep_ << errorString(entry) << _sep_ << guanineRate(entry) << _sep_ << duplicates(entry) << _sep_ ;
			
			
			if (!options.reportDuplicateLocations || duplicates(entry) < 1 || options.duplicatesCutoff <= duplicates(entry) ){
				filehandle << "-" ;
			} else {
				for (int d = 0; d<duplicates(entry); ++d){
					filehandle << ttsIDs[getDuplicateAt(entry, d).i1] << ":" << getDuplicateAt(entry, d).i2 << "-" << getDuplicateAt(entry, d).i2+length(entry) << ";";
				}
			}
			filehandle << ::std::endl;
			
			if (options.prettyString)
				filehandle << prettyString(entry) << ::std::endl;
			else
				filehandle << outputString(entry) << ::std::endl;
			
		} 
		++counter;
	}
	
	
	// output single TFO
	template<typename TFile, typename TEntry, typename TMotifNames>
	inline void printTfoEntry(TFile					&filehandle,
							  TEntry				&entry,
							  unsigned				&counter,
							  TMotifNames const		&tfoIDs,	// duplex names (read from Fasta file, currently unused)
							  Options const			&options
							  ){
				
		char _sep_ = '\t';
		if (options.outputFormat == 0){ // Triplex Format
			filehandle << tfoIDs[getSequenceNo(entry)] << _sep_ << beginPosition(entry) << _sep_ <<  endPosition(entry) << _sep_ ;
			filehandle << score(entry) << _sep_ << entry.motif << _sep_ << ::std::setprecision(2) << (1.0-score(entry)/(endPosition(entry)-beginPosition(entry))) << _sep_ << errorString(entry) << _sep_ << guanineRate(entry) << _sep_ << duplicates(entry) << _sep_ ;
 			
			if (options.prettyString)
				filehandle << prettyString(entry) << _sep_ ;
			else
				filehandle << outputString(entry) << _sep_ ;
			
			if (!options.reportDuplicateLocations || duplicates(entry) < 1 || options.duplicatesCutoff <= duplicates(entry) ){
				filehandle << "-" ;
			} else {
				for (int d = 0; d<duplicates(entry); ++d){
					filehandle << tfoIDs[getDuplicateAt(entry, d).i1] << ":" << getDuplicateAt(entry, d).i2 << "-" << getDuplicateAt(entry, d).i2+length(entry) << ";";
				}
			}
			filehandle << ::std::endl;

			
		} else if (options.outputFormat == 1){ // FASTA
			filehandle << ">" ;
			filehandle  << tfoIDs[getSequenceNo(entry)] << "_" << counter << _sep_;
			filehandle << beginPosition(entry) << "-" <<  endPosition(entry) << " " << entry.motif << _sep_ << score(entry) << _sep_ << errorString(entry) << _sep_ << duplicates(entry) << _sep_ << guanineRate(entry) << _sep_ ;
			
			if (!options.reportDuplicateLocations || duplicates(entry) < 1 || options.duplicatesCutoff <= duplicates(entry) ){
				filehandle << "-" ;
			} else {
				for (int d = 0; d<duplicates(entry); ++d){
					filehandle << tfoIDs[getDuplicateAt(entry, d).i1] << ":" << getDuplicateAt(entry, d).i2 << "-" << getDuplicateAt(entry, d).i2+length(entry) << ";";
				}
			}
			filehandle << ::std::endl;
			
			if (options.prettyString)
				filehandle << prettyString(entry) << ::std::endl;
			else
				filehandle << outputString(entry) << ::std::endl;
			
		}
		++counter;
	}
	
	// output TTS hits
	template <
	typename TFile,
	typename TDuplexSet,
	typename TDuplexNames,
	typename TPotentials
	>
	void dumpTtsMatches(TFile					&filehandle,
						TDuplexSet				&ttsSet,
						TDuplexNames const		&ttsIDs,	// duplex names (read from Fasta file, currently unused)
						TPotentials				&tpots,
						Options					&options)
	{	
		typedef typename Value<TDuplexSet>::Type				TModDuplex;
		typedef typename Iterator<TDuplexSet, Standard>::Type 	TIter;
		typedef typename Iterator<TPotentials>::Type			TPotIter;
		typedef typename Value<TPotentials>::Type				TPot;
		
		// nothing found
		if (length(ttsSet)==0)
			return;
		
		TIter itEnd = end(ttsSet, Standard());
		unsigned counter;
		
		switch (options.outputFormat){
			case 0:	// brief Triplex Format
				counter = 1;
				for(TIter it = begin(ttsSet, Standard()); it != itEnd; ++it){
					printTtsEntry(filehandle, *it, counter, ttsIDs, options);
				}
				break;
				
			case 1:	// extended Triplex Format
				counter = 1;
				for(TIter it = begin(ttsSet, Standard()); it != itEnd; ++it){
					printTtsEntry(filehandle, *it, counter, ttsIDs, options);
				}
				break;
			default:
				break;
		}
		
		// summary file
		char _sep_ = '\t';
		
		for(TPotIter it = begin(tpots, Standard()); it != end(tpots, Standard()); ++it){
			TPot tpot = *it;
			options.summaryFileHandle << ttsIDs[getKey(tpot)] << _sep_ << getCounts(tpot) << _sep_ << ::std::setprecision(3) << (getCounts(tpot)/getNorm(tpot)) << ::std::endl;	
		}
		options.summaryFileHandle.flush();
	}	
		
	//////////////////////////////////////////////////////////////////////////////
	// Output matches
	template <
		typename TFile,
		typename TMotifSet,
		typename TSeqNames,
		typename TPotentials
	>
	void dumpTfoMatches(TFile				&filehandle,
						TMotifSet			&tfoMotifSet,
						TSeqNames const		&tfoIDs,		// tfo names (read from Fasta file, currently unused)
						TPotentials			&tpots,			// holds the matches for all three motifs (Y,R,M) as a Triple
						Options				&options)
	{
		typedef typename Value<TMotifSet>::Type					TMotif;
		typedef typename Iterator<TMotifSet, Standard>::Type	TIter;
		typedef typename Iterator<TPotentials>::Type			TPotIter;
		typedef typename Value<TPotentials>::Type				TPot;
		
		// nothing found?
		if (length(tfoMotifSet)==0){
			return;
		}
			
		SEQAN_PROTIMESTART(dump_time);

		TIter itEnd = end(tfoMotifSet, Standard());
		unsigned counter;
		switch (options.outputFormat)
		{
			case 0:	// brief Triplex Format
				counter = 1;
				for(TIter it = begin(tfoMotifSet, Standard()); it != itEnd; ++it){
					if ((*it).motif == '-')
						continue;
					printTfoEntry(filehandle, *it, counter, tfoIDs, options);
				}
				break;

			case 1:	// Bed format
				counter = 1;
				for(TIter it = begin(tfoMotifSet, Standard()); it != itEnd; ++it){
					printTfoEntry(filehandle, *it, counter, tfoIDs, options);
				}
				break;
			default:
				break;
		}
		
		// summary
		char _sep_ = '\t';
			
		for(TPotIter it = begin(tpots, Standard()); it != end(tpots, Standard()); ++it){
			TPot tpot = *it;
			options.summaryFileHandle << tfoIDs[getKey(tpot)] << _sep_ << getCounts(tpot) << _sep_ << ::std::setprecision(3) << (getCounts(tpot)/getNorm(tpot)) << _sep_;
			options.summaryFileHandle << getCount(tpot,'R') << _sep_ << ::std::setprecision(3) << (getCount(tpot,'R')/getNorm(tpot)) << _sep_;
			options.summaryFileHandle << getCount(tpot,'Y') << _sep_ << ::std::setprecision(3) << (getCount(tpot,'Y')/getNorm(tpot)) << _sep_;
			options.summaryFileHandle << getCount(tpot,'M') << _sep_ << ::std::setprecision(3) << (getCount(tpot,'M')/getNorm(tpot)) << _sep_ << ::std::endl;
			
		}
		options.summaryFileHandle.flush();
		options.timeDumpResults += SEQAN_PROTIMEDIFF(dump_time);
	}

}

#endif

