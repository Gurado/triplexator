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

#ifndef FBUSKE_APPS_TRIPLEXATOR_TRIPLEXATOR_H_
#define FBUSKE_APPS_TRIPLEXATOR_TRIPLEXATOR_H_

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/modifier/modifier_view.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include "triplex.h"

using namespace std;

namespace SEQAN_NAMESPACE_MAIN
{

	// entry function to find TFO/TTS pairs
	template <typename TOligoSet, typename TMotifSet> int mapTriplexes(Options &options);

	// entry function to find TTS 
	template <typename TTargetSet> int investigateTTS(Options &options);

	// entry function to find TFO
	template <typename TOligoSet, typename TMotifSet> int investigateTFO(Options &options);
		
	// application entry function
	int main(int argc, char const ** argv);
	
	//
	// local methods
	//
	
	// setup triplexator command line parser
	void _setupCommandLineParser(CommandLineParser & parser, Options & options);
	
	// parse triplexator command line input
	int _parseCommandLineAndCheck(Options & options,
								 CommandLineParser & parser,
								 int argc,
								 char const ** argv);
	
	// write header to log file
	void _populateLogFile(int argc, const char *argv[], Options	&options);
	
	// find TFO/TTS pairs
	template <
	typename TMotifSet,
	typename TFile,
	typename TShape>
	int _findTriplex(TMotifSet						&tfoMotifSet,
					 StringSet<CharString> const	&tfoNames,
					 TFile							&outputfile,
					 Options						&options,
					 TShape const					&shape);
		
	
	// find TTSs (serial mode)
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
											 Options &options);
	
	// find TTSs (parallel mode)
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
											Options &options);
	// find TFOs
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
								Options &options);
	
	// main helper function
	int _mainWithOptions(int argc, char const ** argv, Options &options);
		
	
}
#endif  // #ifndef FBUSKE_APPS_TRIPLEXATOR_TRIPLEXATOR_H_
