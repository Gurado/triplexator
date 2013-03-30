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

#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_GARDENER_H
#define FBUSKE_APPS_TRIPLEXATOR_HEADER_GARDENER_H

#include <limits>
#include "find_index_qgrams.h"
#include "triplex_alphabet.h"
#include "helper.h"
#include <seqan/seeds2.h>  // Include module under test.
#include <seqan/sequence/adapt_std_list.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>
#include <seqan/misc/misc_dequeue.h>

#if SEQAN_ENABLE_PARALLELISM
#include <seqan/parallel.h>
#endif  // #if SEQAN_ENABLE_PARALLELISM

#ifndef SEQAN_PRAGMA_IF_PARALLEL
#if SEQAN_ENABLE_PARALLELISM
#define STRINGIFY(a) #a
#define SEQAN_PRAGMA_IF_PARALLEL(code) \
_Pragma(STRINGIFY(code))
#else // SEQAN_ENABLE_PARALLELISM
#define SEQAN_PRAGMA_IF_PARALLEL(code)
#endif // SEQAN_ENABLE_PARALLELISM
#endif // SEQAN_PRAGMA_IF_PARALLEL

using namespace seqan;
namespace SEQAN_NAMESPACE_MAIN
{    
	// ============================================================================
	// Tags, Classes, Enums
	// ============================================================================
	
	struct _MULTIPLE_WORKER;
	typedef Tag<_MULTIPLE_WORKER> MULTIPLE_WORKER; // tag for parallel execution
	
	struct _SINGLE_WORKER;
	typedef Tag<_SINGLE_WORKER> SINGLE_WORKER; // tag for serial execution
	
	/**
	 * Class that holds a hit between a needle and a haystack
	 *
	 **/
	template <typename TSpec, typename TPos, typename TId>
	class GardenerHit_
	{
	public:
		TId			hstId;			// haystack sequence id 
		TId			ndlSeqNo;		// needle sequence number
		TPos		hstkPos;		// begin in haystack 
		TPos		ndlPos;			// begin position of hit in needle
		TPos		diag;			// the diagonal
		TPos		score;			// the score
		TPos		hitLength;		// length of the hit
		
		bool operator==(const GardenerHit_<TSpec, TPos, TId>& b) const;
		bool operator!=(const GardenerHit_<TSpec, TPos, TId>& b) const;
		bool operator<(const GardenerHit_<TSpec, TPos, TId>& b) const;
		bool operator>(const GardenerHit_<TSpec, TPos, TId>& b) const;
		
		GardenerHit_(){}
		
		GardenerHit_(TId _hstId, TId _ndlSeqNo, TPos _hstkPos, TPos _ndlPos, TPos _diag, TPos _score, TPos _hitLength):
		hstId(_hstId),
		ndlSeqNo(_ndlSeqNo),
		hstkPos(_hstkPos),
		ndlPos(_ndlPos),
		diag(_diag),
		score(_score),
		hitLength(_hitLength)
		{}
		
		GardenerHit_(GardenerHit_ const &orig): 
		hstId(orig.hstId),
		ndlSeqNo(orig.ndlSeqNo),
		hstkPos(orig.hstkPos),
		ndlPos(orig.ndlPos),
		diag(orig.diag),
		score(orig.score),
		hitLength(orig.hitLength)
		{};
		
		GardenerHit_ & operator = (GardenerHit_ const &orig){
			hstId = orig.hstId;
			ndlSeqNo = orig.ndlSeqNo;
			hstkPos = orig.hstkPos;
			ndlPos = orig.ndlPos;
			diag = orig.diag;
			score = orig.score;
			hitLength = orig.hitLength;
			return *this;
		}
		
		~GardenerHit_(){}
		
		inline TId getHstId(){
			return hstId;
		}
		
		inline TId getNdlSeqNo(){
			return ndlSeqNo;
		}
		
		inline TPos getHstkPos(){
			return hstkPos;
		}
		
		inline TPos getNdlPos(){
			return ndlPos;
		}
		
		inline TPos getDiag(){
			return diag;
		}
		
		inline TPos getScore(){
			return score;
		}
		
		inline TPos getHitLength(){
			return hitLength;
		}
	};
	
	//____________________________________________________________________________
	
	/**
	 * specializations
	 **/
	struct GardenerUngapped_;
	typedef Tag<GardenerUngapped_> GardenerUngapped;
	
	struct GardenerUngappedSegmented_;
	typedef Tag<GardenerUngappedSegmented_> GardenerUngappedSegmented;
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	class Gardener
	{
	public:
		typedef __int64											TPos;
		typedef GardenerHit_<TGardenerSpec, TPos, TId>			TGardenerHit;
		typedef Map<TGardenerHit, Skiplist<> >					THitSet;
		typedef Pair<TId, THitSet* >							THitMapPair;
		typedef Map<THitMapPair, Skiplist< > >					THitMap;	
		typedef typename Iterator<THitSet, Standard>::Type		THitIterator;
		
		THitMap hits; // containing for each  duplex sequence (TId) the list of detected
		
		Gardener<TId, TGardenerSpec>(){}
		
		Gardener<TId, TGardenerSpec>(Gardener<TId, TGardenerSpec> const &orig): hits(orig.hits) {};
		
		Gardener<TId, TGardenerSpec> & operator = (Gardener<TId, TGardenerSpec> const &orig){
			hits = orig.hits;
			return *this;
		}
		
		~Gardener<TId, TGardenerSpec>()
		{
		}
	};

	// seed comparator
	// ATTENTION assumes same diagonal and size of seeds
	template <typename TSeed>
	struct LessRSeed : public ::std::binary_function < TSeed, TSeed, bool >
	{
		inline bool operator() (TSeed const &a, TSeed const &b) const 
		{
			if (getBeginDim0(a) > getBeginDim0(b)) return true;
			return false;
		}
	};
	
	// indicates if two seeds overlap
	// ATTENTION assumes same diagonal and size of seeds
	template <typename TSeed>
	bool isOverlapping(TSeed const & a, TSeed const & b){
		if (getBeginDim0(a) < getBeginDim0(b) and getBeginDim0(b) < getEndDim0(a)) return true;
		if (getBeginDim0(b) < getBeginDim0(a) and getBeginDim0(a) < getEndDim0(b)) return true;
		return false;
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//											Meta Functions		                                                  //
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///.Metafunction.Spec.param.T.type:Class.Gardener
	
	template <typename TId, typename TGardenerSpec>
	struct Spec<Gardener<TId, TGardenerSpec> >
	{
		typedef TGardenerSpec Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Iterator<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::THitIterator Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Cargo<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::TGardenerHit Type;
	};
	
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Id<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::TId Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Value<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::THitMap Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Position<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::TPos Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	inline void reset(Gardener<TId, TGardenerSpec> & me)
	{
		eraseAll(me);
	}
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	inline typename Value<Gardener<TId, TGardenerSpec> >::Type & hits(Gardener<TId, TGardenerSpec> & me)
	{
		return me.hits;
	}
	
	//____________________________________________________________________________
	
	template <typename TSpec, typename TPos, typename TId>
	bool GardenerHit_<TSpec, TPos, TId>::operator==(const GardenerHit_<TSpec, TPos, TId>& b) const {
		if (hstId != b.hstId) return false;
		if (ndlSeqNo != b.ndlSeqNo) return false;
		if (hitLength != b.hitLength) return false;		
		if (ndlPos != b.ndlPos) return false;
		if (hstkPos != b.hstkPos) return false;
		return true;
	}
	
	//____________________________________________________________________________
	
	template <typename TSpec, typename TPos, typename TId>
	bool GardenerHit_<TSpec, TPos, TId>::operator!=(const GardenerHit_<TSpec, TPos, TId>& b) const {
		return !(*this == b);
	}
	
	//____________________________________________________________________________
	
	template <typename TSpec, typename TPos, typename TId>
	bool GardenerHit_<TSpec, TPos, TId>::operator<(const GardenerHit_<TSpec, TPos, TId>& b) const {
		if (hstId < b.hstId) return true;
		if (hstId > b.hstId) return false;
		if (ndlSeqNo < b.ndlSeqNo) return true;
		if (ndlSeqNo > b.ndlSeqNo) return false;
		if (hstkPos < b.hstkPos) return true;
		if (hstkPos > b.hstkPos) return false;
		if (ndlPos < b.ndlPos) return true;
		if (ndlPos > b.ndlPos) return false;
		if (hitLength < b.hitLength) return true;
		return false;
	}
	
	//____________________________________________________________________________
	
	template <typename TSpec, typename TPos, typename TId>
	bool GardenerHit_<TSpec, TPos, TId>::operator>(const GardenerHit_<TSpec, TPos, TId>& b) const {
		return b<*this;
	}


	//____________________________________________________________________________
	
	/** 
	 * returns an iterator to the begin of the results for a given query id
	 */
	template<typename TSpec, typename TId >
	inline bool
	hasAnyHit(Gardener<TId, TSpec> &gardener, 
			  TId &queryid
			  ){
		if (hasKey(gardener.hits, queryid))
			return true;
		else
			return false;
	}
	
	/** 
	 * returns an iterator to the begin of the results for a given query id
	 */
	template<typename TSpec, typename TId >
	inline typename Iterator<typename Gardener<TId, TSpec>::THitSet>::Type 
	harvestBegin(Gardener<TId, TSpec> &gardener, 
				 TId &queryid
				 ){
		if (hasKey(gardener.hits, queryid))
			return begin(*cargo(gardener.hits, queryid));
		else
			return NULL;
	}
	
	//____________________________________________________________________________
	
	/** 
	 * returns an iterator to the end of the results for a given query id
	 */
	template<typename TSpec, typename TId >
	inline typename Iterator<typename Gardener<TId, TSpec>::THitSet>::Type 
	harvestEnd(Gardener<TId, TSpec> &gardener, 
			   TId &queryid
			   ){
		if (hasKey(gardener.hits, queryid))
			return end(*cargo(gardener.hits, queryid));
		else
			return NULL;
	}
	
	//____________________________________________________________________________
	
	/** 
	 * returns the size of the result container
	 */
	template<typename TSpec, typename TId >
	inline typename Size<typename Gardener<TId, TSpec>::THitSet>::Type 
	length(Gardener<TId, TSpec> &gardener, 
		   TId &queryid
		   ){
		if (hasKey(gardener.hits, queryid))
			return length(*cargo(gardener.hits, queryid));
		else
			return 0;
	}
	
	//____________________________________________________________________________
	
	/** 
	 * indicates if a specific query id is registered in the result array 
	 */
	template<typename TSpec, typename TId >
	inline bool hasHit(Gardener<TId, TSpec> &gardener, 
					   TId &queryid, 
					   typename Gardener<TId, TSpec>::TGardenerHit hit
					   ){
		if (hasKey(gardener.hits, queryid) && hasKey(*cargo(gardener.hits,queryid), hit))
			return true;
		else {
			return false;
		}
	}
	

	
	//____________________________________________________________________________
	
	/** 
	 * erase all hits
	 */
	template<typename TSpec, typename TId >
	inline void eraseAll(Gardener<TId, TSpec>  &gardener){
		typedef typename Value<Gardener<TId, TSpec> >::Type	THitMap;
		typedef typename Iterator<THitMap>::Type				THitMapIter;
		typedef typename Value<THitMap>::Type					THitMapPair;
		
		THitMap hitmap = hits(gardener);
		for (THitMapIter it = begin(hitmap);it != end(hitmap);++it){
			THitMapPair hm = *it;
			delete(hm.i2);
		}
		clear(gardener.hits);
	}	

	/**
	 * create runs of mismatches
	 * Triple contains <seed offset mismatch begin, seed offset mismatch end, total errors in sequence>
	 */
	template<typename THaystack, typename TNeedle, typename TSeed, typename TTripleSet>
	inline void	_fillRuns(THaystack const	&haystack,
						  TNeedle const		&needle,
						  TSeed				&seed,
						  TTripleSet		&runs
						  )
	{
		SEQAN_CHECKPOINT
		typedef typename Value<TTripleSet>::Type	TRun;
		typedef typename Value<TRun, 1>::Type		TPos;
		TPos totalErrors = 0;
		TPos mismatchBegin = 0;
		TPos i = mismatchBegin;
		TPos seedLength = getEndDim0(seed)-getBeginDim0(seed);
		
		// append mismatch run starting at beginPosition
		while(i < seedLength && !isMatch(haystack[getBeginDim0(seed)+i], needle[getBeginDim1(seed)+i])) {
			++i;
			++totalErrors;
		}
		appendValue(runs, TRun(mismatchBegin, i, totalErrors));		
		
		// iterate over alignment and append mismatch runs
		while (i < seedLength) {
			// skip matches
			while(i < seedLength && isMatch(haystack[getBeginDim0(seed)+i], needle[getBeginDim1(seed)+i])) {
				++i;
			}
			mismatchBegin = i;
			// skip and count mismatches
			while(i < seedLength && !isMatch(haystack[getBeginDim0(seed)+i], needle[getBeginDim1(seed)+i])) {
				++i;
				++totalErrors;
			}
			appendValue(runs, TRun(mismatchBegin, i, totalErrors));			
		}
	#ifdef TRIPLEX_DEBUG	
		std::cout << infix(needle, getBeginDim1(seed), getEndDim1(seed)) << "\n";
		std::cout << infix(haystack, getBeginDim0(seed), getEndDim0(seed)) << "\n";
		for(unsigned l = 0; l < length(runs); ++l) {
			std::cout << getValue(runs,l).i1 << "  " << getValue(runs,l).i2 << "  " << getValue(runs,l).i3 << std::endl;
		}
	#endif		
	}


	//____________________________________________________________________________


	// checks the error rate of the fragment between end of left and start of right
	template<typename TPos, typename TError>
	inline bool
	_isEpsMatch(Triple<TPos, TPos, TPos> const &left,
				Triple<TPos, TPos, TPos> const &right,
				TError &errorRate
				) {
		SEQAN_CHECKPOINT
		// compute mismatches and length
		TPos errors = right.i3 - left.i3 - (right.i2 - right.i1);
		TPos length = right.i1 - left.i2;
		
		// check error rate
		return errors/(TError)(length) <= errorRate;
	}

	//____________________________________________________________________________


	// counts the number of matches
	template<typename TPos>
	inline TPos
	_countMatches(Triple<TPos, TPos, TPos> const &left,
				  Triple<TPos, TPos, TPos> const &right
				  ) {
		SEQAN_CHECKPOINT
		// compute mismatches and length
		TPos errors = right.i3 - left.i3 - (right.i2 - right.i1);
		TPos length = right.i1 - left.i2;
		
		// check error rate
		return length-errors;
	}


	//____________________________________________________________________________


	/**
	 * Identifies the longest epsilon match and adjusts the seed accordingly
	 */
	template<typename THaystack, typename TNeedle, typename TSeed, typename TPos, typename TError>
	inline TPos _longestEpsMatch(THaystack const	&haystack,
								 TNeedle const		&needle,
								 TSeed				&seed,
								 TPos const			&matchMinLength,
								 TError	const		&errorRate
								 ){
		SEQAN_CHECKPOINT
		// Preprocessing: compute and store mismatch and lengths
		// A run is a triple of mismatch begin position, mismatch end position, 
		// and total number of errors in sequence from begin to end position of this mismatch.
		typedef std::vector<Triple<TPos, TPos, TPos> >	TRuns;
		typedef typename Iterator<TRuns, Rooted>::Type	TRunsIter;
		
		TRuns runs;
		_fillRuns(haystack, needle, seed, runs);
		
		// Identify longest eps match by iterating over combinations of left and right positions
		TRunsIter leftIt = begin(runs, Rooted());
		TRunsIter rightIt = end(runs, Rooted());
		--rightIt;
		
		TPos beginOffset = 0;
		TPos endOffset = 0;
		TPos minLength = matchMinLength - 1;
		TPos matches = 0;
		
		// border up to which the begin of the seed is able to fulfil the minimum length constraint
		TPos seedBeginBorder = max(TPos(0),(*rightIt).i1 - matchMinLength);
		TPos seedEndBorder = (*leftIt).i2 + matchMinLength;
		
		while ((*leftIt).i2 + minLength < (*rightIt).i1 && (*leftIt).i2 <= seedBeginBorder) {
			while ((*leftIt).i2 + minLength < (*rightIt).i1 && (*rightIt).i2 >= seedEndBorder) {
				if(_isEpsMatch(*leftIt, *rightIt, errorRate)) {
					beginOffset = (*leftIt).i2;
					endOffset = (*rightIt).i1; 
					minLength = endOffset - beginOffset;
					matches = _countMatches(*leftIt, *rightIt);
					break;
				}
				--rightIt;
			}
			rightIt = end(runs);
			--rightIt;
			++leftIt;
		}
		
		// adjust the seed borders
		setEndDim0(seed, getBeginDim0(seed) + endOffset);					
		setEndDim1(seed, getBeginDim1(seed) + endOffset);		
		setBeginDim0(seed, getBeginDim0(seed) + beginOffset);
		setBeginDim1(seed, getBeginDim1(seed) + beginOffset);
		return matches;
	}

	//____________________________________________________________________________
	
	/**	
	 * Clips the ends of seeds and add to hitlist
	 */
	template <
	typename THitSet,
	typename TId,
	typename TSeqNo,
	typename TDiag,
	typename THaystack,
	typename TNeedle,
	typename TSeedSet,
	typename TError,
	typename TSize
	>
	inline void _clipSeedEndsAndAdd(THitSet		&hitSet,
									TId			const &queryid,
									TSeqNo		const &seqno,
									TDiag		const &diag,
									THaystack	const &haystack,
									TNeedle		const &needle,
									TSeedSet	&seedset,
									TError		&errorRate,
									TSize		&minLength,
									SINGLE_WORKER const & 
									){
		typedef typename Iterator<TSeedSet, Rooted>::Type	TIter;
		typedef typename Value<TSeedSet>::Type				TSeed;
		typedef typename Position<TSeed>::Type				TPos;
		typedef typename Value<THitSet>::Type				THit;
		
		// cycle through each seedset at a time
		TIter it = begin(seedset, Rooted());
		TIter itEnd = end(seedset, Rooted());
		
		while(it != itEnd){		
			// cut ends to obtain longest epsilon-match that contains the seed alignment
			TPos matches = _longestEpsMatch(haystack, needle, *it, minLength, errorRate);
			
			// size constaint fulfilled ?
			TPos seedLength = getEndDim0(*it)-getBeginDim0(*it);	
			if (seedLength  >= (TPos)minLength){
				
				// create a new hit and append it to the gardeners hit list
				THit hit(queryid,
						 seqno,					// needle seq. number            
						 getBeginDim0(*it),      // begin in haystack      
						 getBeginDim1(*it),		// needle position
						 diag,					// the diagonal
						 matches,
						 seedLength
						 );
				
				// append the hit to the finders hit list if not already contained
				if (!hasKey(hitSet, hit)){
					add(hitSet, hit);
#ifdef TRIPLEX_DEBUG
					::std::cout << "new hit :" << ::std::endl;
					_printHit(hit);
#endif			
				}
			}
			++it;
		}	
	}
	

	//____________________________________________________________________________

	/**	
	 * Clips the ends of seeds and add to hitlist
	 */
	template <
	typename THitSet,
	typename TId,
	typename TSeqNo,
	typename TDiag,
	typename THaystack,
	typename TNeedle,
	typename TSeedSet,
	typename TError,
	typename TSize
	>
	inline void _clipSeedEndsAndAdd(THitSet		&hitSet,
									TId			const &queryid,
									TSeqNo		const &seqno,
									TDiag		const &diag,
									THaystack	const &haystack,
									TNeedle		const &needle,
									TSeedSet	&seedset,
									TError		&errorRate,
									TSize		&minLength,
									MULTIPLE_WORKER const & 
									){
		typedef typename Iterator<TSeedSet, Rooted>::Type	TIter;
		typedef typename Value<TSeedSet>::Type				TSeed;
		typedef typename Position<TSeed>::Type				TPos;
		typedef typename Value<THitSet>::Type				THit;
		
		// cycle through each seed at a time
		TIter it = begin(seedset, Rooted());
		TIter itEnd = end(seedset, Rooted());
		
		while(it != itEnd){		
			// cut ends to obtain longest epsilon-match that contains the seed alignment
			TPos matches = _longestEpsMatch(haystack, needle, *it, minLength, errorRate);
			
			// size constaint fulfilled ?
			TPos seedLength = getEndDim0(*it)-getBeginDim0(*it);	
			if (seedLength  >= (TPos)minLength){
				
				// create a new hit and append it to the gardeners hit list
				THit hit(queryid,
						 seqno,					// needle seq. number            
						 getBeginDim0(*it),      // begin in haystack      
						 getBeginDim1(*it),		// needle position
						 diag,					// the diagonal
						 matches,
						 seedLength
						 );
				
				// append the hit 
				appendValue(hitSet, hit);
			}
			++it;
		}	
	}
	

	//____________________________________________________________________________
	
	/**	
	 * Extends all seeds according to the scoring schema and a X-dropoff 
	 */
	template <
	typename THitSet,
	typename THaystack,
	typename TIndex,
	typename TSpec,
	typename TMap,
	typename TSize,
	typename TScore,
	typename TId
	>
	inline void _extendSeedlings(THitSet					&hitSet,
								 Finder<THaystack, TSpec>	&finder,
								 Pattern<TIndex, TSpec> const &pattern,
								 TMap						&seqmap,
								 Score<TScore, Simple> const &scoreMatrix, 
								 TSize const				&minLength,
								 TScore	const				&scoreDropOff, 
								 TId						&queryid
								 ){
		typedef typename Iterator<TMap>::Type				TMapIter;
		typedef Finder<THaystack, TSpec>					TFinder;
		typedef typename Value<TMap>::Type					TMapPair;
		typedef typename Value<TMapPair,1>::Type			TSeqNo;
		typedef typename Value<TMapPair,2>::Type			TDiagMapPointer;
		typedef typename Value<TDiagMapPointer>::Type		TDiagMap;
		typedef typename Iterator<TDiagMap>::Type			TDiagMapIter;
		typedef typename Value<TDiagMap>::Type				TDiagMapPair;		
		typedef typename Value<TDiagMapPair,1>::Type		TDiag;
		typedef typename Value<TDiagMapPair,2>::Type		TSeedSetPointer;
		typedef typename Value<TSeedSetPointer>::Type		TSeedSet;
		typedef typename Value<TSeedSet>::Type				TSeed;
		typedef typename Position<TSeedSet>::Type			TPos;
		typedef Dequeue<TSeed>								TSeedList;// @TODO known memory leak in seqan string http://trac.mi.fu-berlin.de/seqan/ticket/364
		typedef typename Host<TFinder>::Type				THost;
		typedef typename Value<THitSet>::Type				THit;
		
		// process all needles with entries in the map
		TMapIter itsme = end(seqmap);
		for (TMapIter itsmb = begin(seqmap); itsmb != itsme; ++itsmb){
			TSeqNo &seqno = key(*itsmb);
			// check that seqno is valid
			if (seqno < (TId)countSequences(needle(pattern))){
				TDiagMapPointer diagmapP = cargo(*itsmb);
				TDiagMap diagmap = *diagmapP;
				// process each diagonal at a time
				TDiagMapIter itdme = end(diagmap) ;
				for (TDiagMapIter itdmb = begin(diagmap); itdmb != itdme; ++itdmb){
					TDiag &diag = key(*itdmb);
					// get the seed set
					TSeedSetPointer seedsetP = cargo(*itdmb);
					TSeedSet seedset = *seedsetP;
					THost tmp = host(finder);
					// extend all seeds by first overlapping with succeeding seeds
					TSeedList newset;
					while (!empty(seedset)){
						TSeed seed = front(seedset);
						popFront(seedset);
#ifdef TRIPLEX_DEBUG
						::std::cout << "seed_1:" << getBeginDim0(seed) << ":" << getEndDim0(seed) << " " << getBeginDim1(seed) << ":" << getEndDim1(seed) << ::std::endl;
#endif
						
						while (!empty(seedset) and isOverlapping(seed, front(seedset))){
							TSeed seed2 = front(seedset);
							setBeginDim0(seed, min(getBeginDim0(seed),getBeginDim0(seed2)));
							setBeginDim1(seed, min(getBeginDim1(seed),getBeginDim1(seed2)));
							setEndDim0(seed, max(getEndDim0(seed),getEndDim0(seed2)));
							setEndDim1(seed, max(getEndDim1(seed),getEndDim1(seed2)));
							popFront(seedset);
						}
						
#ifdef TRIPLEX_DEBUG
						::std::cout << "seed_2:" << getBeginDim0(seed) << ":" << getEndDim0(seed) << " " << getBeginDim1(seed) << ":" << getEndDim1(seed) << ::std::endl;
#endif
						// extend seed to both sides as far as possible
						extendSeed(seed, tmp, getSequenceByNo(seqno,needle(pattern)), EXTEND_BOTH, scoreMatrix, scoreDropOff, UnGappedXDrop());

#ifdef TRIPLEX_DEBUG
						::std::cout << "seed_3:" << getBeginDim0(seed) << ":" << getEndDim0(seed) << " " << getBeginDim1(seed) << ":" << getEndDim1(seed) << ::std::endl;
#endif
						
						// merge overlapping windows
						if (getEndDim0(seed)-getBeginDim0(seed) >=  (TPos)minLength){
							if (!empty(newset) && isOverlapping(seed, back(newset))){
								setBeginDim0(seed, min(getBeginDim0(seed),getBeginDim0(back(newset))));
								setBeginDim1(seed, min(getBeginDim1(seed),getBeginDim1(back(newset))));
								setEndDim0(seed, max(getEndDim0(seed),getEndDim0(back(newset))));
								setEndDim1(seed, max(getEndDim1(seed),getEndDim1(back(newset))));
							} else {
								pushBack(newset,seed);
							}
						}
					}
					
					// add all now non-overlapping windows to hitlist
					while(!empty(newset)){
						// create a new hit and append it to the gardeners hit list
						THit hit(queryid,
								 seqno,					// needle seq. number            
								 getBeginDim0(front(newset)),      // begin in haystack      
								 getBeginDim1(front(newset)),		// needle position
								 diag,					// the diagonal
								 0,
								 getEndDim0(front(newset))-getBeginDim0(front(newset))
								 );

#ifdef TRIPLEX_DEBUG
						::std::cout << "newset:" << length(newset) << " hitSet:" << length(hitSet) << ::std::endl;
#endif
						
						// append the hit 
						add(hitSet, hit);
						popFront(newset);
#ifdef TRIPLEX_DEBUG
						::std::cout << "newset:" << length(newset) << " hitSet:" << length(hitSet) << ::std::endl;
#endif
					}
					
//					// find the longest match conform with error rate
//					_clipSeedEndsAndAdd(hitSet, queryid, seqno, diag, host(finder), getSequenceByNo(seqno,needle(pattern)), newset, errorRate, minLength, SINGLE_WORKER() );
					
				} //diagmap
			} else {
				::std::cerr << "Sequence no " << seqno << " exceeds index " << ::std::endl;
			}
		}
	}
	
	//____________________________________________________________________________

	/**
	 * adjust the seed position with respect to the host sequence
	 */
	template <
	typename TSeedSet,
	typename TSequence,
	typename TSpec,
	typename TSize
	>
	inline void _adjustSeeds(TSeedSet					&seedset, 
							 Segment<TSequence, TSpec>	&sequence, 
							 TSize						dim
							 ){
		typedef typename Iterator<TSeedSet>::Type					TSeedSetIter;
		typedef typename Position<Segment<TSequence, TSpec> >::Type	TPos;
		
		TSeedSetIter itsE = end(seedset);
		TPos offset = beginPosition(sequence); 
		for (TSeedSetIter it = begin(seedset); it != itsE; ++it){
			if (dim){
				setBeginDim1(*it, getBeginDim1(*it)+offset);
				setEndDim1(*it, getEndDim1(*it)+offset);
			} else {
				setBeginDim1(*it, getBeginDim0(*it)+offset);
				setEndDim0(*it, getEndDim0(*it)+offset);
			}
		}
	}

	//____________________________________________________________________________
	/**
	 * copy seed passing the qgram lemma
	 */
	template <
	typename TMap,
	typename TId,
	typename TDiag,
	typename TSet
	>
	inline void _putSeedsInMap(TMap			&seqmap,
							   TId	const	&seqNo,
							   TDiag const	&diag,
							   TSet		&posSet
							   ){
		typedef typename Value<TMap>::Type				TMapPair;
		typedef typename Cargo<TMapPair>::Type			TDiagMapPointer;
		typedef typename Value<TDiagMapPointer>::Type 	TDiagMap;
		typedef typename Value<TDiagMap>::Type			TDiagMapPair;
		typedef typename Cargo<TDiagMapPair>::Type		TSeedSetPointer;
		typedef typename Value<TSeedSetPointer>::Type	TSeedSet;
		typedef typename Value<TSeedSet>::Type			TSeed;
		typedef typename Iterator<TSet>::Type			TPosIter;
		typedef typename Value<TSet>::Type				TPos;
		
#ifdef TRIPLEX_DEBUG			
		::std::cout << "add new window:" << length(posSet) << "-" << front(posSet) << " " << back(posSet) << ::std::endl;
#endif	
		// new needle sequence that has no entries yet -- add new needle, and seedset corresponding to diagonal
		if ( ! hasKey(seqmap, seqNo)){		
			// create new needle map
			{
				TDiagMapPointer diagMapPointer = new TDiagMap;
				insert(seqmap, seqNo, diagMapPointer);
			}
			// create new seedset
			{
				TDiagMapPointer tmp_diagMapPointer = cargo(seqmap, seqNo);
				TSeedSetPointer seedSetPointer = new TSeedSet;
				TSeed seed(diag+front(posSet), front(posSet), back(posSet)-front(posSet)+1);
				pushBack(*seedSetPointer,seed);
				insert(*tmp_diagMapPointer, diag, seedSetPointer);
			}
		}
		// needle sequence is known
		else { 
			TDiagMap* diagMapPointer = cargo(seqmap, seqNo);
			// but diagonal index is new
			if (!hasKey(*diagMapPointer, diag)){
				// create new seedset for this diagonal
				TSeedSet* seedSetPointer = new TSeedSet;
				TSeed seed2(diag+front(posSet), front(posSet), back(posSet)-front(posSet)+1);
				pushBack(*seedSetPointer,seed2);
				insert(*diagMapPointer, diag, seedSetPointer);
			} 
			else { // diagonal index is known -- push seed on heap
				TSeedSet* seedSetPointer = cargo(*diagMapPointer, diag);
				TSeed seed3(diag+front(posSet), front(posSet), back(posSet)-front(posSet)+1);
				// check if window simply extends previous one
#ifdef TRIPLEX_DEBUG			
				::std::cout << "extend? " << getBeginDim0(back(*seedSetPointer)) << "=" << getBeginDim0(seed3) << " " << getEndDim0(back(*seedSetPointer)) << " " << getEndDim0(seed3) << ::std::endl;
#endif	
				if (getBeginDim0(back(*seedSetPointer)) <= getBeginDim0(seed3) && getBeginDim0(seed3) <= getEndDim0(back(*seedSetPointer))){
#ifdef TRIPLEX_DEBUG			
					::std::cout << "was: " << getBeginDim0(back(*seedSetPointer)) << "-" << getEndDim0(back(*seedSetPointer)) << " " << getBeginDim1(back(*seedSetPointer)) << "-" << getEndDim1(back(*seedSetPointer)) <<  ::std::endl;
#endif	
					setEndDim0(back(*seedSetPointer), getEndDim0(seed3));
					setEndDim1(back(*seedSetPointer), getEndDim1(seed3));					
#ifdef TRIPLEX_DEBUG			
					::std::cout << "now: " << getBeginDim0(back(*seedSetPointer)) << "-" << getEndDim0(back(*seedSetPointer)) << " " << getBeginDim1(back(*seedSetPointer)) << "-" << getEndDim1(back(*seedSetPointer)) <<  ::std::endl;
#endif	
					
				} else { // add
					pushBack(*seedSetPointer,seed3);
				}
			}
		}	
	}
	
	//____________________________________________________________________________
	/**
	 * get all the hits between needles and haystack
	 * TMap seqmap will be filled using the new operator 
	 * freeing the memory is up to the caller
	 */
	template <
	typename THaystack,
	typename TSpec,
	typename TIndex,
	typename TPos,
	typename TMap
	>
	inline void _collectSeeds(Finder<THaystack, QGramsLookup<TSpec> >		&finder,
							  Pattern<TIndex,  QGramsLookup<TSpec> > const	&pattern,
							  TPos const									&seedsThreshold,
							  TPos const									&minLength,
							  TMap											&seqmap
							  ){
		typedef typename Value<TMap>::Type				TMapPair;
		typedef typename Key<TMapPair>::Type			TId;
		typedef typename Cargo<TMapPair>::Type			TDMPointer;
		typedef typename Value<TDMPointer>::Type		TDM;
		typedef typename Value<TDM>::Type				TDMPair;
		typedef typename Key<TDMPair>::Type				TDiag;
		
		// tmpSeqMap structure constains list of needle position rather than seeds
		typedef TPos									TSeed;
		typedef Dequeue<TSeed>							TSeedSet;
		typedef Pair<TDiag, TSeedSet*>					TDiagMapPair;
		typedef Map<TDiagMapPair, Skiplist< > >			TDiagMap;
		typedef Pair<TId, TDiagMap*>					TSeqMapPair;
		typedef Map<TSeqMapPair, Skiplist< > >			TSeqMap;
		typedef typename Iterator<TDiagMap>::Type		TDiagMapIter;
		typedef typename Iterator<TSeqMap>::Type		TIterM;
		typedef typename Iterator<TDiagMap>::Type		TIterD;
		typedef typename Iterator<TSeedSet>::Type		TIterS;		
		
		TSeqMap tmpSeqmap;
		while (find(finder, pattern)) {
#ifdef TRIPLEX_DEBUG			
			::std::cout << "Q:" << infix(finder) << ::std::endl;
			::std::cout << "T:" << infix(pattern, *finder.curHit) << ::std::endl;
			::std::cout << "H:" << (*finder.curHit).hstkPos << "-N" << (*finder.curHit).ndlSeqNo << ":P" << (*finder.curHit).ndlPos << ":D" << (*finder.curHit).diag << ::std::endl;
			::std::cout << "tmpseqmanSize:" << length(tmpSeqmap) << " ,KeyKnown:" << hasKey(tmpSeqmap, (*finder.curHit).ndlSeqNo) << ::std::endl;
			::std::cout << "seqmanSize:" << length(seqmap) << " ,KeyKnown:" << hasKey(seqmap, (*finder.curHit).ndlSeqNo) << ::std::endl;
#endif			
			// new needle sequence that has no entries yet -- add new needle, and seedset corresponding to diagonal
			if ( ! hasKey(tmpSeqmap, (*finder.curHit).ndlSeqNo)){		
				// create new needle map
				{
					TDiagMap* diagMapPointer = new TDiagMap;
					insert(tmpSeqmap,(*finder.curHit).ndlSeqNo, diagMapPointer);
				}
				// create new seedset
				{
					TDiagMap* tmp_diagMapPointer = cargo(tmpSeqmap, (*finder.curHit).ndlSeqNo);
					TSeedSet* seedSetPointer = new TSeedSet;
					pushBack(*seedSetPointer,(*finder.curHit).ndlPos);
					if (length(*seedSetPointer) >= seedsThreshold){
						_putSeedsInMap(seqmap, (*finder.curHit).ndlSeqNo, (*finder.curHit).diag, *seedSetPointer);
					}
					
					insert(*tmp_diagMapPointer, (*finder.curHit).diag, seedSetPointer);
				}
			}
			// needle sequence is known
			else { 
				TDiagMap* diagMapPointer = cargo(tmpSeqmap, (*finder.curHit).ndlSeqNo);
				// but diagonal index is new
				if (!hasKey(*diagMapPointer, (*finder.curHit).diag)){
					// create new seedset for this diagonal
					TSeedSet* seedSetPointer = new TSeedSet;
					pushBack(*seedSetPointer,(*finder.curHit).ndlPos);
					if (length(*seedSetPointer) >= seedsThreshold){
						_putSeedsInMap(seqmap, (*finder.curHit).ndlSeqNo, (*finder.curHit).diag, *seedSetPointer);
					}					
					insert(*diagMapPointer, (*finder.curHit).diag, seedSetPointer);
					// diagonal index is known -- push seed on heap
				} 
			 	else { 
					TSeedSet* seedSetPointer = cargo(*diagMapPointer, (*finder.curHit).diag);
					// remove positions outside the window
					while( !empty(*seedSetPointer) && front(*seedSetPointer)+minLength < (*finder.curHit).ndlPos+weight(pattern.shape)){
#ifdef TRIPLEX_DEBUG
						::std::cout << "Poping " << (*finder.curHit).ndlSeqNo << " 1st:" << front(*seedSetPointer) << " cur:" << (*finder.curHit).ndlPos << " diag:" << (*finder.curHit).diag << " " << (*finder.curHit).ndlPos+weight(pattern.shape)-minLength << " l" << length(*seedSetPointer) <<::std::endl;
#endif
						popFront(*seedSetPointer);
					}
					pushBack(*seedSetPointer,(*finder.curHit).ndlPos);
					if (length(*seedSetPointer) >= seedsThreshold){
						_putSeedsInMap(seqmap, (*finder.curHit).ndlSeqNo, (*finder.curHit).diag, *seedSetPointer);
					}
				}
			}
		}
		
		// housekeeping
		// empty tmp seqmap properly
		for (TIterM sit=begin(tmpSeqmap); sit != end(tmpSeqmap); ++sit){
			TDiagMap* diagmapPointer = (*sit).i2;
			for (TIterD dit=begin(*diagmapPointer); dit != end(*diagmapPointer); ++dit){
				delete (*dit).i2;
			}
			delete diagmapPointer;
		}
	}

		
	//____________________________________________________________________________
	
	/**	
	 * Performs a q-gram search of the needle set and the haystack and merges the resulting seeds 
	 * subject to a scoring scheme, and an error rate.
	 * Reports hits subject to a minimum length constraint. 
	 * 
	 * Important note 1: The gardener is an ungapped search algorithm
	 * Important note 2: The algorithm is designed so it can be called using multiple haystacks 
	 * but the same set of needles (needles shared in threads)
	 * Important note 3: Gardener should be used with relatively short needles and haystacks only, since its
	 * memory consumption is O(|haystack| * sum(|needle|)) in worst case
	 */
	template <
	typename THitSet,
	typename THaystack,
	typename TSpec,
	typename TIndex,
	typename TError,
	typename TPos,
	typename TDrop,
	typename TId
	>
	inline bool _find(THitSet							&hitSet,
					  Finder<THaystack,  TSpec >		&finder,
					  Pattern<TIndex,  TSpec > const	&pattern,
					  TError const						&errorRate,
					  TPos const						&minLength,
					  TPos const						&seedsThreshold,
					  TDrop const						&xDrop,
					  TId								&queriyid
					  ){
		// used datastructure: Map < NeedleSeqNo, < Map < diagonal, SeedSet > > 
		// for each needle all q-gram hits are stored according to the diagonal they reside in
		typedef int											TScore;
		typedef Seed<Simple, DefaultSeedConfig>				TSeed;
		//		typedef Seed<TPos, SimpleSeed>						TSeed;		
		//		typedef LessRSeed<TSeed>							TLess;
		//		typedef PriorityType<TSeed, TLess, PriorityHeap>	TSeedSet;
		typedef Dequeue<TSeed>								TSeedSet;
		typedef typename MakeSigned_<TPos>::Type			TDiag;
		typedef Pair<TDiag, TSeedSet*>						TDiagMapPair;
		typedef Map<TDiagMapPair, Skiplist< > >				TDiagMap;
		typedef Pair<TId, TDiagMap*>						TSeqMapPair;
		typedef Map<TSeqMapPair, Skiplist< > >				TSeqMap;
		typedef typename Iterator<TSeqMap>::Type			TSeqMapIter;
		typedef typename Iterator<TDiagMap>::Type			TDiagMapIter;
		
		// run gardener on first call
		if (empty(finder) ){
			TSeqMap seqmap;
			// get all maxed seeds for any needle in the haystack (flanked by mismatches)
			_collectSeeds(finder, pattern, seedsThreshold, minLength, seqmap);
			
			// define a scoring scheme
			TScore match = 1;
			TScore mismatch = (TScore)_max((TScore) (-1.0/(errorRate+0.00000001)) + 1, -(TScore)length(host(finder)));
			Score<TScore> scoreMatrix(match, mismatch, std::numeric_limits<int>::max());
			TScore scoreDropOff = (TScore) _max((TScore) xDrop * (-mismatch), minValue<TScore>()+1);
			// extend seeds 
			_extendSeedlings(hitSet, finder, pattern, seqmap, scoreMatrix, minLength, scoreDropOff, queriyid);
			
			// housekeeping
			// free memory from gram-hitlist
			clear(finder.hits);
			// empty seqmap properly
			for (TSeqMapIter sit=begin(seqmap); sit != end(seqmap); ++sit){
				TDiagMap* diagmapPointer = (*sit).i2;
				for (TDiagMapIter dit=begin(*diagmapPointer); dit != end(*diagmapPointer); ++dit){
					delete (*dit).i2;
				}
				delete diagmapPointer;
			}
			
			return true;
		} else 
			return false;
	}
	
	//____________________________________________________________________________

	/** 
	 * prints all hits for a query
	 */
	template< 
	typename TSpec, 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TQuerySet,		// query set (needle)
	typename TId 
	>	
	inline void _printHits(Gardener<TId, TSpec>					&gardener, 
						   Pattern<TIndex, TPatternSpec> const	&pattern,
						   TQuerySet							&queries,
						   TId									&queryid
						   ){
		typedef Gardener<TId, TSpec>						TGardener;
		typedef typename Iterator<TGardener >::Type			TIter;
		typedef typename Value<TGardener>::Type				THitMap;
		typedef typename Cargo<TGardener>::Type				THit;
		
		if (hasAnyHit(gardener, queryid)){
			for (TIter it = harvestBegin(gardener,queryid); it != harvestEnd(gardener, queryid); ++it){
				THit hit = *it;
				::std::cout << hit.ndlSeqNo << "("<< hit.ndlPos <<":" << (hit.ndlPos+hit.hitLength) << ")\t";
				::std::cout << hit.hstId << "("<< hit.hstkPos <<":" << (hit.hstkPos+hit.hitLength) << ")\t";
				::std::cout << "d:" << hit.diag << "\t";
				std::cout << infix(getSequenceByNo(hit.ndlSeqNo, needle(pattern)), hit.ndlPos, hit.ndlPos+hit.hitLength) << "\t";
				::std::cout << hit.ndlSeqNo << "("<< beginPosition(getSequenceByNo(hit.ndlSeqNo, needle(pattern)))+hit.ndlPos;
				::std::cout <<":" << beginPosition(getSequenceByNo(hit.ndlSeqNo, needle(pattern)))+(hit.ndlPos+hit.hitLength) << ")\t";
				std::cout << infix(queries[hit.hstId], hit.hstkPos, hit.hstkPos+hit.hitLength) << "\t";
				::std::cout << hit.hstId << "("<< beginPosition(queries[hit.hstId])+hit.hstkPos;
				::std::cout <<":" << beginPosition(queries[hit.hstId])+(hit.hstkPos+hit.hitLength) << ")\t";
				::std::cout << ::std::endl;
			}
		} else{
			::std::cout << "No Hits found" << ::std::endl;
		}
	}
	
	template< 
	typename THit 
	>	
	inline void _printHit(THit	&hit){
		::std::cout << hit.ndlSeqNo << "("<< hit.ndlPos <<":" << (hit.ndlPos+hit.hitLength) << ")\t";
		::std::cout << hit.hstId << "("<< hit.hstkPos <<":" << (hit.hstkPos+hit.hitLength) << ")\t";
		::std::cout << "d:" << hit.diag << "\t";
		::std::cout << ::std::endl;	
	}
	//////////////////////////////////////////////////////////////////////////////


	//____________________________________________________________________________
	// entry functions

	
	/** 
	 * start gardening by planting from scratch preparing patterns (needles) first
	 */
	template< 
	typename TMotif,		// target set (pattern)
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TShape,		// shape
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TWORKER		// Worker Tag
	>
	void plant(Gardener<TId, TSpec>		&gardener,
			   StringSet<TMotif> const	&targets,
			   TQuerySet				&queries,
			   TError const				&errorRate,
			   TSize const				&minLength,
			   TDrop const				&xDrop,
			   TShape const				&shape,
			   TWORKER
			   ){
		typedef StringSet<TMotif>													TTargetSet;
		typedef Index<TTargetSet, IndexQGram<TShape, OpenAddressing> >				TQGramIndex;
		// configure Pattern
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup> >	TFinder;
		typedef typename Infix<typename GetSequenceByNo<TQGramIndex const>::Type >::Type	TInfix;
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		
		typedef typename Fibre<TQGramIndex, QGramCounts>::Type						TSA;
		typedef typename Iterator<TSA, Standard>::Type								TSAIter;
		typedef typename Position<TSequence>::Type									TPosition;
		typedef typename Gardener<TId, TSpec>::THitSet								THits;
		
		// create index
		TQGramIndex index_qgram(targets);
		resize(indexShape(index_qgram), weight(shape));
		
		// create pattern
		Pattern<TQGramIndex, QGramsLookup< TShape, Standard_QGramsLookup> > pattern(index_qgram,shape);
		plant(gardener, pattern, queries, errorRate, minLength, xDrop, TWORKER() );
	}
		
	/** 
	 * start gardening by planting
	 * prepeare pattern beforehand for reuse
	 * all putative TTSs per duplex will be processed in parallel mode
	 * (for short duplex sequences with presumably few TTSs)
	 */
	template< 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TShape,		// shape
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TWorker
	>
	void plant(Gardener<TId, TSpec>	&gardener,
			   Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> > const	&pattern,
			   TQuerySet			&queries,
			   TError const			&errorRate,
			   TSize const			&minLength,
			   TDrop const			&xDrop,
			   TWorker
			   ){
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		typedef typename Value<Gardener<TId, TSpec> >::Type							THitMap;
		typedef typename Value<THitMap>::Type										THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type								THitSetPointer;
		typedef typename Value<THitSetPointer>::Type								THitSet;
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup > >	TFinder;
		typedef typename Position<TFinder>::Type									TPos;
		
		// q-gram lemma
		// w+1-(k+1)q | w=minimum length, k=errors, q=weight(q-grams)
		TPos minSeedsThreshold = static_cast<TPos>(minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape));
#ifdef TRIPLEX_DEBUG
		::std::cout << "minLength:" << minLength << " errorRate:" << errorRate << " qgram:" << weight(pattern.shape) << ::std::endl;
		::std::cout << (ceil(errorRate*minLength)+1) << " " << ((ceil(errorRate*minLength)+1)*weight(pattern.shape)) << " " << minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape) << ::std::endl;
#endif						
		
		// serial processing
		TId querylen = (TId)length(queries);
		for (TId queryid=0; queryid<querylen; ++queryid){
			THitSetPointer hitsPointer = new THitSet;
			TFinder finder(queries[queryid]); 
			_find(*hitsPointer, finder, pattern, errorRate, (TPos) minLength, minSeedsThreshold, xDrop, queryid );	
			insert(gardener.hits, queryid, hitsPointer);
#ifdef TRIPLEX_DEBUG
			::std::cout << "TTS " << queryid << " : " << queries[queryid] << ::std::endl;
			_printHits(gardener, pattern, queries, queryid);
#endif
		}		
	}
		
	/** 
	 * start gardening by planting
	 * prepeare pattern beforehand for reuse
	 * all putative TTSs per duplex will be processed in parallel mode
	 * (for short duplex sequences with presumably few TTSs)
	 */
	template< 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TShape,		// shape
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TRepeat,		// repeat minimum length
	typename TWorker
	>
	void plant(Gardener<TId, TSpec>	&gardener,
			   Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> > const	&pattern,
			   TQuerySet			&queries,
			   TError const			&errorRate,
			   TSize const			&minLength,
			   TDrop const			&xDrop,
			   TRepeat const		&minRepeatLength,
			   TRepeat const		&maxRepeatPeriod,
			   TWorker
			   ){
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		typedef typename Value<Gardener<TId, TSpec> >::Type							THitMap;
		typedef typename Value<THitMap>::Type										THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type								THitSetPointer;
		typedef typename Value<THitSetPointer>::Type								THitSet;
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup > >	TFinder;
		typedef typename Position<TFinder>::Type									TPos;
		
		// q-gram lemma
		// w+1−(k+1)q | w=minimum length, k=errors, q=weight(q-grams)
		TPos minSeedsThreshold = static_cast<TPos>(minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape));
#ifdef TRIPLEX_DEBUG
		::std::cout << "minLength:" << minLength << " errorRate:" << errorRate << " qgram:" << weight(pattern.shape) << ::std::endl;
		::std::cout << (ceil(errorRate*minLength)+1) << " " << ((ceil(errorRate*minLength)+1)*weight(pattern.shape)) << " " << minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape) << ::std::endl;
#endif								
		// serial processing
		TId querylen = (TId)length(queries);
		for (TId queryid=0; queryid<querylen; ++queryid){
			THitSetPointer hitsPointer = new THitSet;
			TFinder finder(queries[queryid], minRepeatLength, maxRepeatPeriod); 
			_find(*hitsPointer, finder, pattern, errorRate, (TPos) minLength, minSeedsThreshold, xDrop, queryid );	
			insert(gardener.hits, queryid, hitsPointer);
#ifdef TRIPLEX_DEBUG
			::std::cout << "TTS " << queryid << " : " << queries[queryid] << ::std::endl;
			_printHits(gardener, pattern, queries, queryid);
#endif
		}		
	}
		
#if SEQAN_ENABLE_PARALLELISM	
	/** 
	 * start gardening by planting
	 * prepeare pattern beforehand for reuse
	 * all putative TTSs per duplex will be processed in parallel mode
	 * (for long duplex sequences with presumably many TTSs)
	 */
	template< 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TShape,		// shape
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId			// sequence id
	>
	void plant(Gardener<TId, TSpec>	&gardener,
			   Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> > const &pattern,
			   TQuerySet			&queries,
			   TError const			&errorRate,
			   TSize const			&minLength,
			   TDrop const			&xDrop,
			   MULTIPLE_WORKER
			   ){
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		typedef typename Value<Gardener<TId, TSpec> >::Type							THitMap;
		typedef typename Value<THitMap>::Type										THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type								THitSetPointer;
		typedef typename Value<THitSetPointer>::Type								THitSet;
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup > >	TFinder;
		typedef typename Position<TFinder>::Type									TPos;
		
		TId querylen = (TId)length(queries);

		// q-gram lemma
		// w+1−(k+1)q | w=minimum length, k=errors, q=weight(q-grams)
		TPos minSeedsThreshold = static_cast<TPos>(minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape));
#ifdef TRIPLEX_DEBUG
		::std::cout << "minLength:" << minLength << " errorRate:" << errorRate << " qgram:" << weight(pattern.shape) << ::std::endl;
		::std::cout << (ceil(errorRate*minLength)+1) << " " << ((ceil(errorRate*minLength)+1)*weight(pattern.shape)) << " " << minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape) << ::std::endl;
#endif		
		
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel) 
		{
			SEQAN_PRAGMA_IF_PARALLEL(omp for schedule(dynamic) )
			for (TId queryid=0; queryid<querylen; ++queryid){
				THitSetPointer hitsPointer = new THitSet;
				TFinder finder(queries[queryid]);
				_find(*hitsPointer, finder, pattern, errorRate, (TPos) minLength, minSeedsThreshold, xDrop, queryid );	
				
				SEQAN_PRAGMA_IF_PARALLEL(omp critical(addhitmap)  )
				insert(gardener.hits, queryid, hitsPointer);
			}
		}
	}
	
	/** 
	 * start gardening by planting
	 * prepeare pattern beforehand for reuse
	 * all putative TTSs per duplex will be processed in parallel mode
	 * (for long duplex sequences with presumably many TTSs)
	 */
	template< 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TShape,		// shape
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TRepeat		// repeat minimum length
	>
	void plant(Gardener<TId, TSpec>	&gardener,
			   Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> > const &pattern,
			   TQuerySet			&queries,
			   TError const			&errorRate,
			   TSize const			&minLength,
			   TDrop const			&xDrop,
			   TRepeat const		&minRepeatLength,
			   TRepeat const		&maxRepeatPeriod,
			   MULTIPLE_WORKER
			   ){
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		typedef typename Value<Gardener<TId, TSpec> >::Type							THitMap;
		typedef typename Value<THitMap>::Type										THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type								THitSetPointer;
		typedef typename Value<THitSetPointer>::Type								THitSet;
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup > >	TFinder;
		typedef typename Position<TFinder>::Type									TPos;
		
		TId querylen = (TId)length(queries);
		
		// q-gram lemma
		// w+1−(k+1)q | w=minimum length, k=errors, q=weight(q-grams)
		TPos minSeedsThreshold = static_cast<TPos>(minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape));
#ifdef TRIPLEX_DEBUG
		::std::cout << "minLength:" << minLength << " errorRate:" << errorRate << " qgram:" << weight(pattern.shape) << ::std::endl;
		::std::cout << (ceil(errorRate*minLength)+1) << " " << ((ceil(errorRate*minLength)+1)*weight(pattern.shape)) << " " << minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape) << ::std::endl;
#endif	
		
		// create hits pointer for all query ids
		::std::vector<THitSetPointer> tmpPointerList;
		
		for (TId queryid=0; queryid<querylen; ++queryid){
			THitSetPointer hitsPointer = new THitSet;
			appendValue(tmpPointerList, hitsPointer);
		}
		
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel)
		{
			SEQAN_PRAGMA_IF_PARALLEL(omp for schedule(dynamic))
			for (TId queryid=0; queryid<querylen; ++queryid){
				TFinder finder(queries[queryid], minRepeatLength, maxRepeatPeriod); 
				_find(*tmpPointerList[queryid], finder, pattern, errorRate, (TPos) minLength, minSeedsThreshold, xDrop, queryid );	
			}
		}
		
		// copy all hits to the gardener
		for (TId queryid=0; queryid<querylen; ++queryid){
			insert(gardener.hits, queryid, tmpPointerList[queryid]);
		}
			
	}
		
#endif  // SEQAN_ENABLE_PARALLELISM
	
} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_GARDENER_H
