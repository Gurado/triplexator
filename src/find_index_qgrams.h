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

#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_INDEX_QGRAMS_FIND_H
#define FBUSKE_APPS_TRIPLEXATOR_HEADER_INDEX_QGRAMS_FIND_H

#include <seqan/find.h>
#include <seqan/index.h>

namespace SEQAN_NAMESPACE_MAIN
{
	
	template <typename TSpec, typename THstkPos, typename TDiag>
	struct _QGramHit 
	{
		THstkPos	hstkPos;			// begin in haystack 
		unsigned	ndlSeqNo;			// needle sequence number
		THstkPos	ndlPos;				// begin position of hit in needle
		TDiag		diag;				// the diagonal
//		unsigned	length;				// length of the hit
	};
		
	/**
	 .Tag.Index Find Algorithm
	 ..tag.QGram_FIND_Lookup:q-gram search.
	 Finds q-grams in a @Spec.Index_QGrams@ index using the hash table.
	 */

	struct Finder_QGramsLookup_; //Finder that simply looks up the q-gram in the hash table
	typedef Tag<Finder_QGramsLookup_> const Standard_QGramsLookup;
	
	template <typename TShape, typename TSpec = Standard_QGramsLookup>
	struct QGramsLookup;
	
	//////////////////////////////////////////////////////////////////////////////
	// QGram finders
	
	template <typename THaystack, typename TShape, typename TSpec>
	class Finder<THaystack, QGramsLookup<TShape, TSpec> >
	{
	public:
		typedef typename Iterator<THaystack, Rooted>::Type		TIterator;
		typedef typename Position<THaystack>::Type				THstkPos;
		typedef typename MakeSigned_<THstkPos>::Type			TDiag;
		typedef _QGramHit<TSpec, THstkPos, TDiag>				TQGramHit;
		typedef std::list<TQGramHit>							THitString; //@TODO workaround for memory leak in seqan string
		typedef typename Iterator<THitString, Rooted>::Type		THitIterator;
		typedef typename SAValue<THaystack>::Type				TSAValue;
		typedef Repeat<TSAValue, unsigned>						TRepeat;
		typedef String<TRepeat>									TRepeatString;
		typedef typename Iterator<TRepeatString, Rooted>::Type	TRepeatIterator;
		
		TIterator		data_iterator;
		TIterator		haystackEnd;
		bool			_needReinit;	// if true, the Pattern needs to be reinitialized
		THitString		hits;
		THitIterator	curHit, endHit;
		THstkPos		startPos, curPos, endPos;
		THstkPos		dotPos, dotPos2;
		TRepeatString	data_repeats;
		TRepeatIterator	curRepeat, endRepeat;
		TShape			shape;	// shape needs to be saved in Finder since it contains the next hashvalue and pattern needs to stay const 
		bool			hasShape;
		int				maxHitThreshold;
		
		Finder():
		_needReinit(true), hasShape(false) { }
		
		Finder(THaystack &haystack):
		data_iterator(begin(haystack, Rooted())),
		_needReinit(true),
		hasShape(false),
		maxHitThreshold(0){ 
			TShape shape;
		}
		
		
		Finder(THaystack &haystack, int _maxHitThreshold):
		data_iterator(begin(haystack, Rooted())),
		_needReinit(true),
		hasShape(false),
		maxHitThreshold(_maxHitThreshold)
		{ 
			TShape shape;
		}
		
		template <typename TRepeatSize, typename TPeriodSize>
		Finder(THaystack &haystack, TRepeatSize minRepeatLen, TPeriodSize maxPeriod):
		data_iterator(begin(haystack, Rooted())),
		_needReinit(true) ,
		hasShape(false),
		maxHitThreshold(0)
		{
			findRepeats(data_repeats, haystack, minRepeatLen, maxPeriod);
		}
		
		template <typename TRepeatSize, typename TPeriodSize>
		Finder(THaystack &haystack, TRepeatSize minRepeatLen, TPeriodSize maxPeriod, int _maxHitThreshold):
		data_iterator(begin(haystack, Rooted())),
		_needReinit(true) ,
		hasShape(false),
		maxHitThreshold(_maxHitThreshold)
		{
			findRepeats(data_repeats, haystack, minRepeatLen, maxPeriod);
		}
		
		Finder(TIterator &iter):
		data_iterator(iter),
		_needReinit(true),
		hasShape(false){ }
		
		Finder(TIterator const &iter):
		data_iterator(iter),
		_needReinit(true),
		hasShape(false){ }
		
		Finder(Finder const &orig):
		data_iterator(orig.data_iterator),
		haystackEnd(orig.haystackEnd),
		_needReinit(orig._needReinit),
		hits(orig.hits),
		startPos(orig.startPos),
		curPos(orig.curPos),
		endPos(orig.endPos),
		dotPos(orig.dotPos),
		dotPos2(orig.dotPos2),
		data_repeats(orig.data_repeats),
		hasShape(orig.hasShape),
		shape(orig.shape),
		maxHitThreshold(orig.maxHitThreshold)
		{
			curHit = begin(hits, Rooted()) + (orig.curHit - begin(orig.hits, Rooted()));
			endHit = end(hits, Rooted());
			curRepeat = begin(data_repeats, Rooted()) + (orig.curRepeat - begin(orig.data_repeats, Rooted()));
			endRepeat = end(data_repeats, Rooted());
		};
		
		inline typename Reference<TIterator>::Type 
		operator* () { return value(hostIterator(*this)); }
		
		inline typename Reference<TIterator const>::Type 
		operator* () const { return value(hostIterator(*this)); }
		
		operator TIterator () const	{ return data_iterator;	}
        
        Finder & operator = (Finder const &orig) 
        {
            data_iterator = orig.data_iterator;
            haystackEnd = orig.haystackEnd;
            _needReinit = orig._needReinit;
            hits = orig.hits;
            startPos = orig.startPos;
            curPos = orig.curPos;
            endPos = orig.endPos;
            dotPos = orig.dotPos;
            dotPos2 = orig.dotPos2;
            data_repeats = orig.data_repeats;
            curHit = begin(hits, Rooted()) + (orig.curHit - begin(orig.hits, Rooted()));
            endHit = end(hits, Rooted());
            curRepeat = begin(data_repeats, Rooted()) + (orig.curRepeat - begin(orig.data_repeats, Rooted()));
            endRepeat = end(data_repeats, Rooted());
			hasShape = orig.hasShape;
			shape = orig.shape;
			maxHitThreshold = orig.maxHitThreshold;
            return *this;
        }
		
    };
	
	//____________________________________________________________________________
	
	
	template <typename TIndex, typename TShape, typename TSpec>
	class Pattern<TIndex,  QGramsLookup<TShape, TSpec> >
	{
	public:
		typedef typename Size<TIndex>::Type								TSize;
		typedef unsigned												TShortSize;
		typedef typename Fibre<TIndex, Tag<FibreSA> const >::Type		TSA;
		typedef typename Iterator<TSA const, Standard>::Type			TIterator;
		
		Holder<TIndex>	data_host;
		TShape	const	shape;	// cannot be used to compute hash (not thread save) For reference purpose only

		Pattern(TIndex &_index, TShape const &shape): data_host(_index), shape(shape) {
			indexRequire(_index, QGramSADir());
		}
		Pattern(TIndex const &_index, TShape const &shape): data_host(_index), shape(shape){
			indexRequire(_index, QGramSADir());
		}
		
	};	
	

	
	//____________________________________________________________________________
	
	
	template <typename THaystack, typename TSpec>
	inline bool
	atEnd(Finder<THaystack,  QGramsLookup<TSpec> > & me)
	{
		return hostIterator(hostIterator(me)) == hostIterator(me.haystackEnd);
	}
	
	template <typename THaystack, typename TSpec>
	inline void
	goEnd(Finder<THaystack,  QGramsLookup<TSpec>  > & me)
	{
		hostIterator(me) = me.haystackEnd;
	}

	
	//____________________________________________________________________________
	
	
	template <typename TIndex, typename TSpec, typename TSeqNo>
	inline int
	_qgramLemma(Pattern<TIndex, QGramsLookup<TSpec> > const & pattern, 
				TSeqNo seqNo, 
				int errors
	){
		// q-gram lemma: How many conserved q-grams we see at least?
		// each error destroys at most <weight> many (gapped) q-grams
		return 
		sequenceLength(seqNo, host(pattern)) - length(indexShape(host(pattern))) + 1 
		- errors * weight(indexShape(host(pattern)));
	}
	
	template <typename TFinder, typename TIndex, typename TSpec>
	inline bool 
	_nextNonRepeatRange(TFinder &finder,
						Pattern<TIndex,  QGramsLookup<TSpec>  > const &pattern
	){
		typedef typename TFinder::TRepeat		TRepeat;
		typedef typename Value<TRepeat>::Type	TPos;
		
		if (finder.curRepeat == finder.endRepeat) return false;
		
		do 
		{
			finder.startPos = (*finder.curRepeat).endPosition;
			if (++finder.curRepeat == finder.endRepeat) 
			{
				finder.endPos = length(host(finder));
				if (finder.startPos + length(pattern.shape) > finder.endPos)
					return false;
				else
					break;
			} else
				finder.endPos = (*finder.curRepeat).beginPosition;
			// repeat until the shape fits in non-repeat range
		} while (finder.startPos + length(pattern.shape) > finder.endPos);
		
		finder.curPos = finder.startPos;
		hostIterator(finder) = begin(host(finder)) + finder.startPos;
		finder.haystackEnd = begin(host(finder)) + (finder.endPos - length(pattern.shape) + 1);
		
		return true;
	}
	
	template <typename TFinder, typename TIndex, typename TSpec>
	inline bool 
	_firstNonRepeatRange(TFinder &finder,
						 Pattern<TIndex, QGramsLookup<TSpec> > const &pattern
	){
		typedef typename TFinder::TRepeat		TRepeat;
		typedef typename Value<TRepeat>::Type	TPos;
		
		finder.curRepeat = begin(finder.data_repeats, Rooted());
		finder.endRepeat = end(finder.data_repeats, Rooted());
		
		if (finder.curRepeat == finder.endRepeat)
			finder.endPos = length(host(finder));
		else
			finder.endPos = (*finder.curRepeat).beginPosition;
		
		if (length(pattern.shape) > finder.endPos)
			return _nextNonRepeatRange(finder, pattern);
		
		finder.curPos = finder.startPos = 0;
		hostIterator(finder) = begin(host(finder));
		finder.haystackEnd = begin(host(finder)) + (finder.endPos - length(pattern.shape) + 1);
		
		return true;
	}
	
	//____________________________________________________________________________
	

	template< typename THaystack, typename TIndex, typename TShape, typename TSpec>
	inline void setPattern(Finder<THaystack, QGramsLookup<TShape, TSpec> >		&finder,
						   Pattern<TIndex,  QGramsLookup<TShape, TSpec> > const &pattern
	){
		finder.hasShape = true;
		finder.shape = pattern.shape;
	}
	
	template <typename THaystack, typename TIndex, typename TShape, typename TSpec>
	inline bool 
	find(
		 Finder<THaystack,  QGramsLookup<TShape, TSpec> >		&finder,
		 Pattern<TIndex,  QGramsLookup<TShape, TSpec> > const	&pattern
	){
		typedef	typename Value<TShape>::Type				THashValue;
		
		if (empty(finder)){
			// init pattern
			setPattern(finder, pattern);
			
			// init finder
			_finderSetNonEmpty(finder);
			
			finder.dotPos = 100000;
			finder.dotPos2 = 10 * finder.dotPos;
			
			if (!_firstNonRepeatRange(finder, pattern)) return false;
			if (_seedMultiProcessQGram(finder, pattern, hash(finder.shape, hostIterator(hostIterator(finder))))) {
				return true;
			}
			
		} else {
			if (++finder.curHit != finder.endHit) {
				return true;
			}
		}
		
		// all previous matches reported -> search new ones
		clear(finder.hits);
		
		// are we at the end of the text?
		if (atEnd(finder) && finder.curRepeat == finder.endRepeat){
			finder.curHit = finder.endHit;
			return false;
		}
		
		// find next hit
		do{
			if (atEnd(++finder)){
				if (!_nextNonRepeatRange(finder, pattern)){
					return false;
				}
				hash(finder.shape, hostIterator(hostIterator(finder)));
			} else {
				++finder.curPos;
				hashNext(finder.shape, hostIterator(hostIterator(finder)));
			}
			
			if (_seedMultiProcessQGram(finder, pattern, value(finder.shape))){
				return true;
			}
			
		} while (true);
	}
	
	
	template <typename TQGramHit, typename TText, typename TShape>
	inline typename Infix<TText>::Type
	hitInfix(TQGramHit const &hit, TText &text, TShape &shape)
	{
		__int64 hitBegin = hit.hstkPos;
		__int64 hitEnd = hit.hstkPos + weight(shape);
		__int64 textEnd = length(text);
		
		if (hitBegin < 0) hitBegin = 0;
		if (hitEnd > textEnd) hitEnd = textEnd;
		return infix(text, hitBegin, hitEnd);
	}
	
	template <typename THaystack, typename TSpec>
	inline typename Infix<THaystack>::Type
	infix(Finder<THaystack, QGramsLookup<TSpec> > &finder)
	{
		return hitInfix(*finder.curHit, haystack(finder), finder.shape);
	}

	template <typename THaystack, typename TSpec, typename TText>
	inline typename Infix<TText>::Type
	infix(Finder<THaystack, QGramsLookup<TSpec> > &finder, TText &text)
	{
		return hitInfix(*finder.curHit, text);
	}

	
	//____________________________________________________________________________
	
	template <typename TIndex, typename TSpec, typename TQGramHit, typename TText>
	inline typename Infix<TText>::Type
	infix(Pattern<TIndex, QGramsLookup<TSpec> > const & pattern, TQGramHit const &hit, TText &text)
	{
		__int64 hitBegin = hit.ndlPos;
		__int64 hitEnd = hit.ndlPos + weight(pattern.shape);
		__int64 textLength = sequenceLength(hit.ndlSeqNo, needle(pattern));
		
		if (hitEnd > textLength) hitEnd = textLength;
		if (hitBegin < 0) hitBegin = 0;
		
		return infix(text, hitBegin, hitEnd);
	}
	
	template <typename TIndex, typename TSpec, typename TQGramHit>
	inline typename Infix< typename GetSequenceByNo< TIndex const >::Type >::Type
	infix(Pattern<TIndex, QGramsLookup<TSpec> > const & pattern, TQGramHit const &hit)
	{
		return infix(pattern, hit, getSequenceByNo(hit.ndlSeqNo, needle(pattern)));
	}
		
	//////////////////////////////////////////////////////////////////////
	// 
	template <
		typename TFinder,
		typename TIndex,
		typename THashValue,
		typename TSpec
	>
	inline bool _seedMultiProcessQGram(TFinder & finder,
									   Pattern<TIndex, QGramsLookup<TSpec> > const & pattern,
									   THashValue hash)
	{
		typedef Pattern<TIndex,  QGramsLookup<TSpec>  >				TPattern;
		typedef typename Size<TIndex>::Type							TSize;
		typedef typename Fibre<TIndex, QGramSA>::Type				TSA;
		typedef typename Iterator<TSA, Standard>::Type				TSAIter;
		typedef typename TFinder::TQGramHit							THit;
		
		TIndex const &index = host(pattern);
		
		// create an iterator over the positions of the q-gram occurences in pattern
		TSAIter saBegin = begin(indexSA(index), Standard());
		TSAIter occ = saBegin + indexDir(index)[getBucket(index.bucketMap, hash)];
		TSAIter occEnd = saBegin + indexDir(index)[getBucket(index.bucketMap, hash) + 1];
		Pair<unsigned> ndlPos;

		// iterate over all q-gram occurences and do the processing
		for(; occ != occEnd; ++occ)
		{
			posLocalize(ndlPos, *occ, stringSetLimits(index)); // get pair of SeqNo and Pos in needle
			// begin position of the diagonal of q-gram occurence in haystack (possibly negative)
			__int64 diag = finder.curPos;
			diag -= getSeqOffset(ndlPos);
			
			// create a new hit and append it to the finders hit list
			THit hit = {                //                              
				finder.curPos,          // begin in haystack      
				getSeqNo(ndlPos),       // needle seq. number            
				getSeqOffset(ndlPos),	// needle position
				diag					// the diagonal
//				diag,					// the diagonal
//				weight(pattern.shape)
			};  
			
			// append the hit to the finders hit list
			appendValue(finder.hits, hit);
			
		}
		
		finder.curHit = begin(finder.hits, Rooted());
		finder.endHit = end(finder.hits, Rooted());
		
		return !empty(finder.hits);
	}
	
	//////////////////////////////////////////////////////////////////////////////
	
} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_INDEX_QGRAMS_FIND_H
