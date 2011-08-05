// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_REPEAT_BASE_H
#define SEQAN_HEADER_REPEAT_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Class.Repeat
..summary:Store information about a repeat.
..cat:Index
..signature:Repeat<TPos, TPeriod>
..param.TPos:Type to use for storing positions.
...metafunction:Metafunction.Value
..param.TPeriod:Type to use for storing the repeat period.
...metafunction:Metafunction.Size
..include:seqan/index.h
..see:Function.findRepeats

.Memvar.Repeat#beginPosition
..summary:The begin position of the repeat of type $TPos$.
..class:Class.Repeat

.Memvar.Repeat#endPosition
..summary:The end position of the repeat of type $TPos$.
..class:Class.Repeat

.Memvar.Repeat#period
..summary:The period of the repeat of type $TSize$.
..class:Class.Repeat
 */

	template <typename TPos, typename TPeriod>
	struct Repeat {
		TPos		beginPosition;
		TPos		endPosition;
		TPeriod		period;
	};

	template <typename TPos, typename TPeriod>
	struct Value< Repeat<TPos, TPeriod> > {
		typedef TPos Type;
	};

	template <typename TPos, typename TPeriod>
	struct Size< Repeat<TPos, TPeriod> > {
		typedef TPeriod Type;
	};


	template <typename TSize>
	struct RepeatFinderParams {
		TSize minRepeatLen;
		TSize maxPeriod;
	};

	// custom TSpec for our customized wotd-Index
	struct TRepeatFinder;

	template <typename TText>
	struct Cargo<Index<TText, IndexWotd<TRepeatFinder> > > 
	{
		typedef Index<TText, IndexWotd<TRepeatFinder> >	TIndex;
		typedef typename Size<TIndex>::Type					TSize;
		typedef RepeatFinderParams<TSize>					Type;
	};


	// node predicate
	template <typename TText, typename TSpec>
	bool nodePredicate(Iter<Index<TText, IndexWotd<TRepeatFinder> >, TSpec> &it) 
	{
//		return countOccurrences(it) * nodeDepth(it) >= cargo(container(it)).minRepeatLen;
		return countOccurrences(it) * repLength(it) >= cargo(container(it)).minRepeatLen;
	}

	// monotonic hull
	template <typename TText, typename TSpec>
	bool nodeHullPredicate(Iter<Index<TText, IndexWotd<TRepeatFinder> >, TSpec> &it) 
	{
//		return nodeDepth(it) <= cargo(container(it)).maxPeriod;
		return repLength(it) <= cargo(container(it)).maxPeriod;
	}

	template <typename TPos>
	struct RepeatLess_ : public ::std::binary_function<TPos, TPos, bool>
	{
		// key less
		inline bool operator() (TPos const &a, TPos const &b) {
			return posLess(a, b);
		}
	};

	template <typename TValue>
	inline bool _repeatMaskValue(TValue) 
	{
        // TODO(holtgrew): Maybe use unknownValue<TValue>() instead of specializing for all alphabets, especially since we have Rna5 now and might want Rna5Q later.
		return false;
	}

	template <>
	inline bool _repeatMaskValue(Dna5 val) 
	{
		return val.value == 4; // 'N'
	}

	template <>
	inline bool _repeatMaskValue(Dna5Q val) 
	{
        static const Dna5Q n = 'N';
		return val.value == n.value; // 'N'
	}

	template <>
	inline bool _repeatMaskValue(Iupac val) 
	{
		static const Iupac n = 'N';
		return val.value == n.value;
	}
/*
	template <>
	inline bool _repeatMaskValue(AminoAcid val) 
	{
		return val == 'X';
	}
*/
/**
.Function.findRepeats
..summary:Search for repeats in a text.
..cat:Index
..signature:findRepeats(repeatString, text, minRepeatLength[, maxPeriod])
..param.repeatString:A @Class.String@ of @Class.Repeat@ objects.
..param.text:The text to search repeats in.
...type:Class.String
...type:Class.StringSet
..param.minRepeatLength:The minimum length each reported repeat must have.
..param.maxPeriod:Optionally, the maximal period that reported repeats can have.
..remarks:Subsequences of undefined values/$N$s will always be reported.
..include:seqan/index.h
..see:Class.Repeat
 */

	// period-1 optimization
	template <typename TRepeatStore, typename TString, typename TRepeatSize>
	inline void findRepeats(TRepeatStore &repString, TString const &text, TRepeatSize minRepeatLen) 
	{
		typedef typename Value<TRepeatStore>::Type	TRepeat;
		typedef typename Iterator<TString const>::Type	TIterator;
		typedef typename Value<TString>::Type		TValue;
		typedef typename Size<TString>::Type		TSize;

#ifdef SEQAN_PARALLEL_FIXME
        if (length(text) > (TSize)(omp_get_thread_num() * 2 * minRepeatLen)) {
            // Parallel case.

            // NOTE(holtgrew): The minimum text length check above makes it impossible that more than two chunks are
            // required to form an otherwise too short repeat.

            // TODO(holtgrew): Load balancing? Probably not worth it.
            String<TSize> splitters;
            String<TRepeatStore> threadLocalStores;

            // Each threads finds repeats on its chunk in parallel.
            #pragma omp parallel 
            {
                // We have to determine the number of available threads at this point.  We will use the number of thread
                // local stores to determin the number of available threads later on.
                #pragma omp master
                {
                    computeSplitters(splitters, length(text), omp_get_num_threads());
                    resize(threadLocalStores, omp_get_num_threads());
                }
                #pragma omp barrier

                int const t = omp_get_thread_num();
                TRepeatStore & store = threadLocalStores[t];

                TRepeat rep;
                rep.beginPosition = 0;
                rep.endPosition = 0;
                rep.period = 1;

                // Flags used for force-adding repeats for the chunks that have a
                // left/right neighbour.
                bool forceFirst = t > 0;
                bool forceLast = (t + 1) < omp_get_num_threads();

                TIterator it = iter(text, splitters[t], Standard());
                TIterator itEnd = iter(text, splitters[t + 1], Standard());
                if (it != itEnd)
                {
                    TValue last = *it;
                    TSize repLeft = 0;
                    TSize repRight = 1;

                    for (++it; it != itEnd; ++it, ++repRight) 
                    {
                        if (*it != last)
                        {
                            if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen || forceFirst)
                            {
                                forceFirst = false;
                                // insert repeat
                                rep.beginPosition = splitters[t] + repLeft;
                                rep.endPosition = splitters[t] + repRight;
                                appendValue(store, rep);
                            }
                            repLeft = repRight;
                            last = *it;
                        }
                    }
                    if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen || forceLast)
                    {
                        // Insert repeat but only if it is not already in there.
                        if (rep.beginPosition != repLeft && rep.endPosition != repRight) {
                            rep.beginPosition = splitters[t] + repLeft;
                            rep.endPosition = splitters[t] + repRight;
                            appendValue(store, rep);
                        }
                    }
                }
            }

            // Mend the splice points.
            //
            String<Pair<TSize, TSize> > fromPositions;
            resize(fromPositions, length(threadLocalStores));
            fromPositions[0].i1 = 0;
            fromPositions[0].i2 = length(threadLocalStores[0]);
            String<typename Size<TRepeatStore>::Type> outSplitters;
            appendValue(outSplitters, 0);
            int lastNonEmpty = 0;
            for (int i = 0; (unsigned)i < length(threadLocalStores) - 1; ++i) {
                fromPositions[i + 1].i1 = 0;
                fromPositions[i + 1].i2 = length(threadLocalStores[i + 1]);

                if (fromPositions[i].i1 != fromPositions[i].i2)
                    lastNonEmpty = i;

                // We merge the repeats at each split point if they are adjacent and their characters are equal and
                // their length's sums if greater than minRepeatLen.  Otherwise, we might have to remove left or right
                // repeat if the length is smaller than minRepeatLen.
                bool const adjacent = back(threadLocalStores[lastNonEmpty]).endPosition == front(threadLocalStores[i + 1]).beginPosition;
                bool const charsEqual = text[back(threadLocalStores[lastNonEmpty]).beginPosition] == text[front(threadLocalStores[i + 1]).beginPosition];
                bool const sumAboveThreshold = ((TRepeatSize)(front(threadLocalStores[i + 1]).endPosition - back(threadLocalStores[lastNonEmpty]).beginPosition) > minRepeatLen);
                bool const merge = adjacent && charsEqual && sumAboveThreshold;
                if (merge) {
                    fromPositions[i + 1].i1 += 1;
                    back(threadLocalStores[lastNonEmpty]).endPosition = front(threadLocalStores[i + 1]).endPosition;
                } else {
                    // Possibly remove left.
                    if ((TRepeatSize)(back(threadLocalStores[i]).endPosition - back(threadLocalStores[i]).beginPosition) <= minRepeatLen)
                        fromPositions[i].i2 -= 1;
                    // Possibly remove right.
                    if ((TRepeatSize)(front(threadLocalStores[i + 1]).endPosition - front(threadLocalStores[i + 1]).beginPosition) <= minRepeatLen)
                        fromPositions[i + 1].i1 += 1;
                }

                appendValue(outSplitters, back(outSplitters) + fromPositions[i].i2 - fromPositions[i].i1);
            }
            appendValue(outSplitters, back(outSplitters) + back(fromPositions).i2 - back(fromPositions).i1);

            // Allocate memory.
            clear(repString);
            resize(repString, back(outSplitters));

            // Copy back the repeats in parallel.
            unsigned nt = length(threadLocalStores);
            (void) nt;  // Otherwise, GCC 4.6 warns, does not see it used in pragma clause below.
            #pragma omp parallel num_threads(nt)
            {
                int const t = omp_get_thread_num();
                arrayCopy(iter(threadLocalStores[t], fromPositions[t].i1, Standard()),
                          iter(threadLocalStores[t], fromPositions[t].i2, Standard()),
                          iter(repString, outSplitters[t], Standard()));
            }
        } else {
#else  // #ifdef SEQAN_PARALLEL
            // Sequential case.
            TRepeat rep;
            rep.period = 1;
            clear(repString);

            TIterator it = begin(text, Standard());
            TIterator itEnd = end(text, Standard());
            if (it == itEnd) return;

            TValue last = *it;
            TSize repLeft = 0;
            TSize repRight = 1;

            for (++it; it != itEnd; ++it, ++repRight) 
            {
                if (*it != last)
                {
                    if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
                    {
                        // insert repeat
                        rep.beginPosition = repLeft;
                        rep.endPosition = repRight;
                        //					::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
                        appendValue(repString, rep);
                    }
                    repLeft = repRight;
                    last = *it;
                }
            }
            if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
            {
                // insert repeat
                rep.beginPosition = repLeft;
                rep.endPosition = repRight;
                //			::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
                appendValue(repString, rep);
            }
#endif  // #ifdef SEQAN_PARALLEL
#ifdef SEQAN_PARALLEL_FIXME
        }
#endif  // #ifdef SEQAN_PARALLEL
        // #pragma omp critical
        // {
        //     std::cerr << "thread #" << omp_get_thread_num() << " REPEATS:";
        //     for (unsigned i = 0; i < length(repString); ++i) {
        //         std::cerr << " (" << repString[i].beginPosition << ", " << repString[i].endPosition << ", " << repString[i].period << ")";
        //     }
        //     std::cerr << std::endl;
        // }
	}

    // TODO(holtgrew): Why for TString const and StringSet<> const?
	template <typename TRepeatStore, typename TString, typename TSpec, typename TRepeatSize>
	inline void findRepeats(TRepeatStore &repString, StringSet<TString, TSpec> const &text, TRepeatSize minRepeatLen) 
	{
		typedef typename Value<TRepeatStore>::Type	TRepeat;
		typedef typename Iterator<TString>::Type	TIterator;
		typedef typename Value<TString>::Type		TValue;
		typedef typename Size<TString>::Type		TSize;

		TRepeat rep;
		rep.period = 1;
		clear(repString);

		for (unsigned i = 0; i < length(text); ++i)
		{
			TIterator it = begin(text[i], Standard());
			TIterator itEnd = end(text[i], Standard());
			if (it == itEnd) continue;

			TValue last = *it;
			TSize repLeft = 0;
			TSize repRight = 1;
			rep.beginPosition.i1 = i;
			rep.endPosition.i1 = i;

			for (++it; it != itEnd; ++it, ++repRight) 
			{
				if (last != *it)
				{
					if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
					{
						// insert repeat
						rep.beginPosition.i2 = repLeft;
						rep.endPosition.i2 = repRight;
//						::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
						appendValue(repString, rep);
					}
					repLeft = repRight;
					last = *it;
				}
			}
			if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
			{
				// insert repeat
				rep.beginPosition.i2 = repLeft;
				rep.endPosition.i2 = repRight;
//				::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
				appendValue(repString, rep);
			}
		}
	}

	// main function
	template <typename TRepeatStore, typename TText, typename TRepeatSize, typename TPeriodSize>
	void findRepeats(TRepeatStore &repString, TText const &text, TRepeatSize minRepeatLen, TPeriodSize maxPeriod) 
	{
		typedef Index<TText, IndexWotd<TRepeatFinder> >					TIndex;
		typedef typename Size<TIndex>::Type									TSize;
		typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type	TNodeIterator;
		typedef typename Fibre<TIndex, FibreSA>::Type const				TSA;
		typedef typename Infix<TSA>::Type									TOccString;
		typedef typename Iterator<TOccString>::Type							TOccIterator;

		typedef typename Value<TRepeatStore>::Type							TRepeat;
		typedef typename Value<TOccString>::Type							TOcc;

		typedef ::std::map<TOcc,TRepeat,RepeatLess_<TOcc> >					TRepeatList;

		if (maxPeriod < 1) return;
		if (maxPeriod == 1) 
		{
			findRepeats(repString, text, minRepeatLen);
			return;
		}

		TIndex		index(text);
		TRepeatList list;

		// set repeat finder parameters
		cargo(index).minRepeatLen = minRepeatLen;
		cargo(index).maxPeriod = maxPeriod;

		TNodeIterator nodeIt(index);
		TOccIterator itA, itB, itRepBegin, itEnd;
		TRepeat rep;
		for (; !atEnd(nodeIt); goNext(nodeIt))
		{
			if (isRoot(nodeIt)) continue;

			// get occurrences
			TOccString occ = getOccurrences(nodeIt);
			itA = begin(occ, Standard());
			itEnd = end(occ, Standard());
			itRepBegin = itB = itA;

			TSize repLen = repLength(nodeIt);		// representative length
			if ((TSize)minRepeatLen <= repLen) continue;

			TSize diff, period = 0;					// period of current repeat
			TSize repeatLen = 0;					// overall length of current repeat
			TSize minLen = minRepeatLen - repLen;	// minimum repeat length minus length of representative

			for (++itB; itB != itEnd; ++itB)
			{
				diff = posSub(*itB, *itA);
				if (diff != period || getSeqNo(*itA) != getSeqNo(*itB))
				{
					// is the repeat long enough?
					if (repeatLen >= minLen)
						// is the repeat self overlapping or connected?
						if (parentRepLength(nodeIt) < period && period <= repLen)
						{
							// insert repeat
							rep.beginPosition = *itRepBegin;
							rep.endPosition = posAdd(*itA, period);
							rep.period = period;
//							::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
							list.insert(::std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
						}
					itRepBegin = itA;
					period = diff;
					repeatLen = 0;
				}
				repeatLen += period;
				itA = itB;
			}

			// is the last repeat long enough?
			if (repeatLen >= minLen)
				// is the repeat self overlapping or connected?
				if (parentRepLength(nodeIt) < period && period <= repLen)
				{
					// insert repeat
					rep.beginPosition = *itRepBegin;
					rep.endPosition = posAdd(*itA, period);
					rep.period = period;
//					::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
					list.insert(::std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
				}
		}

		// copy low-complex regions to result string
        clear(repString);
		reserve(repString, list.size(), Exact());
		typename TRepeatList::const_iterator lit = list.begin();
		typename TRepeatList::const_iterator litEnd = list.end();
		for (TSize i = 0; lit != litEnd; ++lit, ++i)
			appendValue(repString, (*lit).second);
	}


}	// namespace seqan

#endif
