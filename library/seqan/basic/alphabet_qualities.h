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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Definitions for piggybacking qualities in free bits of bytes.
// ==========================================================================

#ifndef SEQAN_BASIC_ALPHABET_QUALITIES_H_
#define SEQAN_BASIC_ALPHABET_QUALITIES_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Concept.Alphabet With Qualities
..summary:An alphabet where qualities can be attached to the characters.

.Metafunction.QualityValueSize.concept:Concept.Alphabet With Qualities
.Metafunction.HasQualities.concept:Concept.Alphabet With Qualities
.Function.getQualityValue.concept:Concept.Alphabet With Qualities
.Function.convertQuality.concept:Concept.Alphabet With Qualities
.Function.assignQualityValue.concept:Concept.Alphabet With Qualities
.Function.assignQualities.concept:Concept.Alphabet With Qualities
 */

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction QualityValueSize
// ----------------------------------------------------------------------------

/**
.Metafunction.QualityValueSize
..cat:Alphabets
..summary:Return the number of quality values in characters from alphabet with qualities.
..signature:QualityValueSize<TAlphabet>::VALUE
..param.TAlphabet:The alphabet to query for its value size.
..returns:The cardinality of the set of qualities.
..include:seqan/basic.h
 */

// TODO(holtgrew): Do we want a default specialization?
template <typename TValue>
struct QualityValueSize
{
    enum { VALUE = ValueSize<TValue>::VALUE };
};

template <typename TValue>
struct QualityValueSize<TValue const> : QualityValueSize<TValue> {};

/**
.Metafunction.HasQualities
..cat:Alphabets
..summary:Return whether the given type stores qualities besides the alphabet.
..signature:HasQualities<TAlphabet>::VALUE
..signature:HasQualities<TAlphabet>::Type
..param.TAlphabet:The alphabe to query.
..returns:$true$, $false$, $True$, or $False$.
..include:seqan/basic.h
*/

template <typename TValue>
struct HasQualities
{
    enum { VALUE = false };
    typedef False Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getQualityValue()
// ----------------------------------------------------------------------------

/**
.Function.getQualityValue
..cat:Alphabets
..signature:getQualityValue(c)
..summary:Returns the quality of a character from an alphabet with integrated quality.
..param.c:Character to retrieve the quality from.
..returns:Quality value of $c$.
...type:nolink:int
...remarks:The quality value is an integral value between 0 and 62 (inclusive).
..see:Function.assignQualityValue
..see:Function.convertQuality
..include:seqan/basic.h
 */

// ----------------------------------------------------------------------------
// Function convertQuality()
// ----------------------------------------------------------------------------

// TODO(holtgrew): This could use some thought, what about other scales?

/**
.Function.convertQuality
..cat:Alphabets
..signature:convertQuality(c, q)
..summary:Convert an integer quality value into its ASCII representation for FASTQ (Phred scale).
..param.c:Character to store the quality in.
...type:nolink:$char$
..param.q:Value of the quality to convert.
...remarks:The quality value is an integral value between 0 and 62 (inclusive).
...type:nolink:$int$
..see:Function.getQualityValue
..include:seqan/basic.h
 */

inline 
void convertQuality(Ascii & c, int q) 
{
    c = '!' + Ascii(q);
}

// ----------------------------------------------------------------------------
// Function assignQualityValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): What about different quality types? Guess scaling? Look at how other packages do this.
/**
.Function.assignQualityValue
..cat:Alphabets
..signature:assignQualityValue(c, q)
..summary:Assign quality to a character from an alphabet with integrated quality.
..param.c:Target character to assign quality to.
..param.q:Quality to assign to the character.
...type:nolink:int
...type:nolink:char
...remarks:The quality value is an integral value between 0 and 62 (inclusive).
..remarks:If $q$ is a $char$ then $'!'$ is subtracted from $q$. This is useful for ASCII encoded PHRED scores.
..see:Function.getQualityValue
..see:Function.convertQuality
..include:seqan/basic.h
 */

// ----------------------------------------------------------------------------
// Function assignQualities()
// ----------------------------------------------------------------------------

/**
.Function.assignQualities
..cat:Alphabets
..summary:Assign quality values between strings.
..signature:assignQualities(target, source)
..param.target:Target string
...type:nolink:@Class.String@ of any alphabet with qualities, e.g. @Spec.DnaQ@, @Spec.Dna5Q@
..param.source:Source string.
...type:nolink:@Class.String@ of $int$ or $char$.
..remarks:This funciton calls @Function.assignQualityValue@ for all entries of $target$ and $source$, look at the documentation of @Function.assignQualityValue@ on how the values of $source$ are interpreted.
..see:Function.assignQualityValue
..include:seqan/basic.h
 */
template <typename TDest, typename TSource>
void assignQualities(TDest &dst, TSource const &src)
{
    typedef typename Iterator<TDest>::Type TDestIter;
    typedef typename Iterator<TSource>::Type TSourceIter;

    TDestIter itDst = begin(dst, Standard());
    TDestIter itDstEnd = end(dst, Standard());
    TSourceIter itSrcEnd = end(src, Standard());
    
    for (TSourceIter itSrc = begin(src, Standard()); itDst != itDstEnd && itSrc != itSrcEnd; ++itDst, ++itSrc)
        assignQualityValue(*itDst, *itSrc);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_ALPHABET_QUALITIES_H_
