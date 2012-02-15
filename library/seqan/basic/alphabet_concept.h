// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Concept definitions for alphabets.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef CORE_INCLUDE_SEQAN_BASIC_ALPHABET_CONCEPT_H_
#define CORE_INCLUDE_SEQAN_BASIC_ALPHABET_CONCEPT_H_

namespace seqan {

// ============================================================================
// Concepts for generic alphabets
// ============================================================================

/**
.Concept.AlphabetConcept
..summary:Natural container value.
..remarks:Stream output operators are not shown in the function list below, but required.
..remarks:Comparison operators are not shown in the function list below, but required.

.Function.clear.concept:Concept.Aggregate
.Function.value.concept:Concept.Aggregate
.Function.assign.concept:Concept.Aggregate

..Metafunction..concept:Concept.Aggregate
..Metafunction.BitsPerValue.concept:Concept.AlphabetConcept
 */

// minimal requirements for the alphabet of a String class
SEQAN_CONCEPT_REFINE(AlphabetConcept, (TValue), (Assignable)(DefaultConstructible)(CopyConstructible))
{
    TValue val, val2;

    SEQAN_CONCEPT_USAGE(AlphabetConcept)
    {
        // assign must be available as an equivalent to '='
        assign(val, val2);
//      swap(val, val2);
    }
};


// a totally strict ordered alphabet
SEQAN_CONCEPT_REFINE(OrderedAlphabetConcept, (TValue), (AlphabetConcept)(Comparable))
{
    TValue val;

    SEQAN_CONCEPT_USAGE(OrderedAlphabetConcept)
    {
        // type consistency checks
        sameType(minValue(val), val);
        sameType(minValue<TValue>, val);
        sameType(MinValue<TValue>::VALUE, val);
        sameType(maxValue(val), val);
        sameType(maxValue<TValue>, val);
        sameType(MaxValue<TValue>::VALUE, val);

        // sanity checks
        SEQAN_STATIC_ASSERT_MSG(MinValue<TValue>::VALUE <= MaxValue<TValue>::VALUE, "Minimal alphabet value must be less or equal to the maximal value.");
        
        // 0 must be an element of the alphabet, as we want to be able
        // to initialize a TValue variable to omit uninitialized warnings.
        SEQAN_STATIC_ASSERT_MSG(MinValue<TValue>::VALUE <= static_cast<TValue>(0), "0 must be convertible to a valid alphabet value.");
        SEQAN_STATIC_ASSERT_MSG(static_cast<TValue>(0) <= MaxValue<TValue>::VALUE, "0 must be convertible to a valid alphabet value.");
    }
};

// a finite totally strict ordered alphabet
SEQAN_CONCEPT_REFINE(FiniteOrderedAlphabetConcept, (TValue), (OrderedAlphabetConcept))
{
    typedef typename Size<TValue>::Type TSize;

    TValue  val;
    TSize   size;

    SEQAN_CONCEPT_ASSERT((UnsignedIntegerConcept<TSize>));

    SEQAN_CONCEPT_USAGE(FiniteOrderedAlphabetConcept)
    {
        // a finite alphabet must be countable
        sameType(ordValue(val), size);
//      sameType(valueSize<T>(), size);
        sameType(ValueSize<TValue>::VALUE, size);

        // alphabet must be non-empty
        SEQAN_STATIC_ASSERT_MSG(static_cast<TSize>(0) < ValueSize<TValue>::VALUE, "Alphabet size be greater than zero.");
        
        // convert integer to alphabet value
        val = 0;
        val = size;
        TValue val2(0);
        TValue val3(size);
    }
};

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BASIC_ALPHABET_CONCEPT_H_
