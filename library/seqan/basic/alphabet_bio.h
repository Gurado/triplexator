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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Alphabet concepts stemming from biological applications.
// ==========================================================================

#ifndef SEQAN_BASIC_ALPHABET_BIO_H_
#define SEQAN_BASIC_ALPHABET_BIO_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Concept.Alphabet With Gaps
..summary:An alphabet that includes a specific gap character.

.Function.gapValue.concept:Concept.Alphabet With Gaps
.Function.gapValueImpl.concept:Concept.Alphabet With Gaps
 */

/**
.Concept.Alphabet With Unknown Value
..summary:An alphabet which includes a specific "unknown" character.

.Function.unknownValue.concept:Concept.Alphabet With Unknown Value
.Function.unknownValueImpl.concept:Concept.Alphabet With Unknown Value
 */

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function gapValueImpl
// ----------------------------------------------------------------------------

/**
.Function.gapValueImpl
..hidefromindex
..cat:Alphabets
..cat:Alignments
..summary:Implements @Function.gapValue@.
..signature:gapValueImpl(valuePointerTag)
..param.valuePointerTag:A pointer that is used as a tag to specify the value type.
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A gap character.
..see:Function.gapValue
..remarks.text:This function implements @Function.gapValue@.
It is recommended to use @Function.gapValue@ rather than $gapValueImpl$.
..include:seqan/basic.h
*/

template <typename T>
inline T
gapValueImpl(T *)
{
    SEQAN_CHECKPOINT;
    static T const _gap = T();
    return _gap;
}

// ----------------------------------------------------------------------------
// Function gapValue
// ----------------------------------------------------------------------------

/**
.Function.gapValue
..cat:Alphabets
..cat:Alignments
..summary:Return the "gap" value from an alphabet.
..signature:gapValue<T>()
..param.T:The alphabet type to query the "gap" value from.
...type:Concept.Alphabet With Gaps
..returns:The gap character.
..remarks.text:The function is implemented in @Function.gapValueImpl@.
Do not specialize $gapValue$, specialize @Function.gapValueImpl@ instead!
..see:Function.gapValueImpl
..include:seqan/basic.h
 */

template <typename T>
inline T
gapValue()
{
    SEQAN_CHECKPOINT;
    static T * _tag = 0;
    return gapValueImpl(_tag);
}

// ----------------------------------------------------------------------------
// Function unknownValueImpl
// ----------------------------------------------------------------------------

/**
.Function.unknownValueImpl
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.unknownValue@.
..signature:gapValueImpl(valuePointerTag)
..param.valuePointerTag:A pointer that is used as a tag to specify the value type.
...remarks:The pointer needs not to point to a valid object, so it is possible to use a null pointer here.
..returns:A "unknown" character.
..see:Function.unknownValue
..remarks.text:This function implements @Function.unknownValue@.
It is recommended to use @Function.gapValue@ rather than $gapValueImpl$.
..include:seqan/basic.h
*/

template <typename T>
inline T
unknownValueImpl(T *)
{
    SEQAN_CHECKPOINT;
    return 'N';
}

// ----------------------------------------------------------------------------
// Function unknownValue
// ----------------------------------------------------------------------------

/**
.Function.unknownValue
..cat:Alphabets
..summary:Return the "unknown" value from an alphabet.
..signature:unknownValue<T>()
..param.T:The alphabet type to query the "unknown" value from.
...type:Concept.Alphabet With Unknown Value
..returns:The "unknown" value.
 */

template <typename T>
inline T
unknownValue()
{
    SEQAN_CHECKPOINT;
    static T * _tag = 0;
    return unknownValueImpl(_tag);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_ALPHABET_BIO_H_
